#include "Logger.h"
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <mutex>

#if defined(__EMSCRIPTEN__)

#include <emscripten/em_asm.h>

#endif

// Deferral exists only to work around EM_ASM evaluating in the calling
// pthread's own scope, which is a problem exclusive to pthreads-enabled wasm.
// Native and single-threaded wasm compile the whole mechanism out.
#if defined(__EMSCRIPTEN__) && defined(__EMSCRIPTEN_PTHREADS__)
#define DEFER_LOGS 1
#include <emscripten/threading.h>
#endif

namespace {

enum class Level { INFO, WARNING, ERROR };

/** Format buffer size. Matches the previous per-function local. */
constexpr size_t MESSAGE_LIMIT = 1024;

/**
 * Emit one already-formatted message to the host.
 *
 * @param level Which sink to route to.
 * @param message The formatted text.
 */
void emitNow( Level level, const char* message ) {

#if defined(__EMSCRIPTEN__)

    // Three separate EM_ASM bodies rather than one parameterised by sink
    // name: EM_ASM takes a literal, and threading the name through would cost
    // a second UTF8ToString per call on the hot path.
    switch ( level ) {

    case Level::INFO:
        EM_ASM(
            {
                let globalScope =
                    (typeof globalThis !== 'undefined') ? globalThis :
                        ((typeof window !== 'undefined') ? window : global);
                if (typeof globalScope['logInfo'] === 'function') {
                    globalScope['logInfo'](UTF8ToString($0));
                }
            },
            message);
        break;

    case Level::WARNING:
        EM_ASM(
            {
                let globalScope =
                    (typeof globalThis !== 'undefined') ? globalThis :
                        ((typeof window !== 'undefined') ? window : global);
                if (typeof globalScope['logWarning'] === 'function') {
                    globalScope['logWarning'](UTF8ToString($0));
                }
            },
            message);
        break;

    case Level::ERROR:
        EM_ASM(
            {
                let globalScope =
                    (typeof globalThis !== 'undefined') ? globalThis :
                        ((typeof window !== 'undefined') ? window : global);
                if (typeof globalScope['logError'] === 'function') {
                    globalScope['logError'](UTF8ToString($0));
                }
            },
            message);
        break;
    }

#else

    switch ( level ) {

    case Level::INFO:    printf("info   \t%s\n", message); break;
    case Level::WARNING: printf("warning\t%s\n", message); break;
    case Level::ERROR:   printf("error  \t%s\n", message); break;
    }

#endif
}

#if defined(DEFER_LOGS)

/**
 * Deferred-mode state. ALL of it - the depth included - is guarded by
 * bufferLock and nothing else.
 *
 * deferDepth is a plain int, deliberately. Do not add an unsynchronised
 * `deferDepth == 0` fast path to route() to save a lock: worker threads read
 * it on every log call while the owning thread writes it at scope open and
 * close, so an unguarded read is a data race, and the lock is what makes
 * "test the depth and store the message" one atomic step.
 *
 * The buffer is a FIXED, preallocated array rather than a vector of strings.
 * That matters more than it looks: Logger::logError is called from the catch
 * handlers inside FinalizeStagedFaces's parallel_for, so it runs precisely
 * when something has already gone wrong - including wasm OOM, which surfaces
 * as std::bad_alloc under this build's -s ABORTING_MALLOC=0. A logger that
 * allocates could throw out of a catch handler, which on a pool thread is
 * std::terminate and on the calling thread escapes parallel_for while workers
 * still hold its function pointer. Storing into fixed storage cannot throw,
 * cannot allocate, and cannot grow without bound on a degenerate model.
 */
constexpr size_t MAX_DEFERRED = 128;

/** Truncation length for a deferred message, including the terminator. */
constexpr size_t DEFERRED_LIMIT = 256;

struct DeferredMessage {

    Level level;
    char  text[ DEFERRED_LIMIT ];
};

std::mutex      bufferLock;
DeferredMessage buffered[ MAX_DEFERRED ];
size_t          bufferedCount = 0;
size_t          droppedCount  = 0;
int             deferDepth    = 0;

/**
 * Copy one message into a slot, marking it if it did not fit.
 *
 * Callers format into 1024 bytes but slots hold 256, and several of the
 * messages that land here put coordinates ahead of the interesting part -
 * manifold_utils.h and mesh_utils.h format eight %f values before the
 * trailing exception text - so a silent cut removes exactly the words worth
 * reading. The marker at least says the text is short of the original.
 *
 * @param slot Destination.
 * @param level Severity.
 * @param message Formatted text, NUL-terminated.
 */
void fill( DeferredMessage& slot, Level level, const char* message ) {

    constexpr const char* ELLIPSIS = "...";
    constexpr size_t      ELLIPSIS_LENGTH = 3;

    slot.level = level;
    std::strncpy( slot.text, message, DEFERRED_LIMIT - 1 );
    slot.text[ DEFERRED_LIMIT - 1 ] = '\0';

    if ( std::strlen( message ) > DEFERRED_LIMIT - 1 ) {
        std::strcpy( slot.text + ( DEFERRED_LIMIT - 1 - ELLIPSIS_LENGTH ), ELLIPSIS );
    }
}


/**
 * Store a message, or account for it if the buffer is full.
 *
 * Overflow is NOT first-come-first-served, because that loses the wrong
 * message. A degenerate model raises a per-face warning for thousands of
 * faces, which would fill the buffer long before FinalizeStagedFaces's catch
 * handler reports the exception that actually stopped the work - the single
 * message this whole change exists to deliver. An ERROR arriving into a full
 * buffer therefore displaces the oldest non-error entry.
 *
 * Caller must hold bufferLock.
 *
 * @param level Severity.
 * @param message Formatted text, NUL-terminated.
 */
void store( Level level, const char* message ) {

    if ( bufferedCount < MAX_DEFERRED ) {

        fill( buffered[ bufferedCount++ ], level, message );
        return;
    }

    if ( level == Level::ERROR ) {

        for ( size_t where = 0; where < MAX_DEFERRED; ++where ) {

            if ( buffered[ where ].level != Level::ERROR ) {

                fill( buffered[ where ], level, message );
                ++droppedCount;
                return;
            }
        }
    }

    // Count rather than grow. The overflow is reported on flush, so a
    // truncated diagnostic still says it was truncated - silently dropping
    // the tail is the failure mode this whole change exists to remove.
    ++droppedCount;
}

#endif


/**
 * Route one formatted message: store it if a DeferredScope is open,
 * otherwise emit it immediately.
 *
 * @param level Which sink to route to.
 * @param message The formatted text.
 */
void route( Level level, const char* message ) {

#if defined(DEFER_LOGS)

    // Only a message with nowhere to go is buffered. The main runtime thread
    // owns the JS sink, so its EM_ASM reaches a listener whether or not a
    // scope is open — and buffering it would be a REGRESSION, not a fix:
    // parallel_for runs serially below minParallelBatch and the calling
    // thread joins the work, so those iterations used to print immediately
    // and would now vanish on a wasm trap mid-tessellation, which is exactly
    // when a diagnostic matters most.
    //
    // It also drops main-thread traffic out of the fixed buffer entirely,
    // leaving the capacity to the messages that actually need it.
    if ( !emscripten_is_main_runtime_thread() ) {

        std::lock_guard< std::mutex > guard( bufferLock );

        // Tested under the lock, not before it, so that "is a scope open" and
        // "store into its buffer" are one atomic step.
        if ( deferDepth > 0 ) {

            store( level, message );
            return;
        }

        // No scope open: fall through and emit, which on a worker means the
        // message is lost. That is the pre-existing behaviour this change
        // narrows rather than something it introduces — a parallel region
        // that can log needs a DeferredScope around it.
    }

#endif

    emitNow( level, message );
}

}  // namespace


void Logger::logInfo(const char* format, ...) {

    char buffer[MESSAGE_LIMIT];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    route( Level::INFO, buffer );
}

void Logger::logWarning(const char* format, ...) {

    char buffer[MESSAGE_LIMIT];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    route( Level::WARNING, buffer );
}

void Logger::logError(const char* format, ...) {

    char buffer[MESSAGE_LIMIT];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    route( Level::ERROR, buffer );
}

Logger::DeferredScope::DeferredScope() {

#if defined(DEFER_LOGS)

    // A real runtime guard, not an assert. genie.lua applies -DNDEBUG from an
    // unfiltered `configuration {"gmake"}` block, so it reaches the debug
    // emscripten config too - an assert here would be compiled out of every
    // build that defers, leaving the contract unenforced everywhere while the
    // header claimed otherwise.
    //
    // Refusing to open, rather than opening anyway, is what keeps the failure
    // safe: a scope opened on a pool thread would make whichever scope closes
    // last the flusher, and a worker flushing pushes the entire buffer -
    // main-thread messages included - through EM_ASM in a scope with no sink.
    // Inert here costs only the worker messages that a wrong scope could
    // never have delivered anyway.
    if ( !emscripten_is_main_runtime_thread() ) {

        emitNow(
            Level::ERROR,
            "Logger::DeferredScope opened off the main runtime thread; "
            "deferral disabled for it. Worker logs inside it will be lost." );
        return;
    }

    std::lock_guard< std::mutex > guard( bufferLock );

    ++deferDepth;
    opened_ = true;

#endif
}

Logger::DeferredScope::~DeferredScope() {

#if defined(DEFER_LOGS)

    if ( !opened_ ) {
        return;
    }

    DeferredMessage pending[ MAX_DEFERRED ];
    size_t          pendingCount = 0;
    size_t          pendingDrops = 0;

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        if ( --deferDepth > 0 ) {

            // An inner scope closing: leave the buffer for the outer one.
            return;
        }

        // Claimed in the SAME critical section that drops the depth to zero.
        // Two steps would leave a window where the scope reads as closed but
        // its buffer is still pending, so a message raised in that window
        // would take the immediate path and land ahead of messages logged
        // before it.
        pendingCount = bufferedCount;
        pendingDrops = droppedCount;

        for ( size_t where = 0; where < pendingCount; ++where ) {

            pending[ where ] = buffered[ where ];
        }

        bufferedCount = 0;
        droppedCount  = 0;
    }

    // Emitted outside the lock: a JS boundary crossing per message would
    // otherwise block every worker for the duration.
    for ( size_t where = 0; where < pendingCount; ++where ) {

        emitNow( pending[ where ].level, pending[ where ].text );
    }

    if ( pendingDrops > 0 ) {

        char overflow[ MESSAGE_LIMIT ];

        snprintf(
          overflow,
          sizeof( overflow ),
          "%zu further deferred message(s) were dropped (buffer holds %zu)",
          pendingDrops,
          MAX_DEFERRED );

        emitNow( Level::WARNING, overflow );
    }

#endif
}
