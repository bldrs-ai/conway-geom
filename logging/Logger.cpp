#include "Logger.h"
#include <cassert>
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

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        // Tested under the lock, not before it, so that "is a scope open" and
        // "store into its buffer" are one atomic step.
        if ( deferDepth > 0 ) {

            if ( bufferedCount < MAX_DEFERRED ) {

                DeferredMessage& slot = buffered[ bufferedCount++ ];

                slot.level = level;
                std::strncpy( slot.text, message, DEFERRED_LIMIT - 1 );
                slot.text[ DEFERRED_LIMIT - 1 ] = '\0';

            } else {

                // Count rather than grow. The overflow is reported on flush,
                // so a truncated diagnostic still says it was truncated -
                // silently dropping the tail is the failure mode this whole
                // change exists to remove.
                ++droppedCount;
            }

            return;
        }
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

    // The contract is that the flush happens somewhere a JS sink exists, and
    // under emscripten that is knowable rather than a matter of convention.
    // Checking the runtime thread catches a scope opened on a pool thread,
    // which comparing against a previously-recorded owner cannot: at depth 0
    // that comparison is against the value just written, so it is a tautology
    // and passes for any thread.
    assert( emscripten_is_main_runtime_thread() &&
            "Logger::DeferredScope must be opened on the main runtime thread" );

    std::lock_guard< std::mutex > guard( bufferLock );

    ++deferDepth;

#endif
}

Logger::DeferredScope::~DeferredScope() {

#if defined(DEFER_LOGS)

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
