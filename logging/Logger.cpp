#include "Logger.h"
#include <cassert>
#include <cstdarg>
#include <cstdio>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__EMSCRIPTEN__)

#include <emscripten/em_asm.h>

#endif

// Deferral only exists to work around EM_ASM evaluating in the calling
// pthread's scope. Where that cannot happen - native builds, and wasm built
// without pthreads - buffering would be pure downside: an abort or hang
// inside the deferred region (wasm OOM, a stack overflow in a triangulator,
// which is exactly when the per-face logError matters most) never runs the
// destructor, so the whole buffer is discarded where previously those lines
// had already been printed. Keep those builds on immediate emission.
#if defined(__EMSCRIPTEN__) && defined(__EMSCRIPTEN_PTHREADS__)
#define DEFER_LOGS 1
#endif

namespace {

enum class Level { INFO, WARNING, ERROR };

/**
 * Deferred-mode state. ALL of it - the depth included - is guarded by
 * bufferLock and nothing else.
 *
 * In particular deferDepth is a plain int, deliberately. Do not add an
 * unsynchronised `deferDepth == 0` fast path to route() to save a lock:
 * worker threads read it on every log call while the owning thread writes it
 * at scope open and close, so an unguarded read is a data race, and the lock
 * is what makes "test the depth and buffer the message" one atomic step.
 */
std::mutex bufferLock;
std::vector< std::pair< Level, std::string > > buffered;
int deferDepth = 0;

/**
 * The thread that opened the outermost live DeferredScope. Only meaningful
 * while deferDepth > 0.
 */
std::thread::id deferOwner;

/**
 * Emit one already-formatted message to the host.
 *
 * @param level Which sink to route to.
 * @param message The formatted text.
 */
void emitNow( Level level, const char* message ) {

#if defined(__EMSCRIPTEN__)

    // Three separate EM_ASM bodies rather than one parameterised by a JS
    // string: EM_ASM takes a literal, and passing the sink name through would
    // mean a second UTF8ToString per call on the hot path.
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

/**
 * Route one formatted message: buffer it if a DeferredScope is open,
 * otherwise emit it immediately.
 *
 * @param level Which sink to route to.
 * @param message The formatted text.
 */
void route( Level level, const char* message ) {

#if defined(DEFER_LOGS)

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        // Tested under the lock, not before it, so that "is a scope open"
        // and "append to its buffer" are one atomic step. See the note on
        // deferDepth above before adding a lock-free fast path here.
        if ( deferDepth > 0 ) {

            buffered.emplace_back( level, message );
            return;
        }
    }

#endif

    emitNow( level, message );
}

/** Format buffer size. Matches the previous per-function local. */
constexpr size_t MESSAGE_LIMIT = 1024;

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

void Logger::flushDeferred() {

    // Swap out under the lock and emit outside it. Emitting while holding it
    // would block every worker for the duration of a JS boundary crossing per
    // message, and on the native path for the duration of the printfs.
    std::vector< std::pair< Level, std::string > > pending;

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        pending.swap( buffered );
    }

    for ( const auto& entry : pending ) {

        emitNow( entry.first, entry.second.c_str() );
    }
}

Logger::DeferredScope::DeferredScope() {

#if defined(DEFER_LOGS)

    std::lock_guard< std::mutex > guard( bufferLock );

    if ( deferDepth++ == 0 ) {

        deferOwner = std::this_thread::get_id();
    }

    // Nesting is only meaningful on the owning thread. A scope opened on a
    // pool thread while another is live would make the LAST scope to close
    // the flusher, and if that is the worker it would push the whole buffer -
    // main-thread messages included - through EM_ASM in a scope with no sink,
    // which is the bug this class exists to fix.
    assert( deferOwner == std::this_thread::get_id() &&
            "Logger::DeferredScope must not be opened on a pool thread" );

#endif
}

Logger::DeferredScope::~DeferredScope() {

#if defined(DEFER_LOGS)

    std::vector< std::pair< Level, std::string > > pending;

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        if ( --deferDepth > 0 ) {

            // An inner scope closing: leave the buffer for the outer one.
            return;
        }

        // Claimed in the SAME critical section that drops the depth to zero.
        // Doing it in two steps would leave a window where the scope reads as
        // closed but its buffer is still pending, so a message raised in that
        // window would take the immediate path and land ahead of messages
        // logged before it.
        pending.swap( buffered );
    }

    for ( const auto& entry : pending ) {

        emitNow( entry.first, entry.second.c_str() );
    }

#endif
}
