#include "Logger.h"
#include <cstdarg>
#include <cstdio>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#if defined(__EMSCRIPTEN__)

#include <emscripten/em_asm.h>

#endif

namespace {

enum class Level { INFO, WARNING, ERROR };

/**
 * Deferred-mode state.
 *
 * `depth` is atomic because worker threads read it on every log call while
 * the owning thread writes it at scope open/close; the buffer itself is
 * guarded by `lock`. A non-atomic depth would be a data race even though
 * only one thread ever writes it.
 */
std::mutex bufferLock;
std::vector< std::pair< Level, std::string > > buffered;
int deferDepth = 0;

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

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        // Tested under the lock, not before it: a worker must not read
        // "not deferring" while the owning thread is mid-flush and then
        // emit into a dead scope.
        if ( deferDepth > 0 ) {

            buffered.emplace_back( level, message );
            return;
        }
    }

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

    std::lock_guard< std::mutex > guard( bufferLock );

    ++deferDepth;
}

Logger::DeferredScope::~DeferredScope() {

    {
        std::lock_guard< std::mutex > guard( bufferLock );

        if ( --deferDepth > 0 ) {

            // An inner scope closing: leave the buffer for the outer one, so
            // nesting does not flush from a thread that may not be the owner.
            return;
        }
    }

    flushDeferred();
}
