#pragma once
#ifndef LOGGER_H
#define LOGGER_H

/**
 * Emits diagnostics to the host.
 *
 * Under emscripten these calls cross into JS via EM_ASM, which evaluates in
 * *the calling pthread's own scope*. conway installs the sink
 * (globalThis.logWarning / logError / logInfo) on the main thread only, so a
 * log raised on a worker thread lands in a scope nobody is listening to and
 * disappears without a trace.
 *
 * DeferredScope is the way around that. Inside one, every log call - from any
 * thread - is formatted and buffered instead of emitted, and the whole buffer
 * is flushed from the scope's own thread when it closes. Wrap any region that
 * dispatches through ThreadPool in one and its diagnostics survive.
 *
 * The obvious alternative, MAIN_THREAD_EM_ASM, deadlocks: it proxies
 * synchronously to the main thread, and the main thread is blocked inside
 * parallel_for at exactly the moment a worker wants to log.
 * See https://github.com/bldrs-ai/conway/issues/466.
 */
class Logger {
public:
    static void logInfo(const char* format, ...);
    static void logWarning(const char* format, ...);
    static void logError(const char* format, ...);

    /**
     * Buffers every log raised while it is alive, from any thread, and emits
     * them on destruction from the thread that constructed it.
     *
     * Ordering: messages emit in the order the buffer received them, so
     * main-thread messages raised inside the scope are held back with the
     * rest rather than jumping ahead of them. Ordering *between* worker
     * threads is whatever the race produced - inherent to concurrent
     * tessellation, not something this could fix.
     *
     * Nesting is safe; only the outermost scope flushes.
     */
    class DeferredScope {
    public:
        DeferredScope();
        ~DeferredScope();

        DeferredScope(const DeferredScope&) = delete;
        DeferredScope& operator=(const DeferredScope&) = delete;
    };

    /**
     * Emit and clear anything buffered so far, without closing the scope.
     * Call from the thread that owns the scope, at a point where no worker
     * can be logging - between parallel regions, not inside one.
     */
    static void flushDeferred();
};

#endif // LOGGER_H
