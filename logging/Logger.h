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
 * is flushed when it closes. Wrap any region that dispatches through
 * ThreadPool in one, ON THE THREAD THAT DISPATCHES, and its diagnostics
 * survive.
 *
 * Deferral is compiled in only for pthreads-enabled wasm, the one
 * configuration where the problem exists. Elsewhere the scope is inert and
 * logs emit immediately, so an abort inside the region cannot swallow a
 * buffer that was never filled.
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
     * them when the outermost scope closes.
     *
     * MUST be constructed on the thread that will dispatch the parallel work
     * - i.e. one that has a live JS sink - and never on a pool thread. Scopes
     * nest, but only on that same thread; a scope opened on a worker while
     * another is live would make whichever closes last the flusher, and a
     * worker flushing pushes the whole buffer through EM_ASM in a scope with
     * no sink, losing main-thread messages too. Asserted in debug builds.
     *
     * Ordering: messages emit in the order the buffer received them, so
     * main-thread messages raised inside the scope are held back with the
     * rest rather than jumping ahead of them. Ordering *between* worker
     * threads is whatever the race produced - inherent to concurrent
     * tessellation, not something this could fix.
     *
     * The buffer is fixed-size and never allocates, because logError is
     * called from catch handlers that run when memory has already failed.
     * Messages past its capacity are counted and reported as a single
     * "N further messages were dropped" line rather than lost silently, and
     * each stored message is truncated. A deferred log is a diagnostic, not
     * a transcript.
     */
    class DeferredScope {
    public:
        DeferredScope();
        ~DeferredScope();

        DeferredScope(const DeferredScope&) = delete;
        DeferredScope& operator=(const DeferredScope&) = delete;
    };
};

#endif // LOGGER_H
