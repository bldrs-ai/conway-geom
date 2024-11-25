#include "Logger.h"
#include <cstdarg>
#include <cstdio>

#if defined(__EMSCRIPTEN__)

#include <emscripten/em_asm.h>

#endif

void Logger::logInfo(const char* format, ...) {


    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

#if defined(__EMSCRIPTEN__)

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logInfo'] === 'function') {
                globalScope['logInfo'](UTF8ToString($0));
            }
        },
        buffer);

#else 

    printf("info   \t%s\n", buffer);

#endif
}

void Logger::logWarning(const char* format, ...) {

    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

#if defined(__EMSCRIPTEN__)

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logWarning'] === 'function') {
                globalScope['logWarning'](UTF8ToString($0));
            }
        },
        buffer);

#else 

    printf("warning\t%s\n", buffer);

#endif
}

void Logger::logError(const char* format, ...) {

    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

#if defined(__EMSCRIPTEN__)

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logError'] === 'function') {
                globalScope['logError'](UTF8ToString($0));
            }
        },
        buffer);

#else 

    printf("error  \t%s\n", buffer);

#endif
}
