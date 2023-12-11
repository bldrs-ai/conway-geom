#include "Logger.h"
#include <cstdarg>
#include <cstdio>
#include <emscripten/em_asm.h>

void Logger::logInfo(const char* format, ...) {

    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logInfo'] === 'function') {
                globalScope['logInfo'](UTF8ToString($0));
            }
        },
        buffer);
}

void Logger::logWarning(const char* format, ...) {

    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logWarning'] === 'function') {
                globalScope['logWarning'](UTF8ToString($0));
            }
        },
        buffer);
}

void Logger::logError(const char* format, ...) {

    char buffer[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);

    EM_ASM(
        {
            let globalScope = (typeof window !== 'undefined') ? window : global;
            if (typeof globalScope['logError'] === 'function') {
                globalScope['logError'](UTF8ToString($0));
            }
        },
        buffer);
}
