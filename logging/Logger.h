#ifndef LOGGER_H
#define LOGGER_H

class Logger {
public:
    static void logInfo(const char* format, ...);
    static void logWarning(const char* format, ...);
    static void logError(const char* format, ...);
};

#endif // LOGGER_H