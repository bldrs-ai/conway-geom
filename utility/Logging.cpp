/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#include "Logging.h"

#include <iostream>
#include <string>

namespace webifc::utility {

void setLogLevel(const int level) { LOG_LEVEL = static_cast<LogLevel>(level); }
void setLogLevel(const LogLevel level) { LOG_LEVEL = level; }

namespace log {
void debug(const std::string& msg) { log(msg, LogLevel::LOG_LEVEL_DEBUG); }
void info(const std::string& msg) { log(msg, LogLevel::LOG_LEVEL_INFO); }
void warn(const std::string& msg) { log(msg, LogLevel::LOG_LEVEL_WARN); }
void error(const std::string& msg) { log(msg, LogLevel::LOG_LEVEL_ERROR); }

void log(const std::string& msg, const LogLevel& level) {
  if (level >= LOG_LEVEL) {
    std::string fullMsg = msg;
    switch (level) {
      case LogLevel::LOG_LEVEL_DEBUG:
        fullMsg = "DEBUG: " + msg;
        break;
      case LogLevel::LOG_LEVEL_INFO:
        fullMsg = "INFO: " + msg;
        break;
      case LogLevel::LOG_LEVEL_WARN:
        fullMsg = "WARN: " + msg;
        break;
      case LogLevel::LOG_LEVEL_ERROR:
        fullMsg = "ERROR: " + msg;
        break;
      case LogLevel::LOG_LEVEL_OFF:
        return;
      default:
        return;
    }
    std::cout << fullMsg << std::endl;
  }
}
}  // namespace log

}  // namespace webifc::utility
