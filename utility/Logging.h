/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#pragma once

#include <string>

namespace webifc::utility {
enum class LogLevel : int {
  LOG_LEVEL_DEBUG = 0,
  LOG_LEVEL_INFO,
  LOG_LEVEL_WARN,
  LOG_LEVEL_ERROR,
  LOG_LEVEL_OFF
};

static LogLevel LOG_LEVEL = LogLevel::LOG_LEVEL_ERROR;

void setLogLevel(const int level);
void setLogLevel(const LogLevel level);

namespace log {
void debug(const std::string& msg);
void info(const std::string& msg);
void warn(const std::string& msg);
void error(const std::string& msg);
void log(const std::string& msg, const LogLevel& level);
}  // namespace log
}  // namespace webifc::utility
