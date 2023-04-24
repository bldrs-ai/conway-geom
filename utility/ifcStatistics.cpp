#include <stdarg.h>
#include <stdio.h>

#include <cstdint>
#include <map>

namespace webifc::statistics {

bool exportObjs = false;
bool collectStats = false;
bool verboseStats = false;
bool shouldPrintTypeInfo = false;

std::map<uint32_t, uint32_t> uniqueTypeDefs;
int printTypeInfo(const char *format, ...) {
  va_list args;
  va_start(args, format);

  if (shouldPrintTypeInfo) vprintf(format, args);

  va_end(args);
}

void collectStatistics(uint32_t ifcLineType) {
  if (collectStats) {
    std::map<uint32_t, uint32_t>::iterator uniqueTypeDefIterator;
    uniqueTypeDefIterator = uniqueTypeDefs.find(ifcLineType);

    if (uniqueTypeDefIterator != uniqueTypeDefs.end()) {
      uniqueTypeDefIterator->second++;
    } else {
      uniqueTypeDefs.insert(std::pair<uint32_t, uint32_t>(ifcLineType, 1));
    }
  }
}
}  // namespace webifc::statistics