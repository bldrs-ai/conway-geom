/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#pragma once

#include <string>
#include <vector>

#include "Logging.h"

namespace webifc::utility {

enum class LoaderErrorType {
  UNSPECIFIED,
  PARSING,
  BOOL_ERROR,
  UNSUPPORTED_TYPE
};

class LoaderError {
 public:
  LoaderErrorType type;
  std::string message;
  uint32_t expressID;
  uint32_t ifcType;

  LoaderError(const LoaderErrorType t = LoaderErrorType::UNSPECIFIED,
              const std::string m = "", const uint32_t e = 0,
              const uint32_t type = 0);
};

class LoaderErrorHandler {
 public:
  void ReportError(const LoaderErrorType t = LoaderErrorType::UNSPECIFIED,
                   const std::string m = "", const uint32_t e = 0,
                   const uint32_t type = 0);
  void ClearErrors();
  const std::vector<LoaderError> &GetErrors() const;

 private:
  std::vector<LoaderError> _errors;
};

}  // namespace webifc::utility
