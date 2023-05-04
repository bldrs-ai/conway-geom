/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#include "LoaderError.h"

#include <string>
#include <vector>

namespace webifc::utility {

void LoaderErrorHandler::ReportError(const LoaderErrorType t,
                                     const std::string m, const uint32_t e,
                                     const uint32_t type) {
  log::error(m);
  _errors.emplace_back(t, m, e, type);
}

void LoaderErrorHandler::ClearErrors() { _errors.clear(); }

const std::vector<LoaderError>& LoaderErrorHandler::GetErrors() const {
  std::vector<LoaderError> output(_errors);
  return _errors;
}

LoaderError::LoaderError(const LoaderErrorType t, const std::string m,
                         const uint32_t e, const uint32_t type)
    : type(t), message(m), expressID(e), ifcType(type) {}
}  // namespace webifc::utility
