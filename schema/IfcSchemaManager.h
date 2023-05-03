/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "ifc-schema.h"

namespace webifc::schema {
class IfcSchemaManager {
 public:
  IfcSchemaManager();
  const std::vector<IFC_SCHEMA> GetAvailableSchemas() const;
  std::string GetSchemaName(IFC_SCHEMA schema) const;
  uint32_t IfcTypeToTypeCode(const std::string name) const;
  uint32_t IfcTypeToTypeCode(const std::string_view name) const;
  std::string IfcTypeCodeToType(const uint32_t typeCode) const;
  bool IsIfcElement(const uint32_t typeCode) const;
  const std::unordered_set<uint32_t>& GetIfcElementList() const;

 private:
  std::vector<uint32_t> _crcTable;
  std::unordered_set<uint32_t> _ifcElements;
  std::vector<IFC_SCHEMA> _schemas;
  std::vector<std::string> _schemaNames;
  void initSchemaData();
  uint32_t IfcTypeToTypeCode(const void* name, const size_t len) const;
};
}  // namespace webifc::schema
