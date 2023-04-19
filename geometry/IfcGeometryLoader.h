/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#include <unordered_map>
#include <vector>
#include <optional>
#include <glm/glm.hpp>

#include "../parsing/IfcLoader.h"
#include "../schema/IfcSchemaManager.h"

#include "representation/geometry.h"
#include "representation/IfcGeometry.h"
#include "representation/IfcCurve.h"

 // This class takes care of loading raw unprocessed geometry data from the IFC loader

namespace webifc::geometry
{

  class IfcGeometryLoader 
  {
  public:
    IfcGeometryLoader(const webifc::parsing::IfcLoader &loader, webifc::utility::LoaderErrorHandler &errorHandler,const webifc::schema::IfcSchemaManager &schemaManager,uint16_t circleSegments);
    std::array<glm::dvec3,2> GetAxis1Placement(const uint32_t expressID) const;
    glm::dmat3 GetAxis2Placement2D(const uint32_t expressID) const;
    glm::dmat4 GetLocalPlacement(const uint32_t expressID) const;
    glm::dvec3 GetCartesianPoint3D(const uint32_t expressID) const;
    glm::dvec2 GetCartesianPoint2D(const uint32_t expressID) const;
    glm::dvec3 GetVector(const uint32_t expressID) const;
    IfcProfile GetProfile(uint32_t expressID) const;
    IfcProfile3D GetProfile3D(uint32_t expressID) const;
    IfcCurve GetCurve(uint32_t expressID,uint8_t dimensions, bool edge=false) const;
    std::vector<glm::dvec3> ReadIfcCartesianPointList3D(const uint32_t expressID) const;
    std::vector<glm::dvec2> ReadIfcCartesianPointList2D(const uint32_t expressID) const;
    IfcCurve GetOrientedEdge(uint32_t expressID) const;
    IfcCurve GetEdge(uint32_t expressID) const;
    IfcBound3D GetBound(const uint32_t expressID) const;
    IfcCurve GetLoop(const uint32_t expressID) const;
    std::optional<glm::dvec4> GetColor(uint32_t expressID) const;
    IfcAlignment GetAlignment(uint32_t expressID, IfcAlignment alignment = IfcAlignment(), glm::dmat4 transform = glm::dmat4(1)) const;
    bool GetColor(const uint32_t expressID, const glm::dvec4 &outputColor) const; 
    const std::unordered_map<uint32_t, std::vector<uint32_t>> &GetRelVoids() const;
    const std::unordered_map<uint32_t, std::vector<uint32_t>> &GetRelVoidRels() const;
    const std::unordered_map<uint32_t, std::vector<uint32_t>> &GetRelAggregates() const;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> &GetStyledItems() const;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> &GetRelMaterials() const;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> &GetMaterialDefinitions() const;
    double GetLinearScalingFactor() const;
  private:
    IfcCurve GetAlignmentCurve(uint32_t expressID) const;
    IfcProfile GetProfileByLine(uint32_t lineID) const;
    IfcTrimmingSelect GetTrimSelect(uint32_t DIM, std::vector<uint32_t> &tapeOffsets) const;
    void ComputeCurve(uint32_t expressID, IfcCurve &curve, uint8_t dimensions, bool edge, int sameSense = -1, int trimSense = -1, IfcTrimmingArguments trim = {}) const;
    std::vector<IfcSegmentIndexSelect> ReadCurveIndices() const;
    const webifc::parsing::IfcLoader &_loader;
    webifc::utility::LoaderErrorHandler &_errorHandler;
    const webifc::schema::IfcSchemaManager &_schemaManager;
    const std::unordered_map<uint32_t, std::vector<uint32_t>> _relVoidRel;
    const std::unordered_map<uint32_t, std::vector<uint32_t>> _relVoids;
    const std::unordered_map<uint32_t, std::vector<uint32_t>> _relAggregates;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> _styledItems;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> _relMaterials;
    const std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> _materialDefinitions;
    double _linearScalingFactor = 1;
    double _squaredScalingFactor = 1;
    double _cubicScalingFactor = 1;
    double _angularScalingFactor = 1;
    uint16_t _circleSegments;
    std::unordered_map<uint32_t, std::vector<uint32_t>> PopulateRelVoidsMap();
    std::unordered_map<uint32_t, std::vector<uint32_t>> PopulateRelVoidsRelMap();
    std::unordered_map<uint32_t, std::vector<uint32_t>> PopulateRelAggregatesMap();
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> PopulateStyledItemMap();
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> PopulateRelMaterialsMap();
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t, uint32_t>>> PopulateMaterialDefinitionsMap();
    void ReadLinearScalingFactor();
    double ConvertPrefix(const std::string &prefix);
  };
  
}