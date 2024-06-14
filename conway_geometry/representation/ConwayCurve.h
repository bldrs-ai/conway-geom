
/* 
 * Decoupling: https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref: https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcCurve.h
 * */

// Curve Implementation of a Curve

#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <string>

namespace conway::geometry {

struct IfcCurve {
  std::vector<glm::dvec3> points;
  std::vector<uint16_t> indices;
  void Add3d(glm::dvec3 pt);
  void Add2d(glm::dvec2 pt);
  size_t GetPointsSize() const;
  glm::dvec2 Get2d(size_t i) const;
  glm::dvec3 Get3d(size_t i) const;
  void Invert();
  bool IsCCW() const;
  glm::dmat4 getPlacementAtDistance(double length);

  std::string DumpToSVG( const glm::dvec2& size, const glm::dvec2& offset ) const;

  std::string DumpToOBJ() const;

 private:
  static constexpr double EPS_TINY = 1e-9;
};

}  // namespace conway::geometry
