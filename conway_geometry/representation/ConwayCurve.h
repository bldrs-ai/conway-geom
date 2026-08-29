
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

  /**
   * This curve is a closed polygon whose last point does NOT repeat its first.
   *
   * Set by GetLoop's point-list branch, where closure is a fact about the
   * entity (an IFCPOLYLOOP is a closed polygon by definition) rather than
   * something to be re-derived from the geometry. The edge-loop branch leaves
   * it false because those arrive with the head repeated as the tail, which is
   * self-describing.
   *
   * Consumers that need to know whether a bound closes must not measure it:
   * proximity provably cannot separate a closed loop from a coarsely sampled
   * open arc, because the closing gap is fixed by geometry while the sampling
   * scale is not. See reSolveClosedTrimHead and bldrs-ai/conway#655.
   */
  bool closedByConstruction = false;

  bool Add3d( const glm::dvec3& pt);
  bool Add2d( const glm::dvec2& pt);
  size_t GetPointsSize() const;
  glm::dvec2 Get2d(size_t i) const;
  glm::dvec3 Get3d(size_t i) const;
  void Invert();
  bool IsCCW() const;
  IfcCurve Clone() const;

  glm::dmat4 getPlacementAtDistance(double length);

  std::string DumpToSVG( const glm::dvec2& size, const glm::dvec2& offset ) const;

  std::string DumpToOBJ( const std::string& preamble ) const;

 private:
  static constexpr double EPS_TINY = 1e-9;
};

/**
 * Does a point-list loop enclose anything, i.e. should it be treated as a
 * closed polygon whose head is not repeated?
 *
 * The rule GetLoop applies, factored out so it can be tested without linking
 * ConwayGeometryProcessor (whose header reaches for draco and the GLTF SDK, so
 * the standalone tests cannot include it).
 *
 * Fewer than three points cannot enclose anything - a VERTEX_LOOP is one point
 * by design, the degenerate loop at a sphere pole or cone apex.
 */
inline bool pointListLoopIsClosedPolygon(
    const std::vector< glm::dvec3 >& points ) {

  return points.size() >= 3;
}

}  // namespace conway::geometry
