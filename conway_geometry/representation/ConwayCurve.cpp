
/* 
 * Ref: https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcCurve.cpp
 * */

// Curve Implementation of a Curve

#include "ConwayCurve.h"

#include <glm/glm.hpp>
#include <vector>

#include "../operations/geometryutils.h"

namespace conway::geometry {

glm::dvec2 IfcCurve::Get2d(size_t i) const {
  glm::dvec2 ret;
  ret.x = points.at(i).x;
  ret.y = points.at(i).y;
  return ret;
}

glm::dvec3 IfcCurve::Get3d(size_t i) const { return points.at(i); }

void IfcCurve::Add(glm::dvec3 pt) {
  if (points.empty())
    points.push_back(pt);
  else if (!equals(pt, points.back(), EPS_TINY))
    points.push_back(pt);
}

void IfcCurve::Add(glm::dvec2 pt) {
  glm::dvec3 point;
  point.x = pt.x;
  point.y = pt.y;
  point.z = 0;
  Add(point);
}

void IfcCurve::Invert() { std::reverse(points.begin(), points.end()); }

bool IfcCurve::IsCCW() const {
  double sum = 0;

  for (size_t i = 0; i < points.size(); i++) {
    glm::dvec3 pt1 = points.at((i - 1) % points.size());
    glm::dvec3 pt2 = points.at(i);

    sum += (pt2.x - pt1.x) * (pt2.y + pt1.y);
  }

  return sum < 0;
}
}  // namespace conway::geometry