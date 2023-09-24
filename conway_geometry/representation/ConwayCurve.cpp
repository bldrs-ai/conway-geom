
/*
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcCurve.cpp
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

size_t IfcCurve::GetPointsSize() const {
  if (points.empty()) {
    return 0;
  } else {
    return points.size();
  }
}

glm::dvec3 IfcCurve::Get3d(size_t i) const { return points.at(i); }

void IfcCurve::Add3d(glm::dvec3 pt) {
  if (points.empty())
    points.push_back(pt);
  else if (!equals(pt, points.back(), EPS_TINY))
    points.push_back(pt);
}

void IfcCurve::Add2d(glm::dvec2 pt) {
  glm::dvec3 point;
  point.x = pt.x;
  point.y = pt.y;
  point.z = 0;
  Add3d(point);
}

void IfcCurve::Invert() { std::reverse(points.begin(), points.end()); }

bool IfcCurve::IsCCW() const {
  double sum = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    const glm::dvec3& p1 = points[i];
    const glm::dvec3& p2 =
        points[(i + 1) % points.size()];  // Next point (wrapping around)
    sum += (p2.x - p1.x) * (p2.y + p1.y);
  }
  //printf("sum: %.6f\n", sum);
  return sum < 0;
}

/*bool IfcCurve::IsCCW() const
	{
		double sum = 0;

		for (size_t i = 0; i < points.size(); i++)
		{
			glm::dvec3 pt1 = points.at((i + 1) % points.size());
			glm::dvec3 pt2 = points.at(i);
			sum += (pt2.x - pt1.x) * (pt2.y + pt1.y);
		}

    printf("sum test: %.6f\n", sum);

		return sum < 0;
	}*/
}  // namespace conway::geometry
