
/*
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcCurve.cpp
 * */

// Curve Implementation of a Curve

#include "ConwayCurve.h"

#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <sstream>

#include "SVGContext.h"

//#include "../operations/geometryutils.h"

namespace conway::geometry {

std::string IfcCurve::DumpToOBJ() const {
  
  std::stringstream obj;

  for ( const glm::dvec3& t : points ) {
    
    obj << "v " << t.x << " " << t.y << " " << t.z << "\n";
  }

  if ( indices.size() > 1 ) {

      obj << "l";

      for ( uint16_t indice : indices ) {

          obj << " " << indice + 1;
      }

      obj << "\n";
  
  } else if ( points.size() > 1 ) {

      obj << "l";

      for ( size_t where = 1, end = points.size(); where <= end; ++where ) {

          obj << " " << where;
      }
      
      obj << "\n";
  }

  return obj.str();
}

std::string IfcCurve::DumpToSVG( const glm::dvec2& size, const glm::dvec2& offset ) const {
  
  glm::dvec2 min( FLT_MAX );
  glm::dvec2 max( -FLT_MAX );

  for ( const glm::dvec3& point : this->points ) {

    min = glm::min( min, glm::dvec2( point ) );
    max = glm::max( max, glm::dvec2( point ) );
  }

  SVGContext svg( size, offset, min, max );

  svg.header();

  if ( points.size() == 1 ) {

    svg.point( points[ 0 ] );
  
  } else if ( !indices.empty() ) {

    if ( indices.size() == 1 ) {
      
      svg.point( points[ indices[ 0 ] ] );
    } else {

      for ( size_t where = 0, end = indices.size() - 1; where < end; ++where ) {

        svg.line( points[ indices[ where ] ], points[ indices[ where + 1 ] ] );
      }
    }
  
  } else {

      for ( size_t where = 0, end = points.size() - 1; where < end; ++where ) {

        svg.line( points[ where ], points[ where + 1 ] );
      }
  }

  svg.trailer();

  return svg.str();
}

IfcCurve IfcCurve::Clone() const {
  return *this;
}

inline bool equals(glm::dvec3 A, glm::dvec3 B, double eps = 0) {
  return std::fabs(A.x - B.x) <= eps && std::fabs(A.y - B.y) <= eps &&
         std::fabs(A.z - B.z) <= eps;
}

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

bool IfcCurve::IsCCW() const
	{
		double sum = 0;
		auto n = points.size();

		for (size_t i = 0; i < n; i++)
		{
			glm::dvec3 pt1 = points.at((i + n - 1) % n);
			glm::dvec3 pt2 = points.at(i);

			sum += (pt2.x - pt1.x) * (pt2.y + pt1.y);
		}

		return sum < 0;
	}


  glm::dmat4 IfcCurve::getPlacementAtDistance(double length)
	{
		double totalDistance = 0;
		glm::dvec3 pos;
		glm::dvec3 vx = glm::dvec3(1, 0, 0);
		glm::dvec3 vy = glm::dvec3(0, 1, 0);
		glm::dvec3 vz = glm::dvec3(0, 0, 1);
		if (points.size() > 1)
		{
			for (uint32_t i = 0; i < points.size() - 1; i++)
			{
				double distance = glm::distance(points[i], points[i + 1]);
				totalDistance += distance;
				if (totalDistance >= length)
				{
					double factor = (totalDistance - length) / distance;
					pos = points[i] * factor + points[i + 1] * (1 - factor);
					glm::dvec3 tan = points[i + 1] - points[i];
					vx = glm::cross(tan, vz);
					vy = glm::cross(vx, tan);
					vz = tan;
					vx = glm::normalize(vx);
					vy = glm::normalize(vy);
					vz = glm::normalize(vz);

					break;
				}
			}
		}
		return glm::dmat4(
			glm::dvec4(vx, 0),
			glm::dvec4(vy, 0),
			glm::dvec4(vz, 0),
			glm::dvec4(pos, 1));
	}
}  // namespace conway::geometry
