#pragma once

#include <glm/glm.hpp>
#include "structures/winged_edge.h"

namespace conway::geometry {
  
  struct box3 {

    glm::dvec3 min = { DBL_MAX, DBL_MAX, DBL_MAX };
    glm::dvec3 max = { -DBL_MAX, -DBL_MAX, -DBL_MAX }; 

    void merge( const glm::dvec3& with ) {

      min = glm::min( with, min );
      max = glm::max( with, max );
    }

    void merge( const box3& with ) {

      min = glm::min( with.min, min );
      max = glm::max( with.max, max );
    }
  };

  inline bool overlaps( const box3& left, const box3& right ) {

    return
      left.max.x >= right.min.x &&
      left.min.x <= right.max.x &&
      left.max.y >= right.min.y &&
      left.min.y <= right.max.y &&
      left.max.z >= right.min.z &&
      left.min.z <= right.max.z;
  }

  
  inline bool overlaps( const box3& left, const glm::dvec3& right ) {

    return
      left.max.x >= right.x &&
      left.min.x <= right.x &&
      left.max.y >= right.y &&
      left.min.y <= right.y &&
      left.max.z >= right.z &&
      left.min.z <= right.z;
  }

  inline box3 merge( const box3& left, const box3& right ) {

    return {
      glm::min( left.min, right.min ),
      glm::max( left.max, right.max )
    };
  }

  inline box3 make_box( const WingedEdgeMesh< glm::dvec3 >& mesh, size_t triangleIndex, double tolerance = 2.0 * DBL_EPSILON ) {

    const ConnectedTriangle& triangle = mesh.triangles[ triangleIndex ];

    const glm::dvec3& v0 = mesh.vertices[ triangle.vertices[ 0 ] ];

    box3 result = { v0, v0 };

    result.merge( mesh.vertices[ triangle.vertices[ 1 ] ] );
    result.merge( mesh.vertices[ triangle.vertices[ 2 ] ] );

    box3 tolerance = glm::dvec3( tolerance ); 

    result.min -= tolerance;
    result.max += tolerance;

    return result;
  }
}