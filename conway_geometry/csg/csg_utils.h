#pragma once

#include <assert.h>
#include <glm/glm.hpp>
#include <math.h>
#include <vector>
//#include "structures/winged_edge.h"
#include "operations/math_utils.h"

#include "representation/Geometry.h"

namespace conway::geometry {

  inline glm::dvec3 plane_line_segment_intersection(
    const glm::dvec3& t0,
    const glm::dvec3& t1,
    const glm::dvec3& t2,
    const glm::dvec3& l0,
    const glm::dvec3& l1 ) {

    // // all of these predicates should agree in sign from the previous tests
    // // against the triangle, guaranteeing the signs and that the point is within
    // // the triangle.
    // double e0 = predicates::adaptive::orient3d( &t0.x, &t1.x, &l0.x, &l1.x );
    // double e1 = predicates::adaptive::orient3d( &t1.x, &t2.x, &l0.x, &l1.x );
    // double e3 = predicates::adaptive::orient3d( &t2.x, &t0.x, &l0.x, &l1.x );

    // // As all the signs are guaranteed to be correct,
    // // this normalization factor is guaranteed to have the right sign.
    // double invDetSum = 1.0 / ( e0 + e1 + e3 );

    // double a = e0 * invDetSum;
    // double b = e1 * invDetSum;
    // double c = 1.0 - std::min( a + b, 1.0 );

    // return t2 * a + t0 * b + t1 * c; 

    glm::dvec3 direction = l1 - l0;
    
    glm::dvec3 e0     = t1 - t0;
    glm::dvec3 e1     = t2 - t0;
    glm::dvec3 origin = l0 - t0;
    glm::dvec3 normal = glm::cross( e0, e1 );

    double t = dot( origin, normal ) / -dot( direction, normal );

    return l0 + direction * t;
  }

  /** Will get the best 2D projection for a triangle that simply involves truncating an axis
   *  As long as the triangle is non-zero area, given that orient2D is exact, it should
   *  give us the truncated axis projection with the biggest area.
   */
  inline AxisPair best_truncated_projection( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2 ) {

    double bestValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Y ) );

    AxisPair result  = AxisPair::X_Y;

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Z ) ); candidateValue > bestValue ) {

      result    = AxisPair::X_Z;
      bestValue = candidateValue;
    }

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::Y_Z ) ); candidateValue > bestValue ) {

      bestValue = candidateValue;
      result = AxisPair::Y_Z;
    }

    if ( bestValue == 0.0 ) {
      result = AxisPair::NONE;
    }

    return result;
  }

  /** Will get the best 2D projection for a triangle that simply involves truncating an axis
   *  As long as the triangle is non-zero area, given that orient2D is exact, it should
   *  give us the truncated axis projection with the biggest area.
   */
  inline AxisPair best_truncated_projection( const glm::dvec3 (&vertices)[3] ) {

    return best_truncated_projection( vertices[ 0 ], vertices[ 1 ], vertices[ 2 ] );
  }

  inline void extract_vertices(
    const Geometry& first,
    const Geometry& second,
    const uint32_t (&vertexIndices)[3],
    glm::dvec3 (&to)[3] ) {

    const std::vector< glm::dvec3 >& firstVertices  = first.vertices;
    const std::vector< glm::dvec3 >& secondVertices = second.vertices;

    uint32_t partition = static_cast< uint32_t >( firstVertices.size() );

    for ( size_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {
      
      uint32_t vIndex = vertexIndices[ vertInTriangle ];

      to[ vertInTriangle ] = vIndex < partition ? firstVertices[ vIndex ] : secondVertices[ vIndex - partition ];
    }
  }

  inline void extract_vertices(
    const Geometry& mesh,
    const Triangle& triangle,
    glm::dvec3 (&to)[3] ) {

    const std::vector< glm::dvec3 >& vertices = mesh.vertices;

    for ( size_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {
      
      to[ vertInTriangle ] = vertices[ triangle.vertices[ vertInTriangle ] ];
    }
  }

    /** Not exact, but it will compute a centroid deterministically without regards to order or winding */
  inline glm::dvec3 centroid( const std::vector< glm::dvec3 >& vertices, const Triangle& triangle ) {

    uint32_t v0 = triangle.vertices[ 0 ];
    uint32_t v1 = triangle.vertices[ 1 ];
    uint32_t v2 = triangle.vertices[ 2 ];

    return ( ( vertices[ v0 ] + ( vertices[ v1 ] + vertices[ v2 ] ) ) ) * ONE_THIRD; 
  }

  /** Calcualte the centroid for a collection of triangles */
  inline glm::dvec3 centroid( const std::vector< glm::dvec3 >& vertices, const std::vector< Triangle >& triangles ) {

    uint32_t vertexCount = 0;
    glm::dvec3 accumulator( 0 );

    for ( const Triangle& triangle : triangles ) {

      uint32_t v0 = triangle.vertices[ 0 ];
      uint32_t v1 = triangle.vertices[ 1 ];
      uint32_t v2 = triangle.vertices[ 2 ];

      vertexCount += 3;
      accumulator += vertices[ v0 ] + ( vertices[ v1 ] + vertices[ v2 ] ); 
    }

    return vertexCount > 0 ? ( accumulator / double( vertexCount ) ) : accumulator;
  }

  inline void reorder_to_lowest_vertex(
    uint32_t (&triangleIndicesInputOutput)[ 3 ] 
  ) {

    uint32_t lowest      = 0;
    uint32_t lowestValue = triangleIndicesInputOutput[ 0 ];

    if ( triangleIndicesInputOutput[ 1 ] < lowestValue ) {

      lowest      = 1;
      lowestValue = triangleIndicesInputOutput[ 1 ];
    }
    
    if ( triangleIndicesInputOutput[ 2 ] < lowestValue ) {

      lowest      = 2;
    }

    uint32_t v0 = triangleIndicesInputOutput[ lowest ];
    uint32_t v1 = triangleIndicesInputOutput[ ( lowest + 1 ) % 3  ];
    uint32_t v2 = triangleIndicesInputOutput[ ( lowest + 2 ) % 3 ];;

    triangleIndicesInputOutput[ 0 ] = v0;
    triangleIndicesInputOutput[ 1 ] = v1;
    triangleIndicesInputOutput[ 2 ] = v2;
  }

  inline bool lowest_vertex_ordered_parity( const uint32_t (&triangleIndices)[ 3 ] ) {

    return triangleIndices[ 1 ] < triangleIndices[ 2 ];
  }

  inline bool less_lowest_vertex_parity( const uint32_t (&left)[ 3 ], const uint32_t (&right)[ 3 ] ) {

    if ( left[ 0 ] < right[ 0 ]) {
      return true;
    }

    if ( left[ 0 ] > right[ 0 ] ) {
      return false;
    }

    bool leftParity  = lowest_vertex_ordered_parity( left );
    bool rightParity = lowest_vertex_ordered_parity( right );

    uint32_t left1  = left[ 1 ];
    uint32_t left2  = left[ 2 ];
    uint32_t right1 = right[ 1 ];
    uint32_t right2 = right[ 2 ];

    if ( !leftParity ) {

      std::swap( left1, left2 );
    }

    if ( !rightParity ) {

      std::swap( right1, right2 );
    }

    if ( left1 < right1 ) {
      return true;
    }

    if ( left1 > right1 ) {
      return false;
    }

    if ( left2 < right2 ) {
      return true;
    }
    
    // equal or greater-than case.
    return false;
  }



}