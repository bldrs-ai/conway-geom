#pragma once

#include <assert.h>
#include <glm/glm.hpp>
#include <math.h>
#include <vector>
#include "structures/winged_edge.h"
#include <CDT/predicates.h>

namespace conway::geometry {

  constexpr size_t SECOND_AXIS_SHIFT = 2;
  constexpr size_t AXIS_MASK         = ( 1 << SECOND_AXIS_SHIFT ) - 1;
  constexpr size_t X_AXIS_INDEX      = 0;
  constexpr size_t Y_AXIS_INDEX      = 1;
  constexpr size_t Z_AXIS_INDEX      = 2;


  constexpr inline size_t make_axis_pair( size_t first, size_t second ) {

    return ( first ) | ( second << SECOND_AXIS_SHIFT );
  }
  
  constexpr inline  size_t first_axis( AxisPair from ) {
  
    return static_cast< size_t >( from ) & AXIS_MASK;
  }

  constexpr inline  size_t second_axis( AxisPair from ) {
  
    return static_cast< size_t >( from ) & AXIS_MASK;
  }


  enum class AxisPair : size_t {

    X_Y   = make_axis_pair( X_AXIS_INDEX, Y_AXIS_INDEX ),
    X_Z   = make_axis_pair( X_AXIS_INDEX, Z_AXIS_INDEX ),
    Y_Z   = make_axis_pair( X_AXIS_INDEX, Z_AXIS_INDEX )

  };

  inline glm::dvec2 extract( const glm::dvec3& from, AxisPair axes ) {

    return glm::dvec2( from[ first_axis( axes ) ], from[ second_axis( axes ) ] );
  }

  double orient2D( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, AxisPair axes ) {

    glm::dvec2 v0t = extract( v0, axes );
    glm::dvec2 v1t = extract( v1, axes );
    glm::dvec2 v2t = extract( v2, axes );

    return predicates::adaptive::orient2d(
      &(v0t.x),
      &(v1t.x),
      &(v2t.x) );
  }

  /** Will get the best 2D projection for a triangle that simply involves truncating an axis
   *  As long as the triangle is non-zero area, given that orient2D is exact, it should
   *  give us the truncated axis projection with the biggest area.
   */
  AxisPair best_truncated_projection( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2 ) {

    double bestValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Y ) );
    AxisPair result  = AxisPair::X_Y;

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Z ) ); candidateValue > bestValue ) {

      result    = AxisPair::X_Z;
      bestValue = candidateValue;
    }

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::Y_Z ) ); candidateValue > bestValue ) {

      result = AxisPair::Y_Z;
    }

    return result;
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

      result = AxisPair::Y_Z;
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
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
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
    const WingedEdgeMesh< glm::dvec3 >& mesh,
    const ConnectedTriangle& triangle,
    glm::dvec3 (&to)[3] ) {

    const std::vector< glm::dvec3 >& vertices = mesh.vertices;

    for ( size_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {
      
      to[ vertInTriangle ] = vertices[ triangle.vertices[ vertInTriangle ] ];
    }
  }


  inline double length2( const glm::dvec3& v ) {

    return glm::dot( v, v );
  }

  inline bool same_point( const glm::dvec3& v0, const glm::dvec3& v1, double tolerance = 0 ) {

    glm::dvec3 comparison = glm::abs( v0 - v1 );

    return tolerance >= comparison.x && tolerance >= comparison.y && tolerance >= comparison.z;
  }
 
  /**
   * Determinant for 3 vectors (matrix essentially)
   */
  inline double determinant3x3( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2 ) {
    
    return 
      ( v0.x * ( ( v1.y * v2.z ) - ( v1.z * v2.y ) ) ) -
      ( v1.x * ( ( v0.y * v2.z ) - ( v0.z * v2.y ) ) ) +
      ( v2.x * ( ( v0.y * v1.z ) - ( v0.z * v1.y ) ) );
  }

  /** Calculate the signed solid angle of a triangle relative a point */
  inline double triangle_solid_angle( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& from ) {

    glm::dvec3 r0 = v0 - from;
    glm::dvec3 r1 = v1 - from;
    glm::dvec3 r2 = v2 - from;

    double r0l = glm::length( r0 );
    double r1l = glm::length( r1 );
    double r2l = glm::length( r2 );

    double x =
      r0l * r1l * r2l + 
      glm::dot( r0, r1 ) * r2l +
      glm::dot( r1, r2 ) * r0l +
      glm::dot( r2, r0 ) * r1l;
    
    // Note, this is not the accurate version provided by orient3D, but it isn't needed for the solid angle purpose.
    double y = determinant3x3( r0, r1, r2 ); 

    return 2.0 * std::atan2( y, x );
  }

  /**
   * Simple determininant based comparison given 4 points that will calculate 3 vectors relative to the first and then
   * work out the sign of the determinant within tolerance (-1 for negative, 0 for within tolerance of 0 and 1 for positive)
   */
  inline int orient3D( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, double tolerance = 0 ) {

    double result = predicates::adaptive::orient3d< double >( &v0.x, &v1.x, &v2.x, &v3.x );

    if ( result > tolerance ) {

      return 1;
    }

    if ( result < tolerance ) {

      return -1;
    }

    return 0;
  }

  /** Not exact, but it will compute a centroid deterministically without regards to order or winding */
  inline glm::dvec3 centroid( const std::vector< glm::dvec3 >& vertices, const Triangle& triangle ) {

    uint32_t v0 = triangle.vertices[ 0 ];
    uint32_t v1 = triangle.vertices[ 1 ];
    uint32_t v2 = triangle.vertices[ 2 ];

    // Optimal sorting network here, 3 swaps.
    if ( v0 > v2 ) {

      std::swap( v0, v2 );
    }

    if ( v0 > v1 ) {

      std::swap( v0, v1 );
    }

    if ( v1 > v2 ) {

      std::swap( v1, v2 );
    }

    // Deterministic order for centroid calculation.
    return ( ( vertices[ v0 ] + ( vertices[ v1 ] + vertices[ v2 ] ) ) ) * ( 1.0 / 3.0 ); 
  }

  enum class VertexColinearTriangle {

    OUTSIDE = 0,
    INSIDE  = 1,
    EDGE0   = 2,
    EDGE1   = 3,
    EDGE2   = 4
  };

  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    AxisPair axes,
    double tolerance = 0 ) {
    
    int32_t sign = 1;

    // Optimal sorting network here, 5 potential parity swaps.
    if ( v0 > v2 ) {

      std::swap( v0, v2 );
      sign = -1;
    }

    if ( v0 > v1 ) {

      std::swap( v0, v1 );
      sign *= -1;
    }

    if ( v1 > v2 ) {

      std::swap( v1, v2 );
      sign *= -1;
    }

    const std::vector< glm::dvec3 >& firstVertiecs = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    uint32_t partition = static_cast< uint32_t >( firstVertiecs.size() );

    const glm::dvec3& v0a = v0 < partition ? firstVertiecs[ v0 ] : secondVertiecs[ v0 - partition ];
    const glm::dvec3& v1a = v1 < partition ? firstVertiecs[ v1 ] : secondVertiecs[ v1 - partition ];
    const glm::dvec3& v2a = v2 < partition ? firstVertiecs[ v2 ] : secondVertiecs[ v2 - partition ];

    double resultValue = orient2D( v0a, v1a, v2a, axes );

    int result = 0;

    if ( resultValue > tolerance ) {

      return sign;
    }

    if ( resultValue < tolerance ) {

      return -sign;
    }

    return 0;
  }

  
  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const uint32_t (&triangleIndices)[ 3 ],
    uint32_t edge,
    uint32_t opposingVertex,
    AxisPair axes,
    double tolerance = 0 ) {

    return orient2D( 
      first,
      second,
      triangleIndices[ edge ],
      triangleIndices[ ( edge + 1 ) % 3 ],
      opposingVertex,
      axes,
      tolerance );
  }

   inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const uint32_t (&triangleIndices)[ 3 ],
    AxisPair axes,
    double tolerance = 0 ) {

    return orient2D( first, second, triangleIndices[ 0 ], triangleIndices[ 1 ], triangleIndices[ 2 ], axes, tolerance );
  }


  /** Given 2 meshes and 4 vertices, where the vertices use an indexing convention of 0 
   * being the first vertice of the first mesh and first.vertices.size() being the first vertex
   * of the second mesh.
   *
   * This will consistently order the vertices, correcting for the parity flips in the determinant.
   */
  inline int32_t orient3D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    uint32_t v3,
    double tolerance = 0 ) {

    int32_t sign = 1;

    // Optimal sorting network here, 5 potential parity swaps.
    if ( v0 > v2 ) {

      std::swap( v0, v2 );
      sign = -1;
    }

    if ( v1 > v3 ) {

      std::swap( v1, v3 );
      sign *= -1;
    }

    if ( v0 > v1 ) {

      std::swap( v0, v1 );
      sign *= -1;
    }

    if ( v2 > v3 ) {

      std::swap( v2, v3 );
      sign *= -1;
    }

    if ( v1 > v2 ) {

      std::swap( v1, v2 );
      sign *= -1;
    }

    const std::vector< glm::dvec3 >& firstVertiecs = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    uint32_t partition = static_cast< uint32_t >( firstVertiecs.size() );

    const glm::dvec3& v0a = v0 < partition ? firstVertiecs[ v0 ] : secondVertiecs[ v0 - partition ];
    const glm::dvec3& v1a = v1 < partition ? firstVertiecs[ v1 ] : secondVertiecs[ v1 - partition ];
    const glm::dvec3& v2a = v2 < partition ? firstVertiecs[ v2 ] : secondVertiecs[ v2 - partition ];
    const glm::dvec3& v3a = v3 < partition ? firstVertiecs[ v3 ] : secondVertiecs[ v3 - partition ];

    return orient3D( v0a, v1a, v2a, v3a, tolerance ) * sign;
  }
}