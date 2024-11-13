#pragma once

#include <assert.h>
#include <glm/glm.hpp>
#include <math.h>
#include <vector>
#include "structures/winged_edge.h"

#if defined (_MSC_VER)

#pragma warning( push )
#pragma warning( disable : 26812 )

#endif

#include "predicates.h"

#if defined (_MSC_VER)

#pragma warning( pop )

#endif


namespace conway::geometry {

  constexpr size_t SECOND_AXIS_SHIFT = 2;
  constexpr size_t AXIS_MASK         = ( 1 << SECOND_AXIS_SHIFT ) - 1;
  constexpr size_t X_AXIS_INDEX      = 0;
  constexpr size_t Y_AXIS_INDEX      = 1;
  constexpr size_t Z_AXIS_INDEX      = 2;

  constexpr inline size_t make_axis_pair( size_t first, size_t second ) {

    return ( first ) | ( second << SECOND_AXIS_SHIFT );
  }

  enum class AxisPair : size_t {

    X_Y = make_axis_pair( X_AXIS_INDEX, Y_AXIS_INDEX ),
    X_Z = make_axis_pair( X_AXIS_INDEX, Z_AXIS_INDEX ),
    Y_Z = make_axis_pair( Y_AXIS_INDEX, Z_AXIS_INDEX )

  };

  constexpr inline size_t first_axis( AxisPair from ) {
  
    return static_cast< size_t >( from ) & AXIS_MASK;
  }

  constexpr inline size_t second_axis( AxisPair from ) {
  
    return ( static_cast< size_t >( from ) >> SECOND_AXIS_SHIFT );
  }

  inline glm::dvec3 plane_line_segment_intersection(
    const glm::dvec3& t0,
    const glm::dvec3& t1,
    const glm::dvec3& t2,
    const glm::dvec3& l0,
    const glm::dvec3& l1 ) {

    glm::dvec3 direction = l1 - l0;
    
    glm::dvec3 e0     = t1 - t0;
    glm::dvec3 e1     = t2 - t0;
    glm::dvec3 origin = l0 - t0;
    glm::dvec3 normal = glm::cross( e0, e1 );

    double t = dot( origin, normal ) / -dot( direction, normal );

    return l0 + direction * t;
  }

  /** Assuming 2 intersecting line segments, this gets the intersection point */
  inline glm::dvec3 line_segment_line_segment_intersection(
    const glm::dvec3& a0,
    const glm::dvec3& a1,
    const glm::dvec3& b0,
    const glm::dvec3& b1 ) {

    glm::dvec3 e0        = a1 - a0;
    glm::dvec3 direction = b1 - b0;

    // intersecting non-colinear line segments are 
    // coplanar, and have a well defined normal.
    glm::dvec3 e1     = glm::cross( e0, direction );
    glm::dvec3 origin = b0 - a0;
    glm::dvec3 normal = glm::cross( e0, e1 );

    double t = dot( origin, normal ) / -dot( direction, normal );

    return b0 + direction * t;
  }

  inline glm::dvec2 extract( const glm::dvec3& from, AxisPair axes ) {

    return glm::dvec2( from[ first_axis( axes ) ], from[ second_axis( axes ) ] );
  }

  inline std::pair< double, double > extract_pair( const glm::dvec3& from, AxisPair axes ) {

    return std::make_pair( from[ first_axis( axes ) ], from[ second_axis( axes ) ] );
  }

  inline double orient2D( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, AxisPair axes ) {

    glm::dvec2 v0t = extract( v0, axes );
    glm::dvec2 v1t = extract( v1, axes );
    glm::dvec2 v2t = extract( v2, axes );

    return predicates::adaptive::orient2d(
      &(v0t.x),
      &(v1t.x),
      &(v2t.x) );
  }

  inline int32_t orient2D(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, AxisPair axes, double tolerance ) {

    double orientation = orient2D( v0, v1, v2, axes, tolerance );
    
    if ( orientation > tolerance ) {

      return 1;
    }

    if ( orientation < tolerance ) {

      return -1;
    }

    return 0;
  }

  /** Will get the best 2D projection for a triangle that simply involves truncating an axis
   *  As long as the triangle is non-zero area, given that orient2D is exact, it should
   *  give us the truncated axis projection with the biggest area.
   */
  inline AxisPair best_truncated_projection( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2 ) {

    double bestValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Y ) );

    //printf( "%f x_y\n", bestValue );

    AxisPair result  = AxisPair::X_Y;

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::X_Z ) ); candidateValue > bestValue ) {

      //printf("%f x_z\n", candidateValue );

      result    = AxisPair::X_Z;
      bestValue = candidateValue;
    }

    if ( double candidateValue = fabs( orient2D( v0, v1, v2, AxisPair::Y_Z ) ); candidateValue > bestValue ) {

      //printf("%f y_z\n", candidateValue );

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

    // We use the exactly signed version of the determinant here because we really want the sign to be right,
    // as this gives us our accurate face winding up close for inside vs outside.
    // Because the dipoles in the BVH are far away, we don't need the precision, but here we do.
    double y = predicates::adaptive::orient3d< double >( &v0.x, &v1.x, &v2.x, &from.x );

    return 2.0 * std::atan2( y, fabs( x ) );
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
    //if ( v0 > v2 ) {

    //  std::swap( v0, v2 );
    //}

    //if ( v0 > v1 ) {

    //  std::swap( v0, v1 );
    //}

    //if ( v1 > v2 ) {

    //  std::swap( v1, v2 );
    //}

    // Deterministic order for centroid calculation.
    return ( ( vertices[ v0 ] + ( vertices[ v1 ] + vertices[ v2 ] ) ) ) * ( 1.0 / 3.0 ); 
  }

  /** Not exact, but it will compute a centroid deterministically without regards to order or winding */
  inline glm::dvec3 centroid(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const std::vector< glm::dvec3 >& novelVertices,
    const Triangle& triangle ) {

    uint32_t v0 = triangle.vertices[ 0 ];
    uint32_t v1 = triangle.vertices[ 1 ];
    uint32_t v2 = triangle.vertices[ 2 ];

    // Optimal sorting network here, 3 swaps.
    //if ( v0 > v2 ) {

    //  std::swap( v0, v2 );
    //}

    //if ( v0 > v1 ) {

    //  std::swap( v0, v1 );
    //}

    //if ( v1 > v2 ) {

    //  std::swap( v1, v2 );
    //}

    const std::vector< glm::dvec3 >& firstVertices  = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    const std::vector< glm::dvec3 >* vertexSet[ 3 ] = {

      &first.vertices,
      &second.vertices,
      &novelVertices
    };

    uint32_t bPartition = static_cast< uint32_t >( firstVertices.size() );
    uint32_t novelPartition = bPartition + static_cast< uint32_t >( secondVertiecs.size() );

    uint32_t vertexOffset[ 3 ] = {
      0,
      bPartition,
      novelPartition
    };

    uint32_t v0set = ( v0 >= bPartition ) + ( v0 >= novelPartition ); 
    uint32_t v1set = ( v1 >= bPartition ) + ( v1 >= novelPartition );
    uint32_t v2set = ( v2 >= bPartition ) + ( v2 >= novelPartition );

    const glm::dvec3& v0a = (*vertexSet[ v0set ])[ v0 - vertexOffset[ v0set ] ];
    const glm::dvec3& v1a = (*vertexSet[ v1set ])[ v1 - vertexOffset[ v1set ] ];
    const glm::dvec3& v2a = (*vertexSet[ v2set ])[ v2 - vertexOffset[ v2set ] ];

    // Deterministic order for centroid calculation.
    return ( ( v0a + ( v1a + v2a ) ) ) * ( 1.0 / 3.0 ); 
  }


  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    AxisPair axes,
    double tolerance = 0 ) {
    
    //int32_t sign = 1;

    // Optimal sorting network here, 5 potential parity swaps.
    //if ( v0 > v2 ) {

    //  std::swap( v0, v2 );
    //  sign = -1;
    //}

    //if ( v0 > v1 ) {

    //  std::swap( v0, v1 );
    //  sign *= -1;
    //}

    //if ( v1 > v2 ) {

    //  std::swap( v1, v2 );
    //  sign *= -1;
    //}

    const std::vector< glm::dvec3 >& firstVertices = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    uint32_t partition = static_cast< uint32_t >( firstVertices.size() );

    const glm::dvec3& v0a = v0 < partition ? firstVertices[ v0 ] : secondVertiecs[ v0 - partition ];
    const glm::dvec3& v1a = v1 < partition ? firstVertices[ v1 ] : secondVertiecs[ v1 - partition ];
    const glm::dvec3& v2a = v2 < partition ? firstVertices[ v2 ] : secondVertiecs[ v2 - partition ];

    double resultValue = orient2D( v0a, v1a, v2a, axes );

    if ( resultValue > tolerance ) {

      return 1;
    }

    if ( resultValue < tolerance ) {

      return -1;
    }

    return 0;
  }


  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const std::vector< glm::dvec3 >& novelVertices,
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    AxisPair axes,
    double tolerance = 0 ) {
    //
    //int32_t sign = 1;

    //// Optimal sorting network here, 5 potential parity swaps.
    //if ( v0 > v2 ) {

    //  std::swap( v0, v2 );
    //  sign = -1;
    //}

    //if ( v0 > v1 ) {

    //  std::swap( v0, v1 );
    //  sign *= -1;
    //}

    //if ( v1 > v2 ) {

    //  std::swap( v1, v2 );
    //  sign *= -1;
    //}

    const std::vector< glm::dvec3 >& firstVertices = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    const std::vector< glm::dvec3 >* vertexSet[ 3 ] = {

      &first.vertices,
      &second.vertices,
      &novelVertices
    };

    uint32_t bPartition = static_cast< uint32_t >( firstVertices.size() );
    uint32_t novelPartition = bPartition + static_cast< uint32_t >( secondVertiecs.size() );

    uint32_t vertexOffset[ 3 ] = {
      0,
      bPartition,
      novelPartition
    };

    uint32_t v0set = ( v0 >= bPartition ) + ( v0 >= novelPartition ); 
    uint32_t v1set = ( v1 >= bPartition ) + ( v1 >= novelPartition );
    uint32_t v2set = ( v2 >= bPartition ) + ( v2 >= novelPartition );

    const glm::dvec3& v0a = (*vertexSet[ v0set ])[ v0 - vertexOffset[ v0set ] ];
    const glm::dvec3& v1a = (*vertexSet[ v1set ])[ v1 - vertexOffset[ v1set ] ];
    const glm::dvec3& v2a = (*vertexSet[ v2set ])[ v2 - vertexOffset[ v2set ] ];

    double resultValue = orient2D( v0a, v1a, v2a, axes );

    int result = 0;

    if ( resultValue > tolerance ) {

      return 1;
    }

    if ( resultValue < tolerance ) {

      return -1;
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

  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const uint32_t (&triangleIndices)[ 3 ],
    AxisPair axes,
    double tolerance = 0 ) {
    return orient2D( first, second, triangleIndices[ 0 ], triangleIndices[ 1 ], triangleIndices[ 2 ], axes, tolerance );
  }

  
  inline int32_t orient2D(
    const WingedEdgeMesh< glm::dvec3 >& first,
    const WingedEdgeMesh< glm::dvec3 >& second,
    const std::vector< glm::dvec3 >& novelVertices,
    const uint32_t (&triangleIndices)[ 3 ],
    AxisPair axes,
    double tolerance = 0 ) {

    return orient2D( first, second, novelVertices, triangleIndices[ 0 ], triangleIndices[ 1 ], triangleIndices[ 2 ], axes, tolerance );
  }

  inline int32_t orient2D( const glm::dvec3(&vertices)[ 3 ], AxisPair axes, double tolerance = 0 ) {

    double determinant = orient2D( vertices[ 0 ], vertices[ 1 ], vertices[ 2 ], axes );

    if (determinant > tolerance) {

      return 1;
    }

    if (determinant < tolerance) {

      return -1;
    }

    return 0;
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

    //int32_t sign = 1;

    //// Optimal sorting network here, 5 potential parity swaps.
    //if ( v0 > v2 ) {

    //  std::swap( v0, v2 );
    //  sign = -1;
    //}

    //if ( v1 > v3 ) {

    //  std::swap( v1, v3 );
    //  sign *= -1;
    //}

    //if ( v0 > v1 ) {

    //  std::swap( v0, v1 );
    //  sign *= -1;
    //}

    //if ( v2 > v3 ) {

    //  std::swap( v2, v3 );
    //  sign *= -1;
    //}

    //if ( v1 > v2 ) {

    //  std::swap( v1, v2 );
    //  sign *= -1;
    //}

    const std::vector< glm::dvec3 >& firstVertiecs = first.vertices;
    const std::vector< glm::dvec3 >& secondVertiecs = second.vertices;

    uint32_t partition = static_cast< uint32_t >( firstVertiecs.size() );

    const glm::dvec3& v0a = v0 < partition ? firstVertiecs[ v0 ] : secondVertiecs[ v0 - partition ];
    const glm::dvec3& v1a = v1 < partition ? firstVertiecs[ v1 ] : secondVertiecs[ v1 - partition ];
    const glm::dvec3& v2a = v2 < partition ? firstVertiecs[ v2 ] : secondVertiecs[ v2 - partition ];
    const glm::dvec3& v3a = v3 < partition ? firstVertiecs[ v3 ] : secondVertiecs[ v3 - partition ];

    return orient3D( v0a, v1a, v2a, v3a, tolerance );// * sign;
  }
}