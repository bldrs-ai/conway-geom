#pragma once

#include "csg_utils.h"

namespace conway::geometry {

  template < size_t N >
  class MultiMeshVertexIndex {
  public:

    MultiMeshVertexIndex(const WingedEdgeDV3* (&input)[N]) {

      for ( size_t where = 0; where < N; ++where ) {
     
        vertices_[ where ] = &(input[ where ]->vertices );

        if ( where > 0 ) {
          partitions_[ where ] =
            partitions_[ where - 1 ] +
            static_cast< uint32_t >( vertices_[ where - 1 ]->size() );
        }
        else {
          partitions_[ where ] = 0;
        }
      }
    }

    uint32_t size() const {

      return partitions_[ N - 1 ] + static_cast< uint32_t >( vertices_[ N - 1 ]->size() );
    }


    MultiMeshVertexIndex( const std::vector< glm::dvec3 >& singular ) {
    
      for (size_t where = 0; where < N; ++where) {

        vertices_[ where ]   = &singular;
        partitions_[ where ] = 0;
      }
    }

    MultiMeshVertexIndex( const std::vector< glm::dvec3 >* (&input)[N] ) {

      for ( size_t where = 0; where < N; ++where ) {

        vertices_[ where ] = input[ where ];

        if ( where > 0 ) {
          partitions_[ where ] =
            partitions_[ where - 1 ] +
            static_cast< uint32_t >( vertices_[ where - 1 ]->size() );
        }
        else {
          partitions_[ where ] = 0;
        }
      }
    }


    /** Not exact, but it will compute a centroid deterministically without regards to order or winding */
    glm::dvec3 centroid( const Triangle& triangle ) {

      uint32_t v0 = triangle.vertices[ 0 ];
      uint32_t v1 = triangle.vertices[ 1 ];
      uint32_t v2 = triangle.vertices[ 2 ];

      const glm::dvec3& v0a = (*this)[ v0 ];
      const glm::dvec3& v1a = (*this)[ v1 ];
      const glm::dvec3& v2a = (*this)[ v2 ];

      // Deterministic order for centroid calculation.
      return ( ( v0a + ( v1a + v2a ) ) ) * (1.0 / 3.0);
    }

    const glm::dvec3& operator[]( uint32_t index ) const {

      static_assert( N <= 3 );

      if constexpr ( N == 2 ) {

        uint32_t meshIndex = index >= partitions_[ 1 ];

        assert( vertices_[ meshIndex ]->size() > ( index - partitions_[ meshIndex ] ) );

        return ( *vertices_[ meshIndex ])[ index - partitions_[ meshIndex ] ];
      }
       
      if constexpr ( N == 3 ) {

        uint32_t meshIndex =
          static_cast< uint32_t >( index >= partitions_[ 1 ] ) + 
          static_cast< uint32_t >( index >= partitions_[ 2 ] );

        return (*vertices_[ meshIndex ])[ index - partitions_[ meshIndex ] ];
      }
    }

    uint32_t operator()( uint32_t mesh, uint32_t index ) const {
      
      return partitions_[ mesh ] + index;
    }

    void extract( uint32_t meshIndex, const uint32_t(&vertexIndices)[ 3 ], glm::dvec3 ( &vertices )[ 3 ] ) const {

      const std::vector< glm::dvec3 >& meshVertices =  *(vertices_[ meshIndex ]);

      vertices[ 0 ] = meshVertices[ vertexIndices[ 0 ] ];
      vertices[ 1 ] = meshVertices[ vertexIndices[ 1 ] ];
      vertices[ 2 ] = meshVertices[ vertexIndices[ 2 ] ];
    }

    bool isZeroAreaTriangle( const uint32_t (&triangleIndices)[ 3 ], double tolerance = 0 ) const {

      return isZeroAreaTriangle( triangleIndices[ 0 ], triangleIndices[ 1 ], triangleIndices[ 2 ] );
    }

    bool isZeroAreaTriangle( uint32_t v0, uint32_t v1, uint32_t v2, double tolerance = 0 ) const {

      const glm::dvec3& v0a = (*this)[ v0 ];
      const glm::dvec3& v1a = (*this)[ v1 ];
      const glm::dvec3& v2a = (*this)[ v2 ];

      return is_zero_area_triangle( v0a, v1a, v2a, tolerance );
    }

    int32_t orient2D( const uint32_t (&triangleIndices)[ 3 ], AxisPair axes, double tolerance = 0 ) const {
    
      return orient2D( triangleIndices[ 0 ], triangleIndices[ 1 ], triangleIndices[ 2 ], axes, tolerance );
    }

    int32_t orient2D( uint32_t v0, uint32_t v1, uint32_t v2, AxisPair axes, double tolerance = 0 ) const {

      const glm::dvec3& v0a = (*this)[ v0 ];
      const glm::dvec3& v1a = (*this)[ v1 ];
      const glm::dvec3& v2a = (*this)[ v2 ];

      return conway::geometry::orient2D( v0a, v1a, v2a, axes, tolerance );
    }

    inline int32_t orient2D(
      const uint32_t(&triangleIndices)[3],
      uint32_t edge,
      uint32_t opposingVertex,
      AxisPair axes,
      double tolerance = 0) const {

      return orient2D(
        triangleIndices[ edge ],
        triangleIndices[ ( edge + 1 ) % 3 ],
        opposingVertex,
        axes,
        tolerance);
    }

    inline int32_t orient3D(
      uint32_t v0,
      uint32_t v1,
      uint32_t v2,
      uint32_t v3,
      double tolerance = 0 ) const {

      const glm::dvec3& v0a = (*this)[ v0 ];
      const glm::dvec3& v1a = (*this)[ v1 ];
      const glm::dvec3& v2a = (*this)[ v2 ];
      const glm::dvec3& v3a = (*this)[ v3 ];

      return conway::geometry::orient3D( v0a, v1a, v2a, v3a, tolerance );
    }

  private:

    const std::vector< glm::dvec3 >* vertices_[ N ] {};

    uint32_t partitions_[ N ] {};
  };

  inline MultiMeshVertexIndex< 2 > multi_mesh_vertex_index( const WingedEdgeDV3& a, const WingedEdgeDV3& b ) {
  
    const WingedEdgeDV3* values[ 2 ] = { &a, &b };

    return MultiMeshVertexIndex< 2 >( values );
  }

  inline MultiMeshVertexIndex< 2 > multi_mesh_vertex_index( const WingedEdgeDV3& a, const std::vector< glm::dvec3 >& novel ) {

    const std::vector< glm::dvec3 >* values[2] = { &a.vertices, &novel };

    return MultiMeshVertexIndex< 2 >( values );
  }

  inline MultiMeshVertexIndex< 3 > multi_mesh_vertex_index( const WingedEdgeDV3& a, const WingedEdgeDV3& b, const std::vector< glm::dvec3 >& novel ) {

    const std::vector< glm::dvec3 >* values[ 3 ] = { &a.vertices, &b.vertices, &novel };

    return MultiMeshVertexIndex< 3 >( values );
  }

}