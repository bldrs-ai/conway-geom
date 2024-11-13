#pragma once

#include "csg_utils.h"

namespace conway::geometry {

  template < size_t N >
  class MultiMeshVertexIndex {
  public:

    MultiMeshVertexIndex(const WingedEdgeDV3* (&input)[N]) {

      size_t offset = 0;

      for ( size_t where = 0; where < N; ++where ) {
     
        vertices_[ where ] = &(input[ where ]->vertices );

        if ( where > 0 ) {
          partitions_[where] = partitions_[ where - 1 ] + vertices_[ where - 1 ]->size();
        }
        else {
          partitions_[ where ] = 0;
        }
      }
    }

    MultiMeshVertexIndex( const std::vector< glm::dvec3 >* (&input)[N]) {

      size_t offset = 0;

      for ( size_t where = 0; where < N; ++where ) {

        vertices_[ where ] = input[ where ];

        if ( where > 0 ) {
          partitions_[ where ] = partitions_[ where - 1 ] + vertices_[ where - 1 ]->size();
        }
        else {
          partitions_[ where ] = 0;
        }
      }
    }

    const glm::dvec3& operator[]( uint32_t index ) const {

      if constexpr ( N == 2 ) {

        uint32_t meshIndex = index >= partitions_[ 1 ];

        return (*vertices_[ meshIndex ])[ index - partitions_[ meshIndex ] ];
      }
       
      if constexpr ( N == 3 ) {

        uint32_t meshIndex =
          static_cast< uint32_t >( index >= partitions_[ 1 ] ) + 
          static_cast< uint32_t >( index >= partitions_[ 2 ] );

        return (*vertices_[ meshIndex ])[ index - partitions_[ meshIndex ] ];
      }
      else {

        for (size_t where = 0; where < ( N - 1 ); ++where) {

          if (index < partitions_[ where + 1 ]) {

            return (*vertices_[ where ])[ index - partitions_[ where ] ];
          }
        }

        return (*vertices_[ N - 1 ])[ index - partitions_[ N - 1 ] ];
      }
    }

    uint32_t operator()(uint32_t mesh, uint32_t index) const {
      
      return partitions_[ mesh ] + index;
    }

    int32_t orient2D( uint32_t v0, uint32_t v1, uint32_t v2, AxisPair axes, double tolerance = 0 ) const {

      int32_t sign = 1;

      // Optimal sorting network here, 5 potential parity swaps.
      if (v0 > v2) {

        std::swap(v0, v2);
        sign = -1;
      }

      if (v0 > v1) {

        std::swap(v0, v1);
        sign *= -1;
      }

      if (v1 > v2) {

        std::swap(v1, v2);
        sign *= -1;
      }

      const glm::dvec3& v0a = *this( v0 );
      const glm::dvec3& v1a = *this( v1 );
      const glm::dvec3& v2a = *this( v2 );

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

    inline int32_t orient3D(
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

      const glm::dvec3& v0a = *this( v0 );
      const glm::dvec3& v1a = *this( v1 );
      const glm::dvec3& v2a = *this( v2 );
      const glm::dvec3& v3a = *this( v3 );

      return orient3D( v0a, v1a, v2a, v3a, tolerance ) * sign;
    }

  private:

    const std::vector< glm::dvec3 >* vertices_[ N ];

    uint32_t partitions_[ N ];
  };

  inline MultiMeshVertexIndex< 2 > multi_mesh_vertex_index( const WingedEdgeDV3& a, const WingedEdgeDV3& b ) {
  
    const WingedEdgeDV3* values[ 2 ] = { &a, &b };

    return MultiMeshVertexIndex< 2 >( values );
  }

  inline MultiMeshVertexIndex< 3 > multi_mesh_vertex_index( const WingedEdgeDV3& a, const WingedEdgeDV3& b, const std::vector< glm::dvec3 >& novel ) {

    const std::vector< glm::dvec3 >* values[ 3 ] = { &a.vertices, &b.vertices, &novel };

    return MultiMeshVertexIndex< 3 >( values );
  }

}