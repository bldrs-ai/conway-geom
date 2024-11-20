#pragma once

#include "hash_functions.h"
#include "winged_edge.h"
#include "csg/csg_utils.h"

namespace conway::geometry {

  struct VertexWelder {

    std::unordered_map< glm::dvec3, uint32_t > unique;
    std::vector< uint32_t >                    remap;
    std::vector< ConnectedTriangle >           old_triangles;

    void weld( WingedEdgeDV3& toWeld ) {

      std::vector< glm::dvec3 >& vertices = toWeld.vertices;

      old_triangles.clear();
      remap.resize( toWeld.vertices.size() );
      unique.clear();

      for (
        uint32_t where = 0, end = static_cast<uint32_t>( toWeld.vertices.size() );
        where < end;
        ++where) {

        glm::dvec3 localVertex = vertices[ where ];

        auto [ to, success ] = unique.try_emplace( localVertex, static_cast< uint32_t >( unique.size() ) );

        remap[ where ] = to->second;
      }

      vertices.resize( unique.size() );

      for ( const auto& [ vertex, index ] : unique ) {

        vertices[ index ] = vertex;

      }

      std::swap( old_triangles, toWeld.triangles );

      toWeld.edges.clear();
      toWeld.edge_map.clear();

      for ( const ConnectedTriangle& triangle : old_triangles ) {

        uint32_t i0 = remap[ triangle.vertices[ 0 ] ];
        uint32_t i1 = remap[ triangle.vertices[ 1 ] ];
        uint32_t i2 = remap[ triangle.vertices[ 2 ] ];

        if ( i0 == i1 || i1 == i2 || i2 == i0 ) {

          continue;
        }

        if ( is_zero_area_triangle(
          vertices[ i0 ],
          vertices[ i1 ],
          vertices[ i2 ] ) ) {

          continue;
        }

        toWeld.makeTriangle( i0, i1, i2 );
      }
    }

  };

}