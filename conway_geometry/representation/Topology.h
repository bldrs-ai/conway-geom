#pragma once

#include <stdint.h>
#include <utility>

namespace conway::geometry {

  constexpr uint32_t EMPTY_INDEX = 0xFFFFFFFF;

  struct Edge {
    
    uint32_t triangles[ 2 ] = { EMPTY_INDEX, EMPTY_INDEX };

    uint32_t vertices[ 2 ] = {};

    uint32_t non_manifold_linked_edge = EMPTY_INDEX;

    uint32_t otherVertex( uint32_t vertexIndex ) const {

      return vertices[ 0 ] == vertexIndex ? vertices[ 1 ] : vertices[ 0 ];
    }
    
    uint32_t otherTriangle( uint32_t triangleIndex ) const {

      return triangles[ 0 ] == triangleIndex ? triangles[ 1 ] : triangles[ 0 ];
    }

    bool hasVertex( uint32_t vertexIndex ) const {
      return vertices[ 0 ] == vertexIndex || vertices[ 1 ] == vertexIndex;
    }

    bool border() const { return triangles[ 0 ] == EMPTY_INDEX || triangles[ 1 ] == EMPTY_INDEX; }
  };

  struct Triangle {

    uint32_t vertices[ 3 ] = { EMPTY_INDEX, EMPTY_INDEX, EMPTY_INDEX };

    bool operator== ( const Triangle& against ) const {

      return
        against.vertices[ 0 ] == vertices[ 0 ] &&
        against.vertices[ 1 ] == vertices[ 1 ] && 
        against.vertices[ 2 ] == vertices[ 1 ];
    }

    bool operator!= ( const Triangle& against ) const {

      return
        against.vertices[0] == vertices[0] &&
        against.vertices[1] == vertices[1] &&
        against.vertices[2] == vertices[1];
    }

    uint32_t otherVertex( const Edge& edge ) const;

    uint32_t vertexInTriangle( uint32_t vertexIndex ) const;

  };
    
  struct TriangleEdges {

    uint32_t edges[ 3 ] = { EMPTY_INDEX, EMPTY_INDEX, EMPTY_INDEX };

  };

  struct ConnectedTriangle : public Triangle, public TriangleEdges {};

  inline uint32_t Triangle::vertexInTriangle( uint32_t vertexIndex ) const {

    for ( uint32_t where = 0; where < 3; ++where ) {

      if ( vertexIndex == vertices[ where ] ) {

        return where;
      }
    }

    assert( false );

    return EMPTY_INDEX;
  }

  inline uint32_t Triangle::otherVertex( const Edge& edge ) const {

   uint32_t ev0 = edge.vertices[0];
   uint32_t ev1 = edge.vertices[1];

   for ( uint32_t localVertex = 0; localVertex < 3; ++localVertex ) {

     uint32_t vertex = vertices[ localVertex ];

     if (vertex != ev0 && vertex != ev1) {

       return vertex;
     }
   }

   // we should never get here.
   return EMPTY_INDEX;
  }

  inline uint64_t edgeCompoundID( uint32_t vertex1, uint32_t vertex2 ) {

    if ( vertex1 > vertex2 ) {
      std::swap( vertex1, vertex2 );
    }

    return ( uint64_t( vertex1 ) << uint64_t( 32 ) ) | uint64_t( vertex2 );
  }

}