#pragma once

#include <stdint.h>
#include <vector>
#include <unordered_map>
#include "representation/IfcGeometry.h"
#include <optional>
#include "structures/aabb_tree.h"

namespace conway::geometry {

  constexpr uint32_t EMPTY_INDEX = 0xFFFFFFFF;

  struct Edge {
    
    uint32_t triangles[ 2 ] = { EMPTY_INDEX, EMPTY_INDEX };

    uint32_t vertices[ 2 ] = {};

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

    uint32_t vertices[ 3 ] = {};

  };

  struct ConnectedTriangle : public Triangle {

    uint32_t edges[ 3 ] = {};

    uint32_t otherVertex( const Edge& edge ) const {

      uint32_t ev0 = edge.vertices[ 0 ];
      uint32_t ev1 = edge.vertices[ 1 ];

      for ( uint32_t localVertex = 0; localVertex < 3; ++localVertex ) {

        uint32_t vertex = vertices[ localVertex ];

        if ( vertex != ev0 && vertex != ev1 ) {

          return vertex;
        }
      }

      // we should never get here.
      return EMPTY_INDEX;
    }    

  };

  inline uint64_t edgeCompoundID( uint32_t vertex1, uint32_t vertex2 ) {

    if ( vertex1 > vertex2 ) {
      std::swap( vertex1, vertex2 );
    }

    return ( uint64_t( vertex1 ) << uint64_t( 32 ) ) | uint64_t( vertex2 );
  }

  template < typename VertexType >
  struct WingedEdgeMesh {

    std::vector< ConnectedTriangle > triangles;
    
    std::vector< Edge > edges;

    std::vector< VertexType > vertices;

    std::unordered_map< uint64_t, uint32_t > edge_map;

    std::optional< AABBTree > bvh;

    void makeTriangle( uint32_t a, uint32_t b, uint32_t c, uint32_t index ) {

      ConnectedTriangle& triangle = triangles[ index ];

      triangle.vertices[ 0 ]  = a;
      triangle.vertices[ 1 ]  = b;
      triangle.vertices[ 2 ]  = c;

      triangle.edges[ 0 ] = makeEdge( a, b, index );
      triangle.edges[ 1 ] = makeEdge( b, c, index );
      triangle.edges[ 2 ] = makeEdge( c, a, index );
    }

    void clear() {
      bvh.reset();
      edge_map.clear();
      vertices.clear();
      edges.clear();
      triangles.clear();
    }

    /** Construct a BVH for this. */
    void makeBVH() {

      if ( !bvh.has_value() ) {
        bvh.emplace( *this );
      }
    }

    std::optional< uint32_t > getEdge( uint32_t v0, uint32_t v1 ) const {

      auto mapIterator = edge_map.find( edgeCompoundID( v0, v1 ) );

      if ( mapIterator == edge_map.end() ) {
        return std::nullopt;
      }

      return std::optional< uint32_t >( mapIterator->second );
    }

    void deleteTriangle( uint32_t index ) {

      const ConnectedTriangle& toDelete = triangles[ index ];

      for ( uint32_t localEdge = 0; localEdge < 3; ++localEdge ) {

        Edge& edge = edges[ toDelete.edges[ localEdge ] ];

        if ( edge.triangles[ 0 ] == index ) {

          edge.triangles[ 0 ] = edge.triangles[ 1 ];

        } 

        edge.triangles[ 1 ] = EMPTY_INDEX;
      }

      uint32_t backIndex = static_cast< uint32_t >( triangles.size() - 1 );
      if ( index != backIndex ) {

        const ConnectedTriangle& back = triangles.back();
      
        for ( uint32_t localEdge = 0; localEdge < 3; ++localEdge ) {

          Edge& edge = edges[ back.edges[ localEdge ] ];

          for ( uint32_t onEdge = 0; onEdge < 2; ++onEdge ) {

            if ( edge.triangles[ onEdge ] == backIndex ) {
              edge.triangles[ onEdge ] = index;
              break;
            }
          }
        }

        triangles[ index ] = back;
      } 

      triangles.pop_back();
    }

    uint32_t makeVertex( const VertexType& value ) {

      uint32_t index = static_cast< uint32_t >( vertices.size() );

      vertices.push_back( value );

      return index;
    }

    uint32_t makeTriangle( uint32_t a, uint32_t b, uint32_t c ) {

      uint32_t index =
        static_cast< uint32_t >( triangles.size() );

      triangles.push_back( ConnectedTriangle {} );

      makeTriangle( a, b, c, index );

      return index;
    }

    uint32_t makeEdge( uint32_t v1, uint32_t v2, uint32_t triangleIndex ) {
      uint64_t edgeIdentifier = edgeCompoundID( v1, v2 );

      auto [ mapIterator, emplaced ] = edge_map.try_emplace( edgeIdentifier, static_cast< uint32_t >( edges.size() ) );

      if ( emplaced ) {
        edges.push_back( Edge { { triangleIndex, EMPTY_INDEX }, { v1, v2 } } );               
      }
      else {
        Edge& currentEdge = edges[ mapIterator->second ];

        currentEdge.triangles[ currentEdge.triangles[ 0 ] == EMPTY_INDEX ? 0 : 1 ] = triangleIndex;
      }

      return mapIterator->second;
    }
  };
}