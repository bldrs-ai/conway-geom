#pragma once

#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <optional>
#include <string>
#include <glm/glm.hpp>
#include "structures/aabb_tree.h"
#include <stdio.h>

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

  };

  struct ConnectedTriangle : public Triangle {

    uint32_t edges[ 3 ] = { EMPTY_INDEX, EMPTY_INDEX, EMPTY_INDEX };

    uint32_t otherVertex( const Edge& edge ) const;
  };

  // Specialise this to extract different vertex types of OBJ dump
  template < typename VertexType >
  struct VertexExtractor {

      glm::dvec3 operator()( const VertexType& value ) const {

        return glm::dvec3( 0 );
      }
  };

  // Specialised for the glm::dvec3 case.
  template <>
  struct VertexExtractor< glm::dvec3 > {

      const glm::dvec3& operator()( const glm::dvec3& value ) const {

          return value;
      }
  };

  template < typename VertexType >
  struct WingedEdgeMesh {

    std::vector< ConnectedTriangle > triangles;
    
    std::vector< Edge > edges;

    std::vector< VertexType > vertices;

    std::unordered_map< uint64_t, uint32_t > edge_map;

    std::optional< AABBTree > bvh;

    void makeTriangle( uint32_t a, uint32_t b, uint32_t c, uint32_t index );

    void clear();

    /** Construct a BVH for this. */
    AABBTree& makeBVH() {

      if ( !bvh.has_value() ) {
        bvh.emplace( *this );
      }

      return *bvh;
    }

    std::optional< uint32_t > getEdge( uint32_t v0, uint32_t v1 ) const;

    void deleteTriangle( uint32_t index );

    uint32_t makeVertex( const VertexType& value ) ;

    uint32_t makeTriangle( uint32_t a, uint32_t b, uint32_t c );

    uint32_t makeEdge( uint32_t v1, uint32_t v2, uint32_t triangleIndex );

    std::string dumpToOBJ( const std::string& preamble = "" ) const;

    void append( const WingedEdgeMesh< VertexType >& other );
  };

  using WingedEdgeDV3 = WingedEdgeMesh< glm::dvec3 >;

  inline uint64_t edgeCompoundID( uint32_t vertex1, uint32_t vertex2 ) {

    if ( vertex1 > vertex2 ) {
      std::swap( vertex1, vertex2 );
    }

    return ( uint64_t( vertex1 ) << uint64_t( 32 ) ) | uint64_t( vertex2 );
  }

  inline uint32_t ConnectedTriangle::otherVertex( const Edge& edge ) const {

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


  template < typename VertexType >
  inline void WingedEdgeMesh< VertexType >::makeTriangle(uint32_t a, uint32_t b, uint32_t c, uint32_t index) {

    ConnectedTriangle& triangle = triangles[index];
    
    assert( a < vertices.size() );
    assert( b < vertices.size() );
    assert( c < vertices.size() );

    triangle.vertices[ 0 ] = a;
    triangle.vertices[ 1 ] = b;
    triangle.vertices[ 2 ] = c;

    triangle.edges[ 0 ] = makeEdge( a, b, index );
    triangle.edges[ 1 ] = makeEdge( b, c, index );
    triangle.edges[ 2 ] = makeEdge( c, a, index );
  }

  template < typename VertexType >
  inline void WingedEdgeMesh< VertexType >::clear() {
    bvh.reset();
    edge_map.clear();
    vertices.clear();
    edges.clear();
    triangles.clear();
  }

  template < typename VertexType >
  inline std::optional< uint32_t > WingedEdgeMesh< VertexType >::getEdge(uint32_t v0, uint32_t v1) const {

    auto mapIterator =
      edge_map.find( edgeCompoundID( v0, v1 ) );

    if ( mapIterator == edge_map.end() ) {
      return std::nullopt;
    }

    return std::optional< uint32_t >( mapIterator->second );
  }

  template < typename VertexType >
  inline void WingedEdgeMesh< VertexType >::deleteTriangle( uint32_t index ) {

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

  template < typename VertexType >
  inline uint32_t WingedEdgeMesh< VertexType >::makeVertex( const VertexType& value ) {

    uint32_t index = static_cast<uint32_t>( vertices.size() );

    vertices.push_back( value );

    return index;
  }

  template < typename VertexType >
  inline uint32_t WingedEdgeMesh< VertexType >::makeTriangle( uint32_t a, uint32_t b, uint32_t c ) {

    uint32_t index =
      static_cast<uint32_t>(triangles.size());

    triangles.push_back( ConnectedTriangle{} );

    makeTriangle( a, b, c, index );

    return index;
  }

  template < typename VertexType >
  inline uint32_t WingedEdgeMesh< VertexType >::makeEdge( uint32_t v1, uint32_t v2, uint32_t triangleIndex ) {

    uint64_t edgeIdentifier = edgeCompoundID( v1, v2 );

    while ( true ) {

      auto [ mapIterator, emplaced ] = edge_map.try_emplace( edgeIdentifier, static_cast<uint32_t>( edges.size() ) );

      if ( emplaced ) {
        
        edges.push_back( Edge{ { triangleIndex, EMPTY_INDEX }, { v1, v2 } });

      }
      else {

        Edge& currentEdge = edges[ mapIterator->second ];

        // Note, this allows non manifold edges to exist,
        // but edges will not point back to triangles
        if ( currentEdge.triangles[ 0 ] == EMPTY_INDEX || currentEdge.triangles[ 1 ] == EMPTY_INDEX ) {

          currentEdge.triangles[ currentEdge.triangles[ 0 ] == EMPTY_INDEX ? 0 : 1 ] = triangleIndex;
        }
        else {

          currentEdge.non_manifold_linked_edge = static_cast< uint32_t >( edges.size() );

          // This is a work around for non-manifold cases.
          edge_map.erase( mapIterator );
          continue;

        }
      }

      return mapIterator->second;
    }
  }

  template < typename VertexType >
  inline std::string WingedEdgeMesh< VertexType >::dumpToOBJ( const std::string& preamble ) const {

    std::string obj;

    size_t numPoints = vertices.size();
    size_t numFaces = triangles.size();

    obj.reserve( numPoints * 32 + numFaces * 32 ); // preallocate memory

    const char* vFormat = "v %.6f %.6f %.6f\n";
    const char* fFormat = "f %zu// %zu// %zu//\n";

    obj.append( preamble );

    VertexExtractor< VertexType > vertexExtractor;

    for ( size_t i = 0; i < numPoints; ++i ) {

      glm::dvec3 t = vertexExtractor( vertices[ i ] );
 
      char vBuffer[ 128 ];
      
      snprintf( vBuffer, sizeof( vBuffer ), vFormat, t.x, t.y, t.z );
      
      obj.append(vBuffer);
    }

    for ( size_t i = 0; i < numFaces; ++i ) {
      
      const ConnectedTriangle& triangle = triangles[ i ];

      size_t f1 = triangle.vertices[ 0 ] + 1;
      size_t f2 = triangle.vertices[ 1 ] + 1;
      size_t f3 = triangle.vertices[ 2 ] + 1;

      char fBuffer[128];
    
      snprintf( fBuffer, sizeof( fBuffer ), fFormat, f1, f2, f3 );
      
      obj.append( fBuffer );
    
    }

    return obj;
  }

  template < typename VertexType >
  void WingedEdgeMesh< VertexType >::append( const WingedEdgeMesh< VertexType >& other ) {

    uint32_t vertexOffset = static_cast< uint32_t >( vertices.size() );
    
    bool hasBVH     = bvh.has_value();
    bool hasDipoles = hasBVH && bvh->hasDipoles();

    bvh.reset();

    vertices.insert( vertices.end(), other.vertices.begin(), other.vertices.end() );

    for ( const ConnectedTriangle& triangle : other.triangles ) {

      makeTriangle(
        triangle.vertices[ 0 ] + vertexOffset,
        triangle.vertices[ 1 ] + vertexOffset,
        triangle.vertices[ 2 ] + vertexOffset );
    }

    if ( hasBVH ) {

      AABBTree& localBVH = makeBVH();

      if ( hasDipoles ) {

        localBVH.dipoles( *this );
      }
    }
  }
}
