#pragma once

#include <stdint.h>

namespace conway::geometry {

  enum class ContactRegion : uint8_t {

    E0   = 0,
    E1   = 1,
    E2   = 2,
    FACE = 3,
    V0   = 4,
    V1   = 5,
    V2   = 6
  };
 
  constexpr uint32_t edgeIndex( ContactRegion a ) {
  
    return static_cast< uint32_t >( a );
  }

  constexpr uint32_t vertexIndex( ContactRegion a ) {
  
    return static_cast< uint32_t >( a ) - static_cast< uint32_t >( ContactRegion::V0 );
  }

  constexpr ContactRegion vertex0( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V0;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          static_cast< uint32_t >( a ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr ContactRegion vertex1( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V1;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          ( ( static_cast< uint32_t >( a ) + 1 ) % 3 ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr ContactRegion vertex2( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V1;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          ( ( static_cast< uint32_t >( a ) + 2 ) % 3 ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr bool isEdge( ContactRegion a ) {

    return a < ContactRegion::FACE;
  }

  constexpr bool isVertex( ContactRegion a ) {

    return a > ContactRegion::FACE;
  }

  constexpr bool isFace( ContactRegion a ) {

    return a == ContactRegion::FACE;
  }

  constexpr bool isSameEdge( ContactRegion a, ContactRegion b ) {

    if ( a == ContactRegion::FACE || b == ContactRegion::FACE ) {

      return false;
    }

    // if they are both edges, they must be the equal to be the same edge
    if ( a < ContactRegion::FACE && b < ContactRegion::FACE ) {
      
      return a == b;
    }

    // all vertices share an edge in a triangle
    if ( a > ContactRegion::FACE && b > ContactRegion::FACE ) {
    
      return true;
    }

    // this is an edge/vertex or vertex/edge comparison
    // make sure the edge is in a and the vertex in b.
    if ( a < b ) {

      std::swap( a, b );
    }

    return vertex0( a ) == b || vertex1( a ) == b;
  }

  static constexpr ContactRegion vertex( uint32_t vertexInFace ) {

    return static_cast< ContactRegion >( vertexInFace + static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  static constexpr ContactRegion edge( uint32_t edgeInFace ) {

    return static_cast< ContactRegion >( edgeInFace );
  }

  struct ContactPair {

    ContactRegion with : 3    = ContactRegion::FACE;
    ContactRegion against : 3 = ContactRegion::FACE;

    ContactPair() = default;

    ContactPair( const ContactPair& ) = default;

    ContactPair& operator=( const ContactPair& ) = default;
    
    ContactPair( ContactRegion _with, ContactRegion _against ) : with( _with ), against( _against ) {}
  
    bool isEdgeEdge() const {

      return isEdge( with ) && isEdge( against );
    }

    bool isFaceEdge() const {

      return isFace( with ) && isEdge( against );
    }

    bool isEdgeFace() const {

      return isEdge( with ) && isFace( against ) ;
    }

    bool isVertexVertex() const { 

      return isVertex( with ) && isVertex( against );
    }

    bool isFaceVertex() const { 

      return isFace( with ) && isVertex( against );
    }

    bool isVertexFace() const { 

      return isVertex( with ) && isFace( against );
    }
    
    bool isVertexEdge() const { 

      return isVertex( with ) && isEdge( against );
    }

    bool isEdgeVertex() const { 

      return isEdge( with ) && isVertex( against );
    }
    
    bool isVertexFaceOrEdge() const { 

      return isVertex( with ) && !isVertex( against );
    }
    
    bool isFaceOrEdgeVertex() const { 

      return !isVertex( with ) && isVertex( against );
    }
  };

}