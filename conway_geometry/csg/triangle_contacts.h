#pragma once 

#include "contact_pair.h"
#include "structures/fixed_stack.h"
#include "csg_utils.h"

namespace conway::geometry {

  enum class FaceFace : uint8_t {

    NONE                 = 0,
    COLINEAR             = 1,
    STRADDLES_OR_TOUCHES = 2 // faces straddle each other's planes.
  };
  
  struct TriangleContacts {

    FixedStack< ContactPair, 6 > pairs;

    void insertIfNotDuplicate( ContactRegion with, ContactRegion against, ContactType type ) {

      ContactPair candidate( with, against, type );

      std::byte candidateByte = std::bit_cast< std::byte >( candidate );

      for ( ContactPair pair : pairs ) {

        std::byte currentByte = std::bit_cast< std::byte >( pair );

        if ( currentByte == candidateByte ) {

          return;
        }
      }

      pairs.push( candidate );
    }
     
    void faceVertex( uint32_t vertexInOtherTriangle ) {

      insertIfNotDuplicate( 
        ContactRegion::FACE,
        vertex( vertexInOtherTriangle ),
        ContactType::COLINEAR );
    }

    void vertexFace( uint32_t vertexInThisTriangle ) {

      insertIfNotDuplicate( 
        vertex( vertexInThisTriangle ),
        ContactRegion::FACE,
        ContactType::COLINEAR );
    }

    void faceEdge( uint32_t edgeInOtherTriangle ) {

      insertIfNotDuplicate( 
        ContactRegion::FACE,
        edge( edgeInOtherTriangle ),
        ContactType::STRADDLING );
    }

    void edgeFace( uint32_t edgeInThisTriangle ) {

      insertIfNotDuplicate( 
        edge( edgeInThisTriangle ),
        ContactRegion::FACE,
        ContactType::STRADDLING );
    }

    void edgeEdge( uint32_t edgeInThisTriangle, uint32_t edgeInOtherTriangle ) {

      insertIfNotDuplicate( 
        edge( edgeInThisTriangle ),
        edge( edgeInOtherTriangle ),
        ContactType::STRADDLING );
    }

    void edgeEdge2D( uint32_t edgeInThisTriangle, uint32_t edgeInOtherTriangle ) {

      insertIfNotDuplicate( 
        edge( edgeInThisTriangle ),
        edge( edgeInOtherTriangle ),
        ContactType::COLINEAR );
    }

    void edgeVertex( uint32_t edgeInThisTriangle, uint32_t vertexInOtherTriangle ) {

      insertIfNotDuplicate( 
        edge( edgeInThisTriangle ),
        vertex( vertexInOtherTriangle ),
        ContactType::COLINEAR );
    }

    void vertexEdge( uint32_t vertexInThisTriangle, uint32_t edgeInOtherTriangle ) {

      insertIfNotDuplicate( 
        vertex( vertexInThisTriangle ),
        edge( edgeInOtherTriangle ),
        ContactType::COLINEAR );
    }

    void vertexVertex( uint32_t vertexInThisTriangle, uint32_t vertexInOtherTriangle ) {

      insertIfNotDuplicate( 
        vertex( vertexInThisTriangle ),
        vertex( vertexInOtherTriangle ),
        ContactType::COLINEAR );
    }

    bool empty() const {
      return pairs.empty();
    }

    FaceFace face_to_face = FaceFace::NONE;

    uint32_t this_triangle_index {};
      
    uint32_t other_triangle_index {};
  };
}