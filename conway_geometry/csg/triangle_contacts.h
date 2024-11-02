#pragma once 

#include "contact_pair.h"
#include "structures/fixed_stack.h"

namespace conway::geometry {

  
    struct TriangleContacts {

      FixedStack< ContactPair, 6 > pairs;

      void insertIfNotDuplicate( ContactRegion with, ContactRegion against ) {

        ContactPair candidate( with, against );
        uint8_t     candidateByte = std::bit_cast< uint8_t >( candidate );

        for ( ContactPair pair : pairs ) {

          uint8_t currentByte = std::bit_cast< uint8_t >( pair );

          if ( currentByte == candidateByte ) {

            return;
          }
        }
      }
     
      void faceVertex( uint32_t vertexInOtherTriangle ) {

        insertIfNotDuplicate( 
          ContactRegion::FACE,
          vertex( vertexInOtherTriangle ) );
      }

      void vertexFace( uint32_t vertexInThisTriangle ) {

        insertIfNotDuplicate( 
          vertex( vertexInThisTriangle ),
          ContactRegion::FACE );
      }

      void faceEdge( uint32_t edgeInOtherTriangle ) {

        insertIfNotDuplicate( 
          ContactRegion::FACE,
          edge( edgeInOtherTriangle ) );
      }

      void edgeFace( uint32_t edgeInThisTriangle ) {

        insertIfNotDuplicate( 
          edge( edgeInThisTriangle ),
          ContactRegion::FACE );
      }

      void edgeEdge( uint32_t edgeInThisTriangle, uint32_t edgeInOtherTriangle ) {

        insertIfNotDuplicate( 
          edge( edgeInThisTriangle ),
          edge( edgeInOtherTriangle ) );
      }
      
      void edgeVertex( uint32_t edgeInThisTriangle, uint32_t vertexInOtherTriangle ) {

        insertIfNotDuplicate( 
          edge( edgeInThisTriangle ),
          vertex( vertexInOtherTriangle ) );
      }

      void vertexEdge( uint32_t vertexInThisTriangle, uint32_t edgeInOtherTriangle ) {

        insertIfNotDuplicate( 
          vertex( vertexInThisTriangle ),
          edge( edgeInOtherTriangle ) );
      }

      void vertexVertex( uint32_t vertexInThisTriangle, uint32_t vertexInOtherTriangle ) {

        insertIfNotDuplicate( 
          vertex( vertexInThisTriangle ),
          vertex( vertexInOtherTriangle ) );
      }


      bool empty() const {
        return pairs.empty();
      }

      FaceFace face_to_face = FaceFace::NONE;

      uint32_t this_triangle_index;
      
      uint32_t other_triangle_index;
    };
}