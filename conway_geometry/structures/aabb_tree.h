#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include "box.h"
#include <bit>

namespace conway::geometry {

  template < typename VertexType >
  struct WingedEdgeMesh;

  class AABBTree {
  public:

    AABBTree( const WingedEdgeMesh< glm::dvec3 >& mesh );

    double gwn(
      const WingedEdgeMesh< glm::dvec3 >& mesh,
      const glm::dvec3& against );

    /** Intersect this AABB tree with another AABB tree and call back intersectionFunction for the
     * leaves of both trees that intersect, with the relative triangle indices.
     */
    template < typename OnIntersectionFunctionType >
    void intersect( const AABBTree& other, OnIntersectionFunctionType intersectionFunction ) const {

      static_assert(
        std::is_invocable_v< OnIntersectionFunctionType, uint32_t, uint32_t >,
        "Intersect requires a callable function that receives 2 uint32_t triangle indices." );

      if ( other.triangles_.empty() || this->triangles_.empty() ) {
        return;
      }

      // This possibly can be lowered in practice - CS
      Cursor stack[ MAX_RECURSION_DEPTH * 4 ];

      stack[ 0 ] = { 0, 1, 0, triangles_.size() };
      stack[ 1 ] = { 0, 1, 0, other.triangles_.size() };

      size_t end = 2;

      while ( end > 0 ) {

        const Cursor& otherCursor = stackThis[ --end ];
        const Cursor& thatCursor  = stackThis[ --end ];

        if ( !overlaps( boxes_[ thisCursor.node ], other.boxes_[ otherCursor.node ] ) ) {

          continue;
        }

        size_t thisSize  = thisCursor.end - thisCursor.start;
        size_t otherSize = otherCursor.end - otherCursor.start;

        if ( thisSize == 1 && otherSize == 1 ) {

          intersectionFunction( triangles_[ thisCursor.start ], other.triangles_[ otherCursor.start ] );
          continue;
        }

        if ( otherSize > thisSize ) {

          size_t partition      = std::bit_floor( otherSize );
          size_t leftChild      = otherCursor.children;
          size_t rightChild     = leftChild + 1;
          size_t leftGChild     = leftChild + 2;
          size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
          size_t partitionPoint = otherCursor.start + partition;
 
          // We skip over erasing thatCursor, so no aliasing issue.
          ++end;
          stack[ end++ ] = { rightChild, rightGChild, partitionPoint, otherCursor.end };
          stack[ end++ ] = thatCursor;
          stack[ end++ ] =  { leftChild, leftGChild, otherCursor.start, partitionPoint };
 
        } else {

          size_t partition      = std::bit_floor( thisSize );
          size_t leftChild      = thisCursor.children;
          size_t rightChild     = leftChild + 1;
          size_t leftGChild     = leftChild + 2;
          size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
          size_t partitionPoint = thisCursor.start + partition;

          // We skip over erasing otherCursor, so no aliasing issue.
          stack[ end++ ] = { leftChild, leftGChild, thisCursor.start, partitionPoint };
          ++end;
          stack[ end++ ] = { rightChild, rightGChild, partitionPoint, thisCursor.end };
          stack[ end++ ] = otherCursor;
        }
      }
    }

    /** Intersect this AABB tree with a single box and return the triangles (as indices) that
     * 
     */
    template < typename OnIntersectionFunctionType >
    void intersect( const box3& with, OnIntersectionFunctionType intersectionFunction ) const {

      static_assert(
        std::is_invocable_v< OnIntersectionFunctionType, uint32_t >,
        "Intersect requires a callable function that receives a uint32_t triangle indice." );

      if ( triangles_.empty() ) {
        return;
      }

      Cursor stack[ MAX_RECURSION_DEPTH ];

      stack[ 0 ] = { 0, 1, 0, triangles_.size() },;
      end        = 1;

      while ( end > 0 ) {

        const Cursor& thisCursor = stack[ --end ];
  
        if ( !overlaps( boxes_[ thisCursor.node ], with ) ) {

          continue;
        }

        size_t size  = cursor.end - cursor.start;

        if ( size == 1 ) {

          intersectionFunction( triangles_[ thisCursor.start ], other.triangles_[ otherCursor.start ] );
          continue;
        }

        size_t partition      = std::bit_floor( size );
        size_t leftChild      = cursor.children;
        size_t rightChild     = leftChild + 1;
        size_t leftGChild     = leftChild + 2;
        size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
        size_t partitionPoint = cursor.start + partition;

        stack_[ end++ ] = { rightChild, rightGChild, partitionPoint, thisCursor.end };
        stack_[ end++ ] = { leftChild, leftGChild, thisCursor.start, partitionPoint };
      }
    }

  private:

    static constexpr size_t MAX_RECURSION_DEPTH = 32;

    /** Cursor for traversing the tree */
    struct Cursor {

      size_t node;     // actual node
      size_t children; // offset of any children of node
      size_t start;    // the start of the span of leaves associated
      size_t end;      // the end of the span of leaves associated

    };

    /** Cursor for traversing the tree */
    struct CursorWithParent {

      size_t node;     // actual node
      size_t children; // offset of any children of node
      size_t start;    // the start of the span of leaves associated
      size_t end;      // the end of the span of leaves associated
      size_t parent;   // the parent node
      bool   walked;   // parent node has been walked.
    };

    /** Generate the dipoles to accelerate gwn calculation */
    void dipoles( 
      const WingedEdgeMesh< glm::dvec3 >& mesh );


    /** Get the size of a subtree given a number of leaves, including the root of the subtree */
    static size_t size( size_t leafCount ) {
      
      size_t remainder = leafCount;

      if ( remainder <= 1 ) {
        return remainder;
      }

      size_t nodeCount = 0;

      while ( true ) {
      
        size_t partitionPoint = std::bit_floor( leafCount );
        
        remainder -= partitionPoint;
        nodeCount += ( partitionPoint * 2 ) - 1;

        if ( remainder == partitionPoint ) {
          nodeCount += ( partitionPoint * 2 );
          break;
        } else {
          // Add a bridging parent for the right subtree + full left subtree.
          ++nodeCount;
        }
      }

      return nodeCount;
    }

    /**
     * Construct a subtree given a span leaves.
     */
    void construct( const WingedEdgeMesh< glm::dvec3 >& mesh, Cursor cursor, size_t depth );

    struct Dipole {

      /* The normal of the diapole is an area weighted sum of the normals below */
      glm::dvec3 normal;

      /* The position is the area weighted average of its children*/
      glm::dvec3 position;

      /** The area of this dipole */
      double area;

      /** Radius squared of this dipole */
      double radius2;
    };
    
    std::vector< uint32_t > triangles_;
    std::vector< box3 >     boxes_;
    std::vector< Dipole >   dipoles_;

    size_t depth_ = 0;
  };

}