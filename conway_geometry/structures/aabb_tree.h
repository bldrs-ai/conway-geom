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

    AABBTree( const WingedEdgeMesh< glm::dvec3 >& mesh, double tolerance = 0 );


    /** This GWN specifically works out  */
    double gwn(
      const WingedEdgeMesh< glm::dvec3 >& mesh,
      const glm::dvec3& against,
      std::optional< uint32_t > triangleInMeshIndex = std::nullopt ) const;

    /** Does this have dipoles */
    bool hasDipoles() const {

      return dipoles_.size() == boxes_.size();
    }

    /** Generate the dipoles for gwn calculation, if not called before GWN,
        GWN will always return 0. */
    void dipoles( const WingedEdgeMesh< glm::dvec3 >& mesh );

    struct Dipole {

      /* The normal of the diapole is an area weighted sum of the normals below */
      glm::dvec3 normal = glm::dvec3( 0.0 );

      /* The position is the area weighted average of its children*/
      glm::dvec3 position = glm::dvec3( 0.0 );

      /** The area of this dipole */
      double area = 0;

      /** Radius squared of this dipole */
      double radius2 = 0;
    };

    /** Get the representative dipole of this. (gwn or dipoles must have been called first).*/
    Dipole dipole() const {

      if (dipoles_.size() == 0) {

        return Dipole();
      }

      return dipoles_[0];
    }

    /** Get the root node box */
    box3 bounds() const {

      if (boxes_.size() == 0) {

        return box3();
      }

      return boxes_[ 0 ];
    }

    /** Intersect this AABB tree with itself, returning triangle pairs with overlapping boxes,
      * where the triangle pairs are not the same.
      */
    template < typename OnIntersectionFunctionType >
    std::invoke_result_t< OnIntersectionFunctionType, uint32_t, uint32_t > intersect(
        OnIntersectionFunctionType intersectionFunction ) const;

    /** Intersect this AABB tree with another AABB tree and call back intersectionFunction for the
     * leaves of both trees that intersect, with the relative triangle indices.
     */
    template < typename OnIntersectionFunctionType >
    std::invoke_result_t< OnIntersectionFunctionType, uint32_t, uint32_t > intersect(
      const AABBTree& other,
      OnIntersectionFunctionType intersectionFunction ) const;

    /** Intersect this AABB tree with a single box and return the triangles (as indices) that
     *
     * If the intersection function returns a boolean, it will terminate when "false" is returned (and return false).
     */
    template < typename OnIntersectionFunctionType >
    std::invoke_result_t< OnIntersectionFunctionType, uint32_t > intersect(
        const box3& with, OnIntersectionFunctionType intersectionFunction ) const;

  private:

    static constexpr size_t MAX_RECURSION_DEPTH = 32;

    /** Cursor for traversing the tree */
    struct Cursor {

      size_t node;     // actual node
      size_t start;    // the start of the span of leaves associated
      size_t end;      // the end of the span of leaves associated

    };

    /** Cursor for traversing the tree */
    struct CursorWithParent {

      size_t node;     // actual node
      size_t start;    // the start of the span of leaves associated
      size_t end;      // the end of the span of leaves associated
      size_t parent;   // the parent node
      bool   walked;   // parent node has been walked.
    };

    /** Dual cursor structure for mutual recursion of two trees (or this tree twice) */
    struct DualCursor {

      Cursor left;
      Cursor right;
    };

    /** Get the size of a subtree given a number of leaves, including the root of the subtree */
    static size_t size( size_t leafCount ) {
      
      size_t remainder = leafCount;

      if ( remainder <= 1 ) {
        return remainder;
      }

      size_t nodeCount = 0;

      while ( remainder > 1 ) {
      
        size_t partitionPoint = std::bit_floor( remainder - 1 );
        
        remainder -= partitionPoint;
        nodeCount += ( partitionPoint * 2 ) - 1;

        if ( remainder == partitionPoint ) {
          nodeCount += ( partitionPoint * 2 );
          remainder  = 0;
          break;
        } else {
          // Add a bridging parent for the right subtree + full left subtree.
          ++nodeCount;
        }
      }

      return nodeCount + remainder * 2;
    }

    /**
     * Construct a subtree given a span leaves.
     */
    void construct( const WingedEdgeMesh< glm::dvec3 >& mesh, double tolerance = 0 );
        
    std::vector< uint32_t > triangles_;
    std::vector< box3 >     boxes_;
    std::vector< Dipole >   dipoles_;
  };
}


/** Intersect this AABB tree with another AABB tree and call back intersectionFunction for the
 * leaves of both trees that intersect, with the relative triangle indices.
 */
template < typename OnIntersectionFunctionType >
std::invoke_result_t< OnIntersectionFunctionType, uint32_t, uint32_t >
conway::geometry::AABBTree::intersect(
    const AABBTree& other,
    OnIntersectionFunctionType intersectionFunction ) const {

  static_assert(
    std::is_invocable_v< OnIntersectionFunctionType, uint32_t, uint32_t >,
    "Intersect requires a callable function that receives 2 uint32_t triangle indices." );

  if ( other.triangles_.empty() || this->triangles_.empty() ) {
    return;
  }

  // This possibly can be lowered in practice
  // the "largest first" strategy should
  // put a bound on how far you can go down one tree and then the other,
  // as opposed to the worst case "all my nodes, then all your nodes".
  DualCursor stack[ MAX_RECURSION_DEPTH * 2 ];

  stack[ 0 ] = { { 0, 0, triangles_.size() }, { 0, 0, other.triangles_.size() } };

  size_t end = 1;

  while ( end > 0 ) {

    DualCursor&   cursor      = stack[ --end ];
    const Cursor& thisCursor  = cursor.left;
    const Cursor& otherCursor = cursor.right;

    if ( !overlaps( boxes_[ thisCursor.node ], other.boxes_[ otherCursor.node ] ) ) {

      continue;
    }

    size_t thisSize  = thisCursor.end - thisCursor.start;
    size_t otherSize = otherCursor.end - otherCursor.start;

    if ( thisSize == 1 && otherSize == 1 ) {

      if constexpr ( std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t, uint32_t >::value ) {
        
        if ( !intersectionFunction( triangles_[ thisCursor.start ], other.triangles_[ otherCursor.start ] ) ) {

          return false;
        }
      }
      else {

        intersectionFunction( triangles_[ thisCursor.start ], other.triangles_[ otherCursor.start ] );
      }

      continue;
    }

    if ( otherSize > thisSize ) {

      size_t partition      = std::bit_floor( otherSize - 1 );
      size_t leftChild      = otherCursor.node + 1;
      size_t rightChild     = otherCursor.node + ( partition * size_t( 2 ) ); // see note above 
      size_t partitionPoint = otherCursor.start + partition;
      size_t cursorStart    = otherCursor.start; 
      size_t cursorEnd      = otherCursor.end; 

      // We skip over erasing thisCursor, so no aliasing issue.
      cursor.right = { rightChild, partitionPoint, cursorEnd } ;

      ++end;

      stack[ end ].right = { leftChild, cursorStart, partitionPoint };
      stack[ end ].left  = thisCursor;

      ++end;

    } else {

      size_t partition      = std::bit_floor( thisSize - 1 );
      size_t leftChild      = thisCursor.node + 1;
      size_t rightChild     = thisCursor.node + ( partition * size_t( 2 ) ); // see note above 
      size_t partitionPoint = thisCursor.start + partition;
      size_t cursorStart    = thisCursor.start; 
      size_t cursorEnd      = thisCursor.end;

      // We skip over erasing otherCursor, so no aliasing issue.
      cursor.left = { rightChild, partitionPoint, cursorEnd };
      
      ++end;
      
      stack[ end ].left  = { leftChild, cursorStart, partitionPoint };
      stack[ end ].right = otherCursor;

      ++end;
    }
  }

  if constexpr ( std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t >::value ) {

    return true;
  }
}



/** Intersect this AABB tree with another AABB tree and call back intersectionFunction for the
 * leaves of both trees that intersect, with the relative triangle indices.
 */
template < typename OnIntersectionFunctionType >
std::invoke_result_t< OnIntersectionFunctionType, uint32_t, uint32_t >
conway::geometry::AABBTree::intersect(
    OnIntersectionFunctionType intersectionFunction ) const {

  static_assert(
    std::is_invocable_v< OnIntersectionFunctionType, uint32_t, uint32_t >,
    "Intersect requires a callable function that receives 2 uint32_t triangle indices." );

  if ( triangles_.empty() ) {
    return;
  }

  // This possibly can be lowered in practice
  // the "largest first" strategy should
  // put a bound on how far you can go down one tree and then the other,
  // as opposed to the worst case "all my nodes, then all your nodes".
  DualCursor stack[ MAX_RECURSION_DEPTH * 2 ];

  stack[ 0 ] = { { 0, 0, triangles_.size() }, { 0, 0, triangles_.size() } };

  size_t end = 1;

  while ( end > 0 ) {

    DualCursor&   cursor      = stack[ --end ];
    const Cursor& thisCursor  = cursor.left;
    const Cursor& otherCursor = cursor.right;

    if (
      otherCursor.end <= thisCursor.start ||
      !overlaps( boxes_[ thisCursor.node ], boxes_[ otherCursor.node ] ) ) {

      continue;
    }

    size_t thisSize  = thisCursor.end - thisCursor.start;
    size_t otherSize = otherCursor.end - otherCursor.start;

    if ( thisSize == 1 && otherSize == 1 ) {

      if ( thisCursor.node == otherCursor.node ) {
        continue;
      }

      if constexpr ( std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t, uint32_t >::value ) {
        
        if ( !intersectionFunction( triangles_[ thisCursor.start ], triangles_[ otherCursor.start ] ) ) {

          return false;
        }
      }
      else {

        intersectionFunction( triangles_[ thisCursor.start ], triangles_[ otherCursor.start ] );
      }

      continue;
    }

    if ( otherSize > thisSize ) {

      size_t partition      = std::bit_floor( otherSize - 1 );
      size_t leftChild      = otherCursor.node + 1;
      size_t rightChild     = otherCursor.node + ( partition * size_t( 2 ) ); // see note above 
      size_t partitionPoint = otherCursor.start + partition;
      size_t cursorStart    = otherCursor.start; 
      size_t cursorEnd      = otherCursor.end; 

      // We skip over erasing thisCursor, so no aliasing issue.
      cursor.right = { rightChild, partitionPoint, cursorEnd } ;

      ++end;

      stack[ end ].right = { leftChild, cursorStart, partitionPoint };
      stack[ end ].left  = thisCursor;

      ++end;

    } else {

      size_t partition      = std::bit_floor( thisSize - 1 );
      size_t leftChild      = thisCursor.node + 1;
      size_t rightChild     = thisCursor.node + ( partition * size_t( 2 ) ); // see note above 
      size_t partitionPoint = thisCursor.start + partition;
      size_t cursorStart    = thisCursor.start; 
      size_t cursorEnd      = thisCursor.end;

      // We skip over erasing otherCursor, so no aliasing issue.
      cursor.left = { rightChild, partitionPoint, cursorEnd };
      
      ++end;
      
      stack[ end ].left  = { leftChild, cursorStart, partitionPoint };
      stack[ end ].right = otherCursor;

      ++end;
    }
  }

  if constexpr ( std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t >::value ) {

    return true;
  }
}


template < typename OnIntersectionFunctionType >
std::invoke_result_t< OnIntersectionFunctionType, uint32_t >
conway::geometry::AABBTree::intersect(
    const box3& with, OnIntersectionFunctionType intersectionFunction ) const {

  static_assert(
    std::is_invocable_v< OnIntersectionFunctionType, uint32_t >,
    "Intersect requires a callable function that receives a uint32_t triangle indice." );

  if ( triangles_.empty() ) {
    return;
  }

  Cursor stack[ MAX_RECURSION_DEPTH ];

  stack[ 0 ] = { 0, 0, triangles_.size() };

  size_t end = 1;

  while ( end > 0 ) {

    const Cursor& cursor = stack[ --end ];

    if ( !overlaps( boxes_[ cursor.node ], with ) ) {

      continue;
    }

    size_t size  = cursor.end - cursor.start;

    if ( size == 1 ) {

      if constexpr (std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t >::value) {
      
        if ( !intersectionFunction( triangles_[ cursor.start ] ) ) {

            return false;
        }
      }
      else {

        intersectionFunction( triangles_[ cursor.start ] );
      }
      continue;
    }

    size_t partition      = std::bit_floor( size - 1 );
    size_t leftChild      = cursor.node + 1;
    size_t rightChild     = cursor.node + ( partition * size_t( 2 ) ); // see note above
    // copy because this cursor will be overwritten.
    size_t cursorStart    = cursor.start;
    size_t cursorEnd      = cursor.end;
    size_t partitionPoint = cursor.start + partition;

    stack[ end++ ] = { rightChild, partitionPoint, cursorEnd };
    stack[ end++ ] = { leftChild, cursorStart, partitionPoint };
  }

  if constexpr ( std::is_invocable_r< bool, OnIntersectionFunctionType, uint32_t >::value ) {

    return true;
  }
}

