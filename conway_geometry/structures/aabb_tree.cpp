#include "csg/csg_utils.h"
#include "aabb_tree.h"
#include "winged_edge.h"

constexpr double GWN_PRECISION              = 3.0;
constexpr double INVERSE_AREA_GWN_PRECISION = 1.0 / GWN_PRECISION;

namespace conway::geometry {

  box3 make_box( const WingedEdgeMesh< glm::dvec3 >& mesh, size_t triangleIndex, double tolerance = 0 ) {

    const ConnectedTriangle& triangle = mesh.triangles[ triangleIndex ];

    const glm::dvec3& v0 = mesh.vertices[ triangle.vertices[ 0 ] ];

    box3 result = { v0, v0 };

    result.merge( mesh.vertices[ triangle.vertices[ 1 ] ] );
    result.merge( mesh.vertices[ triangle.vertices[ 2 ] ] );

    glm::dvec3 toleranceBox( tolerance ); 

    result.min -= toleranceBox;
    result.max += toleranceBox;

    return result;
  }
}

    /** Generate the dipoles for this */
void conway::geometry::AABBTree::dipoles(
  const WingedEdgeMesh< glm::dvec3 >& mesh ) {

  const std::vector< ConnectedTriangle >& triangles = mesh.triangles;
  const std::vector< glm::dvec3 >&        vertices  = mesh.vertices;

  double result = 0;

  dipoles_.clear();
  dipoles_.resize( boxes_.size(), {} );

  CursorWithParent stack[ MAX_RECURSION_DEPTH * 2 ];

  stack[ 0 ] = { 0, 0, triangles_.size(), 0 };

  size_t end = 1;

  while ( end > 0 ) {

    CursorWithParent& cursor = stack[ end - 1 ];

    if ( cursor.walked ) {

      if ( cursor.parent != cursor.node ) {

        const Dipole& dipole = dipoles_[ cursor.node ];
        Dipole& parentDipole = dipoles_[ cursor.parent ];

        parentDipole.normal   += dipole.normal;
        parentDipole.area     += dipole.area;
        parentDipole.position += dipole.position;
      }

      --end;
      continue;
    }

    cursor.walked = true;

    size_t size = cursor.end - cursor.start;

    if ( size == 1 ) {
  
      Dipole& dipole = dipoles_[ cursor.node ];

      const ConnectedTriangle& triangle = mesh.triangles[ cursor.start ];

      const glm::dvec3& v0 = mesh.vertices[ triangle.vertices[ 0 ] ];
      const glm::dvec3& v1 = mesh.vertices[ triangle.vertices[ 1 ] ];
      const glm::dvec3& v2 = mesh.vertices[ triangle.vertices[ 2 ] ];

      dipole.normal   = 0.5 * glm::cross( v1 - v0, v2 - v0 );
      dipole.area     = glm::length( dipole.normal );
      dipole.position = ( ( v0 + v1 + v2 ) / 3.0 ) * dipole.area;

      continue;
    }

    size_t partition      = std::bit_floor( size - 1 );
    size_t leftChild      = cursor.node + 1;
    size_t rightChild     = cursor.node + ( partition * size_t(2) );
    size_t partitionPoint = cursor.start + partition;
    
    stack[ end++ ] = { leftChild, cursor.start, partitionPoint, cursor.node, false };
    stack[ end++ ] = { rightChild, partitionPoint, cursor.end, cursor.node, false };
  }

  for (
    size_t node = 0, nodeEnd = boxes_.size(); 
    node < nodeEnd;
    ++node ) {
  
    Dipole&     dipole = dipoles_[ node ];
    const box3& box    = boxes_[ node ];

    dipole.position /= dipole.area;

    dipole.radius2 = length2( glm::max( box.max - dipole.position, dipole.position - box.min ) );

    dipole.radius2 *= GWN_PRECISION;
  }
}

conway::geometry::AABBTree::AABBTree( const WingedEdgeMesh< glm::dvec3 >& mesh ) {

  size_t leafCount = mesh.triangles.size();

  if ( leafCount == 0 ) {
    return;
  }

  triangles_.resize( leafCount );
  
  std::iota( triangles_.begin(), triangles_.end(), 0 );

  size_t nodeCount = size( leafCount );

  boxes_.resize( nodeCount );

  construct( mesh );
}


double conway::geometry::AABBTree::gwn( 
  const WingedEdgeMesh< glm::dvec3 >& mesh,
  const glm::dvec3& against ) const {

  const std::vector< ConnectedTriangle >& triangles = mesh.triangles;
  const std::vector< glm::dvec3 >&        vertices  = mesh.vertices;

  double result = 0;

  Cursor stack[ MAX_RECURSION_DEPTH ];

  stack[ 0 ] = { 0, 0, triangles_.size() };

  size_t end = 1;

  while ( end > 0 ) {

    const Cursor& cursor = stack[ --end ];

    size_t size = cursor.end - cursor.start;
    
    const Dipole& dipole = dipoles_[ cursor.node ];

    glm::dvec3 againstDipole = dipole.position - against;
    double     aD2           = length2( againstDipole );
    bool calculateSubArea    = false;

    if ( aD2 > dipole.radius2 ) {

      result += glm::dot( againstDipole, dipole.normal ) / ( aD2 * sqrt( aD2 ) );
      continue; 
    }

    if ( size == 1 ) {

      const ConnectedTriangle& triangle = triangles[ triangles_[ cursor.start ] ];

      result += 
        triangle_solid_angle(
          vertices[ triangle.vertices[ 0 ] ],
          vertices[ triangle.vertices[ 1 ] ],
          vertices[ triangle.vertices[ 2 ] ],
          against );

      continue;
    }

    assert( size > 0 );

    size_t partition      = std::bit_floor( size - 1 );
    size_t leftChild      = cursor.node + 1;
    size_t rightChild     = cursor.node + ( partition * size_t( 2 ) );
    // copy because this cursor will be overwritten.
    size_t cursorStart    = cursor.start;
    size_t cursorEnd      = cursor.end;
    size_t partitionPoint = cursor.start + partition;
    
    stack[ end++ ] = { leftChild, cursorStart, partitionPoint };
    stack[ end++ ] = { rightChild, partitionPoint, cursorEnd };
  }

  return result / ( 4.0 * M_PI );
}

void conway::geometry::AABBTree::construct( const WingedEdgeMesh< glm::dvec3 >& mesh ) {

  const std::vector< ConnectedTriangle >& triangles = mesh.triangles;
  const std::vector< glm::dvec3 >&        vertices  = mesh.vertices;

  Cursor stack[ MAX_RECURSION_DEPTH ];

  stack[ 0 ] = { 0, 0, triangles_.size() };

  size_t end    = 1;
  size_t visits = 0;

  while ( end > 0 ) {

    const Cursor& cursor = stack[ --end ];
    size_t        size   = cursor.end - cursor.start;

    ++visits;

    assert( visits <= boxes_.size() );

    if (size == 1) {

      box3& nodeBox = boxes_[ cursor.node ] = make_box( mesh, cursor.start );

      nodeBox.min -= glm::dvec3( DBL_EPSILON );
      nodeBox.max += glm::dvec3( DBL_EPSILON );
      continue;
    }

    // Here instead of getting the "mid" point we partition to make it so the left sub-tree
    // will be the highest power of 2 lower than the size of the span.
    // This means the left subtree will always be binary-complete (i.e. a perfect completely full balanced
    // binary tree), and the right subtree will be the "left-overs".
    // Because we only ever need to skip over the *left* subtree and left subtrees are always full,
    // it means on the left we have to skip over ( partition point * 2 ) - 1 nodes on the left side.
    size_t partition      = std::bit_floor( size - 1 );
    size_t leftChild      = cursor.node + 1;
    size_t rightChild     = cursor.node + ( partition * size_t( 2 ) );
    size_t partitionPoint = cursor.start + partition;

    glm::dvec3 intervalIntegral = glm::dvec3( 0 );

    box3& nodeBox = boxes_[ cursor.node ];

    // We could recursively merge boxes and it would strictly be logically better
    // but one of the problems with fuzzy bools is the alternating axis split heuristic
    // which can lead to pathological cases where parents and children have similar node boxes.
    // Which make it much more likely to degrade to the N^2 everything overlaps everything
    // worst case.
    // Here we construct the box ahead of time using the full run of triangles,
    // because we want to use that information to work out which axis to partition by.
    // It's also significantly cheaper than estimating surface area heuristic,
    // which while it would be better again, is not *that* much better.
    for (
      auto where = triangles_.begin() + cursor.start,
      trianglesEnd = triangles_.begin() + cursor.end;
      where < trianglesEnd;
      ++where) {

      box3 triangleBox = make_box( mesh, *where );

      nodeBox.merge( triangleBox );

      intervalIntegral += triangleBox.max - triangleBox.min;
    }

    // The actual heuristic is based around the size of the average interval
    // of the box of the triangle in an axis vs the width of the box in that axis.
    // where the box is large vs the triangles in it, likely it will 
    // partition well... if intervals are similar, it degrades to longest axis.
    glm::dvec3 axisHeuristic = ( nodeBox.max - nodeBox.min ) / intervalIntegral;

    double   bestValue = 0;
    uint32_t axis      = 0;

    for ( uint32_t candidateAxis = 0; candidateAxis < 3; ++candidateAxis ) {

      if ( intervalIntegral[ candidateAxis ] > 0 && axisHeuristic[ candidateAxis ] > bestValue ) {

        bestValue = axisHeuristic[ candidateAxis ];
        axis      = candidateAxis;
      }
    }

    nodeBox.min -= glm::dvec3( DBL_EPSILON );
    nodeBox.max += glm::dvec3( DBL_EPSILON );

    if ( size > 2 ) {

      // This will cluster by centroid on "axis"
      // by using the partitionth element to separate the items.
      std::nth_element(
        triangles_.begin() + cursor.start,
        triangles_.begin() + ( partitionPoint ),
        triangles_.begin() + cursor.end,
        [ &, axis ]( uint32_t left, uint32_t right ) {

          const ConnectedTriangle& leftTriangle  = triangles[ left ];
          const ConnectedTriangle& rightTriangle = triangles[ right ];

          double leftA = 
            vertices[ leftTriangle.vertices[ 0 ] ][ axis ] +
            vertices[ leftTriangle.vertices[ 1 ] ][ axis ] +
            vertices[ leftTriangle.vertices[ 2 ] ][ axis ];

          double rightA =
            vertices[ rightTriangle.vertices[ 0 ] ][ axis ] +
            vertices[ rightTriangle.vertices[ 1 ] ][ axis ] +
            vertices[ rightTriangle.vertices[ 2 ] ][ axis ];

          return leftA < rightA;
        });
    }

    size_t cursorStart = cursor.start;
    size_t cursorEnd   = cursor.end;

    stack[ end++ ] = { rightChild, partitionPoint, cursorEnd };
    stack[ end++ ] = { leftChild, cursorStart, partitionPoint };
  }

  assert( visits == boxes_.size() );
}