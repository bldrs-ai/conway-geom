#include "csg/csg_utils.h"
#include "aabb_tree.h"
#include "winged_edge.h"

constexpr double GWN_PRECISION              = 1.5;
constexpr double INVERSE_AREA_GWN_PRECISION = 1.0 / GWN_PRECISION;

    /** Generate the dipoles for this */
void conway::geometry::AABBTree::dipoles(
  const WingedEdgeMesh< glm::dvec3 >& mesh ) {

  const std::vector< ConnectedTriangle >& triangles = mesh.triangles;
  const std::vector< glm::dvec3 >&        vertices  = mesh.vertices;

  double result = 0;

  CursorWithParent stack[ MAX_RECURSION_DEPTH ];

  stack[ 0 ] = { 0, 1, 0, triangles_.size() };

  size_t end = 1;

  while ( end > 0 ) {

    CursorWithParent& cursor = stack[ end ];

    if ( cursor.walked ) {

      if ( cursor.parent != cursor.node ) {

          const Dipole& dipole = dipoles_[ cursor.node ];
          Dipole& parentDipole = dipoles_[ cursor.node ];

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

      dipole.normal   = glm::cross( v1 - v0, v2 - v0 );
      dipole.area     = 0.5 * glm::length( dipole.normal );
      dipole.position = ( ( v0 + v1 + v2 ) / 3.0 ) * dipole.area;

      continue;
    }

    size_t partition      = std::bit_floor( size );
    size_t leftChild      = cursor.children;
    size_t rightChild     = leftChild + 1;
    size_t leftGChild     = leftChild + 2;
    size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
    size_t partitionPoint = cursor.start + partition;
    
    stack[ end++ ] = { leftChild, leftGChild, cursor.start, partitionPoint, cursor.node, false };
    stack[ end++ ] = { rightChild, leftGChild, partitionPoint, cursor.end, cursor.node, false };
  }



  for (
    size_t node = 0, nodeEnd = boxes_.size(); 
    node < end;
    ++node ) {
  
    Dipole&     dipole = dipoles_[ node ];
    const box3& box    = boxes_[ node ];

    dipole.position /= dipole.area;

    double radius2 = length2( glm::max( box.max - dipole.position, dipole.position - box.min ) );

    // Not in the original paper, but there is a case
    // where a large high curvature piece of surface has
    // a dipole centred in the box, giving it the smallest
    // radius possible, while actually making it a worse
    // approximation up close.
    // In this case we up the radius^2 to the area (approximating the disk)
    // or effectively square the precision factor if that's smaller.
    if ( dipole.area > radius2 ) {

      double cappedArea = GWN_PRECISION * radius2;

      dipole.radius2 = cappedArea < dipole.area ? cappedArea : dipole.area;
    }

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

  construct( mesh, { 0, 1, 0, leafCount }, 0 );
}

/** Inner intersection function that co-traverses two AABB trees */
double conway::geometry::AABBTree::gwn( 
  const WingedEdgeMesh< glm::dvec3 >& mesh,
  const glm::dvec3& against ) {

  const std::vector< ConnectedTriangle >& triangles = mesh.triangles;
  const std::vector< glm::dvec3 >&        vertices  = mesh.vertices;

  double result = 0;

  Cursor stack[ MAX_RECURSION_DEPTH ];

  stack[ 0 ] = { 0, 1, 0, triangles_.size() };

  size_t end = 1;

  while ( end > 0 ) {

    const Cursor& cursor = stack[ --end ];

    size_t size = cursor.end - cursor.start;
    
    if ( size == 1 ) {

      continue;
    }

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

    size_t partition      = std::bit_floor( size );
    size_t leftChild      = cursor.children;
    size_t rightChild     = leftChild + 1;
    size_t leftGChild     = leftChild + 2;
    size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
    size_t partitionPoint = cursor.start + partition;
    
    stack[ end++ ] = { leftChild, leftGChild, cursor.start, partitionPoint };
    stack[ end++ ] = { rightChild, leftGChild, partitionPoint, cursor.end };
  }
}

void conway::geometry::AABBTree::construct( const WingedEdgeMesh< glm::dvec3 >& mesh, Cursor cursor, size_t depth ) {
  
  size_t size = cursor.end - cursor.start;

  ++depth;

  if ( depth > depth_ ) {

    depth_ = depth;
  }

  if ( size == 1 ) {

    boxes_[ cursor.node ] = make_box( mesh, cursor.start );
    return;
  }
  
  // Here instead of getting the "mid" point we partition to make it so the left sub-tree
  // will be the highest power of 2 lower than the size of the span.
  // This means the left subtree will always be binary-complete (i.e. a perfect completely full balanced
  // binary tree), and the right subtree will be the "left-overs".
  // Because we only ever need to skip over the *left* subtree and left subtrees are always full,
  // it means on the left we have to skip over ( partition point * 2 ) - 1 nodes on the left side,
  // but because we keep the right child local, we have to skip that too.
  size_t partition = std::bit_floor( size );

  size_t leftChild      = cursor.children;
  size_t rightChild     = leftChild + 1;
  size_t leftGChild     = leftChild + 2;
  size_t rightGChild    = leftChild + ( partition * size_t( 2 ) ); // see note above 
  size_t partitionPoint = cursor.start + partition;

  glm::dvec3 intervalIntegral = glm::dvec3( 0 );

  box3 nodeBox;

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
         end   = triangles_.begin() + cursor.end;
    where < end; 
    ++where ) {

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

  boxes_[ cursor.node ] = nodeBox;

  // This will cluster by by centroid on "axis"
  // by using the partitionth element to separate the items.
  // Consider choosing an axis via heuristic here first instead. - CS
  std::nth_element( 
    triangles_.begin() + cursor.start,
    triangles_.begin() + partitionPoint,
    triangles_.begin() + cursor.end, 
    [=]( uint32_t left, uint32_t right ) {

      const ConnectedTriangle& leftTriangle  = mesh.triangles[ left ];
      const ConnectedTriangle& rightTriangle = mesh.triangles[ right ];

      double leftA = mesh.vertices[ leftTriangle.vertices[ 0 ] ][ axis ] + 
        mesh.vertices[ leftTriangle.vertices[ 1 ] ][ axis ] + 
        mesh.vertices[ leftTriangle.vertices[ 2 ] ][ axis ];
        
      double rightA = mesh.vertices[ rightTriangle.vertices[ 0 ] ][ axis ] + 
        mesh.vertices[ rightTriangle.vertices[ 1 ] ][ axis ] + 
        mesh.vertices[ rightTriangle.vertices[ 2 ] ][ axis ];

        return leftA < rightA;
    });

  construct( mesh, { leftChild, leftGChild, cursor.start, partitionPoint }, depth );
  construct( mesh, { rightChild, rightGChild, partitionPoint, cursor.end }, depth );
}