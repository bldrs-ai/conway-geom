#include "csg_mesher.h"
#include "structures/winged_edge.h"
#include "multi_mesh_vertex_index.h"

#include <execution>
#include <algorithm>
#include <format>
#include <sstream>

#include "representation/SVGContext.h"
#include "structures/thread_pool.h"

constexpr bool OUTPUT_A                        = true;
constexpr bool OUTPUT_B                        = true;
constexpr bool OUTPUT_BOUNDARY                 = true;
constexpr bool OUTPUT_ORIGINAL                 = true;
constexpr bool OUTPUT_SHARED_FACES             = true;
constexpr bool DUMP_SVGS_ON_TRIANGLE_ERROR     = false;
constexpr size_t MAXIMUM_TRIANGULATION_RETRIES = 16;

constexpr double WELD_EDGE_TOLERANCE = ( 1.0 / ( double ( 1 << 24 ) ) ) * ( 1.0 / ( double ( 1 << 24 ) ) ) * 0.5;

constexpr double GWN_TOLERANCE          = 0.5 + DBL_EPSILON;

constexpr bool OUTPUT_A_BOUNDARY = OUTPUT_A && OUTPUT_BOUNDARY;
constexpr bool OUTPUT_B_BOUNDARY = OUTPUT_B && OUTPUT_BOUNDARY;
constexpr bool OUTPUT_A_ORIGINAL = OUTPUT_A && OUTPUT_ORIGINAL;
constexpr bool OUTPUT_B_ORIGINAL = OUTPUT_B && OUTPUT_ORIGINAL;
constexpr bool OUTPUT_A_SHARED   = OUTPUT_A && OUTPUT_SHARED_FACES;

namespace conway::geometry {

  namespace {
  
    inline bool lexicographical_less( const glm::dvec3& a, const glm::dvec3& b ) {

      return a.x < b.x || ( a.x == b.x && a.y < b. y);
    }
  }
}


void conway::geometry::CSGMesher::process(
  Geometry& a, // these are not const because we lazily generate dipoles.
  const std::vector< TriangleTriangleContactPair >& contacts,
  const PrefixSumMap& aContactMap,
  const std::vector< bool >& boundarySet,
  Geometry& output ) {

  reset();

  unifiedPlanes_.allocate( static_cast< uint32_t >( a.triangles.size() ) ); 

  // We don't need to know the number of novel vertices yet, because it's the last in the partition.
  MultiMeshVertexIndex< 2 > vertices = multi_mesh_vertex_index( a, novelVertices_ );

  unifiedVertices_.allocate( static_cast< uint32_t >( a.vertices.size() ) );

  // unify vertex pairs.
  for( const TriangleTriangleContactPair& trianglePair : contacts ) {

    const TriangleContacts& triangleContacts = trianglePair.triangles[ 0 ];

    if ( triangleContacts.face_to_face == FaceFace::COLINEAR ) {
    
      unifiedPlanes_.merge(
        triangleContacts.this_triangle_index,
        triangleContacts.other_triangle_index );
    }

    for ( ContactPair contact : triangleContacts.pairs ) {
      
      if ( contact.isVertexVertex() ) {

        const Triangle& aTriangle = a.triangles[ triangleContacts.this_triangle_index ];
        const Triangle& bTriangle = a.triangles[ triangleContacts.other_triangle_index ];

        uint32_t vertexA = aTriangle.vertices[ vertexIndex( contact.with ) ];
        uint32_t vertexB = bTriangle.vertices[ vertexIndex( contact.against ) ];

        unifiedVertices_.merge( vertexA, vertexB );
      }
    }
  }

  // Optimise here, because this will make the mappings above
  // single hop, but the novel vertices don't need to be, because they self reference
  unifiedVertices_.optimize();
  unifiedPlanes_.optimize();

  const std::vector< glm::dvec3 >& aVertices = a.vertices;

  for (
    uint32_t aTriangleIndex = 0,
             aTriangleEnd   = static_cast< uint32_t >( a.triangles.size() );
    aTriangleIndex < aTriangleEnd;
    ++aTriangleIndex ) {

    std::span< const uint32_t > trianglePairs = aContactMap.get( aTriangleIndex );

    if ( !boundarySet[ aTriangleIndex ] || trianglePairs.empty() ) {
      continue;
    }

    const Triangle& aTriangle           = a.triangles[ aTriangleIndex ];
    const TriangleEdges& aTriangleEdges = a.triangle_edges[ aTriangleIndex ];

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {
      insertLocalVertex( aTriangle.vertices[ vertexInTriangle ] );
    }

    uint32_t aPlaneIndex = unifiedPlanes_.find( aTriangleIndex );

    for ( uint32_t pairIndex : trianglePairs ) {

      uint32_t triangleInPair = contacts[ pairIndex ].triangles[ 0 ].this_triangle_index == aTriangleIndex ? 0 : 1;

      const TriangleContacts&  trianglePair = contacts[ pairIndex ].triangles[ triangleInPair ];
      const Triangle& bTriangle             = a.triangles[ trianglePair.other_triangle_index ];
      const TriangleEdges& bTriangleEdges   = a.triangle_edges[ trianglePair.other_triangle_index ];

      assert( aTriangleIndex == trianglePair.this_triangle_index );

      std::span< const ContactPair > contactPairs = trianglePair.pairs.values();

      FixedStack< uint32_t, 6 > additionalVertices;

      uint32_t bPlaneIndex = unifiedPlanes_.find( trianglePair.other_triangle_index );

      for ( ContactPair contact : contactPairs ) {

        if ( contact.isEdgeEdge() ) {

          uint32_t edgeInTriangleA = edgeIndex( contact.with );
          uint32_t edge0           = aTriangleEdges.edges[ edgeInTriangleA ];
          uint32_t edge1           = bTriangleEdges.edges[ edgeIndex( contact.against ) ];

          uint32_t canidateIndex = vertices( 1, static_cast<uint32_t>( novelVertices_.size() ) );

          if ( edge0 > edge1 ) {
          
            std::swap( edge0, edge1 );
          }

          const auto [ to, success ] =
            edgeEdgeVertices_.try_emplace(
              std::make_pair( edge0, edge1 ),
              canidateIndex );

          if ( success ) {

            glm::dvec3 e0v0 = aVertices[ a.edges[ edge0 ].vertices[ 0 ] ];
            glm::dvec3 e0v1 = aVertices[ a.edges[ edge0 ].vertices[ 1 ] ];
            glm::dvec3 e1v0 = aVertices[ a.edges[ edge1 ].vertices[ 0 ] ];
            glm::dvec3 e1v1 = aVertices[ a.edges[ edge1 ].vertices[ 1 ] ];

            // while we can't guarantee exact results for this operation,
            // knowing there are some duplicate-vertices-abutting-faces cases
            // means that we can at least make it deterministic by ordering 
            if ( lexicographical_less( e0v1, e0v0 ) ) {

              std::swap( e0v1, e0v0 );
            }
            
            if ( lexicographical_less( e1v1, e1v0 ) ) {

              std::swap( e1v1, e1v0 );
            }

            glm::dvec3 novelVertex =
              line_segment_line_segment_intersection( e0v0, e0v1, e1v0, e1v1 );

            novelVertices_.push_back( novelVertex );
            unifiedVertices_.allocate();
          }

          additionalVertices.push(
            insertLocalVertexOnEdge(
              to->second,
              edgeInTriangleA ) );

        } else if ( contact.isFaceEdge() ) {              

          uint32_t edge          = bTriangleEdges.edges[ edgeIndex( contact.against ) ];
          uint32_t canidateIndex = vertices( 1, static_cast<uint32_t>( novelVertices_.size() ) );

          const auto [ to, success ] =
            faceEdgeVertices_[ 0 ].try_emplace(
              std::make_pair( aPlaneIndex, edge ),
              canidateIndex );

          if ( success ) {            

            const glm::dvec3& t0 = aVertices[ aTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = aVertices[ aTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = aVertices[ aTriangle.vertices[ 2 ] ];

            glm::dvec3 ev0 = aVertices[ a.edges[ edge ].vertices[ 0 ] ];
            glm::dvec3 ev1 = aVertices[ a.edges[ edge ].vertices[ 1 ] ];

            // Because these vertices might be on a different edge without 
            // the same identity in some cases, we sort to make this deterministic...
            // although, the triangle itself should be deterministic
            if ( lexicographical_less( ev1, ev0 ) ) {

              std::swap( ev0, ev1 );
            }

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            novelVertices_.push_back( novelVertex ); 

            unifiedVertices_.allocate();
          }

          additionalVertices.push( insertLocalVertex( to->second ) );

        } else if ( contact.isEdgeFace() ) {

          uint32_t edgeInTriangle = edgeIndex( contact.with );
          uint32_t edge           = aTriangleEdges.edges[ edgeInTriangle ];

          uint32_t canidateIndex = vertices( 1, static_cast<uint32_t>( novelVertices_.size() ) );

          const auto [ to, success ] =
            faceEdgeVertices_[ 0 ].try_emplace(
              std::make_pair( bPlaneIndex, edge ),
              canidateIndex );

          if ( success ) {

            const glm::dvec3& t0 = aVertices[ bTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = aVertices[ bTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = aVertices[ bTriangle.vertices[ 2 ] ];

            glm::dvec3 ev0 = aVertices[ a.edges[ edge ].vertices[ 0 ] ];
            glm::dvec3 ev1 = aVertices[ a.edges[ edge ].vertices[ 1 ] ];

            // Because these vertices might be on a different edge without 
            // the same identity in some cases, we sort to make this deterministic...
            // although, the triangle itself should be deterministic
            if ( lexicographical_less( ev1, ev0 ) ) {

              std::swap( ev0, ev1 );
            }

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            novelVertices_.push_back( novelVertex );
            unifiedVertices_.allocate();
          }

          additionalVertices.push(
            insertLocalVertexOnEdge( 
              to->second,
              edgeInTriangle ) );

        } else if ( contact.isFaceVertex() ) {

          additionalVertices.push( insertLocalVertex( bTriangle.vertices[ vertexIndex( contact.against ) ] ) );

        }
        else if ( contact.isEdgeVertex() ) {

          additionalVertices.push(
            insertLocalVertexOnEdge(
              bTriangle.vertices[ vertexIndex( contact.against ) ],
              edgeIndex( contact.with ) ) );

        } else {

          // vertex face or vertex edge
          additionalVertices.push( insertLocalVertex( aTriangle.vertices[ vertexIndex( contact.with ) ] ) );

        }
      }

      addEdges( contactPairs, additionalVertices );
    }

    triangulate( a, vertices, aTriangleIndex, false, 0 );
  }

  // now we do a winding-insensitive sort of the triangles.

  std::sort(
    initialChartTriangles_[ 0 ].begin(),
    initialChartTriangles_[ 0 ].end(),
    [] ( const std::pair< Triangle, uint32_t >& left, const std::pair< Triangle, uint32_t >& right ) {

      return less_lowest_vertex_parity( left.first.vertices, right.first.vertices );
    });


  outside_[ 0 ].clear();
  outside_[ 0 ].resize( initialChartTriangles_[ 0 ].size(), 0 );

  auto aWhere = initialChartTriangles_[ 0 ].begin();
  auto aEnd   = initialChartTriangles_[ 0 ].end();

  AABBTree& aBVH = *a.bvh;

  ThreadPool::instance().parallel_for( 0, initialChartTriangles_[ 0 ].size(), [&]( size_t where ) {

    const std::pair< Triangle, uint32_t >& aTriangle = initialChartTriangles_[ 0 ][ where ];

    glm::dvec3 aCentre = vertices.centroid3( aTriangle.first );

    double gwn = aBVH.gwn( a, aCentre, 3.0, aTriangle.second );

    outside_[ 0 ][ where ] = fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );
  });

//  std::transform( std::execution::par, aWhere, aEnd, outside_[ 0 ].begin(), [&]( const std::pair< Triangle, uint32_t >& aTriangle ) {

//    glm::dvec3 aCentre = vertices.centroid3( aTriangle.first );

//    double gwn = aBVH.gwn( a, aCentre, 3.0, aTriangle.second );

//    return fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );
// });

  auto aGWN = outside_[ 0 ].begin();

  if ( initialChartTriangles_[ 0 ].size() > 1 ) {

    while ( ( aWhere + 1 ) < aEnd ) {

      const Triangle& aTriangle = aWhere->first;
      const Triangle& bTriangle = (aWhere + 1)->first;

      if ( less_lowest_vertex_parity( aTriangle.vertices, bTriangle.vertices ) ) {

        bool outside = ( *aGWN ) == 1;

        if ( outside ) {

          outputTriangleStream_.push_back( aTriangle );
        }

        ++aGWN;
        ++aWhere;
        continue;
      }

      bool windingParity =
        lowest_vertex_ordered_parity( aTriangle.vertices ) ==
        lowest_vertex_ordered_parity( bTriangle.vertices );

      // if these triangles are wound the same, keep A (cos we flip the winding on B
      // for subtraction, we can always keep A and throw away B if the winding is the same
      // and throw away both if the winding is different)
      if ( ( *aGWN ) == 1 && windingParity ) {

        outputTriangleStream_.push_back( aTriangle );

      }

      ++aWhere;
      ++aWhere;
      ++aGWN;
      ++aGWN;

      if ( aWhere == aEnd ) {

        break;
      }
    }
  }

  while ( aWhere < aEnd ) {

   bool outside = (*aGWN) == 1;

    if ( outside ) {

      const Triangle& aTriangle = aWhere->first;

      outputTriangleStream_.push_back( aTriangle );
   }

    ++aGWN;
    ++aWhere;
  }

  for (
    uint32_t aTriangleIndex = 0,
              end = static_cast< uint32_t >( a.triangles.size() );
    aTriangleIndex < end;
    ++aTriangleIndex ) {
  
    if ( !boundarySet[ aTriangleIndex ] ) {
    
      outputTriangleStream_.push_back( a.triangles[ aTriangleIndex ] );
    }
  }

  vertexUsed_.clear();
  vertexUsed_.resize( unifiedVertices_.size(), false );

  globalVertexMap_.clear();
  globalVertexMap_.resize( unifiedVertices_.size(), EMPTY_INDEX );

  // Remap and compact all the vertices as we go, outputting the triangle stream.
  for ( Triangle& triangle : outputTriangleStream_ ) {

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle) {
      
      uint32_t originalVertexIndex = triangle.vertices[ vertexInTriangle ];
      uint32_t unifiedVertexIndex  = unifiedVertices_.find( originalVertexIndex );
      uint32_t mappedVertex;

      if( !vertexUsed_[ unifiedVertexIndex ] ) {

        mappedVertex = output.MakeVertex( vertices[ unifiedVertexIndex ] );

        globalVertexMap_[ unifiedVertexIndex ] = mappedVertex;
        vertexUsed_[ unifiedVertexIndex ]      = true;

      } else {

        mappedVertex = globalVertexMap_[ unifiedVertexIndex ];
      }

      triangle.vertices[ vertexInTriangle ] = mappedVertex;
    }
    
    output.MakeTriangle( triangle.vertices[ 0 ], triangle.vertices[ 1 ], triangle.vertices[ 2 ] );
  }
}

void conway::geometry::CSGMesher::process(
  Geometry& a, // these are not const because we lazily generate dipoles.
  Geometry& b,
  const std::vector< TriangleTriangleContactPair >& contacts,
  const PrefixSumMap& aContactMap,
  const PrefixSumMap& bContactMap,
  const std::vector< bool > (&boundarySet)[ 2 ],
  bool aOutside,
  bool bOutside,
  bool flipBWinding,
  Geometry& output ) {
  
  reset();

  assert( a.triangles.size() == a.triangle_edges.size() );
  assert( b.triangles.size() == b.triangle_edges.size() );

  unifiedPlanes_.allocate( static_cast< uint32_t >( a.triangles.size() + b.triangles.size() ) ); 

  // We don't need to know the number of novel vertices yet, because it's the last in the partition.
  MultiMeshVertexIndex< 3 > vertices = multi_mesh_vertex_index(a, b, novelVertices_);

  unifiedVertices_.allocate( vertices( 2, 0 ) );

  uint32_t bTriangleOffset = static_cast< uint32_t >( a.triangles.size() );

  // unify vertex pairs.
  for( const TriangleTriangleContactPair& trianglePair : contacts ) {

    const TriangleContacts& triangleContacts = trianglePair.triangles[ 0 ];

    if ( triangleContacts.face_to_face == FaceFace::COLINEAR ) {
    
      unifiedPlanes_.merge(
        triangleContacts.this_triangle_index,
        bTriangleOffset + triangleContacts.other_triangle_index );
    }

    for ( ContactPair contact : triangleContacts.pairs ) {
      
      if ( contact.isVertexVertex() ) {

        const Triangle& aTriangle = a.triangles[ triangleContacts.this_triangle_index ];
        const Triangle& bTriangle = b.triangles[ triangleContacts.other_triangle_index ];

        uint32_t vertexA = aTriangle.vertices[ vertexIndex( contact.with ) ];
        uint32_t vertexB = vertices( 1, bTriangle.vertices[ vertexIndex( contact.against ) ] );

        unifiedVertices_.merge( vertexA, vertexB );
      }
    }
  }

  // Optimise here, because this will make the mappings above
  // single hop, but the novel vertices don't need to be, because they self reference
  unifiedVertices_.optimize();
  unifiedPlanes_.optimize();

  const std::vector< glm::dvec3 >& aVertices = a.vertices;
  const std::vector< glm::dvec3 >& bVertices = b.vertices;

  for (
    uint32_t aTriangleIndex = 0,
             aTriangleEnd   = static_cast< uint32_t >( a.triangles.size() );
    aTriangleIndex < aTriangleEnd;
    ++aTriangleIndex ) {

    std::span< const uint32_t > trianglePairs = aContactMap.get( aTriangleIndex );

    if ( !boundarySet[ 0 ][ aTriangleIndex ] ||trianglePairs.empty() ) {
      continue;
    }

    const Triangle&      aTriangle      = a.triangles[ aTriangleIndex ];
    const TriangleEdges& aTriangleEdges = a.triangle_edges[ aTriangleIndex ];

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {
      insertLocalVertex( aTriangle.vertices[ vertexInTriangle ] );
    }

    uint32_t aPlaneIndex = unifiedPlanes_.find( aTriangleIndex );

    for ( uint32_t pairIndex : trianglePairs ) {

      const TriangleContacts&  trianglePair   = contacts[ pairIndex].triangles[ 0 ];
      const Triangle&          bTriangle      = b.triangles[ trianglePair.other_triangle_index ];
      const TriangleEdges&     bTriangleEdges = b.triangle_edges[ trianglePair.other_triangle_index ];

      assert( aTriangleIndex == trianglePair.this_triangle_index );

      std::span< const ContactPair > contactPairs = trianglePair.pairs.values();

      FixedStack< uint32_t, 6 > additionalVertices;

      uint32_t bPlaneIndex = unifiedPlanes_.find( trianglePair.other_triangle_index + bTriangleOffset );

      for ( ContactPair contact : contactPairs ) {

        if ( contact.isEdgeEdge() ) {

          uint32_t edgeInTriangleA = edgeIndex( contact.with );
          uint32_t edge0           = aTriangleEdges.edges[ edgeInTriangleA ];
          uint32_t edge1           = bTriangleEdges.edges[ edgeIndex( contact.against ) ];

          uint32_t canidateIndex = vertices( 2, static_cast<uint32_t>( novelVertices_.size() ) );

          const auto [ to, success ] =
            edgeEdgeVertices_.try_emplace(
              std::make_pair( edge0, edge1 ),
              canidateIndex );

          if ( success ) {

            glm::dvec3 e0v0 = aVertices[ a.edges[ edge0 ].vertices[ 0 ] ];
            glm::dvec3 e0v1 = aVertices[ a.edges[ edge0 ].vertices[ 1 ] ];
            glm::dvec3 e1v0 = bVertices[ b.edges[ edge1 ].vertices[ 0 ] ];
            glm::dvec3 e1v1 = bVertices[ b.edges[ edge1 ].vertices[ 1 ] ];

            // while we can't guarantee exact results for this operation,
            // knowing there are some duplicate-vertices-abutting-faces cases
            // means that we can at least make it deterministic by ordering 
            if ( lexicographical_less( e0v1, e0v0 ) ) {

              std::swap( e0v1, e0v0 );
            }
            
            if ( lexicographical_less( e1v1, e1v0 ) ) {

              std::swap( e1v1, e1v0 );
            }

            glm::dvec3 novelVertex =
              line_segment_line_segment_intersection( e0v0, e0v1, e1v0, e1v1 );

            //const auto [ existing, novelSuccess ] = 
            //  duplicateVertexMap_.try_emplace( novelVertex, canidateIndex );

            //if ( novelSuccess ) {

              novelVertices_.push_back( novelVertex ); // todo, add correct vertex generation.
              unifiedVertices_.allocate();
         /*   }
            else {

              to->second = existing->second;
            }*/
          }

          additionalVertices.push(
            insertLocalVertexOnEdge(
              to->second,
              edgeInTriangleA ) );

        } else if ( contact.isFaceEdge() ) {              

          uint32_t edge          = bTriangleEdges.edges[ edgeIndex( contact.against ) ];
          uint32_t canidateIndex = vertices( 2, static_cast<uint32_t>( novelVertices_.size() ) );

          const auto [ to, success ] =
            faceEdgeVertices_[ 0 ].try_emplace(
              std::make_pair( aPlaneIndex, edge ),
              canidateIndex );

          if ( success ) {            

            const glm::dvec3& t0 = aVertices[ aTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = aVertices[ aTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = aVertices[ aTriangle.vertices[ 2 ] ];

            glm::dvec3 ev0 = bVertices[ b.edges[ edge ].vertices[ 0 ] ];
            glm::dvec3 ev1 = bVertices[ b.edges[ edge ].vertices[ 1 ] ];

            // Because these vertices might be on a different edge without 
            // the same identity in some cases, we sort to make this deterministic...
            // although, the triangle itself should be deterministic
            if ( lexicographical_less( ev1, ev0 ) ) {

              std::swap( ev0, ev1 );
            }

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            //const auto [ existing, novelSuccess ] =
            //  duplicateVertexMap_.try_emplace( novelVertex, canidateIndex );

            //if ( novelSuccess ) {

              novelVertices_.push_back( novelVertex ); 
              unifiedVertices_.allocate();
         /*   }
            else {

              to->second = existing->second;
            }*/
          }

          additionalVertices.push( insertLocalVertex( to->second ) );

        } else if ( contact.isEdgeFace() ) {

          uint32_t edgeInTriangle = edgeIndex( contact.with );
          uint32_t edge           = aTriangleEdges.edges[ edgeInTriangle ];

          uint32_t canidateIndex = vertices( 2, static_cast<uint32_t>( novelVertices_.size() ) );

          const auto [ to, success ] =
            faceEdgeVertices_[ 1 ].try_emplace(
              std::make_pair( bPlaneIndex, edge ),
              canidateIndex );

          if ( success ) {

            const glm::dvec3& t0 = bVertices[ bTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = bVertices[ bTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = bVertices[ bTriangle.vertices[ 2 ] ];

            glm::dvec3 ev0 = aVertices[ a.edges[ edge ].vertices[ 0 ] ];
            glm::dvec3 ev1 = aVertices[ a.edges[ edge ].vertices[ 1 ] ];

            // Because these vertices might be on a different edge without 
            // the same identity in some cases, we sort to make this deterministic...
            // although, the triangle itself should be deterministic
            if ( lexicographical_less( ev1, ev0 ) ) {

              std::swap( ev0, ev1 );
            }
            

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            //const auto [ existing, novelSuccess ] =
            //  duplicateVertexMap_.try_emplace( novelVertex, canidateIndex );

            //if ( novelSuccess ) {

              novelVertices_.push_back( novelVertex ); // todo, add correct vertex generation.
              unifiedVertices_.allocate();
         /*   }
            else {

              to->second = existing->second;
            }*/
          }

          additionalVertices.push(
            insertLocalVertexOnEdge( 
              to->second,
              edgeInTriangle ) );

        } else if ( contact.isFaceVertex() ) {

          additionalVertices.push( insertLocalVertex( vertices( 1, bTriangle.vertices[ vertexIndex( contact.against ) ] ) ) );

        }
        else if ( contact.isEdgeVertex() ) {

          additionalVertices.push(
            insertLocalVertexOnEdge(
              vertices( 1, bTriangle.vertices[ vertexIndex( contact.against ) ] ),
              edgeIndex( contact.with ) ) );

        } else {

          // vertex face or vertex edge
          additionalVertices.push( insertLocalVertex( aTriangle.vertices[ vertexIndex( contact.with ) ] ) );

        }
      }

      addEdges( contactPairs, additionalVertices );
    }

    triangulate( a, vertices, aTriangleIndex, false, 0 );
  }

  for (
    uint32_t bTriangleIndex = 0, bTriangleEnd = static_cast< uint32_t >( b.triangles.size() );
    bTriangleIndex < bTriangleEnd;
    ++bTriangleIndex ) {

    std::span< const uint32_t > trianglePairs = bContactMap.get( bTriangleIndex );

    if ( !boundarySet[ 1 ][ bTriangleIndex ] || trianglePairs.empty() ) {
      continue;
    }

    const Triangle&      bTriangle      = b.triangles[ bTriangleIndex ];
    const TriangleEdges& bTriangleEdges = b.triangle_edges[ bTriangleIndex ];

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {

      insertLocalVertex( vertices( 1, bTriangle.vertices[ vertexInTriangle ] ) );
    }

    uint32_t bPlaneIndex = unifiedPlanes_.find( bTriangleIndex + bTriangleOffset );

    for ( uint32_t pairIndex : trianglePairs ) {

      const TriangleContacts&  trianglePair   = contacts[ pairIndex ].triangles[ 1 ];
      const Triangle&          aTriangle      = a.triangles[ trianglePair.other_triangle_index ];
      const TriangleEdges&     aTriangleEdges = a.triangle_edges[ trianglePair.other_triangle_index ];

      assert( bTriangleIndex == trianglePair.this_triangle_index );

      std::span< const ContactPair > contactPairs = trianglePair.pairs.values();

      FixedStack< uint32_t, 6 > additionalVertices;

      uint32_t aPlaneIndex = unifiedPlanes_.find( trianglePair.other_triangle_index );

      for ( ContactPair contact : contactPairs ) {

        if ( contact.isEdgeEdge() ) {

          uint32_t edgeInTriangleB = edgeIndex( contact.with );
          uint32_t edge0           = bTriangleEdges.edges[ edgeInTriangleB ];
          uint32_t edge1           = aTriangleEdges.edges[ edgeIndex( contact.against ) ];

          auto to = edgeEdgeVertices_.find( std::make_pair( edge1, edge0 ) );

          assert( to != edgeEdgeVertices_.end() );

          additionalVertices.push(
            insertLocalVertexOnEdge( to->second, edgeInTriangleB ) );

        }
        else if ( contact.isFaceEdge() ) {

          uint32_t edge = aTriangleEdges.edges[ edgeIndex( contact.against ) ];

          const auto to =
            faceEdgeVertices_[1].find(
              std::make_pair( bPlaneIndex, edge ) );

          assert( to != faceEdgeVertices_[ 1 ].end() );

          additionalVertices.push(
            insertLocalVertex( to->second ) );

        }
        else if ( contact.isEdgeFace() ) {

          uint32_t edgeInTriangle = edgeIndex( contact.with );
          uint32_t edge           = bTriangleEdges.edges[ edgeInTriangle ];

          const auto to =
            faceEdgeVertices_[ 0 ].find(
              std::make_pair( aPlaneIndex, edge ) );

          assert( to != faceEdgeVertices_[ 0 ].end() );

          additionalVertices.push(
            insertLocalVertexOnEdge(
              to->second,
              edgeInTriangle ) );

        } else if ( contact.isFaceVertex() ) {

          additionalVertices.push( insertLocalVertex( aTriangle.vertices[ vertexIndex( contact.against ) ] ) );

        }
        else if ( contact.isEdgeVertex() ) {

          additionalVertices.push(
            insertLocalVertexOnEdge(
              aTriangle.vertices[ vertexIndex( contact.against ) ],
              edgeIndex( contact.with ) ) );

        } else {

          // vertex face or vertex edge
          additionalVertices.push( insertLocalVertex( vertices( 1, bTriangle.vertices[ vertexIndex( contact.with ) ] ) ) );
        }
      }

      addEdges( contactPairs, additionalVertices );
    }
    
    triangulate( b, vertices, bTriangleIndex, flipBWinding, 1 );
  }

  // now we do a winding-insensitive sort of the triangles.
  for ( uint32_t triangleSet = 0; triangleSet < 2; ++triangleSet ) {

    std::sort(
      initialChartTriangles_[ triangleSet ].begin(),
      initialChartTriangles_[ triangleSet ].end(),
      [] ( const std::pair< Triangle, uint32_t >& left, const std::pair< Triangle, uint32_t >& right ) {

        return less_lowest_vertex_parity( left.first.vertices, right.first.vertices );
      });
  }

  outside_[ 0 ].clear();
  outside_[ 1 ].clear();

  outside_[ 0 ].resize( initialChartTriangles_[ 0 ].size(), 0 );
  outside_[ 1 ].resize( initialChartTriangles_[ 1 ].size(), 0 );

  auto aWhere = initialChartTriangles_[ 0 ].begin();
  auto aEnd   = initialChartTriangles_[ 0 ].end();
  auto bWhere = initialChartTriangles_[ 1 ].begin();
  auto bEnd   = initialChartTriangles_[ 1 ].end();

  AABBTree& aBVH = *a.bvh;
  AABBTree& bBVH = *b.bvh;

// #if defined(__EMSCRIPTEN__)
//   std::transform( aWhere, aEnd, outside_[0].begin(), [&]( const Triangle& aTriangle ) {
// #else
//  std::transform( std::execution::par, aWhere, aEnd, outside_[0].begin(), [&]( const std::pair< Triangle, uint32_t >& aTriangle ) {
// //#endif
//     glm::dvec3 aCentre = vertices.centroid3( aTriangle.first );

//     double gwn = bBVH.gwn( b, aCentre, 3.0 );

//     return fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );
//   });

  ThreadPool::instance().parallel_for( 0, initialChartTriangles_[ 0 ].size(), [&]( size_t where ) {

    const std::pair< Triangle, uint32_t >& aTriangle = initialChartTriangles_[ 0 ][ where ];
  
    glm::dvec3 aCentre = vertices.centroid3( aTriangle.first );

    double gwn = bBVH.gwn( b, aCentre, 3.0 );

    outside_[ 0 ][ where ] = fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );

  });


// #if defined(__EMSCRIPTEN__)
//   std::transform( bWhere, bEnd, outside_[1].begin(), [&]( const Triangle& bTriangle ) {
// #else
//   std::transform( std::execution::par, bWhere, bEnd, outside_[1].begin(), [&]( const std::pair< Triangle, uint32_t >& bTriangle ) {
// //#endif
//     glm::dvec3 bCentre = vertices.centroid3( bTriangle.first );

//     double gwn = aBVH.gwn( a, bCentre, 3.0 );

//     return fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );
//   });

ThreadPool::instance().parallel_for( 0, initialChartTriangles_[ 1 ].size(), [&]( size_t where ) {

    const std::pair< Triangle, uint32_t >& bTriangle = initialChartTriangles_[ 1 ][ where ];

    glm::dvec3 bCentre = vertices.centroid3( bTriangle.first );

    double gwn = aBVH.gwn( a, bCentre, 3.0 );

    outside_[ 1 ][ where ] = fabs( gwn ) < GWN_TOLERANCE ? uint8_t( 1 ) : uint8_t( 0 );

  });

  auto aGWN = outside_[ 0 ].begin();
  auto bGWN = outside_[ 1 ].begin();

  while ( aWhere < aEnd && bWhere < bEnd ) {

    const Triangle& aTriangle = aWhere->first;
    const Triangle& bTriangle = bWhere->first;

    if ( less_lowest_vertex_parity( aTriangle.vertices, bTriangle.vertices ) ) {

      bool outside = ( *aGWN ) == 1;

      if ( outside == aOutside && OUTPUT_A_BOUNDARY ) {

        outputTriangleStream_.push_back( aTriangle );
      }

      ++aGWN;
      ++aWhere;
      continue;
    }

    if ( less_lowest_vertex_parity( bTriangle.vertices, aTriangle.vertices ) ) {

      bool outside = ( *bGWN ) == 1;

      if ( outside == bOutside && OUTPUT_B_BOUNDARY ) {

        outputTriangleStream_.push_back( bTriangle );
      }

      ++bWhere;
      ++bGWN;
      continue;
    }

    bool windingParity =
      lowest_vertex_ordered_parity( aTriangle.vertices ) ==
      lowest_vertex_ordered_parity( bTriangle.vertices );

    // if these triangles are wound the same, keep A (cos we flip the winding on B
    // for subtraction, we can always keep A and throw away B if the winding is the same
    // and throw away both if the winding is different)
    if ( windingParity && OUTPUT_A_SHARED ) {

      outputTriangleStream_.push_back( aTriangle );
    }

    ++aWhere;
    ++bWhere;

    ++aGWN;
    ++bGWN;
  }
    
  while ( aWhere < aEnd && OUTPUT_A_BOUNDARY ) {

    bool outside = (*aGWN) == 1;

    if ( outside == aOutside ) {

      const Triangle& aTriangle = aWhere->first;

      outputTriangleStream_.push_back( aTriangle );
    }

    ++aGWN;
    ++aWhere;
  }

  while ( bWhere < bEnd && OUTPUT_B_BOUNDARY ) {

    bool outside = (*bGWN) == 1;

    if ( outside == bOutside ) {

      const Triangle& bTriangle = bWhere->first;

      outputTriangleStream_.push_back( bTriangle );
    }

    ++bGWN;
    ++bWhere;
  }

  if constexpr ( OUTPUT_A_ORIGINAL ) {
    walkAndInsertNonBoundary( aOutside, boundarySet[ 0 ], a, b, false, 0 );
  }

  if constexpr ( OUTPUT_B_ORIGINAL ) {
    walkAndInsertNonBoundary( bOutside, boundarySet[ 1 ], b, a, flipBWinding, vertices( 1, 0 ) );
  }

  vertexUsed_.clear();
  vertexUsed_.resize( unifiedVertices_.size(), false );

  globalVertexMap_.clear();
  globalVertexMap_.resize( unifiedVertices_.size(), EMPTY_INDEX );

  assert( unifiedVertices_.size() == ( novelVertices_.size() + a.vertices.size() + b.vertices.size() ) );

  // Remap and compact all the vertices as we go, outputting the triangle stream.
  for ( Triangle& triangle : outputTriangleStream_ ) {

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle) {
      
      uint32_t originalVertexIndex = triangle.vertices[ vertexInTriangle ];
      uint32_t unifiedVertexIndex  = unifiedVertices_.find( originalVertexIndex );

      uint32_t mappedVertex;

      if( !vertexUsed_[ unifiedVertexIndex ] ) {

        mappedVertex = output.MakeVertex( vertices[ unifiedVertexIndex ] );

        globalVertexMap_[ unifiedVertexIndex ] = mappedVertex;
        vertexUsed_[ unifiedVertexIndex ]      = true;

      } else {

        mappedVertex = globalVertexMap_[ unifiedVertexIndex ];
      }

      triangle.vertices[ vertexInTriangle ] = mappedVertex;
    }
    
    output.MakeTriangle( triangle.vertices[ 0 ], triangle.vertices[ 1 ], triangle.vertices[ 2 ] );
  }
}


void conway::geometry::CSGMesher::reset() {

  unifiedVertices_.reset();
  novelVertices_.clear();
  edgeEdgeVertices_.clear();
  faceEdgeVertices_[ 0 ].clear();
  faceEdgeVertices_[ 1 ].clear();
  vertexUsed_.clear();
}

void conway::geometry::CSGMesher::walkAndInsertNonBoundary(
  bool outside,
  const std::vector< bool >& boundarySet,
  const Geometry& mesh,
  const Geometry& otherMesh,
  bool flippedWinding,
  uint32_t vertexOffset ) {

  walked_.clear();

  const AABBTree& bvh = *otherMesh.bvh;
  
  assert( boundarySet.size() == mesh.triangles.size() );

  walked_.reserve( boundarySet.size() );
  walked_.insert( walked_.begin(), boundarySet.begin(), boundarySet.end() );

  uint32_t walkedCursor = 0;
  uint32_t walkedEnd    = static_cast< uint32_t >( boundarySet.size() );

  for ( ; walkedCursor < walkedEnd; ++walkedCursor ) {

    if ( walked_[ walkedCursor ] ) {
      continue;
    }

    const Triangle& initialTriangle = mesh.triangles[ walkedCursor ];
    glm::dvec3      centre          = centroid( mesh.vertices, initialTriangle );
    bool            triangleOutside = fabs( bvh.gwn( otherMesh, centre ) ) < GWN_TOLERANCE;

    triangleStack_.push_back(  walkedCursor );

    if ( triangleOutside == outside ) {

      while ( !triangleStack_.empty() ) {

        uint32_t nextTriangleIndex = triangleStack_.back();

        triangleStack_.pop_back();

        if ( walked_[ nextTriangleIndex ] ) {

          continue;
        }

        walked_[ nextTriangleIndex ] = true;

        const Triangle&      triangle      = mesh.triangles[ nextTriangleIndex ];
        const TriangleEdges& triangleEdges = mesh.triangle_edges[ nextTriangleIndex ];

        // Copy explicitly instead of slicing to avoid warnings - CS
        Triangle outputTriangle = { 
            triangle.vertices[ 0 ],
            triangle.vertices[ 1 ],
            triangle.vertices[ 2 ] 
        };

        if ( flippedWinding ) {

          std::swap( outputTriangle.vertices[ 1 ], outputTriangle.vertices[ 2 ] );
        } 

        outputTriangle.vertices[ 0 ] += vertexOffset;
        outputTriangle.vertices[ 1 ] += vertexOffset;
        outputTriangle.vertices[ 2 ] += vertexOffset;

        outputTriangleStream_.push_back( outputTriangle );

        for ( uint32_t triangleInEdge = 0; triangleInEdge < 3; ++triangleInEdge ) {
          
          uint32_t edgeIndex = triangleEdges.edges[ triangleInEdge ];

          if ( edgeIndex == EMPTY_INDEX ) {

            continue;
          } 

          uint32_t opposingTriangle =
            mesh.edges[ edgeIndex ].otherTriangle( nextTriangleIndex );

          if ( opposingTriangle != EMPTY_INDEX && !walked_[ opposingTriangle ] ) {

            triangleStack_.push_back( opposingTriangle );
          }
        }
      }
    } else {

      while ( !triangleStack_.empty() ) {

        uint32_t nextTriangleIndex = triangleStack_.back();

        triangleStack_.pop_back();

        if ( walked_[ nextTriangleIndex ] ) {

          continue;
        }

        walked_[ nextTriangleIndex ] = true;

        const TriangleEdges& triangleEdges = mesh.triangle_edges[ nextTriangleIndex ];

        for ( uint32_t triangleInEdge = 0; triangleInEdge < 3; ++triangleInEdge ) {

          uint32_t edgeIndex = triangleEdges.edges[ triangleInEdge ];

          if ( edgeIndex == EMPTY_INDEX ) {

            continue;
          } 

          uint32_t opposingTriangle = mesh.edges[ edgeIndex ].otherTriangle( nextTriangleIndex );

          if ( opposingTriangle != EMPTY_INDEX && !walked_[ opposingTriangle ] ) {

            triangleStack_.push_back( opposingTriangle );
          }
        }
      }
    }
  }
}

void conway::geometry::CSGMesher::addEdges(
  std::span< const ContactPair > contactPairs,
  const FixedStack< uint32_t, 6 >& additionalVertices ) {

  assert( contactPairs.size() == additionalVertices.size() );

  // Intersection is a point.
  if ( contactPairs.size() <= 1 ) {

    return;
  }

  // Intersection is a segment.
  if ( contactPairs.size() == 2 ) {

    ContactPair contact0 = contactPairs[ 0 ];
    ContactPair contact1 = contactPairs[ 1 ];

    uint32_t additionalVertice0 = additionalVertices[ 0 ];
    uint32_t additionalVertice1 = additionalVertices[ 1 ];

    if (
      isSameEdge( contact0.with, contact1.with ) || 
      unifiedVertices_.find( additionalVertice0 ) == unifiedVertices_.find( additionalVertice1 ) ) {

      return;
    }

    edges_.emplace_back( additionalVertice0, additionalVertice1 );
    return;
  }

  size_t innerEnd = contactPairs.size();

  // Intersection is a polygon
  for ( size_t outer = 0, end = contactPairs.size() - 1; outer < end; ++outer ) {

    ContactPair outerContact = contactPairs[ outer ];

    for ( size_t inner = outer + 1; inner < innerEnd; ++inner  ) {

      ContactPair innerContact = contactPairs[ inner ];

      // If this pair of contacts with the same triangle b
      // both are a point on the same edge (including end vertices) of
      // the triangle we are intersecting, but
      // not on the same edge on *this* triangle, it represents an edge constraint.
      if (
        !isSameEdge( outerContact.with, innerContact.with ) &&
        isSameEdge( outerContact.against, innerContact.against ) ) {

        edges_.emplace_back( additionalVertices[ outer ], additionalVertices[ inner ] );
      }
    }
  }
}

template < size_t N >
void conway::geometry::CSGMesher::triangulate(
  const Geometry& mesh,
  const MultiMeshVertexIndex< N >& vertices,
  uint32_t triangleInMeshIndex,
  bool flippedWinding,
  uint32_t outputStreamIndex ) {

  assert( mesh.triangles.size() > triangleInMeshIndex );

  const Triangle& triangle = mesh.triangles[ triangleInMeshIndex ];

#if !defined( __EMSCRIPTEN__ )
  for ( CDT::Edge& edge : edges_ ) {

    glm::dvec3 v1 = vertices[ edge.v1() ];
    glm::dvec3 v2 = vertices[ edge.v2() ];

    contraintEdge_.emplace_back( v1, v2 );
  }
#endif

  localVertexMap_.clear();

  // Find any duplicates (including merged vertices) and make them unique.
  for ( uint32_t localVertex : localVertices_ ) {

    const glm::dvec3& vertex = vertices[ localVertex ];

    if ( isnan( vertex.x ) || isnan( vertex.y ) || isnan( vertex.z ) ) {
      continue;
    }

    localVertexMap_.try_emplace( unifiedVertices_.find( localVertex ), static_cast< uint32_t >( localVertexMap_.size() ) );
  }

  // Cut off the size of what's been removed.
  localVertices_.resize( localVertexMap_.size() );

  // set the vertices back to their new locations.
  for ( auto [ localVertex, index ] : localVertexMap_ ) {

    localVertices_[ index ] = localVertex;
  }

  edges_.erase( 
    std::remove_if(
      edges_.begin(),
      edges_.end(),
      [&]( const CDT::Edge& edge ){

        uint32_t ev1 = unifiedVertices_.find( edge.v1() );
        uint32_t ev2 = unifiedVertices_.find( edge.v2() );

        return
          ev1 == ev2 ||
          !localVertexMap_.contains( ev1 ) ||
          !localVertexMap_.contains( ev2 );
      } ),
      edges_.end() );

  // Remap edges with the new unique vertices.
  for ( CDT::Edge& edge: edges_ ) {

    uint32_t edgeV1 = edge.v1();
    uint32_t edgeV2 = edge.v2();

    edge =
      CDT::Edge(
        localVertexMap_[ unifiedVertices_.find(edgeV1)],
        localVertexMap_[ unifiedVertices_.find(edgeV2)]);
  }

  // Now use the best 2D projection to extract the vertices.
  glm::dvec3 triangleVertices[ 3 ];

  extract_vertices( mesh, triangle, triangleVertices );

  AxisPair axesToExtract = best_truncated_projection( triangleVertices );

  int32_t winding = orient2D( triangleVertices, axesToExtract );

  // 3 vertices are already fully constrained and triangulated.
  if ( localVertices_.size() == 3 ) {

    Triangle localTriangle { { 
        unifiedVertices_.find( localVertices_[ 0 ] ),
        unifiedVertices_.find( localVertices_[ 1 ] ),
        unifiedVertices_.find( localVertices_[ 2 ] )
    }};

    edges_.clear();

    localVertexMap_.clear();
    localVertices_.clear();
    local2DVertices_.clear();
    localVertexEdgeFlags_.clear();

    for ( auto& where : std::span( onEdgeVertices_ ) ) {

      where.clear();
    }

    int32_t localWinding = vertices.orient2D( localTriangle.vertices, axesToExtract );

    if ( localTriangle.vertices[ 0 ] == localTriangle.vertices[ 1 ] ||
          localTriangle.vertices[ 1 ] == localTriangle.vertices[ 2 ] ||
          localTriangle.vertices[ 2 ] == localTriangle.vertices[ 0 ] ||
          localWinding == 0 ) {
      return;
    }

    bool flipWinding = ( ( flippedWinding ? -1 : 1 ) * localWinding * winding ) < 0;

    reorder_to_lowest_vertex( localTriangle.vertices );

    if ( flipWinding ) {
      
      std::swap( localTriangle.vertices[ 1 ], localTriangle.vertices[ 2 ] );

    }

    initialChartTriangles_[ outputStreamIndex ].emplace_back( localTriangle, triangleInMeshIndex );
    return;
  }

  local2DVertices_.clear();
  local2DVertices_.reserve( localVertices_.size() );
  localVertexEdgeFlags_.clear();
  localVertexEdgeFlags_.resize( localVertices_.size(), 0 );

  glm::length_t firstAxis  = first_axis( axesToExtract );
  glm::length_t secondAxis = second_axis( axesToExtract );

  assert( axesToExtract != AxisPair::NONE );

  for ( uint32_t partitionedIndice : localVertices_ )  {

    const glm::dvec3& inputVertex = vertices[ partitionedIndice ];

    local2DVertices_.push_back( CDT::V2d< double >::make( inputVertex[ firstAxis ], inputVertex[ secondAxis ] ) );
  }
  
  for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle ) {

    std::vector< uint32_t >& edgeVertices = onEdgeVertices_[ edgeInTriangle ];

    uint32_t v0Index = edgeInTriangle;
    uint32_t v1Index = ( edgeInTriangle + 1 ) % 3;

    uint8_t edgeFlag = 1 << static_cast< uint8_t >( edgeInTriangle );

    localVertexEdgeFlags_[ v0Index ] |= edgeFlag;
    localVertexEdgeFlags_[ v1Index ] |= edgeFlag;

    for ( uint32_t onEdgeVertex : edgeVertices ) {

      onEdgeVertex = localVertexMap_[ unifiedVertices_.find( onEdgeVertex ) ];

      localVertexEdgeFlags_[ onEdgeVertex ] |= edgeFlag;
    }
  }

  edges_.erase(
    std::remove_if(
      edges_.begin(),
      edges_.end(),
      [](const CDT::Edge& edge) { return edge.v1() == edge.v2(); } ),
    edges_.end() );

  // Despite de-duplicating vertices in 3D, there is a numerical case
  // where the *dropped* axis can be *very* slightly different in (effectively being the same value),
  // but the other two values are the same.
  // This is a degenerate case in general, and these should be merged, and the triangulation requires it,
  // so we'll allow it to detect those cases in 2D and remap.
  CDT::DuplicatesInfo duplicates = CDT::RemoveDuplicatesAndRemapEdges( local2DVertices_, edges_ );

  std::vector< uint32_t > localVertRemapping; 

  localVertRemapping.resize( local2DVertices_.size(), 0 );

  for ( size_t duplicate : duplicates.duplicates ) {

    unifiedVertices_.merge( localVertices_[ duplicate ], localVertices_[ duplicates.mapping[ duplicate ] ] );

    localVertexEdgeFlags_[ duplicates.mapping[ duplicate ] ] |= localVertexEdgeFlags_[ duplicate ];
  }

  for ( size_t where = 0, end = duplicates.mapping.size(); where < end; ++where ) {

    if ( localVertRemapping[ duplicates.mapping[ where ] ] == 0 ) {
      localVertRemapping[ duplicates.mapping[ where ] ] = static_cast< uint32_t >( where );
    }
  }

  CDT::Triangulation< double > triangulation( CDT::VertexInsertionOrder::AsProvided, CDT::IntersectingConstraintEdges::NotAllowed, 0 );
  
  size_t retryOperations = 0;

  while ( true ) {

    try {

      triangulation.insertVertices( local2DVertices_ );
      triangulation.insertEdges( edges_ );
      triangulation.eraseSuperTriangle();

      break;

    } catch ( const CDT::DuplicateVertexError& duplicateVertex ) {

      const glm::dvec3& v1 = vertices[ localVertices_[ localVertRemapping[ duplicateVertex.v1() ] ] ];
      const glm::dvec3& v2 = vertices[ localVertices_[ localVertRemapping[ duplicateVertex.v2() ] ] ];

      triangulation.triangles.clear();

      printf(
        "Duplicate vertex v1: %.20f %.20f %.20f v2: %.20f %.20f %.20f (indices - local + global) %u %u %u %u\n",
        v1.x,
        v1.y,
        v1.z,
        v2.x,
        v2.y,
        v2.z,
        duplicateVertex.v1(),
        duplicateVertex.v2(),
        localVertices_[ duplicateVertex.v1() ],
        localVertices_[ duplicateVertex.v2() ] );

      if constexpr ( DUMP_SVGS_ON_TRIANGLE_ERROR ) {

        printf( "\n\n\n%s\n\n\n", dumpEdgeAndVertsToSVG().c_str() );
      }

      break;

    } catch ( const CDT::IntersectingConstraintsError& constraintError ) {

      if ( ++retryOperations > MAXIMUM_TRIANGULATION_RETRIES ) {
        printf( "Intersecting constraint error couldn't be resolved:\n\t%s\n", constraintError.what() );
        break;
      }

      const CDT::V2d< double >& e00 = local2DVertices_[ constraintError.e1().v1() ];
      const CDT::V2d< double >& e01 = local2DVertices_[ constraintError.e1().v2() ];
      const CDT::V2d< double >& e10 = local2DVertices_[ constraintError.e2().v1() ];
      const CDT::V2d< double >& e11 = local2DVertices_[ constraintError.e2().v2() ];

      double e0v0 = predicates::adaptive::orient2d( &e00.x, &e01.x, &e10.x );
      double e0v1 = predicates::adaptive::orient2d( &e00.x, &e01.x, &e11.x );
      double e1v0 = predicates::adaptive::orient2d( &e10.x, &e01.x, &e00.x );
      double e1v1 = predicates::adaptive::orient2d( &e10.x, &e01.x, &e01.x );

      auto edge1 = std::find( edges_.begin(), edges_.end(), constraintError.e1() );
      auto edge2 = std::find( edges_.begin(), edges_.end(), constraintError.e2() );

      // Try and resolve the constraint by welding
      if ( fabs( e0v0 ) < WELD_EDGE_TOLERANCE && edge1 != edges_.end() ) {

        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e1() ) );

        edges_.emplace_back( constraintError.e1().v1(), constraintError.e2().v1() );
        edges_.emplace_back( constraintError.e2().v1(), constraintError.e1().v2() );

      } else if ( fabs( e0v1 ) < WELD_EDGE_TOLERANCE && edge1 != edges_.end()  ) {

        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e1() ) );

        edges_.emplace_back( constraintError.e1().v1(), constraintError.e2().v2() );
        edges_.emplace_back( constraintError.e2().v2(), constraintError.e1().v2() );

      } else if ( fabs( e1v0 ) < WELD_EDGE_TOLERANCE && edge2 != edges_.end() ) {

        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e2() ) );

        edges_.emplace_back( constraintError.e2().v1(), constraintError.e1().v1() );
        edges_.emplace_back( constraintError.e1().v1(), constraintError.e2().v2() );

      } else if ( fabs( e1v1 ) < WELD_EDGE_TOLERANCE && edge2 != edges_.end() ) {

        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e2() ) );
        edges_.emplace_back( constraintError.e2().v1(), constraintError.e1().v2() );
        edges_.emplace_back( constraintError.e1().v2(), constraintError.e2().v2() );

      } else if ( ( e0v0 + e0v1 < e1v0 + e1v1 || edges_.end() == edge2 ) && edge1 != edges_.end()) {
        
        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e1() ) );

      } else if ( edge2 != edges_.end() ) {

        edges_.erase( std::find( edges_.begin(), edges_.end(), constraintError.e2() ) );

      } else {

        printf( "Intersecting constraint error couldn't be resolved due to missing edge in set:\n\t%s\n", constraintError.what() );
        break;
      }

      triangulation.~Triangulation();

      new (&triangulation) CDT::Triangulation< double >( CDT::VertexInsertionOrder::AsProvided, CDT::IntersectingConstraintEdges::NotAllowed, 0 );

      if constexpr ( DUMP_SVGS_ON_TRIANGLE_ERROR ) {
      
        printf( "\n\n\n%s\n\n\n", dumpEdgeAndVertsToSVG().c_str() );
      }
    }
  }

  int32_t foundWinding = 0;

  for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {

    // This is here because after an intersecting constraint error, there may be a "virtual" vertex
    // created.
    if ( cdtTriangle.vertices[ 0 ] >= localVertRemapping.size() ||
         cdtTriangle.vertices[ 1 ] >= localVertRemapping.size() ||
         cdtTriangle.vertices[ 2 ] >= localVertRemapping.size() ) {

      continue;
    }

    uint32_t localTriangle[ 3 ] = {

      localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 0 ] ] ],
      localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 1 ] ] ],
      localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 2 ] ] ]

    };

    foundWinding = vertices.orient2D( localTriangle, axesToExtract );

    if ( foundWinding != 0 ) {

      break;
    }
  }

  assert( foundWinding != 0 );

  std::vector< std::pair< Triangle, uint32_t > >& outputStream = initialChartTriangles_[ outputStreamIndex ];

  bool flipWinding = ( ( flippedWinding ? -1 : 1 ) * foundWinding * winding ) < 0;

  if ( flipWinding ) {

    for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {

      if ( cdtTriangle.vertices[ 0 ] >= localVertRemapping.size() ||
           cdtTriangle.vertices[ 1 ] >= localVertRemapping.size() ||
           cdtTriangle.vertices[ 2 ] >= localVertRemapping.size() ) {

        continue;
      }

      uint8_t edgeMask =
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 0 ] ] ] &
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 1 ] ] ] &
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 2 ] ] ];

      // Hey CDT, are you disrespecting me with triangles on the same edge,
      // despite that case being full constrained? 
      // My guess is that the triangle exists *prior* to edges being inserted,
      // and although CDT has tolerance for points being "on" edges, the near
      // zero area triangle exists within the boundary and can't be edge-flipped
      // in a compliant way, being technically an overlap.
      // If we were using exact maths an a matching exact CDT, then we
      // wouldn't have to worry, but for this purpose we're in-exact,
      // so we must deal with this particular case, despite constraining the edges.
      if ( edgeMask != 0 ) {
        continue;
      }

      Triangle localTriangle { { 
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 0 ] ] ] ),
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 1 ] ] ] ),
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 2 ] ] ] )
      }};
      
      if ( localTriangle.vertices[ 0 ] == localTriangle.vertices[ 1 ] ||
            localTriangle.vertices[ 1 ] == localTriangle.vertices[ 2 ] ||
            localTriangle.vertices[ 2 ] == localTriangle.vertices[ 0 ] ||
            vertices.orient2D( localTriangle.vertices, axesToExtract ) == 0 ) {
        continue;
      }

      reorder_to_lowest_vertex( localTriangle.vertices );

      std::swap( localTriangle.vertices[ 1 ], localTriangle.vertices[ 2 ] );

      outputStream.emplace_back( localTriangle, triangleInMeshIndex );
    }

  } else {

    for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {

      if ( cdtTriangle.vertices[ 0 ] >= localVertRemapping.size() ||
           cdtTriangle.vertices[ 1 ] >= localVertRemapping.size() ||
           cdtTriangle.vertices[ 2 ] >= localVertRemapping.size() ) {

        continue;
      }

      uint8_t edgeMask =
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 0 ] ] ] &
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 1 ] ] ] &
        localVertexEdgeFlags_[ localVertRemapping[ cdtTriangle.vertices[ 2 ] ] ];

      if ( edgeMask != 0 ) {
        continue;
      }

      Triangle localTriangle { { 
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 0 ] ] ] ),
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 1 ] ] ] ),
        unifiedVertices_.find( localVertices_[ localVertRemapping[ cdtTriangle.vertices[ 2 ] ] ] )
      }};

      if ( localTriangle.vertices[ 0 ] == localTriangle.vertices[ 1 ] ||
            localTriangle.vertices[ 1 ] == localTriangle.vertices[ 2 ] ||
            localTriangle.vertices[ 2 ] == localTriangle.vertices[ 0 ] ||
            vertices.orient2D( localTriangle.vertices, axesToExtract ) == 0 ) {
        continue;
      }

      reorder_to_lowest_vertex( localTriangle.vertices );

      outputStream.emplace_back( localTriangle, triangleInMeshIndex );
    }
  }

  edges_.clear();

  localVertexMap_.clear();
  localVertices_.clear();
  local2DVertices_.clear();
  localVertexEdgeFlags_.clear();

  for ( auto& where : std::span( onEdgeVertices_ ) ) {

    where.clear();
  }
}

#if !defined( __EMSCRIPTEN__ )
std::string conway::geometry::CSGMesher::dumpConstraints( const std::string& preamble ) const {

  std::stringstream obj;

  obj << preamble;

  size_t indice = 0;

  for ( const auto& [ a, b ] : contraintEdge_ ) {

    obj << "v " << a.x << " " << a.y << " " << a.z << "\n";
    obj << "v " << b.x << " " << b.y << " " << b.z << "\n";

    size_t indice1 = ++indice;
    size_t indice2 = ++indice;

    obj << "l " << indice1 << " " << indice2 << "\n";
  }

  return obj.str();
}
#endif

uint32_t conway::geometry::CSGMesher::insertLocalVertex( uint32_t inputVertex ) {

  uint32_t result = /*unifiedVertices_.find*/( inputVertex );

  localVertices_.push_back( result );

  return result;
}

uint32_t conway::geometry::CSGMesher::insertLocalVertexOnEdge( uint32_t inputVertex, uint32_t edgeInTriangle ) {
  
  uint32_t result = /*unifiedVertices_.find*/( inputVertex );

  localVertices_.push_back( result );

  onEdgeVertices_[ edgeInTriangle ].push_back( result );

  return result;
}

std::string conway::geometry::CSGMesher::dumpNovelVertices( const std::string& preamble ) const {

  std::string result = preamble;

  result.append(
    std::format(
      "ply\n"
      "format ascii 1.0\n"
      "element vertex {}\n"
      "property float x\n"
      "property float y\n"
      "property float z\n"
      "property uchar red\n"
      "property uchar green\n"
      "property uchar blue\n"
      "property uchar alpha\n"
      "end_header\n", novelVertices_.size() ) );

  for ( const glm::dvec3& novelVertex : novelVertices_ ) {

    result.append( 
      std::format(
        "{} {} {} 255 0 255 255\n",
        novelVertex.x,
        novelVertex.y,
        novelVertex.z) );
  }
 
  return result;
}

constexpr double SVG_RIGHT_PAD_RATIO = 0.1;

std::string conway::geometry::CSGMesher::dumpEdgeAndVertsToSVG() const {

  glm::dvec2 min( DBL_MAX );
  glm::dvec2 max( -DBL_MAX );

  for (  const CDT::V2d< double >& vert: local2DVertices_ ) {

    glm::dvec2 point( vert.x, vert.y );

    min = glm::min( min, point );
    max = glm::max( max, point );
  }

  SVGContext svg(
    glm::dvec2( 1536, 1536 ),
    glm::dvec2( 50, 50 ),
    min,
    max + ( SVG_RIGHT_PAD_RATIO * ( max - min ) ) );

  svg.header();

  for ( const CDT::Edge& edge: edges_ ) {

    const CDT::V2d< double >& v0 = local2DVertices_[ edge.v1() ];
    const CDT::V2d< double >& v1 = local2DVertices_[ edge.v2() ];

    svg.line( glm::dvec2( v0.x, v0.y ), glm::dvec2( v1.x, v1.y ) );
  }
  
  size_t index = 0; 

  for ( const CDT::V2d< double >& vert: local2DVertices_ ) {

    std::stringstream labelText;

    labelText << index;

    svg.point( labelText.str().c_str(), glm::dvec2( vert.x, vert.y ) );
  
    ++index;
  }

  svg.trailer();

  return svg.str();
}