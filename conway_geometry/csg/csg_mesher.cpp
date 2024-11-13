#include "csg_mesher.h"
#include "structures/winged_edge.h"
#include "multi_mesh_vertex_index.h"

void conway::geometry::CSGMesher::process(
  WingedEdgeDV3& a, // these are not const because we lazily generate dipoles.
  WingedEdgeDV3& b,
  const std::vector< TriangleContacts >& aContacts,
  const PrefixSumMap& aContactMap,
  const std::vector< TriangleContacts >& bContacts,
  const PrefixSumMap& bContactMap,
  const std::vector< bool > (&boundarySet)[ 2 ],
  bool aOutside,
  bool bOutside,
  bool flipBWinding,
  WingedEdgeMesh< glm::dvec3 >& output ) {
  
  reset();

  unifiedVertices_.allocate( a.vertices.size() + b.vertices.size() );

  const std::vector< glm::dvec3 >* abNovel[ 3 ] = { &a.vertices, &b.vertices, &novelVertices_ };

  // We don't need to know the number of novel vertices yet, because it's the last in the partition.
  MultiMeshVertexIndex< 3 > vertices = multi_mesh_vertex_index( a, b, novelVertices_ );

  uint32_t bOffset = static_cast< uint32_t >( b.vertices.size() );

  // unify vertex pairs.
  for( const TriangleContacts& trianglePair : aContacts ) {

    for ( ContactPair contact : trianglePair.pairs ) {
      
      if ( contact.isVertexVertex() ) {

        const ConnectedTriangle& aTriangle = a.triangles[ trianglePair.this_triangle_index ];
        const ConnectedTriangle& bTriangle = b.triangles[ trianglePair.other_triangle_index ];

        uint32_t vertexA = aTriangle.vertices[ vertexIndex( contact.with ) ];
        uint32_t vertexB = vertices( 1, bTriangle.vertices[ vertexIndex( contact.against ) ] );

        unifiedVertices_.merge( vertexA, vertexB );
      }
    }
  }

  // Optimise here, because this will make the mappings above
  // single hop, but the novel vertices don't need to be, because they self reference
  unifiedVertices_.optimize();

  const std::vector< glm::dvec3 >& aVertices = a.vertices;
  const std::vector< glm::dvec3 >& bVertices = b.vertices;

  uint32_t novelOffset = bOffset + bVertices.size();

  for (
    uint32_t aTriangleIndex = 0, aTriangleEnd = static_cast< uint32_t >( a.triangles.size() );
    aTriangleIndex < aTriangleEnd;
    ++aTriangleIndex ) {

    std::span< const uint32_t > trianglePairs = aContactMap.get( aTriangleIndex );

    if ( trianglePairs.empty() ) {
      continue;
    }

    const ConnectedTriangle& aTriangle = a.triangles[ aTriangleIndex ];

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {
      insertLocalVertex( aTriangle.vertices[ vertexInTriangle ] );
    }

    for ( uint32_t pairIndex : trianglePairs ) {

      const TriangleContacts&  trianglePair = aContacts[ pairIndex ];
      const ConnectedTriangle& bTriangle    = b.triangles[ trianglePair.other_triangle_index ];

      assert( aTriangleIndex == trianglePair.this_triangle_index );

      std::span< const ContactPair > contactPairs = trianglePair.pairs.values();

      FixedStack< uint32_t, 6 > additionalVertices;

      for ( ContactPair contact : contactPairs ) {

        if ( contact.isEdgeEdge() ) {

          uint32_t edgeInTriangleA = edgeIndex( contact.with );
          uint32_t edge0           = aTriangle.edges[ edgeInTriangleA ];
          uint32_t edge1           = bTriangle.edges[ edgeIndex( contact.against ) ];

          const auto [ to, success ] =
            edgeEdgeVertices_.try_emplace(
              std::make_pair( edge0, edge1 ),
              static_cast< uint32_t >( novelVertices_.size() ) );

          if ( success ) {

            const glm::dvec3& e0v0 = aVertices[ a.edges[ edge0 ].vertices[ 0 ] ];
            const glm::dvec3& e0v1 = aVertices[ a.edges[ edge0 ].vertices[ 1 ] ];
            const glm::dvec3& e1v0 = bVertices[ b.edges[ edge1 ].vertices[ 0 ] ];
            const glm::dvec3& e1v1 = bVertices[ b.edges[ edge1 ].vertices[ 1 ] ];

            glm::dvec3 novelVertex =
              line_segment_line_segment_intersection( e0v0, e0v1, e1v0, e1v1 );

            novelVertices_.push_back( novelVertex ); // todo, add correct vertex generation.
            unifiedVertices_.allocate();
          }

          additionalVertices.push(
            insertLocalVertexOnEdge(
              vertices( 2, to->second ),
              edgeInTriangleA ) );

        } else if ( contact.isFaceEdge() ) {              

          uint32_t edge = bTriangle.edges[ edgeIndex( contact.against ) ];

          const auto [ to, success ] =
            faceEdgeVertices_[ 0 ].try_emplace(
              std::make_pair( trianglePair.this_triangle_index, edge ),
              static_cast< uint32_t >( novelVertices_.size() ) );

          if ( success ) {            

            const glm::dvec3& t0 = aVertices[ aTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = aVertices[ aTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = aVertices[ aTriangle.vertices[ 2 ] ];

            const glm::dvec3& ev0 = bVertices[ b.edges[ edge ].vertices[ 0 ] ];
            const glm::dvec3& ev1 = bVertices[ b.edges[ edge ].vertices[ 1 ] ];

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            novelVertices_.push_back( novelVertex );
            unifiedVertices_.allocate();
          }

          additionalVertices.push( insertLocalVertex( vertices( 2, to->second ) ) );

        } else if ( contact.isEdgeFace() ) {

          uint32_t edgeInTriangle = edgeIndex( contact.with );
          uint32_t edge           = aTriangle.edges[ edgeInTriangle ];

          const auto [ to, success ] =
            faceEdgeVertices_[ 1 ].try_emplace(
              std::make_pair( trianglePair.other_triangle_index, edge ),
              static_cast< uint32_t >( novelVertices_.size() ) );

          if ( success ) {

            const glm::dvec3& t0 = bVertices[ bTriangle.vertices[ 0 ] ];
            const glm::dvec3& t1 = bVertices[ bTriangle.vertices[ 1 ] ];
            const glm::dvec3& t2 = bVertices[ bTriangle.vertices[ 2 ] ];

            const glm::dvec3& ev0 = aVertices[ a.edges[ edge ].vertices[ 0 ] ];
            const glm::dvec3& ev1 = aVertices[ a.edges[ edge ].vertices[ 1 ] ];

            glm::dvec3 novelVertex =
              plane_line_segment_intersection( t0, t1, t2, ev0, ev1 );

            novelVertices_.push_back( novelVertex );
            unifiedVertices_.allocate();
          }

          additionalVertices.push(
            insertLocalVertexOnEdge( 
              vertices( 2, to->second ),
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

    triangulate( a, a, b, aTriangle, false, 0 );
  }

  for (
    uint32_t bTriangleIndex = 0, bTriangleEnd = static_cast< uint32_t >( b.triangles.size() );
    bTriangleIndex < bTriangleEnd;
    ++bTriangleIndex ) {

    std::span< const uint32_t > trianglePairs = bContactMap.get( bTriangleIndex );

    if ( trianglePairs.empty() ) {
      continue;
    }

    const ConnectedTriangle& bTriangle = b.triangles[ bTriangleIndex ];

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {
      insertLocalVertex( vertices( 1, bTriangle.vertices[ vertexInTriangle ] ) );
    }

    for ( uint32_t pairIndex : trianglePairs ) {

      const TriangleContacts&  trianglePair = bContacts[ pairIndex ];
      const ConnectedTriangle& aTriangle    = a.triangles[ trianglePair.other_triangle_index ];

      assert( bTriangleIndex == trianglePair.this_triangle_index );

      std::span< const ContactPair > contactPairs = trianglePair.pairs.values();

      FixedStack< uint32_t, 6 > additionalVertices;

      for ( ContactPair contact : contactPairs ) {

        if ( contact.isEdgeEdge() ) {

          uint32_t edgeInTriangleB = edgeIndex( contact.with );
          uint32_t edge0           = bTriangle.edges[ edgeInTriangleB ];
          uint32_t edge1           = aTriangle.edges[ edgeIndex( contact.against ) ];

          auto to = edgeEdgeVertices_.find( std::make_pair( edge1, edge0 ) );

          assert( to != edgeEdgeVertices_.end() );

          additionalVertices.push(
            insertLocalVertexOnEdge( vertices( 2, to->second ), edgeInTriangleB ) );

        }
        else if ( contact.isFaceEdge() ) {

          uint32_t edge = aTriangle.edges[ edgeIndex( contact.against ) ];

          const auto to =
            faceEdgeVertices_[1].find(
              std::make_pair( trianglePair.this_triangle_index, edge ) );

          assert( to != faceEdgeVertices_[ 1 ].end() );

          additionalVertices.push(
            insertLocalVertex( vertices( 2, to->second ) ) );

        }
        else if ( contact.isEdgeFace() ) {

          uint32_t edgeInTriangle = edgeIndex( contact.with );
          uint32_t edge           = bTriangle.edges[ edgeInTriangle ];

          const auto to =
            faceEdgeVertices_[0].find(
              std::make_pair( trianglePair.other_triangle_index, edge ) );

          assert( to != faceEdgeVertices_[ 0 ].end() );

          additionalVertices.push(
            insertLocalVertexOnEdge(
              vertices( 2, to->second ),
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
    
    triangulate( b, a, b, bTriangle, flipBWinding, 1 );
  }

  // now we do a winding-insensitive sort of the triangles.
  for ( uint32_t triangleSet = 0; triangleSet < 2; ++triangleSet ) {

    std::sort(
      initialChartTriangles_[ triangleSet ].begin(),
      initialChartTriangles_[ triangleSet ].end(),
      [] ( const Triangle& left, const Triangle& right ) {

        return less_lowest_vertex_parity( left.vertices, right.vertices );
      });
  }

  auto aWhere = initialChartTriangles_[ 0 ].begin();
  auto aEnd   = initialChartTriangles_[ 0 ].end();
  auto bWhere = initialChartTriangles_[ 1 ].begin();
  auto bEnd   = initialChartTriangles_[ 1 ].end();

  AABBTree& aBVH = *a.bvh; 
  AABBTree& bBVH = *b.bvh; 

  aBVH.dipoles( a );
  bBVH.dipoles( b );

  while ( aWhere < aEnd && bWhere < bEnd ) {

    const Triangle& aTriangle = *aWhere;
    const Triangle& bTriangle = *bWhere;

    if ( less_lowest_vertex_parity( aTriangle.vertices, bTriangle.vertices ) ) {

      glm::dvec3 aCentre = centroid( a, b, novelVertices_, aTriangle );

      double gwn = bBVH.gwn( b, aCentre );

      bool outside = fabs( gwn ) < 0.5;

      if ( outside == aOutside ) {

        outputTriangleStream_.push_back( aTriangle );
      }

      ++aWhere;
      continue;
    }

    if ( less_lowest_vertex_parity( bTriangle.vertices, aTriangle.vertices ) ) {

      glm::dvec3 bCentre = centroid( a, b, novelVertices_, bTriangle );

      double gwn = aBVH.gwn( a, bCentre );

      bool outside = fabs( gwn ) < 0.5;

      if ( outside == bOutside ) {

        outputTriangleStream_.push_back( bTriangle );
      }

      ++bWhere;
      continue;
    }

    bool windingParity =
      lowest_vertex_ordered_parity( aTriangle.vertices ) ==
      lowest_vertex_ordered_parity( bTriangle.vertices );

    // if these triangles are wound the same, keep A (cos we flip the winding on B
    // for subtraction, we can always keep A and throw away B if the winding is the same
    // and throw away both if the winding is different)
    if ( windingParity ) {

      outputTriangleStream_.push_back( aTriangle );
    }

    ++aWhere;
    ++bWhere;
  }
    
  while ( aWhere < aEnd ) {

    const Triangle& aTriangle = *aWhere;
    glm::dvec3      aCentre   = centroid( a, b, novelVertices_, aTriangle );

    double gwn = bBVH.gwn( b, aCentre );

    bool outside = fabs( gwn ) < 0.5;

    if ( outside == aOutside ) {

      outputTriangleStream_.push_back( aTriangle );
    }

    ++aWhere;
  }

  while ( bWhere < bEnd ) {

    const Triangle& bTriangle = *bWhere;
    glm::dvec3      bCentre   = centroid( a, b, novelVertices_, bTriangle );

    double gwn = aBVH.gwn( a, bCentre );

    bool outside = fabs( gwn ) < 0.5;

    if ( outside == bOutside ) {

      outputTriangleStream_.push_back( bTriangle );
    }

    ++bWhere;
  }

  walkAndInsertNonBoundary( aOutside, bBVH, boundarySet[ 0 ], a, b, false, 0 );
  walkAndInsertNonBoundary( bOutside, aBVH, boundarySet[ 1 ], b, a, flipBWinding, bOffset );
  
  uint32_t novelPartition = bOffset + b.vertices.size();

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

        mappedVertex = output.makeVertex( vertices[ unifiedVertexIndex ] );

        globalVertexMap_[ unifiedVertexIndex ] = mappedVertex;
        vertexUsed_[ unifiedVertexIndex ]      = true;

      } else {

        mappedVertex = globalVertexMap_[ unifiedVertexIndex ];
      }

      triangle.vertices[ vertexInTriangle ] = mappedVertex;
    }
    
    output.makeTriangle( triangle.vertices[ 0 ], triangle.vertices[ 1 ], triangle.vertices[2 ] );
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
  AABBTree& bvh,
  const std::vector< bool >& boundarySet,
  const WingedEdgeMesh< glm::dvec3 >& mesh,
  const WingedEdgeMesh< glm::dvec3 >& otherMesh,
  bool flippedWinding,
  uint32_t vertexOffset ) {

  walked_.clear();
  
  walked_.reserve( boundarySet.size() );
  walked_.insert( walked_.begin(), boundarySet.begin(), boundarySet.end() );

  uint32_t walkedCursor = 0;
  uint32_t walkedEnd    = static_cast< uint32_t >( boundarySet.size() );

  while ( walkedCursor < walkedEnd ) {

    if ( !walked_[ walkedCursor ] ) {

      const ConnectedTriangle& initialTriangle = mesh.triangles[ walkedCursor ];
      glm::dvec3               centre          = centroid( mesh.vertices, initialTriangle );
      bool                     triangleOutside = fabs( bvh.gwn( otherMesh, centre ) ) < 0.5;

      triangleStack_.push_back(  walkedCursor );

      if ( triangleOutside == outside ) {

        while ( !triangleStack_.empty() ) {

          uint32_t nextTriangleIndex = triangleStack_.back();

          triangleStack_.pop_back();

          if ( walked_[ nextTriangleIndex ] ) {

            continue;
          }

          walked_[ nextTriangleIndex ] = true;

          const ConnectedTriangle& triangle = mesh.triangles[ nextTriangleIndex ];

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
            
            uint32_t opposingTriangle =
              mesh.edges[ triangle.edges[ triangleInEdge ] ].otherTriangle( nextTriangleIndex );

            if ( opposingTriangle != EMPTY_INDEX && !walked_[ opposingTriangle ] ) {

              triangleStack_.push_back( opposingTriangle );
            }
          }
        }
      } else {

        while ( !triangleStack_.empty() ) {

          uint32_t nextTriangleIndex = triangleStack_.back();

          triangleStack_.pop_back();

          if (walked_[nextTriangleIndex]) {

            continue;
          }

          walked_[nextTriangleIndex] = true;

          const ConnectedTriangle& triangle = mesh.triangles[ nextTriangleIndex ];

          for ( uint32_t triangleInEdge = 0; triangleInEdge < 3; ++triangleInEdge ) {
          
            uint32_t opposingTriangle = mesh.edges[ triangle.edges[ triangleInEdge ] ].otherTriangle( nextTriangleIndex );

            if ( opposingTriangle != EMPTY_INDEX && !walked_[ opposingTriangle ] ) {

              triangleStack_.push_back( opposingTriangle );
            }
          }
        }
      }
   }

    ++walkedCursor;
  }
}

void conway::geometry::CSGMesher::addEdges(
  std::span< const ContactPair > contactPairs,
  FixedStack< uint32_t, 6 >& additionalVertices ) {

  size_t innerEnd = contactPairs.size();

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

void conway::geometry::CSGMesher::triangulate(
  const WingedEdgeMesh< glm::dvec3 >& mesh,
  const WingedEdgeMesh< glm::dvec3 >& a,
  const WingedEdgeMesh< glm::dvec3 >& b,
  const ConnectedTriangle& triangle,
  bool flippedWinding,
  uint32_t outputStreamIndex ) {

  // Find any duplicates (including merged vertices) and make them unique.
  for ( uint32_t localVertex : localVertices_ ) {

    uint32_t foundVertice = unifiedVertices_.find( localVertex );

    localVertexMap_.try_emplace( foundVertice, static_cast< uint32_t >( localVertexMap_.size() ) );
  }

  // Cut off the size of what's been removed.
  localVertices_.resize( localVertexMap_.size() );

  // set the vertices back to their new locations.
  for ( auto [ localVertex, index ] : localVertexMap_ ) {

    localVertices_[ index ] = localVertex;
  }

  // Remap edges with the new unique vertices.
  for ( CDT::Edge& edge: edges_ ) {

    edge =
      CDT::Edge(
        localVertexMap_[ unifiedVertices_.find( edge.v1() ) ],
        localVertexMap_[ unifiedVertices_.find( edge.v2() ) ] );
  }

  // Now use the best 2D projection to extract the vertices.
  glm::dvec3 triangleVertices[ 3 ];

  extract_vertices( mesh, triangle, triangleVertices );

  AxisPair axesToExtract = best_truncated_projection(triangleVertices);

  int32_t winding = orient2D( triangleVertices, axesToExtract );

  MultiMeshVertexIndex< 3 > vertices = multi_mesh_vertex_index( a, b, novelVertices_ );

  local2DVertices_.reserve( localVertices_.size() );
  localVertexEdgeFlags_.clear();
  localVertexEdgeFlags_.resize( localVertices_.size(), 0 );

  size_t firstAxis  = first_axis(axesToExtract);
  size_t secondAxis = second_axis(axesToExtract);

  for ( uint32_t partitionedIndice : localVertices_)  {

    const glm::dvec3& inputVertex = vertices[ partitionedIndice ];

    local2DVertices_.push_back( CDT::V2d< double >::make( inputVertex[ firstAxis ], inputVertex[ secondAxis ] ) );
  }

  double edgeErrorFactor = 0;

  for ( size_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle ) {

    std::vector< uint32_t >& edgeVertices = onEdgeVertices_[ edgeInTriangle ];

    uint32_t v0Index = edgeInTriangle;
    uint32_t v1Index = (edgeInTriangle + 1) % 3;

    uint8_t edgeFlag = 1 << static_cast< uint8_t >( edgeInTriangle );

    localVertexEdgeFlags_[ v0Index ] |= edgeFlag;
    localVertexEdgeFlags_[ v1Index ] |= edgeFlag;

    if ( edgeVertices.empty() ) {

      edges_.emplace_back( v0Index, v1Index );
      continue;
    }

    double maxOrient2DError = 0;

    glm::dvec2 v0 = extract( triangleVertices[ v0Index ], axesToExtract );
    glm::dvec2 v1 = extract( triangleVertices[ v1Index ], axesToExtract );

    for ( uint32_t& onEdgeVertex : edgeVertices ) {

      onEdgeVertex = localVertexMap_[ unifiedVertices_.find( onEdgeVertex ) ];

      const CDT::V2d< double >& v2 = local2DVertices_[ onEdgeVertex ];

      localVertexEdgeFlags_[ onEdgeVertex ] |= edgeFlag;

      maxOrient2DError = std::max( fabs( predicates::adaptive::orient2d( &v0.x, &v1.x, &v2.x ) ), maxOrient2DError );
    }

    if ( edgeVertices.size() > 1 ) {

      glm::dvec2 interval    = v1 - v0;
      glm::dvec2 absInterval = glm::abs( interval );

      size_t sortAxis   = first_axis( axesToExtract );
      size_t sortAxis2D = 0;

      if ( absInterval.y > absInterval.x ) {

        sortAxis    = second_axis ( axesToExtract );
        sortAxis2D = 1;
      }

      bool ascendsTowardsV1 = interval[ sortAxis2D ] > 0;

      if ( ascendsTowardsV1 ) {

        std::sort(
          edgeVertices.begin(),
          edgeVertices.end(),
          [ sortAxis, &vertices ]( uint32_t left, uint32_t right ) {

            return vertices[ left ][ sortAxis ] < vertices[ right ][ sortAxis ];
          }
        );
      }
      else {

        std::sort(
          edgeVertices.begin(),
          edgeVertices.end(),
          [ sortAxis, &vertices ]( uint32_t left, uint32_t right ) {

            return vertices[ left ][ sortAxis ] > vertices[ right ][ sortAxis ];
          }
        );
      }

      {
        const CDT::V2d< double >& e0 = local2DVertices_[ v0Index ];
        const CDT::V2d< double >& e1 = local2DVertices_[ edgeVertices[ 1 ] ];
        const CDT::V2d< double >& v2 = local2DVertices_[ edgeVertices[ 0 ] ];

        double orient2DError = fabs(predicates::adaptive::orient2d(&e0.x, &e1.x, &v2.x));

        orient2DError *= 1 + DBL_EPSILON;

        edgeErrorFactor = std::max(edgeErrorFactor, orient2DError / (CDT::distance(e0, e1) * (1 - DBL_EPSILON)));
      }

      {
        const CDT::V2d< double >& e0 = local2DVertices_[ edgeVertices[ edgeVertices.size() - 2 ] ];
        const CDT::V2d< double >& e1 = local2DVertices_[ v1Index ];
        const CDT::V2d< double >& v2 = local2DVertices_[ edgeVertices.back() ];

        double orient2DError = fabs(predicates::adaptive::orient2d(&e0.x, &e1.x, &v2.x));

        orient2DError *= 1 + DBL_EPSILON;

        edgeErrorFactor = std::max(edgeErrorFactor, orient2DError / (CDT::distance(e0, e1) * (1 - DBL_EPSILON)));
      }
    }

    // Add the leading/trailing constraints
    edges_.emplace_back( v0Index, edgeVertices.front() );
    edges_.emplace_back( edgeVertices.back(), v1Index );

    for (
      size_t vertInEdgeIndex = 0, end = edgeVertices.size() - 1;
      vertInEdgeIndex < end;
      ++vertInEdgeIndex ) {
      
      edges_.emplace_back( edgeVertices[ vertInEdgeIndex ], edgeVertices[ vertInEdgeIndex + 1 ] );
    }

    for (
      size_t vertInEdgeIndex = 1, end = edgeVertices.size() - 1;
      vertInEdgeIndex < end;
      ++vertInEdgeIndex) {
      
      const CDT::V2d< double >& e0 = local2DVertices_[ edgeVertices[ vertInEdgeIndex - 1 ] ];
      const CDT::V2d< double >& e1 = local2DVertices_[ edgeVertices[ vertInEdgeIndex + 1 ] ];
      const CDT::V2d< double >& v2 = local2DVertices_[ edgeVertices[ vertInEdgeIndex ] ];

      double orient2DError = fabs( predicates::adaptive::orient2d( &e0.x, &e1.x, &v2.x ) );

      orient2DError *= 1 + DBL_EPSILON;
      
      edgeErrorFactor = std::max( edgeErrorFactor, orient2DError / ( CDT::distance( e0, e1 ) * ( 1 - DBL_EPSILON ) ) );

      edges_.emplace_back( edgeVertices[ vertInEdgeIndex ], edgeVertices[ vertInEdgeIndex + 1 ] );
    }

    // This should effectively give a rounding hedged version of the error
    // in determinining these are on a constraint edge.
    maxOrient2DError *= 1 + DBL_EPSILON;
    edgeErrorFactor   = std::max( edgeErrorFactor, maxOrient2DError / ( glm::distance( v0, v1 ) * ( 1 - DBL_EPSILON ) ) );
  }

  CDT::Triangulation< double > triangulation( CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::DontCheck, 2.0 * edgeErrorFactor );

  triangulation.insertVertices( local2DVertices_ );
  triangulation.insertEdges( edges_ );
  triangulation.eraseOuterTriangles();

  assert( triangulation.triangles.size() > 0 );

  int32_t foundWinding         = 0;
  size_t  windingTriangleIndex = 0;

  for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {

    uint32_t localTriangle[ 3 ] = {

      localVertices_[ cdtTriangle.vertices[ 0 ] ],
      localVertices_[ cdtTriangle.vertices[ 1 ] ],
      localVertices_[ cdtTriangle.vertices[ 2 ] ]

    };

    foundWinding = orient2D( a, b, novelVertices_, localTriangle, axesToExtract );

    if ( foundWinding != 0 ) {

      break;
    }
  }

  std::vector< Triangle >& outputStream = initialChartTriangles_[ outputStreamIndex ];

  bool flipWinding = ( ( flippedWinding ? -1 : 1 ) * foundWinding * winding ) < 0;

  if ( flipWinding ) {

    for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {

      uint8_t edgeMask =
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 0 ] ] &
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 1 ] ] &
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 2 ] ];

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
        localVertices_[ cdtTriangle.vertices[ 0 ] ],
        localVertices_[ cdtTriangle.vertices[ 1 ] ],
        localVertices_[ cdtTriangle.vertices[ 2 ] ]
      }};

      reorder_to_lowest_vertex( localTriangle.vertices );

      std::swap( localTriangle.vertices[ 1 ], localTriangle.vertices[ 2 ] );

      outputStream.push_back( localTriangle );
    }

  } else {

    for ( const CDT::Triangle& cdtTriangle : triangulation.triangles ) {
      
      uint8_t edgeMask =
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 0 ] ] &
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 1 ] ] &
        localVertexEdgeFlags_[ cdtTriangle.vertices[ 2 ] ];

      if ( edgeMask != 0 ) {
        continue;
      }

      Triangle localTriangle { { 
        localVertices_[ cdtTriangle.vertices[ 0 ] ],
        localVertices_[ cdtTriangle.vertices[ 1 ] ],
        localVertices_[ cdtTriangle.vertices[ 2 ] ]
      }};

      reorder_to_lowest_vertex( localTriangle.vertices );

      outputStream.push_back( localTriangle );
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

uint32_t conway::geometry::CSGMesher::insertLocalVertex( uint32_t inputVertex ) {

  localVertices_.push_back( inputVertex );

  return inputVertex;
}

uint32_t conway::geometry::CSGMesher::insertLocalVertexOnEdge( uint32_t inputVertex, uint32_t edgeInTriangle ) {

  localVertices_.push_back( inputVertex );

  onEdgeVertices_[ edgeInTriangle ].push_back( inputVertex );

  return inputVertex;
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