#include "csg.h"

/**
 * Creates an index of the two meshes interactions so queries can be run.
 * 
 * This is called by run internally, and does not need to be called before calling run.
 */
void conway::geometry::CSG::index( WingedEdgeMesh< glm::dvec3 >& a, WingedEdgeMesh< glm::dvec3 >& b, double tolerance = 0 ) {

  double tolerance2 = tolerance * tolerance;
  double tolerance3 = tolerance2 * tolerance;

  reset();

  a.makeBVH();
  b.makeBVH();

  AABBTree& bvhA = *a.bvh;
  AABBTree& bvhB = *b.bvh;

  std::vector< TriangleTriangleClassification >& candidatePairs = candidate_pairs;

  uint32_t bVertexOffset = a.vertices.size();
  uint32_t bEdgeOffset   = a.edges.size();
  uint32_t bFaceOffset   = a.triangles.size();

  uint32_t meshVertexOffsets[ 2 ] = { 0, bVertexOffset };
  uint32_t meshEdgeOffsets[ 2 ]   = { 0, bEdgeOffset };
  uint32_t meshFaceOffsets[ 2 ]   = { 0, bFaceOffset };

  boundary[ 0 ].resize( a.triangles.size(), false );
  boundary[ 1 ].resize( b.triangles.size(), false );

  face_planes.allocate( a.triangles.size() + b.triangles.size() );

  bvhA.intersect( bvhB, [&]( uint32_t left, uint32_t right ) {

    TriangleTriangleClassification candidate;

    const ConnectedTriangle& aTriangle = a.triangles[ left ];
    const ConnectedTriangle& bTriangle = b.triangles[ right ];

    uint32_t vertIndices[ 2 ][ 3 ];
  
    for ( uint32_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {

      vertIndices[ 0 ][ vertInTriangle ] = aTriangle.vertices[ vertInTriangle ];
      vertIndices[ 1 ][ vertInTriangle ] = bTriangle.vertices[ vertInTriangle ] + bVertexOffset;
    }

    // pre-calculate our face to vertex values and face-to-face colinearity.
    for ( uint32_t triangleInFace = 0; triangleInFace < 2; ++triangleInFace ) {

      candidate.face_to_face = 
        vertexFaceSigns(
          a,
          b,
          vertIndices[ triangleInFace ],
          vertIndices[ triangleInFace ^ 1 ],
          candidate.triangles[ triangleInFace ],
          tolerance3 );

      // if a face doesn't straddle the other face's plane, don't add it to the candidate list.
      // even if we end up in some degenerate AABB tree case, this should stop memory exploding,
      // apart from the one weird degenerate case where someone has made a whole bunch of triangles
      // of roughly even sized ad different angles with roughly the same centroid.
      if ( candidate.face_to_face == FaceFace::NONE ) {
        return;
      }

      if ( candidate.face_to_face == FaceFace::COLINEAR ) {

        face_planes.merge( left, right + a.triangles.size() );
      }
    }

    candidate.triangles[ 0 ].triangle_index = left;
    candidate.triangles[ 1 ].triangle_index = right;

    candidatePairs.push_back( candidate );
  } );

  vertices.allocate( a.vertices.size() + b.vertices.size() );

  // Second pass, with all our unified verts, we can now calculate predicates
  // and populate the interactions between the triangles to find intersection verts
  for ( TriangleTriangleClassification& candidate : candidatePairs ) {

    const ConnectedTriangle& aTriangle = a.triangles[ candidate.triangles[ 0 ].triangle_index ];
    const ConnectedTriangle& bTriangle = b.triangles[ candidate.triangles[ 1 ].triangle_index ];

    const ConnectedTriangle* triangles[ 2 ] = { &aTriangle, &bTriangle };

    uint32_t vertIndices[ 2 ][ 3 ];
  
    for ( uint32_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {

      vertIndices[ 0 ][ vertInTriangle ] = aTriangle.vertices[ vertInTriangle ];
      vertIndices[ 1 ][ vertInTriangle ] = bTriangle.vertices[ vertInTriangle ] + bVertexOffset;
    }

    TriangleInteractions& aPredicates = candidate.triangles[ 0 ];
    TriangleInteractions& bPredicates = candidate.triangles[ 1 ];

    uint32_t unityCount = 0;

    uint32_t hasColinear[ 2 ]    = {};
    AxisPair safeProjection[ 2 ] = {AxisPair::X_Y, AxisPair::X_Y };

    for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {

      hasColinear[ 0 ] += bPredicates.vertex_face_signs[ vertexInTriangle ] == 0;
      hasColinear[ 1 ] += aPredicates.vertex_face_signs[ vertexInTriangle ] == 0;
    }

    glm::dvec3 localVertices[ 2 ][ 3 ];

    extract_vertices( a, aTriangle, localVertices[ 0 ] );
    extract_vertices( b, bTriangle, localVertices[ 1 ] );

    // Signs for colinear vertices relative to edges.
    int32_t vertexEdgeSigns[ 2 ][ 3 ][ 3 ] = {};

    std::optional< AxisPair > safeProjections[ 2 ];

    TriangleContacts contacts[ 2 ];

    for ( uint32_t triangleInPair = 0; triangleInPair < 2; ++triangleInPair ) {

      contacts[ triangleInPair ].this_triangle_index  =
        candidate.triangles[ triangleInPair ].triangle_index;
      contacts[ triangleInPair ].other_triangle_index =
        candidate.triangles[ triangleInPair ^ 1 ].triangle_index;
    }

    // If one of the vertices in a face is colinear to the other, make sure we have a safe 2D projection.
    // This grabs the truncated axis orthogonal projection with the largest area.
    // It is more numerically robust than using the normal.
    for ( uint32_t triangleInPair = 0; triangleInPair < 2; ++triangleInPair ) {

      std::optional< AxisPair >& safeProjection = safeProjections[ triangleInPair ];

      for ( uint32_t vertexInOpposingTriangle = 0; vertexInOpposingTriangle < 3; ++vertexInOpposingTriangle ) {

        uint32_t opposingTriangleInPair = triangleInPair ^ 1;
        int8_t   opposingPointSign      = candidate.triangles[ opposingTriangleInPair ].vertex_face_signs[ vertexInOpposingTriangle ];

        if ( opposingPointSign == 0 ) {

          if ( safeProjection == std::nullopt ) {

            safeProjection = best_truncated_projection( localVertices[ triangleInPair ] );
          }

          uint32_t opposingVertexIndex = vertIndices[ opposingTriangleInPair ][ vertexInOpposingTriangle ];

          int      zeroCount     = 0;
          int      positiveCount = 0; 
          uint32_t lastZeroEdge  = 0;

          int32_t (&edgeSigns)[ 3 ] = vertexEdgeSigns[ opposingTriangleInPair ][ vertexInOpposingTriangle ];

          for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle ) {

            int vertexEdgeSign = orient2D(
                a,
                b,
                vertIndices[ triangleInPair ],
                edgeInTriangle,
                opposingVertexIndex,
                *safeProjection,
                tolerance2
              );

            edgeSigns[ edgeInTriangle ] = vertexEdgeSign;

            if ( vertexEdgeSign == 0 ) {

              lastZeroEdge = edgeInTriangle;
              ++zeroCount;
            }

            positiveCount += vertexEdgeSign > 0;
          }

          // this is a degenerate case, we cannot be on all 3 edges unless
          // they are coincident at this vertex, which means a zero area triangle.
          assert( zeroCount != 3 );

          if ( zeroCount == 2 ) {

            // In this case we *have* to be inside the opposing edge.
            // as we're a vertex on the other face.
            assert( positiveCount == 1 );

            uint32_t matchingVertexInFace =
              ( lastZeroEdge + ( edgeSigns[ lastZeroEdge - 1 ] != 0 ) ) % 3;

            // unify these vertices, as this is a colinear vertex also colinear two edges,
            // and therefore equal opposing edge vertex.
            vertices.merge( opposingVertexIndex, vertIndices[ triangleInPair ][ matchingVertexInFace ] );
                
            contacts[ triangleInPair ].vertexVertex( matchingVertexInFace, vertexInOpposingTriangle );
            contacts[ opposingTriangleInPair ].vertexVertex( vertexInOpposingTriangle, matchingVertexInFace );

          } else if ( zeroCount == 1 && positiveCount == 2 ) {

            // This vertex is colinear to an edge, and if it's inside of two edges,
            // as well as colinear the face, it's a vertex edge pair. 
            // Note, unlike the paper, we don't need an edge 1D case
            // by being within the wedge of the edges, we have implicitly 
            // decided this vertex is on an edge, and
            // if the next or previous is on the same edge, that will
            // drop out.
            contacts[ triangleInPair ].edgeEdge( lastZeroEdge, vertexInOpposingTriangle );
            contacts[ opposingTriangleInPair ].edgeEdge( vertexInOpposingTriangle, lastZeroEdge );

          } else if ( positiveCount == 3 ) {

            // This vertex is colinear to the face and within all 3 edges
            // so it is a face vertex pair
            contacts[ triangleInPair ].faceVertex( vertexInOpposingTriangle );
            contacts[ opposingTriangleInPair ].vertexFace( vertexInOpposingTriangle );
          }
        }
      }
    }

    bool foundColinear = false;

    for ( uint32_t edgeInTriangle0 = 0; edgeInTriangle0 < 3; ++edgeInTriangle0 ) {

      uint32_t t0v0 = edgeInTriangle0;
      uint32_t t0v1 = ( edgeInTriangle0 + 1 ) % 3;

      int32_t v0sign = aPredicates.vertex_face_signs[ edgeInTriangle0 ]; 
      int32_t v1sign = aPredicates.vertex_face_signs[ ( edgeInTriangle0 + 1 ) % 3 ];

      bool v0Zero = v0sign == 0;
      bool v1Zero = v1sign == 0;

      if ( v0Zero && v1Zero ) {

        foundColinear = true;

        // Colinear case.
        for ( uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1 ) {

          uint32_t t1v0 = edgeInTriangle1;
          uint32_t t1v1 = ( edgeInTriangle1 + 1 ) % 3;

          int32_t t0v0e1sign = vertexEdgeSigns[ 0 ][ t0v0 ][ edgeInTriangle1 ];
          int32_t t0v1e1sign = vertexEdgeSigns[ 0 ][ t0v1 ][ edgeInTriangle1 ];

          if ( t0v0e1sign * t0v1e1sign != -1 ) {
            continue;
          }

          int32_t opposingv0sign = bPredicates.vertex_face_signs[ edgeInTriangle1 ]; 
          int32_t opposingv1sign = bPredicates.vertex_face_signs[ ( edgeInTriangle1 + 1 ) % 3 ];

          // if signs are the same, we can only straddle if co-linear to the face
          if ( opposingv0sign == opposingv1sign ) {

            if ( opposingv0sign != 0 ) {
              continue; 
            }

            // it's a double face co-linear edge straddle showdown
            int32_t t1v0e0sign = vertexEdgeSigns[ 1 ][ t0v0 ][ edgeInTriangle0 ];
            int32_t t1v1e0sign = vertexEdgeSigns[ 1 ][ t0v1 ][ edgeInTriangle0 ];
        
            if ( t1v0e0sign * t1v1e0sign != -1 ) {
              continue;
            }

            contacts[ 0 ].edgeEdge( edgeInTriangle0, edgeInTriangle1 );
            contacts[ 1 ].edgeEdge( edgeInTriangle1, edgeInTriangle0 );

          } else if ( opposingv0sign * opposingv1sign == -1 ) {

            uint32_t t1v0indice = vertIndices[ 1 ][ t1v0 ];
            uint32_t t1v1indice = vertIndices[ 1 ][ t1v1 ];

            // we have a case where the first edge is colinear to the second face and straddles the line of the second edge,
            // this can only be an edge-edge case, and the second face edge has to straddle the line of the first,
            // but in the 2D projection of the second face.
            int32_t t1v0Sign = orient2D( a, b, vertIndices[ 0 ], edgeInTriangle0, t1v0indice, *safeProjections[ 1 ], tolerance2 );
            int32_t t1v1Sign = orient2D( a, b, vertIndices[ 0 ], edgeInTriangle0, t1v1indice, *safeProjections[ 1 ], tolerance2 );

            if ( t1v0Sign * t1v1Sign != -1 ) {

              continue;
            }
            
            contacts[ 0 ].edgeEdge( edgeInTriangle0, edgeInTriangle1 );
            contacts[ 1 ].edgeEdge( edgeInTriangle1, edgeInTriangle0 );
          }
        }
      } else if ( !v0Zero && !v1Zero && v0sign != v1sign ) {

        for ( uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1 ) {

          int32_t opposingv0sign = bPredicates.vertex_face_signs[ edgeInTriangle1 ]; 
          int32_t opposingv1sign = bPredicates.vertex_face_signs[ ( edgeInTriangle1 + 1 ) % 3 ];

          // if signs are the same, we can only straddle if co-linear to the face
          if ( opposingv0sign != opposingv1sign || opposingv0sign != 0 ) {
            continue;
          }

          foundColinear = true; 

          uint32_t t0v0indice = vertIndices[ 0 ][ t0v0 ];
          uint32_t t0v1indice = vertIndices[ 0 ][ t0v1 ];

          // we have a case where the second edge is colinear to the first face and straddles the line of the first edge,
          // this can only be an edge-edge case, and the first face edge has to straddle the line of the second,
          // but in the 2D projection of the first face.
          int32_t t0v0Sign = orient2D( a, b, vertIndices[ 1 ], edgeInTriangle1, t0v0indice, *safeProjections[ 0 ], tolerance2 );
          int32_t t0v1Sign = orient2D( a, b, vertIndices[ 1 ], edgeInTriangle1, t0v1indice, *safeProjections[ 0 ], tolerance2 );

          // note we can break here, only one edge can be colinear to the face on this triangle
          // if the edge on the other triangle isn't colinear.
          if ( t0v0Sign * t0v1Sign != -1 ) {

            break;
          }

          contacts[ 0 ].edgeEdge( edgeInTriangle0, edgeInTriangle1 );
          contacts[ 1 ].edgeEdge( edgeInTriangle1, edgeInTriangle0 );
          break;
        }
      }
    }
  
    // There can't be any face edge cases if we have a colinear edge.
    if ( foundColinear || candidate.face_to_face == FaceFace::COLINEAR ) {

      add( contacts );
      continue;
    }

    std::optional< int32_t > orient3Dcache[ 3 ][ 3 ];

    for ( uint32_t edgeInTriangle0 = 0; edgeInTriangle0 < 3; ++edgeInTriangle0 ) {

      uint32_t t0v0 = edgeInTriangle0;
      uint32_t t0v1 = ( edgeInTriangle0 + 1 ) % 3;

      int32_t v0sign = aPredicates.vertex_face_signs[ edgeInTriangle0 ]; 
      int32_t v1sign = aPredicates.vertex_face_signs[ ( edgeInTriangle0 + 1 ) % 3 ];

      if ( v0sign * v1sign != -1 ) {
        continue;
      }
      
      int32_t t0v0indice = vertIndices[ 0 ][ edgeInTriangle0 ]; 
      int32_t t0v1indice = vertIndices[ 0 ][ ( edgeInTriangle0 + 1 ) % 3 ];

      int32_t t1v0indice = vertIndices[ 1 ][ 0 ]; 
      int32_t t1v1indice = vertIndices[ 1 ][ 1 ];
      int32_t t1v2indice = vertIndices[ 1 ][ 2 ];

      std::optional< int32_t > (&o)[ 3 ] = orient3Dcache[ edgeInTriangle0 ];

      o[ 0 ] = orient3D( a, b, t0v0indice, t0v1indice, t1v0indice, t1v1indice, tolerance3 );
      o[ 1 ] = orient3D( a, b, t0v0indice, t0v1indice, t1v1indice, t1v2indice, tolerance3 );

      if ( *o[ 0 ] * *o[ 1 ] < 0 ) {

        continue;
      }

      o[ 2 ] = orient3D( a, b, t0v0indice, t0v1indice, t1v2indice, t1v0indice, tolerance3 );

      if ( *o[ 1 ] * *o[ 2 ] < 0  || *o[ 2 ] * *o[ 0 ] < 0 ) {

        continue;
      }

      uint32_t inEdgeCount = 0;

      for ( uint32_t opposingEdge = 0; opposingEdge < 3; ++opposingEdge ) {

        // note, unlike the paper, we've already dealt with edge-vertex cases before here, because
        // they require a vertex->face co-linearity for the vertex and will be dealt with in 2D.
        if ( *o[ opposingEdge ] == 0 ) {

          if ( 
            *o[ ( opposingEdge + 1 ) % 3 ] != 0 && 
            *o[ ( opposingEdge + 2 ) % 3 ] != 0 ) {

            int32_t ov0sign = bPredicates.vertex_face_signs[ opposingEdge ]; 
            int32_t ov1sign = bPredicates.vertex_face_signs[ ( opposingEdge + 1 ) % 3 ];

            // guarantee this is not an edge vertex case
            if ( ov0sign * ov1sign == -1 ) {

              contacts[ 0 ].edgeEdge( edgeInTriangle0, opposingEdge );
              contacts[ 1 ].edgeEdge( opposingEdge, edgeInTriangle0 );
            }
            break;
          }
        } else {

          ++inEdgeCount;
        }
      }

      if ( inEdgeCount == 3 ) {

        contacts[ 0 ].edgeFace( edgeInTriangle0 );
        contacts[ 1 ].faceEdge( edgeInTriangle0 );
      }
    }

    {
      int32_t t0v0indice = vertIndices[ 0 ][ 0 ]; 
      int32_t t0v1indice = vertIndices[ 0 ][ 1 ];
      int32_t t0v2indice = vertIndices[ 0 ][ 2 ];

      for ( uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1 ) {

        int32_t o[ 3 ] = {};

        int32_t t1v0indice = vertIndices[ 1 ][ edgeInTriangle1 ]; 
        int32_t t1v1indice = vertIndices[ 1 ][ ( edgeInTriangle1 + 1 ) % 3 ];

        o[ 0 ] =
          orient3Dcache[ 0 ][ edgeInTriangle1 ].has_value() ? 
            *orient3Dcache[ 0 ][ edgeInTriangle1 ] :
            orient3D( a, b, t0v0indice, t0v1indice, t1v0indice, t1v1indice, tolerance3 );

        o[ 1 ] =
          orient3Dcache[ 1 ][ edgeInTriangle1 ].has_value() ? 
            *orient3Dcache[ 1 ][ edgeInTriangle1 ] :
            orient3D( a, b, t0v1indice, t0v2indice, t1v0indice, t1v1indice, tolerance3 );

        if ( o[ 0 ] * o[ 1 ] <= 0 ) {

          continue;
        }

        o[ 2 ] =
          orient3Dcache[ 2 ][ edgeInTriangle1 ].has_value() ? 
            *orient3Dcache[ 2 ][ edgeInTriangle1 ] :
            orient3D( a, b, t0v2indice, t0v0indice, t1v0indice, t1v1indice, tolerance3 );

        if ( o[ 1 ] * o[ 2 ] <= 0 || o[ 2 ] * o[ 0 ] <= 0 ) {

          continue;
        }

        // we already handled the edge-edge cases in the previous loop with the other triangle,
        // this is strictly this edge stabbing the first face.
        contacts[ 0 ].faceEdge( edgeInTriangle1 );
        contacts[ 1 ].edgeFace( edgeInTriangle1 );
      }
    }

    add( contacts );
  }

  for ( size_t where = 0; where < 2; ++where ) {

    face_contact_map[ where ].construct(
      contacts[ where ],
      static_cast< uint32_t >( a.triangles.size() ),
      []( const TriangleContacts& contact ) { return contact.this_triangle_index; } );
  }
}