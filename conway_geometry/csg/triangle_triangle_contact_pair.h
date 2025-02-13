#pragma once

#include <stdint.h>
#include <utility>
#include "triangle_contacts.h"
#include "multi_mesh_vertex_index.h"

namespace conway::geometry {


  struct TriangleTriangleContactPair {

    TriangleContacts triangles[ 2 ] = {};

    bool empty() const {

      return
        triangles[ 0 ].face_to_face == FaceFace::NONE ||
        triangles[ 1 ].face_to_face == FaceFace::NONE ||
        triangles[ 0 ].pairs.empty() ||
        triangles[ 1 ].pairs.empty();
    }
  };


  inline FaceFace vertex_face_signs(
    const glm::dvec3( &triangleVertices )[ 3 ],
    const glm::dvec3( &opposingTriangleVertices )[ 3 ],
    /*out*/ int32_t(&vertexFaceSigns)[3],
    double tolerance
  ) {

    int32_t signs[ 3 ] = {};

    for ( uint32_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle ) {

      // if a vertex has been unified, orient3D will *always* be 0.
      // and we should propagate that, even if orient3D should 
      // be *exact* in theory and we are paranoid about input
      // ordering to make it as deterministic as possible.
      int32_t sign =
        orient3D(
          opposingTriangleVertices[ 0 ],
          opposingTriangleVertices[ 1 ],
          opposingTriangleVertices[ 2 ],
          triangleVertices[ vertInTriangle ],
          tolerance );

      signs[ sign + 1 ] += 1;

      // Classify this vertex relative a face.
      vertexFaceSigns[ vertInTriangle ] = sign;
    }

    if ( signs[ 1 ] == 3 ) {

      return FaceFace::COLINEAR;
    }

    if ( signs[ 0 ] == 3 || signs[ 2 ] == 3 ) {

      return FaceFace::NONE;
    }

    return FaceFace::STRADDLES_OR_TOUCHES;
  }

  /** Find the combination of intersecting elements for a candidate pair if they intersect */
  inline TriangleTriangleContactPair find_intersections(
    const MultiMeshVertexIndex< 2 >& vertices,
    const uint32_t( &triangleInMeshIndices )[ 2 ],
    const Triangle* ( &triangles )[ 2 ],
    bool selfIntersectionOnly = false,
    bool shareEdge = false,
    double tolerance = 0 ) {

    TriangleTriangleContactPair result;

    result.triangles[ 1 ].other_triangle_index =
    result.triangles[ 0 ].this_triangle_index  = triangleInMeshIndices[ 0 ];
    
    result.triangles[ 0 ].other_triangle_index =
    result.triangles[ 1 ].this_triangle_index  = triangleInMeshIndices[ 1 ];

    const Triangle&      aTriangle      = *triangles[ 0 ];
    const Triangle&      bTriangle      = *triangles[ 1 ];


    glm::dvec3 localVertices[ 2 ][ 3 ] = {};

    vertices.extract( 0, aTriangle.vertices, localVertices[ 0 ] );
    vertices.extract( 1, bTriangle.vertices, localVertices[ 1 ] );

    // These uniquely identify vertices using a partition, and also
    // handle the case where we are doing self intersections.
    uint32_t vertIndices[ 2 ][ 3 ] = {};

    for (uint32_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle) {

      vertIndices[ 0 ][ vertInTriangle ] = vertices( 0, aTriangle.vertices[ vertInTriangle ] );
      vertIndices[ 1 ][ vertInTriangle ] = vertices( 1, bTriangle.vertices[ vertInTriangle ] );
    }

    int32_t vertexFaceSigns[ 2 ][ 3 ] = {};

    int32_t (&aFaceSigns)[ 3 ] = vertexFaceSigns[ 0 ];
    int32_t (&bFaceSigns)[ 3 ] = vertexFaceSigns[ 1 ];

    FaceFace aFaceFace = vertex_face_signs( localVertices[ 0 ], localVertices[ 1 ], aFaceSigns, 0 );
    
    result.triangles[ 0 ].face_to_face = aFaceFace;

    // If we share an edge, we only care about the co-linear case (degenerate)
    if ( aFaceFace == FaceFace::NONE ||
      ( selfIntersectionOnly && shareEdge && aFaceFace != FaceFace::COLINEAR ) ) {

      return result;
    }
    
    FaceFace bFaceFace = vertex_face_signs( localVertices[ 1 ], localVertices[ 0 ], bFaceSigns, 0 );
        
    result.triangles[ 1 ].face_to_face = bFaceFace;

    if ( bFaceFace == FaceFace::NONE ) {

      return result;
    }

    AxisPair safeProjections[2] = {
      best_truncated_projection( localVertices[ 0 ] ),
      best_truncated_projection( localVertices[ 1 ] )
    };

    // This gets us our axis pairs, but also guarantees the triangle is non-zero area.
    if ( safeProjections[ 0 ] == AxisPair::NONE || safeProjections[ 1 ] == AxisPair::NONE ) {

      return result;
    }

    if ( aFaceFace != bFaceFace ) {

      // using exact maths, we should never get here if both of the projections above aren't none
      // however on some very borderline cases we've seen this happen in practice... which 
      // maybe down to some slight bugs in floating point code generation.
      // However, this definitely implies some kind of degenerate configuration where there 
      // is a linear dependency (duplicate vertex or some kind of colinearity).
      return result;
    }

    // Signs for colinear vertices relative to edges.
    int32_t vertexEdgeSigns[ 2 ][ 3 ][ 3 ] = {};

    // If one of the vertices in a face is colinear to the other, make sure we have a safe 2D projection.
    // This grabs the truncated axis orthogonal projection with the largest area.
    // It is more numerically robust than using the normal.
    for (uint32_t triangleInPair = 0; triangleInPair < 2; ++triangleInPair) {

      AxisPair safeProjection = safeProjections[ triangleInPair ];

      for ( uint32_t vertexInOpposingTriangle = 0; vertexInOpposingTriangle < 3; ++vertexInOpposingTriangle ) {

        uint32_t opposingTriangleInPair = triangleInPair ^ 1;
        int32_t  opposingPointSign      = vertexFaceSigns[ opposingTriangleInPair ][ vertexInOpposingTriangle ];

        if ( opposingPointSign == 0 ) {

          uint32_t opposingVertexIndex = vertIndices[ opposingTriangleInPair ][ vertexInOpposingTriangle ];

          int      zeroCount     = 0;
          int      positiveCount = 0;
          int      negativeCount = 0;
          uint32_t lastZeroEdge  = 0;

          int32_t (&edgeSigns)[ 3 ] = vertexEdgeSigns[ opposingTriangleInPair ][ vertexInOpposingTriangle ];

          for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle ) {

            int vertexEdgeSign =
              vertices.orient2D(
                vertIndices[ triangleInPair ],
                edgeInTriangle,
                opposingVertexIndex,
                safeProjection,
                0 );

            edgeSigns[ edgeInTriangle ] = vertexEdgeSign;

            if ( vertexEdgeSign == 0 ) {

              lastZeroEdge = edgeInTriangle;
              ++zeroCount;
            }

            positiveCount += vertexEdgeSign > 0;
            negativeCount += vertexEdgeSign < 0;
          }

          if ( zeroCount == 3 ) {

            result.triangles[ 0 ].pairs.clear();
            result.triangles[ 1 ].pairs.clear();

            // Degenerate case, must be an infinitely small triangle.
            return result;
          }

          if (zeroCount == 2) {

            uint32_t matchingVertexInFace =
              (lastZeroEdge + (edgeSigns[lastZeroEdge - 1] != 0)) % 3;

            assert( localVertices[ triangleInPair ][ matchingVertexInFace ] == localVertices[ opposingTriangleInPair ][ vertexInOpposingTriangle ] );

            if ( !selfIntersectionOnly || vertIndices[ triangleInPair ][ matchingVertexInFace ] != vertIndices[ opposingTriangleInPair ][ vertexInOpposingTriangle ] ) {
              result.triangles[ triangleInPair ].vertexVertex( matchingVertexInFace, vertexInOpposingTriangle );
              result.triangles[ opposingTriangleInPair ].vertexVertex( vertexInOpposingTriangle, matchingVertexInFace );
            }
          }
          else if (zeroCount == 1 && (positiveCount == 2 || negativeCount == 2)) {

            // This vertex is colinear to an edge, and if it's inside of two edges,
            // as well as colinear the face, it's a vertex edge pair. 
            // Note, unlike the paper, we don't need an edge 1D case
            // by being within the wedge of the edges, we have implicitly 
            // decided this vertex is on an edge, and
            // if the next or previous is on the same edge, that will
            // drop out.
            result.triangles[ triangleInPair ].edgeVertex( lastZeroEdge, vertexInOpposingTriangle );
            result.triangles[ opposingTriangleInPair ].vertexEdge( vertexInOpposingTriangle, lastZeroEdge );
          }
          else if (positiveCount == 3 || negativeCount == 3) {

            // This vertex is colinear to the face and within all 3 edges
            // so it is a face vertex pair
            result.triangles[ triangleInPair ].faceVertex( vertexInOpposingTriangle );
            result.triangles[ opposingTriangleInPair ].vertexFace( vertexInOpposingTriangle );
          }
        }
      }
    }

    bool foundColinear = false;

    for (uint32_t edgeInTriangle0 = 0; edgeInTriangle0 < 3; ++edgeInTriangle0) {

      uint32_t t0v0 = edgeInTriangle0;
      uint32_t t0v1 = (edgeInTriangle0 + 1) % 3;

      int32_t v0sign = aFaceSigns[ t0v0 ];
      int32_t v1sign = aFaceSigns[ t0v1 ];

      bool v0Zero = v0sign == 0;
      bool v1Zero = v1sign == 0;

      if (v0Zero && v1Zero) {

        foundColinear = true;

        // Colinear case.
        for (uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1) {

          uint32_t t1v0 = edgeInTriangle1;
          uint32_t t1v1 = (edgeInTriangle1 + 1) % 3;

          int32_t t0v0e1sign = vertexEdgeSigns[ 0 ][ t0v0 ][ edgeInTriangle1 ];
          int32_t t0v1e1sign = vertexEdgeSigns[ 0 ][ t0v1 ][ edgeInTriangle1 ];

          if ( t0v0e1sign * t0v1e1sign != -1 ) {
            continue;
          }

          int32_t opposingv0sign = bFaceSigns[ t1v0 ];
          int32_t opposingv1sign = bFaceSigns[ t1v1 ];

          // if signs are the same, we can only straddle if co-linear to the face
          if ( opposingv0sign == opposingv1sign ) {

            if ( opposingv0sign != 0 ) {
              continue;
            }

            // it's a double face co-linear edge straddle showdown
            int32_t t1v0e0sign = vertexEdgeSigns[ 1 ][ t1v0 ][ edgeInTriangle0 ];
            int32_t t1v1e0sign = vertexEdgeSigns[ 1 ][ t1v1 ][ edgeInTriangle0 ];

            if ( t1v0e0sign * t1v1e0sign != -1 ) {
              continue;
            }

            result.triangles[ 0 ].edgeEdge2D( edgeInTriangle0, edgeInTriangle1 );
            result.triangles[ 1 ].edgeEdge2D( edgeInTriangle1, edgeInTriangle0 );
          }
          else if ( opposingv0sign * opposingv1sign == -1 ) {

            uint32_t t1v0indice = vertIndices[ 1 ][ t1v0 ];
            uint32_t t1v1indice = vertIndices[ 1 ][ t1v1 ];

            //// we have a case where the first edge is colinear to the second face and straddles the line of the second edge,
            //// this can only be an edge-edge case, and the second face edge has to straddle the line of the first,
            //// but in the 2D projection of the second face.
            int32_t t1v0Sign = vertices.orient2D( vertIndices[ 0 ], edgeInTriangle0, t1v0indice, safeProjections[ 1 ], 0 );
            int32_t t1v1Sign = vertices.orient2D( vertIndices[ 0 ], edgeInTriangle0, t1v1indice, safeProjections[ 1 ], 0 );

            if (t1v0Sign * t1v1Sign != -1) {

              continue;
            }

            result.triangles[ 0 ].edgeEdge2D( edgeInTriangle0, edgeInTriangle1 );
            result.triangles[ 1 ].edgeEdge( edgeInTriangle1, edgeInTriangle0 );
          }
        }
      }
      else if ( v0sign * v1sign == -1 ) {

        for ( uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1 ) {

          uint32_t t1v0 = edgeInTriangle1;
          uint32_t t1v1 = ( edgeInTriangle1 + 1 ) % 3;

          int32_t opposingv0sign = bFaceSigns[ t1v0 ];
          int32_t opposingv1sign = bFaceSigns[ t1v1 ];

          // if signs are the same, we can only straddle if co-linear to the face
          // for the 2D case
          if ( opposingv0sign != 0 || opposingv1sign != 0 ) {
            continue;
          }

          foundColinear = true;

          int32_t t1v0e1sign = vertexEdgeSigns[ 1 ][ t1v0 ][ edgeInTriangle0 ];
          int32_t t1v1e1sign = vertexEdgeSigns[ 1 ][ t1v1 ][ edgeInTriangle0 ];

          if ( t1v0e1sign * t1v1e1sign != -1 ) {
            continue;
          }

          uint32_t t0v0indice = vertIndices[ 0 ][ t0v0 ];
          uint32_t t0v1indice = vertIndices[ 0 ][ t0v1 ];

          // we have a case where the second edge is colinear to the first face and straddles the line of the first edge,
          // this can only be an edge-edge case, and the first face edge has to straddle the line of the second,
          // but in the 2D projection of the first face.
          int32_t t0v0Sign = vertices.orient2D( vertIndices[ 1 ], edgeInTriangle1, t0v0indice, safeProjections[ 0 ], 0 );
          int32_t t0v1Sign = vertices.orient2D( vertIndices[ 1 ], edgeInTriangle1, t0v1indice, safeProjections[ 0 ], 0 );

          // note we can break here, only one edge can be colinear to the face on this triangle
          // if the edge on the other triangle isn't colinear.
          if ( t0v0Sign * t0v1Sign != -1 ) {

            break;
          }

          result.triangles[ 0 ].edgeEdge( edgeInTriangle0, edgeInTriangle1 );
          result.triangles[ 1 ].edgeEdge2D( edgeInTriangle1, edgeInTriangle0 );
          break;
        }
      }
    }

    // There can't be any face edge cases if we have a colinear edge.
    if ( foundColinear || aFaceFace == FaceFace::COLINEAR ) {

      return result;
    }

    for (uint32_t edgeInTriangle0 = 0; edgeInTriangle0 < 3; ++edgeInTriangle0) {

      uint32_t t0v0 = edgeInTriangle0;
      uint32_t t0v1 = (edgeInTriangle0 + 1) % 3;

      int32_t v0sign = aFaceSigns[ t0v0 ];
      int32_t v1sign = aFaceSigns[ t0v1 ];

      if ( ( v0sign * v1sign ) != -1 ) {
        continue;
      }

      uint32_t t0v0indice = vertIndices[ 0 ][ t0v0 ];
      uint32_t t0v1indice = vertIndices[ 0 ][ t0v1 ];

      uint32_t t1v0indice = vertIndices[ 1 ][ 0 ];
      uint32_t t1v1indice = vertIndices[ 1 ][ 1 ];
      uint32_t t1v2indice = vertIndices[ 1 ][ 2 ];

      int32_t o[ 3 ] = {};

      o[ 0 ] = vertices.orient3D( t0v0indice, t0v1indice, t1v0indice, t1v1indice, 0 );
      o[ 1 ] = vertices.orient3D( t0v0indice, t0v1indice, t1v1indice, t1v2indice, 0 );

      if ( o[ 0 ] * o[ 1 ] < 0) {

        continue;
      }

      o[ 2 ] = vertices.orient3D( t0v0indice, t0v1indice, t1v2indice, t1v0indice, 0 );

      if ( o[ 1 ] * o[ 2 ] < 0 || o[ 2 ] * o[ 0 ] < 0) {

        continue;
      }

      uint32_t inEdgeCount = 0;

      for (uint32_t opposingEdge = 0; opposingEdge < 3; ++opposingEdge) {

        // note, unlike the paper, we've already dealt with edge-vertex cases before here, because
        // they require a vertex->face co-linearity for the vertex and will be dealt with in 2D.
        if ( o[ opposingEdge ] == 0 ) {

          uint32_t t1v1 = ( opposingEdge + 1 ) % 3;
          uint32_t t1v2 = ( opposingEdge + 2 ) % 3;

          if ( o[ t1v1 ] * o[ t1v2 ] == 1 ) {

            int32_t ov0sign = bFaceSigns[ opposingEdge ];
            int32_t ov1sign = bFaceSigns[ t1v1 ]; // NOTE - changed this because it should be the second vertex int he opposing edge - CS

            // guarantee this is not an edge vertex case
            if ( ov0sign * ov1sign == -1 ) {

              result.triangles[ 0 ].edgeEdge( edgeInTriangle0, opposingEdge );
              result.triangles[ 1 ].edgeEdge( opposingEdge, edgeInTriangle0 );
            }
            break;
          }
        }
        else {

          ++inEdgeCount;
        }
      }

      if ( inEdgeCount == 3 ) {

        result.triangles[ 0 ].edgeFace( edgeInTriangle0 );
        result.triangles[ 1 ].faceEdge( edgeInTriangle0 );
      }
    }

    {
      int32_t t0v0indice = vertIndices[0][0];
      int32_t t0v1indice = vertIndices[0][1];
      int32_t t0v2indice = vertIndices[0][2];

      for (uint32_t edgeInTriangle1 = 0; edgeInTriangle1 < 3; ++edgeInTriangle1) {

        int32_t o[3] = {};

        uint32_t v0sign = bFaceSigns[ edgeInTriangle1 ];
        uint32_t v1sign = bFaceSigns[ (edgeInTriangle1 + 1) % 3 ];

        if ( v0sign * v1sign != -1 ) {
          continue;
        }

        int32_t t1v0indice = vertIndices[1][edgeInTriangle1];
        int32_t t1v1indice = vertIndices[1][(edgeInTriangle1 + 1) % 3];

        o[0] = vertices.orient3D(t1v0indice, t1v1indice, t0v0indice, t0v1indice, 0 );
        o[1] = vertices.orient3D(t1v0indice, t1v1indice, t0v1indice, t0v2indice, 0 );

        if ( o[0] * o[1] <= 0 ) {

          continue;
        }

        o[2] =
          vertices.orient3D(t1v0indice, t1v1indice, t0v2indice, t0v0indice, 0 );

        if (o[1] * o[2] <= 0 || o[2] * o[0] <= 0) {

          continue;
        }

        // we already handled the edge-edge cases in the previous loop with the other triangle,
        // this is strictly this edge stabbing the first face.
        result.triangles[ 0 ].faceEdge( edgeInTriangle1 );
        result.triangles[ 1 ].edgeFace( edgeInTriangle1 );
      }

      return result;
    }
  }
}
