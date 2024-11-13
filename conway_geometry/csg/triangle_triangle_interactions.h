#pragma once

#include "triangle_contacts.h"
#include "structures/winged_edge.h"
#include "csg_utils.h"
#include <optional>

namespace conway::geometry {

  /** State for calculating triangle triangle intersections, this is expected
    * to live on the stack and outlive a and b (which should have no triangles added)
    */
  struct TriangleTriangleInteractions {

    // these are bucket constants
    static constexpr size_t NEGATIVE = 0;
    static constexpr size_t ZERO     = 1;
    static constexpr size_t POSITIVE = 2;

    static constexpr size_t ALL_TRI_VERTICES = 3;
    static constexpr size_t NO_TRI_VERTICES  = 3;

    int32_t vertex_signs[ 2 ][ 3 ]       = {};
    
    uint32_t vertex_sign_counts[ 2 ][ 3 ] = {};

    int32_t vertex_edge_signs[ 2 ][ 3 ][ 3 ] = {};
    
    std::optional< AxisPair > projection_axis[ 2 ];

    const WingedEdgeDV3* meshes[ 2 ];

    const ConnectedTriangle* triangles[ 2 ];

    uint32_t triangle_in_mesh[ 2 ];

    TriangleContacts contacts[ 2 ];

    glm::dvec3 vertices[ 2 ][ 3 ];

    double tolerance;

    TriangleTriangleInteractions( const WingedEdgeDV3& a, const WingedEdgeDV3& b, uint32_t triangleInMeshA, uint32_t triangleInMeshB, double _tolerance ) :
      meshes { &a, &b },
      triangle_in_mesh { triangleInMeshA, triangleInMeshB },
      tolerance( _tolerance ),
      triangles { &a.triangles[ triangleInMeshA ], &b.triangles[ triangleInMeshB ] } {

      extract_vertices( a, a.triangles[ triangleInMeshA ], vertices[ 0 ] );
      extract_vertices( b, b.triangles[ triangleInMeshB ], vertices[ 1 ] );
    }

    bool evaluate();

  private:

    bool evaluateVertexSigns( uint32_t triangleInPair ) {

      uint32_t opposingTriangleInPair = triangleInPair ^ 1;

      const glm::dvec3 (&triangleVertices)[ 3 ] = vertices[ triangleInPair ];
      const glm::dvec3 (&opposingVertices)[ 3 ] = vertices[ opposingTriangleInPair ];

      int32_t (&vertexSigns)[ 3 ]       = vertex_signs[ triangleInPair ];
      uint32_t (&vertexSignCounts)[ 3 ] = vertex_sign_counts[ triangleInPair ];

      for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {

        int32_t sign = 
            orient3D( 
              opposingVertices[ 0 ],
              opposingVertices[ 1 ],
              opposingVertices[ 2 ],
              triangleVertices[ vertexInTriangle ],
              tolerance );

        vertexSigns[ vertexInTriangle ] = sign;
        vertexSignCounts[ sign + 1 ]   += 1;
      }

      // all vertices are either one side or the other
      if (
        vertexSignCounts[ NEGATIVE ] == ALL_TRI_VERTICES ||
        vertexSignCounts[ POSITIVE ] == ALL_TRI_VERTICES ) {

        return false;
      }
    }

    bool evaluateColinearVertices( uint32_t triangleInPair ) {

      uint32_t (&vertexSignCounts)[ 3 ] = vertex_sign_counts[ triangleInPair ];

      if ( vertexSignCounts[ ZERO ] == NO_TRI_VERTICES ) {

        return false;
      }

      bool result = false;

      uint32_t         opposingTriangleInPair       = triangleInPair ^ 1;      
      int32_t          (&vertexSigns)[ 3 ]          = vertex_signs[ triangleInPair ];
      const glm::dvec3 (&opposingVertices)[ 3 ]     = vertices[ opposingTriangleInPair ];
      int32_t          (&vertexEdgeSigns)[ 3 ][ 3 ] = vertex_edge_signs[ triangleInPair ];

      AxisPair projectionAxis = best_truncated_projection( opposingVertices );

      projection_axis[ opposingTriangleInPair ] = projectionAxis;

      const glm::dvec3 (&triangleVertices)[ 3 ] = vertices[ triangleInPair ];

      for ( uint32_t vertexInTriangle = 0; vertexInTriangle < 3; ++vertexInTriangle ) {

        int32_t vertexSign = vertexSigns[ vertexInTriangle ];

        if ( vertexSign != 0 ) {

          continue;
        }
        
        const glm::dvec3& vertex  = triangleVertices[ vertexInTriangle ];
        int32_t (&edgeSigns)[ 3 ] = vertexEdgeSigns[ vertexInTriangle ];

        uint32_t signCounts[ 3 ] = {};

        for ( uint32_t edgeVertex0 = 0; edgeVertex0 < 3; ++edgeVertex0 ) {

          uint32_t edgeVertex1 = ( edgeVertex0 + 1 ) % 3;

          int32_t sign =
            orient2D(
              opposingVertices[ edgeVertex0 ],
              opposingVertices[ edgeVertex1 ],
              vertex,
              projectionAxis,
              tolerance );

          signCounts[ sign + 1 ]  += 1;
          edgeSigns[ edgeVertex0 ] = sign;
        }

        if ( signCounts[ NEGATIVE ] == 3 || signCounts[ POSITIVE ] == 3 ) {

          contacts[ triangleInPair ].vertexFace( vertexInTriangle );
          contacts[ opposingTriangleInPair ].faceVertex( vertexInTriangle );

          result = true;
          continue;
        }

        assert( signCounts[ ZERO ] != 3 );

        if ( signCounts[ ZERO ] == 2 ) {

          // An edge is indexed is the same as the first vertex in it,
          // and in this case we test the preceeding edge (which has its second vertex)
          // followed by the next edge (which has this as its first).
          // If both have zero counts, it means the "next edge" is a vertex
          // vertex pair.
          for ( uint32_t preceedingEdge = 0; preceedingEdge < 3; ++preceedingEdge ) {

            if ( edgeSigns[ preceedingEdge ] != 0 ) {

              continue;
            }

            uint32_t nextEdge = ( preceedingEdge + 1 ) % 3;

            if ( edgeSigns[ nextEdge ] != 0 ) {

              continue;
            }

            contacts[ triangleInPair ].vertexVertex( vertexInTriangle, nextEdge );
            contacts[ opposingTriangleInPair ].vertexVertex( nextEdge, vertexInTriangle );

            result = true;
            break;
          }
        }
        else if ( signCounts[ ZERO ] == 1 ) {

          // In this case, the vertex is colinear to an edge
          // find the edge, and make sure we are on the same side
          // of the other edges.
          for ( uint32_t edge = 0; edge < 3; ++edge ) {
          
            if ( edgeSigns[ edge ] != 0 ) {

              continue;
            }

            if ( edgeSigns[ ( edge + 1 ) % 3 ] != edgeSigns[ ( edge + 2 ) % 3 ] ) {

              break;
            }

            contacts[ triangleInPair ].vertexEdge( vertexInTriangle, edge );
            contacts[ opposingTriangleInPair ].edgeVertex( edge, vertexInTriangle );

            result = true;
            break;
          }
        }
      }

      return result;
    }

    bool evaluateVertexSigns() {

      return evaluateVertexSigns( 0 ) && evaluateVertexSigns( 1 );
    }

    bool evaluateColinearVertices() {

      bool result = evaluateVertexSigns( 0 );

      result |= evaluateVertexSigns( 1 );

      return result;
    }
  };

  inline bool TriangleTriangleInteractions::evaluate() {

    bool hasStraddleOrColinear = evaluateVertexSigns();

    if ( !hasStraddleOrColinear ) {
    
      return false;
    }

    bool result = evaluateColinearVertices()

    return !contacts[ 0 ].pairs.empty() || !contacts[ 1 ].pairs.empty();
  }
}