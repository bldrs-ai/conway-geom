#pragma once

#include <cmath>
#include <glm/glm.hpp>

#include "../representation/winged_edge.h"
#include <queue>
#include "../representation/IfcGeometry.h"

namespace conway::geometry {

  /**
   * A UV parameterized vertex on 2 parameter surface.
   */
  struct ParameterVertex {

    glm::dvec3 point;
    glm::dvec2 uv;
    
  };

  /**
   * Compute the normal of a parameter vertex.
   * 
   * Note we normalize the edge vectors for stability.
   */
  glm::dvec3 computeNormal(
    const ParameterVertex& v0,
    const ParameterVertex& v1,
    const ParameterVertex& v2
  ) {

    glm::dvec3 v01(v1.point - v0.point);
    glm::dvec3 v02(v2.point - v0.point);

    glm::dvec3 norm = 
      glm::cross( glm::normalize( v01 ), glm::normalize( v02 ) );

    return glm::normalize( norm );
  }

  /**
   * Calulate the area of the triangle
   */
  double computeArea(
    const ParameterVertex& v0,
    const ParameterVertex& v1,
    const ParameterVertex& v2
  ) {
    glm::dvec3 v01(v1.point - v0.point);
    glm::dvec3 v02(v2.point - v0.point);

    glm::dvec3 norm = glm::cross(v01, v02);

    double result = glm::length( norm );

    if ( std::isnan( result ) ) {
      result = 0.001;
    }

    return result;
  }

  /**
   * Is a triangle wound counterclockwise in UV space.
   */
  bool isCCW(
    const ParameterVertex& v0,
    const ParameterVertex& v1,
    const ParameterVertex& v2
  ) {
    double a = v1.uv.x * v0.uv.y + v2.uv.x * v1.uv.y + v0.uv.x * v2.uv.y;
    double b = v0.uv.x * v1.uv.y + v1.uv.x * v2.uv.y + v2.uv.x * v0.uv.y;

    return a < b;
  }

  /**
   * A candidate edge for splitting with the parameter vertex.
   */
  struct CandidateEdge {

    double   deflection;
    uint32_t edge;

    ParameterVertex vertex;
  };

  /**
   * Sorting operator for candidate edge priority.
   */
  inline bool operator<( const CandidateEdge& left, const CandidateEdge& right ) {

    return ( left.deflection < right.deflection ) ||
      ( left.deflection == right.deflection && left.edge < right.edge );
  }

  /**
   * Append a winged edge mesh 
   */
  inline void appendMeshToGeometry( WingedEdgeMesh< ParameterVertex >& mesh, IfcGeometry& geometry ) {

    uint32_t baseVertex = geometry.vertexData.size() / 6;

    std::vector< glm::dvec3 > normals( mesh.vertices.size(), glm::dvec3( 0 ) );

    for ( Triangle& triangle : mesh.triangles ) {

      if ( !isCCW(  
        mesh.vertices[ triangle.vertices[ 0 ] ],
        mesh.vertices[ triangle.vertices[ 1 ] ],
        mesh.vertices[ triangle.vertices[ 2 ] ]) ) {

        std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
        std::swap( triangle.edges[ 0 ], triangle.edges[ 2 ] );
      }

      glm::dvec3 triangleAreaNormal = 
        computeNormal(  
          mesh.vertices[ triangle.vertices[ 0 ] ],
          mesh.vertices[ triangle.vertices[ 1 ] ],
          mesh.vertices[ triangle.vertices[ 2 ] ] ) *
        ( computeArea(   
          mesh.vertices[ triangle.vertices[ 0 ] ],
          mesh.vertices[ triangle.vertices[ 1 ] ],
          mesh.vertices[ triangle.vertices[ 2 ] ] ) );

      if ( std::isnan( triangleAreaNormal.x ) ) {
        continue;
      }

      for ( uint32_t localVertex = 0; localVertex < 3; ++localVertex ) {
        normals[ triangle.vertices[ localVertex ] ] += triangleAreaNormal;
      }
    }

    for ( size_t vertexIndex = 0, end = mesh.vertices.size(); vertexIndex < end; ++vertexIndex ) {

      glm::dvec3 vertex = mesh.vertices[ vertexIndex ].point;
      const glm::dvec3& prenormal = normals[ vertexIndex ];
      glm::dvec3 normal = glm::normalize( prenormal );


      // Note, we have to have local versions of vertex and normal
      // cos addpoint isn't const correct.
      geometry.AddPoint( vertex, normal );
    }

    for ( const Triangle& triangle : mesh.triangles ) {   

      geometry.AddFace(
        baseVertex + triangle.vertices[ 0 ],
        baseVertex + triangle.vertices[ 1 ],
        baseVertex + triangle.vertices[ 2 ] );
    }
  }

  /**
   * Given a parameterized surface (UV)->(XYZ),
   * this will take a starting mesh with parameterized vertices and tesselate the internal triangles
   */
  template< typename SurfacePointFunction >
  inline void tesselate(
    WingedEdgeMesh< ParameterVertex >& mesh,
    SurfacePointFunction surface,
    int32_t maximumTriangles,
    double minimumDeflection ) {

    std::priority_queue< CandidateEdge > candidates;

    auto addCandidate = [&]( uint32_t edgeIndex ) {

      if ( edgeIndex == EMPTY_INDEX  ) {
        return;
      }

      const Edge& edge = mesh.edges[ edgeIndex ];

      if ( edge.border() ) {
        return;
      }

      const ParameterVertex& v0   = mesh.vertices[ edge.vertices[ 0 ] ];
      const ParameterVertex& v1   = mesh.vertices[ edge.vertices[ 1 ] ];

      glm::dvec3 averagePoint = ( v0.point + v1.point ) * 0.5;
      glm::dvec2 newUV        = ( v0.uv + v1.uv ) * 0.5;
      glm::dvec3 newPoint     = surface( averagePoint, newUV );

      double deflection = glm::distance( averagePoint, newPoint );

      if ( minimumDeflection > deflection ) {
        return;
      }

      candidates.push( CandidateEdge { 
        deflection * glm::distance( v0.point, v1.point ),
        edgeIndex,
        ParameterVertex { newPoint, newUV } 
        } );
    };

    for (
      uint32_t edgeIndex = 0, end = static_cast< uint32_t >( mesh.edges.size() );
      edgeIndex < end;
      ++edgeIndex ) {

      addCandidate( edgeIndex );
    }

    maximumTriangles -= mesh.triangles.size();

    while ( !candidates.empty() && maximumTriangles > 0 ) {

      const CandidateEdge& candidate    = candidates.top();
      // copy edge because it mutates later
      // as may the references as the vector re-allocates.
      Edge                 edge         = mesh.edges[ candidate.edge ];

      const Triangle&      t0           = mesh.triangles[ edge.triangles[ 0 ] ];
      const Triangle&      t1           = mesh.triangles[ edge.triangles[ 1 ] ];
      uint32_t             otherVertex0 = t0.otherVertex( edge );
      uint32_t             otherVertex1 = t1.otherVertex( edge );
      uint32_t             newVertex    = mesh.makeVertex( candidate.vertex );

      candidates.pop();

      mesh.deleteTriangle( edge.triangles[ 1 ] );
      mesh.deleteTriangle( edge.triangles[ 0 ] );

      mesh.makeTriangle( otherVertex0, edge.vertices[ 0 ], newVertex );
      mesh.makeTriangle( newVertex, edge.vertices[ 1 ], otherVertex0 );
      mesh.makeTriangle( newVertex, edge.vertices[ 0 ], otherVertex1 );
      mesh.makeTriangle( otherVertex1, edge.vertices[ 1 ], newVertex );

      addCandidate( mesh.getEdge( otherVertex0, newVertex ).value_or( EMPTY_INDEX ) );
      addCandidate( mesh.getEdge( otherVertex1, newVertex ).value_or( EMPTY_INDEX ) );
      addCandidate( mesh.getEdge( edge.vertices[ 0 ], newVertex ).value_or( EMPTY_INDEX ) );
      addCandidate( mesh.getEdge( edge.vertices[ 1 ], newVertex ).value_or( EMPTY_INDEX ) );

      maximumTriangles -= 2;
    }
  }
}
