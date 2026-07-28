#pragma once

#include <cmath>
#include <stdexcept>
#include <glm/glm.hpp>

#include "structures/winged_edge.h"
#include "structures/alloc_telemetry.h"
#include <atomic>
#include <queue>
#include "representation/Geometry.h"
#include "representation/IfcGeometryReps.h"
#include "structures/hash_functions.h"
#include "representation/ConwayCurve.h"

#include <iostream>
#include <fstream>

#include <unordered_set>
#include <unordered_map>

#include "representation/Geometry.h"

#if defined(_MSC_VER)

#pragma warning(push)
#pragma warning(disable : 26812)

#endif

#include "CDT.h"

#if defined(_MSC_VER)

#pragma warning(pop)

#endif

namespace conway::geometry
{
  struct LoopEdge
  {

    uint32_t vertex1;
    uint32_t vertex2;
    uint32_t loop;

    LoopEdge(uint32_t v1, uint32_t v2, uint32_t _loop = EMPTY_INDEX) :
      vertex1(v1),
      vertex2(v2),
      loop( _loop )
    {
      if (v2 < v1)
      {
        std::swap(vertex1, vertex2);
      }
    }

    bool operator==(const LoopEdge &other) const
    {
      return ( vertex1 == other.vertex1 && vertex2 == other.vertex2 );
    }
  };

  using VertexMap = std::unordered_map<glm::dvec3, uint32_t>;

  struct LoopEdgeMesh
  {
    VertexMap vertexMap;
    std::vector<glm::dvec3> normalFormVertices;
    std::vector<LoopEdge> loopEdges;
  };


  template <
      typename ParameterizationFunction
      >
  inline void tesselatePlane(
      Geometry& output,
      const std::vector< IfcBound3D > &boundaries,
      ParameterizationFunction parameterization )
  {
    LoopEdgeMesh boundaryMesh;

    VertexMap               &vertexMap = boundaryMesh.vertexMap;
    std::vector<glm::dvec3> &vertices  = output.vertices;
    std::vector< LoopEdge > &loopEdges = boundaryMesh.loopEdges;

    uint32_t loopIndex = 0;
 
     // Build connected edge structure with unique vertices from boundary loops.
    for ( const IfcBound3D &boundary : boundaries )
    {
      const std::vector<glm::dvec3> &points =
          boundary.curve.points;

      bool     firstVertexInLoop   = true;
      uint32_t previousVertexIndex = 0;
      uint32_t firstIndexInLoop    = 0;

      for ( const glm::dvec3 &point : points )
      {
        auto foundVertex = vertexMap.find( point );

        uint32_t vertexIndex = 0;

        if ( foundVertex == vertexMap.end() )
        {
          vertexIndex = static_cast< uint32_t >( vertices.size() );

          vertices.push_back( point );

          vertexMap[ point ] = vertexIndex;
        }
        else
        {
          vertexIndex = foundVertex->second;
        }

        if ( firstVertexInLoop )
        {
          previousVertexIndex = vertexIndex;
          firstIndexInLoop    = vertexIndex;
          firstVertexInLoop   = false;
          continue;
        }

        if ( previousVertexIndex == vertexIndex )
        {
          continue;
        }
        
        loopEdges.emplace_back( previousVertexIndex, vertexIndex, loopIndex );
        
        previousVertexIndex = vertexIndex;
      }

      if ( previousVertexIndex != firstIndexInLoop )
      {
        loopEdges.emplace_back( previousVertexIndex, firstIndexInLoop, loopIndex );
      }

      ++loopIndex;
    }
    
    std::vector< glm::dvec2 > parameterizationVertices;

    size_t vertexCount = vertices.size();

    parameterizationVertices.reserve( vertexCount );

    for ( const glm::dvec3 &vertex : vertices )
    {
      parameterizationVertices.emplace_back( parameterization( vertex ) );
    }

    std::vector< CDT::Edge > cdtEdges;
    std::vector<uint32_t>    vertexRemapping;
    std::vector<uint32_t>    reverseVertexRemapping(
      vertexCount,
      EMPTY_INDEX );
    std::vector< CDT::V2d< double > > cdtVertices;

    vertexRemapping.reserve( vertexCount );
    cdtVertices.reserve( vertexCount );

    // Setup constraints, discard illegal edges and remap vertices for CDT.
    for ( const LoopEdge &loopEdge : boundaryMesh.loopEdges )
    {
      auto [ v1, v2, _ ] = loopEdge;

      const glm::dvec2 &param1 = parameterizationVertices[ v1 ];
      const glm::dvec2 &param2 = parameterizationVertices[ v2 ];

      uint32_t remappedV1 = reverseVertexRemapping[ v1 ];
      uint32_t remappedV2 = reverseVertexRemapping[ v2 ];

      if ( remappedV1 == EMPTY_INDEX )
      {
        remappedV1 = static_cast<uint32_t>( cdtVertices.size() );

        vertexRemapping.push_back( v1 );

        reverseVertexRemapping[ v1 ] = remappedV1;

        cdtVertices.emplace_back( param1.x, param1.y );
      }

      if ( remappedV2 == EMPTY_INDEX )
      {
        remappedV2 = static_cast< uint32_t >( cdtVertices.size() );

        vertexRemapping.push_back( v2 );

        reverseVertexRemapping[ v2 ] = remappedV2;

        cdtVertices.emplace_back( param2.x, param2.y );
      }

      if ( remappedV1 == remappedV2 )
      {
        continue; // Skip edges with same vertices.
      }

      cdtEdges.emplace_back( remappedV1, remappedV2 );
    }

    // atomic: tesselatePlane may run concurrently on the thread pool during
    // staged face finalization.
    static std::atomic< uint32_t > svgIndex    = 0;
    uint32_t                       outputIndex = svgIndex++;

// Output SVG for debugging.
#if (OUTPUT_SVG_DEBUG == 1)

    printf( "Writing planar triangulation SVG %d\n", outputIndex );

    std::ofstream svgFile( "planar_" + std::to_string( outputIndex ) + ".svg" );

    glm::dvec2 minBound(  std::numeric_limits<double>::max(),  std::numeric_limits<double>::max() );
    glm::dvec2 maxBound( -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() );

    for ( const CDT::Edge &edge : cdtEdges )
    {
      const CDT::V2d< double > &v1 = cdtVertices[ edge.v1() ];
      const CDT::V2d< double > &v2 = cdtVertices[ edge.v2() ];

      glm::dvec2 v1d( v1.x, v1.y );
      glm::dvec2 v2d( v2.x, v2.y );

      minBound = glm::min( minBound, glm::min( v1d, v2d ) );
      maxBound = glm::max( maxBound, glm::max( v1d, v2d ) );
    }

    double base   = std::min( minBound.x, minBound.y );
    double extent = std::max( maxBound.x, maxBound.y ) - base;

    auto svgScale = [=]( double value ) {

      return 50 + ( 1024.0 * ( ( value - base ) / extent ) ); 

    };

    svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
              << "width=\"1124\" height=\"1124\">\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"512\" style=\"stroke:rgb(255, 132, 0);stroke-width:2\" fill=\"none\"/>\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"256\" style=\"stroke:rgb(0, 0, 255);stroke-width:2\" fill=\"none\"/>\n";

    for ( const CDT::Edge &edge : cdtEdges )
    {
      const CDT::V2d< double > &v1 = cdtVertices[ edge.v1() ];
      const CDT::V2d< double > &v2 = cdtVertices[ edge.v2() ];

      svgFile << "<line x1=\"" << svgScale( v1.x ) << "\" y1=\"" << svgScale( v1.y )
              << "\" x2=\"" << svgScale( v2.x ) << "\" y2=\"" << svgScale( v2.y )
              << "\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
    }
    
    for ( uint32_t vertexIndex = 0; vertexIndex < cdtVertices.size(); ++vertexIndex )
    {
      const CDT::V2d< double > &vertex       = cdtVertices[ vertexIndex ];

      svgFile << "<circle cx=\"" << svgScale( vertex.x ) << "\" cy=\"" << svgScale( vertex.y )
              << "\" r=\"2\" fill=\"red\" />\n";
    }

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"2\" fill=\"green\" />\n";

    svgFile << "</svg>\n";

    svgFile.close();

#endif

    // No edges to process, return early.
    if ( cdtEdges.empty() )
    {
      return;
    }

    // Non-finite CDT input guard — same rationale as the dual-parameterization
    // guards below. Here the planar projection basis can come out NaN when
    // `GetBasisFromCoplanarPoints` returns nearly-coincident vectors and
    // `glm::normalize` divides a zero-length difference; every projected
    // coordinate is then NaN and CDT balloons hundreds of MiB before failing
    // (measured: supercap2/supercap3 pegged the 4 GiB wasm ceiling from 2 MB
    // sources through this path). The face is dropped either way.
    for ( const CDT::V2d< double > &v : cdtVertices )
    {
      if ( !std::isfinite( v.x ) || !std::isfinite( v.y ) )
      {
        throw std::runtime_error(
          "conway: non-finite planar projection; dropping face" );
      }
    }

    // TryResolve (here and at the other face-tessellation CDT sites in this
    // file and mesh_utils.h): split intersecting constraints at their
    // crossing instead of throwing — the throw dropped the whole face,
    // leaving a hole in the mesh. Real CAD input reaches this with
    // near-degenerate sliver boundaries (Onshape AmazingHand
    // Right_Hand.step); a locally-imperfect triangulation beats a missing
    // face. The csg_mesher sites deliberately keep NotAllowed: the boolean
    // pipeline relies on exact constraint topology.
    CDT::Triangulation< double > triangulation(
      CDT::VertexInsertionOrder::AsProvided,
      CDT::IntersectingConstraintEdges::TryResolve, 0);

    try
    {
      conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
      triangulation.insertVertices( cdtVertices );
      triangulation.insertEdges( cdtEdges );

      triangulation.eraseSuperTriangle();
    }
    catch (const CDT::IntersectingConstraintsError &e) {

      const CDT::V2d< double >& ev1 = cdtVertices[ e.e1().v1() ];
      const CDT::V2d< double >& ev2 = cdtVertices[ e.e1().v2() ];
      const CDT::V2d< double >& ev3 = cdtVertices[ e.e2().v1() ];
      const CDT::V2d< double >& ev4 = cdtVertices[ e.e2().v2() ];

      Logger::logError( "CDT Exception (plane svg: %d) ((%f,%f),(%f,%f)) -> ((%f,%f),(%f,%f)): %s", outputIndex, ev1.x, ev1.y, ev2.x, ev2.y, ev3.x, ev3.y, ev4.x, ev4.y, e.what() );
      return;
    }
    catch (const CDT::Error &e)
    {
      Logger::logError( "CDT Exception (plane svg: %d): %s", outputIndex, e.what() );
      return;
    }

    const CDT::TriangleVec &triangles = triangulation.triangles;

    for ( const CDT::Triangle &triangle : triangles )
    {
      auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

      if (
        cdtv1 == cdtv2 ||
        cdtv2 == cdtv3 ||
        cdtv3 == cdtv1 )
      {
        continue;
      }

      // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
      if (
        cdtv1 >= vertexRemapping.size() ||
        cdtv2 >= vertexRemapping.size() ||
        cdtv3 >= vertexRemapping.size() )
      {
        continue;
      }

      uint32_t v1 = vertexRemapping[ cdtv1 ];
      uint32_t v2 = vertexRemapping[ cdtv2 ];
      uint32_t v3 = vertexRemapping[ cdtv3 ];

      output.MakeTriangle( v1, v2, v3 );
    }
    
    std::unordered_map< uint64_t, uint32_t > loopEdgeSet;

    loopEdgeSet.reserve( boundaryMesh.loopEdges.size() );
    
    output.EnableConnectivity();

    for ( const LoopEdge& loopEdge : loopEdges )
    {
      uint64_t edgeID =
        edgeCompoundID(
          loopEdge.vertex1,
          loopEdge.vertex2 );

      loopEdgeSet.insert( { edgeID, loopEdge.loop } );
    }

    std::vector< bool > toDelete( output.triangles.size(), false );
    std::vector< bool > touched( output.triangles.size(), false );

    std::stack< uint32_t > triangleStack;

    uint64_t outerBoundaryIndex = EMPTY_INDEX;

    // Find triangle that has a free edge that is off the boundary loop.
    for ( uint32_t triangleIndex = 0, 
                   end = static_cast< uint32_t >( output.triangles.size() );
     triangleIndex < end;
     ++triangleIndex )
    {
      const TriangleEdges& triangle = output.triangle_edges[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( edge.otherTriangle( triangleIndex ) != EMPTY_INDEX )
        {
          continue;
        }

        uint64_t edgeID = edge.compoundID();

        auto foundLoopEdge = loopEdgeSet.find( edgeID );

        if ( foundLoopEdge == loopEdgeSet.end() )
        {
          toDelete[ triangleIndex ] = true; // Mark triangle for deletion.
          touched[ triangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( triangleIndex ); 
          break;
        } else if ( outerBoundaryIndex == EMPTY_INDEX ) {
          outerBoundaryIndex = foundLoopEdge->second;
        }
      }
    }

    bool flipMesh = false;

    if ( outerBoundaryIndex != EMPTY_INDEX ) {

      IfcCurve test;

      const IfcBound3D& outerBounds = boundaries[ outerBoundaryIndex ];

      for ( const glm::dvec3& pt : outerBounds.curve.points) {

        glm::dvec2 proj = parameterization( pt );

        test.Add2d( proj );
      }

      // if the outer bound is clockwise under the current projection (v12,v13,n),
      // we invert the projection
      if (!test.IsCCW()) {
        
        flipMesh = true;
      }
    }
    
    // Mark triangles that are connected to triangles outside the outer
    // loop deleted by traversing connectivity.
    while ( !triangleStack.empty() )
    {
      uint32_t triangleIndex = triangleStack.top();

      triangleStack.pop();

      const TriangleEdges& triangle = output.triangle_edges[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {
          continue; // Edge is part of the loop, skip it.
        }

        uint32_t otherTriangleIndex = edge.otherTriangle( triangleIndex );

        if ( otherTriangleIndex != EMPTY_INDEX )
        {
          if ( touched[ otherTriangleIndex ] )
          {
            continue; // Triangle already processed.
          }

          touched[ otherTriangleIndex ]  = true; // Mark triangle as touched.
          toDelete[ otherTriangleIndex ] = true; // Mark triangle for deletion.

          triangleStack.push( otherTriangleIndex ); // Add to stack for processing.
        }
      }
    }
   
    // Find triangle that has a free edge but is on the loop.
    for ( uint32_t triangleIndex = 0, 
                   end = static_cast< uint32_t >( output.triangles.size() );
     triangleIndex < end;
     ++triangleIndex )
    {
      const TriangleEdges& triangle = output.triangle_edges[ triangleIndex ];

      if ( touched[ triangleIndex ] )
      {
        continue; // Triangle already processed.
      }

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( edge.otherTriangle( triangleIndex ) != EMPTY_INDEX )
        {
          continue;
        }

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {        
          touched[ triangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( triangleIndex ); 
          break;
        }
      }
    }

    // Mark triangles that are connected to triangles outside the outer
    // loop deleted by traversing connectivity.
    while ( !triangleStack.empty() )
    {
      uint32_t triangleIndex = triangleStack.top();

      triangleStack.pop();

      const TriangleEdges& triangle = output.triangle_edges[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {
          continue; // Edge is part of the loop, skip it.
        }

        uint32_t otherTriangleIndex = edge.otherTriangle( triangleIndex );

        if ( otherTriangleIndex != EMPTY_INDEX )
        {
          if ( touched[ otherTriangleIndex ] )
          {
            continue; // Triangle already processed.
          }

          touched[ otherTriangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( otherTriangleIndex ); // Add to stack for processing.
        }
      }
    }    

    // Go back through the triangles in reverse order and delete those marked for deletion,
    // as deletion moves the last in vector to the current index.
    for ( uint32_t triangleIndex1 = static_cast< uint32_t >( output.triangles.size() ); 
          triangleIndex1 > 0; 
          --triangleIndex1 )
    {
      uint32_t triangleIndex = triangleIndex1 - 1;

      // If a triangle was marked for deletion, or wasn't touched in the walk,
      // its outside the outer bounds or in a hole.
      if ( toDelete[ triangleIndex ] || !touched[ triangleIndex ] )
      {
        output.DeleteTriangle( triangleIndex );
      }
    }

    if ( flipMesh )
    {
      output.ReverseFaces();
    }
  }  

  template <
      typename ParameterizationFunction,
      typename EquatorSideFunction,
      typename DiscardEdgeFunction>
  inline void tesselateHalfDualParameterization(
      WingedEdgeMesh<glm::dvec2> &outputMesh,
      std::vector<LoopEdge> &equatorBoundaryEdges,
      const LoopEdgeMesh &boundaryMesh,
      uint32_t expectedEquatorSide,
      ParameterizationFunction parameterization,
      EquatorSideFunction equatorSide,
      DiscardEdgeFunction discardEdge)
  {
    // auto&: WingedEdgeMesh::vertices is a std::pmr::vector (AFTP arena-backing).
    auto &parameterizationVertices = outputMesh.vertices;

    size_t vertexCount = boundaryMesh.normalFormVertices.size();

    parameterizationVertices.reserve( vertexCount );

    const std::vector<glm::dvec3> &normalFormVertices  = boundaryMesh.normalFormVertices;
    bool                           vertexOnEquatorSide = false;
    
    for ( const glm::dvec3 &vertex : normalFormVertices )
    {
      if ( equatorSide( vertex ) == expectedEquatorSide )
      {
        vertexOnEquatorSide = true;
        break;
      }
    }

    // Early exit if no vertices are on the expected side of the equator.
    if ( !vertexOnEquatorSide )
    {
      return;
    }

    for ( const glm::dvec3 &vertex : normalFormVertices )
    {
      parameterizationVertices.emplace_back( parameterization( vertex ) );
    }

    std::vector< CDT::Edge >      cdtEdges;
    std::vector<uint32_t>         vertexRemapping;
    std::vector<uint32_t>         reverseVertexRemapping(
      vertexCount,
      EMPTY_INDEX );
    std::vector< CDT::V2d< double > > cdtVertices;

    vertexRemapping.reserve( vertexCount );
    cdtVertices.reserve( vertexCount );

    bool edgeDiscarded = false;

    // Setup constraints, discard illegal edges and remap vertices for CDT.
    for ( const LoopEdge &loopEdge : boundaryMesh.loopEdges )
    {
      auto [ v1, v2, _ ] = loopEdge;

      const glm::dvec3 &normal1 = normalFormVertices[ v1 ];
      const glm::dvec3 &normal2 = normalFormVertices[ v2 ];
      const glm::dvec2 &param1  = parameterizationVertices[ v1 ];
      const glm::dvec2 &param2  = parameterizationVertices[ v2 ];

      if ( discardEdge( normal1, normal2, param1, param2 ) )
      {
        edgeDiscarded = true;
        continue;
      }

      uint32_t remappedV1 = reverseVertexRemapping[ v1 ];
      uint32_t remappedV2 = reverseVertexRemapping[ v2 ];

      if ( remappedV1 == EMPTY_INDEX )
      {
        remappedV1 = static_cast<uint32_t>( cdtVertices.size() );

        vertexRemapping.push_back( v1 );

        reverseVertexRemapping[ v1 ] = remappedV1;

        cdtVertices.emplace_back( param1.x, param1.y );
      }

      if ( remappedV2 == EMPTY_INDEX )
      {
        remappedV2 = static_cast< uint32_t >( cdtVertices.size() );

        vertexRemapping.push_back( v2 );

        reverseVertexRemapping[ v2 ] = remappedV2;

        cdtVertices.emplace_back( param2.x, param2.y );
      }

      cdtEdges.emplace_back( remappedV1, remappedV2 );
    }

    // atomic: this can run concurrently on the thread pool during staged
    // face finalization (same fix as tesselatePlane above).
    static std::atomic< uint32_t > svgIndex    = 0;
    uint32_t                       outputIndex = svgIndex++;

#if (OUTPUT_SVG_DEBUG == 1)

    std::ofstream svgFile( "cdt_" + std::to_string( expectedEquatorSide ) + "_" + std::to_string( outputIndex ) + ".svg" );

    auto svgScale = []( double value ) {

      return 50 + ( 1024.0 * ( value + 2.0 ) / 4.0 ); 

    };

    svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
              << "width=\"1124\" height=\"1124\">\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"512\" style=\"stroke:rgb(255, 132, 0);stroke-width:2\" fill=\"none\"/>\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"256\" style=\"stroke:rgb(0, 0, 255);stroke-width:2\" fill=\"none\"/>\n";

    for ( const CDT::Edge &edge : cdtEdges )
    {
      const CDT::V2d< double > &v1 = cdtVertices[ edge.v1() ];
      const CDT::V2d< double > &v2 = cdtVertices[ edge.v2() ];

      svgFile << "<line x1=\"" << svgScale( v1.x ) << "\" y1=\"" << svgScale( v1.y )
              << "\" x2=\"" << svgScale( v2.x ) << "\" y2=\"" << svgScale( v2.y )
              << "\" style=\"stroke:rgb(0,0,0);stroke-width:2\" />\n";
    }
    
    for ( uint32_t vertexIndex = 0; vertexIndex < cdtVertices.size(); ++vertexIndex )
    {
      const CDT::V2d< double > &vertex       = cdtVertices[ vertexIndex ];
      const glm::dvec3         &normalVertex = normalFormVertices[ vertexRemapping[ vertexIndex ] ];

      svgFile << "<circle cx=\"" << svgScale( vertex.x ) << "\" cy=\"" << svgScale( vertex.y )
              << "\" r=\"4\" fill=\"red\" /> <!-- " << normalVertex.x << ", " << normalVertex.y << ", " << normalVertex.z << " -->\n";
    }

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"4\" fill=\"green\" />\n";

    svgFile << "</svg>\n";

    svgFile.close();

#endif

    // No edges to process, return early.
    if ( cdtEdges.empty() )
    {
      return;
    }

    // Non-finite CDT input guard. A dual-parameterization singularity — a
    // full-loop toroidal/spherical face whose boundary maps onto the
    // projection pole, where `normalize(originalPlanar)` divides a
    // zero-length vector — yields NaN/inf parameter coordinates for every
    // vertex. CDT chokes on that input, ballooning its triangle store by
    // hundreds of MiB per face before it aborts, and the whole face is
    // dropped regardless (it commits no geometry) — on Arty_Z7, 8 such
    // fillet faces ratcheted the grow-only wasm heap to the 4 GiB ceiling.
    // Detect it up front and drop the face the same way, without the
    // transient. Byte-identical: CDT of non-finite input never produced
    // committed geometry. (Correctly tessellating pole-crossing faces is a
    // separate, digest-moving change.)
    for ( const CDT::V2d< double > &v : cdtVertices )
    {
      if ( !std::isfinite( v.x ) || !std::isfinite( v.y ) )
      {
        throw std::runtime_error(
          "conway: non-finite dual-parameterization (projection pole); dropping face" );
      }
    }

    CDT::Triangulation< double > triangulation(
      CDT::VertexInsertionOrder::Auto,
      CDT::IntersectingConstraintEdges::TryResolve, 0);

    try
    {
      conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
      triangulation.insertVertices( cdtVertices );
      triangulation.insertEdges( cdtEdges );

      if ( edgeDiscarded ) {
        triangulation.eraseSuperTriangle();
      } else {
        triangulation.eraseOuterTrianglesAndHoles();
      }
    }
    catch (const CDT::IntersectingConstraintsError &e) {

      const CDT::V2d< double >& ev1 = cdtVertices[ e.e1().v1() ];
      const CDT::V2d< double >& ev2 = cdtVertices[ e.e1().v2() ];
      const CDT::V2d< double >& ev3 = cdtVertices[ e.e2().v1() ];
      const CDT::V2d< double >& ev4 = cdtVertices[ e.e2().v2() ];

      Logger::logError( "CDT Exception (hemisphere: %d, svg: %d) ((%f,%f),(%f,%f)) -> ((%f,%f),(%f,%f)): %s", expectedEquatorSide, outputIndex, ev1.x, ev1.y, ev2.x, ev2.y, ev3.x, ev3.y, ev4.x, ev4.y, e.what() );
      return;
    }
    catch (const CDT::Error &e)
    {
      Logger::logError( "CDT Exception (hemisphere: %d, svg: %d): %s", expectedEquatorSide, outputIndex, e.what() );
      return;
    }

    const CDT::TriangleVec &triangles = triangulation.triangles;

    for ( const CDT::Triangle &triangle : triangles )
    {
      auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

      if (
        cdtv1 == cdtv2 ||
        cdtv2 == cdtv3 ||
        cdtv3 == cdtv1 )
      {
        continue;
      }

      // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
      if (
        cdtv1 >= vertexRemapping.size() ||
        cdtv2 >= vertexRemapping.size() ||
        cdtv3 >= vertexRemapping.size() )
      {
        continue;
      }

      uint32_t v1 = vertexRemapping[ cdtv1 ];
      uint32_t v2 = vertexRemapping[ cdtv2 ];
      uint32_t v3 = vertexRemapping[ cdtv3 ];

      const glm::dvec3 &vertex1 = boundaryMesh.normalFormVertices[ v1 ];
      const glm::dvec3 &vertex2 = boundaryMesh.normalFormVertices[ v2 ];
      const glm::dvec3 &vertex3 = boundaryMesh.normalFormVertices[ v3 ];

      if ( edgeDiscarded ) {
        bool v1equatorSide = equatorSide( vertex1 ) == expectedEquatorSide;
        bool v2equatorSide = equatorSide( vertex2 ) == expectedEquatorSide;
        bool v3equatorSide = equatorSide( vertex3 ) == expectedEquatorSide;
      
        // Add edges to the equator boundary if they are on the expected 
        // side of the equator.
        if ( v1equatorSide &&
            v2equatorSide &&
            !v3equatorSide )
        {
          equatorBoundaryEdges.emplace_back( v1, v2 );
        }

        if ( v2equatorSide &&
            v3equatorSide &&
            !v1equatorSide )
        {
          equatorBoundaryEdges.emplace_back( v2, v3 );
        }

        if ( v3equatorSide &&
            v1equatorSide &&
            !v1equatorSide)
        {
          equatorBoundaryEdges.emplace_back( v3, v1 );
        }

        if ( !v1equatorSide ||
            !v2equatorSide ||
            !v3equatorSide )
        {
          continue; // At least one vertex is on the wrong side of the equator, skip triangle.
        }
      }

      outputMesh.makeTriangle( v1, v2, v3 );
    }
  }

  template <
      typename NormalFormFunction,
      typename Parameterization1Function,
      typename Parameterization2Function,
      typename EquatorSideFunction,
      typename DiscardEdge1Function,
      typename DiscardEdge2Function>
  inline void tesselateDualParametrization(
      WingedEdgeMesh< glm::dvec3 > &output,
      const std::vector<IfcBound3D> &boundaries,
      NormalFormFunction normalForm,
      Parameterization1Function parameterization1,
      Parameterization2Function parameterization2,
      EquatorSideFunction equatorSide,
      DiscardEdge1Function discardEdge1,
      DiscardEdge2Function discardEdge2 )
  {

    LoopEdgeMesh boundaryMesh;

    VertexMap               &vertexMap          = boundaryMesh.vertexMap;
    auto &vertices           = output.vertices;  // pmr::vector (AFTP)
    std::vector<glm::dvec3> &normalFormVertices = boundaryMesh.normalFormVertices;
    std::vector< LoopEdge > &loopEdges          = boundaryMesh.loopEdges;

     // Build connected edge structure with unique vertices from boundary loops.
    for ( const IfcBound3D &boundary : boundaries )
    {
      const std::vector<glm::dvec3> &points =
          boundary.curve.points;

      bool     firstVertexInLoop   = true;
      uint32_t previousVertexIndex = 0;
      uint32_t firstIndexInLoop    = 0;

      for ( const glm::dvec3 &point : points )
      {
        glm::dvec3 normalFormPoint = normalForm( point );

        auto foundVertex = vertexMap.find( normalFormPoint );

        uint32_t vertexIndex = 0;

        if ( foundVertex == vertexMap.end() )
        {
          vertexIndex = static_cast< uint32_t >( vertices.size() );

          vertices.push_back( point );
          normalFormVertices.push_back( normalFormPoint );

          vertexMap[ normalFormPoint ] = vertexIndex;
        }
        else
        {
          vertexIndex = foundVertex->second;
        }

        if ( firstVertexInLoop )
        {
          previousVertexIndex = vertexIndex;
          firstIndexInLoop    = vertexIndex;
          firstVertexInLoop   = false;
          continue;
        }

        if ( previousVertexIndex == vertexIndex )
        {
          continue;
        }

        loopEdges.emplace_back( previousVertexIndex, vertexIndex );

        previousVertexIndex = vertexIndex;
      }

      if ( previousVertexIndex != firstIndexInLoop )
      {
        loopEdges.emplace_back( previousVertexIndex, firstIndexInLoop );
      }
    }

    WingedEdgeMesh<glm::dvec2> side0;
    WingedEdgeMesh<glm::dvec2> side1;

    std::vector<LoopEdge> side0EquatorEdges;
    std::vector<LoopEdge> side1EquatorEdges;

    tesselateHalfDualParameterization(
        side0,
        side0EquatorEdges,
        boundaryMesh,
        0,
        parameterization1,
        equatorSide,
        discardEdge1);

    if ( !side0EquatorEdges.empty() || side0.triangles.empty() )
    {
      tesselateHalfDualParameterization(
          side1,
          side1EquatorEdges,
          boundaryMesh,
          1,
          parameterization2,
          equatorSide,
          discardEdge2);
    }

    std::vector< CDT::Edge > cdtEdges;

    size_t vertexCount = boundaryMesh.normalFormVertices.size();

    std::vector< uint32_t >           vertexRemapping;
    std::vector< uint32_t >           reverseVertexRemapping( vertexCount, EMPTY_INDEX );
    std::vector< CDT::V2d< double > > cdtVertices;

    vertexRemapping.reserve( vertexCount );
    cdtVertices.reserve( vertexCount );

    auto &parameterization1Vertices = side0.vertices;  // pmr::vector (AFTP)

    for ( const Triangle &triangle : side0.triangles )
    {
      auto [ v1, v2, v3 ] = triangle.vertices;

      if ( v1 == v2 || v2 == v3 || v3 == v1 )
      {
        continue;
      }

      output.makeTriangle( v1, v2, v3 );
    }

    for ( const Triangle &triangle : side1.triangles )
    {
      auto [ v1, v2, v3 ] = triangle.vertices;

      if ( v1 == v2 || v2 == v3 || v3 == v1 )
      {
        continue;
      }

      output.makeTriangle( v1, v2, v3 );
    }

    if ( !side0EquatorEdges.empty() ||
         !side1EquatorEdges.empty() )
    {
      // Setup constraints for both sides for the equator seem, and remap vertices for CDT.
      for ( const LoopEdge &loopEdge : side0EquatorEdges )
      {
        auto [ v1, v2, _ ] = loopEdge;

        const glm::dvec2 &param1 = parameterization1Vertices[ v1 ];
        const glm::dvec2 &param2 = parameterization1Vertices[ v2 ];

        uint32_t remappedV1 = reverseVertexRemapping[ v1 ];
        uint32_t remappedV2 = reverseVertexRemapping[ v2 ];

        if ( remappedV1 == EMPTY_INDEX )
        {
          remappedV1 = static_cast<uint32_t>( cdtVertices.size() );

          vertexRemapping.push_back(v1);

          reverseVertexRemapping[ v1 ] = remappedV1;

          cdtVertices.emplace_back( param1.x, param1.y );
        }

        if ( remappedV2 == EMPTY_INDEX )
        {
          remappedV2 = static_cast<uint32_t>( cdtVertices.size() );

          vertexRemapping.push_back( v2 );

          reverseVertexRemapping[ v2 ] = remappedV2;

          cdtVertices.emplace_back( param2.x, param2.y );
        }

        cdtEdges.emplace_back( remappedV1, remappedV2 );
      }

      // Note we use the parameterization from the first side for the second side,
      // because they all need to be in a single frame of reference for the CDT.
      for ( const LoopEdge &loopEdge : side1EquatorEdges )
      {
        auto [ v1, v2, _ ] = loopEdge;

        const glm::dvec2 &param1 = parameterization1Vertices[ v1 ];
        const glm::dvec2 &param2 = parameterization1Vertices[ v2 ];

        uint32_t remappedV1 = reverseVertexRemapping[ v1 ];
        uint32_t remappedV2 = reverseVertexRemapping[ v2 ];

        if ( remappedV1 == EMPTY_INDEX )
        {
          remappedV1 = static_cast<uint32_t>( cdtVertices.size() );

          vertexRemapping.push_back( v1 );

          reverseVertexRemapping[ v1 ] = remappedV1;

          cdtVertices.emplace_back(param1.x, param1.y);
        }

        if ( remappedV2 == EMPTY_INDEX )
        {
          remappedV2 = static_cast<uint32_t>(cdtVertices.size());

          vertexRemapping.push_back( v2 );

          reverseVertexRemapping[ v2 ] = remappedV2;

          cdtVertices.emplace_back(param2.x, param2.y);
        }

        cdtEdges.emplace_back(remappedV1, remappedV2);
      }

      for ( const LoopEdge &loopEdge : loopEdges )
      {
        auto [v1, v2, _] = loopEdge;

        const glm::dvec3 &normal1 = normalFormVertices[ v1 ];
        const glm::dvec3 &normal2 = normalFormVertices[ v2 ];

        const glm::dvec2 &param1 = parameterization1Vertices[ v1 ];
        const glm::dvec2 &param2 = parameterization1Vertices[ v2 ];

        // We only want loop edges that straddle the equator.
        if ( equatorSide( normal1 ) == equatorSide( normal2 ) )
        {
          continue;
        }

        uint32_t remappedV1 = reverseVertexRemapping[ v1 ];
        uint32_t remappedV2 = reverseVertexRemapping[ v2 ];

        if ( remappedV1 == EMPTY_INDEX )
        {
          remappedV1 = static_cast<uint32_t>( cdtVertices.size() );

          vertexRemapping.push_back( v1 );

          reverseVertexRemapping[ v1 ] = remappedV1;

          cdtVertices.emplace_back( param1.x, param1.y );
        }

        if ( remappedV2 == EMPTY_INDEX )
        {
          remappedV2 = static_cast<uint32_t>( cdtVertices.size() );

          vertexRemapping.push_back( v2 );

          reverseVertexRemapping[ v2 ] = remappedV2;

          cdtVertices.emplace_back( param2.x, param2.y );
        }

        cdtEdges.emplace_back( remappedV1, remappedV2 );
      }
    }

    if ( !cdtEdges.empty() )
    {
      // Same non-finite input guard as tesselateHalfDualParameterization
      // above — the equator-stitch CDT reads the side-0 parameterization,
      // which is NaN/inf when the boundary sits on the projection pole.
      for ( const CDT::V2d< double > &v : cdtVertices )
      {
        if ( !std::isfinite( v.x ) || !std::isfinite( v.y ) )
        {
          throw std::runtime_error(
            "conway: non-finite dual-parameterization (projection pole); dropping face" );
        }
      }

      CDT::Triangulation< double > triangulation(
        CDT::VertexInsertionOrder::AsProvided,
        CDT::IntersectingConstraintEdges::TryResolve,
        0);

      try
      {
        conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
        triangulation.insertVertices( cdtVertices );
        triangulation.insertEdges( cdtEdges );
        triangulation.eraseOuterTriangles();
      }
      catch ( const CDT::Error &e )
      {
        Logger::logError( "CDT Exception (hole fill): %s", e.what() );
        return;
      }

      const CDT::TriangleVec &triangles = triangulation.triangles;

      for (const CDT::Triangle &triangle : triangles)
      {
        auto [ v1, v2, v3 ] = triangle.vertices;

        if ( v1 == v2 || v2 == v3 || v3 == v1 )
        {
          continue;
        }

        // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
        if (
          v1 >= reverseVertexRemapping.size() ||
          v2 >= reverseVertexRemapping.size() ||
          v3 >= reverseVertexRemapping.size() )
        {
          continue;
        }

        uint32_t outputV1 = reverseVertexRemapping[ v1 ];
        uint32_t outputV2 = reverseVertexRemapping[ v2 ];
        uint32_t outputV3 = reverseVertexRemapping[ v3 ];

        output.makeTriangle( outputV1, outputV2, outputV3 );
      }
    }

    std::unordered_set< uint64_t > loopEdgeSet;

    loopEdgeSet.reserve( boundaryMesh.loopEdges.size() );

    for ( const LoopEdge& loopEdge : loopEdges )
    {
      uint64_t edgeID =
        edgeCompoundID(
          loopEdge.vertex1,
          loopEdge.vertex2 );

      loopEdgeSet.insert( edgeID );
    }

    std::vector< bool > toDelete( output.triangles.size(), false );
    std::vector< bool > touched( output.triangles.size(), false );

    std::stack< uint32_t > triangleStack;

    // Find triangle that has a free edge that is off the boundary loop.
    for ( uint32_t triangleIndex = 0, 
                   end = static_cast< uint32_t >( output.triangles.size() );
     triangleIndex < end;
     ++triangleIndex )
    {
      const ConnectedTriangle& triangle = output.triangles[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( edge.otherTriangle( triangleIndex ) != EMPTY_INDEX )
        {
          continue;
        }

        if ( loopEdgeSet.count( edge.compoundID() ) == 0 )
        {        
          toDelete[ triangleIndex ] = true; // Mark triangle for deletion.
          touched[ triangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( triangleIndex ); 
          break;
        }
      }
    }

    // Mark triangles that are connected to triangles outside the outer
    // loop deleted by traversing connectivity.
    while ( !triangleStack.empty() )
    {
      uint32_t triangleIndex = triangleStack.top();

      triangleStack.pop();

      const ConnectedTriangle& triangle = output.triangles[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {
          continue; // Edge is part of the loop, skip it.
        }

        uint32_t otherTriangleIndex = edge.otherTriangle( triangleIndex );

        if ( otherTriangleIndex != EMPTY_INDEX )
        {
          if ( touched[ otherTriangleIndex ] )
          {
            continue; // Triangle already processed.
          }

          touched[ otherTriangleIndex ]  = true; // Mark triangle as touched.
          toDelete[ otherTriangleIndex ] = true; // Mark triangle for deletion.

          triangleStack.push( otherTriangleIndex ); // Add to stack for processing.
        }
      }
    }
   
    // Find triangle that has a free edge but is on the loop.
    for ( uint32_t triangleIndex = 0, 
                   end = static_cast< uint32_t >( output.triangles.size() );
     triangleIndex < end;
     ++triangleIndex )
    {
      const ConnectedTriangle& triangle = output.triangles[ triangleIndex ];

      if ( touched[ triangleIndex ] )
      {
        continue; // Triangle already processed.
      }

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( edge.otherTriangle( triangleIndex ) != EMPTY_INDEX )
        {
          continue;
        }

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {        
          touched[ triangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( triangleIndex ); 
          break;
        }
      }
    }

    // Mark triangles that are connected to triangles outside the outer
    // loop deleted by traversing connectivity.
    while ( !triangleStack.empty() )
    {
      uint32_t triangleIndex = triangleStack.top();

      triangleStack.pop();

      const ConnectedTriangle& triangle = output.triangles[ triangleIndex ];

      for ( uint32_t edgeInTriangle = 0; edgeInTriangle < 3; ++edgeInTriangle )
      {
        const Edge& edge = output.edges[ triangle.edges[ edgeInTriangle ] ];

        if ( loopEdgeSet.count( edge.compoundID() ) > 0 )
        {
          continue; // Edge is part of the loop, skip it.
        }

        uint32_t otherTriangleIndex = edge.otherTriangle( triangleIndex );

        if ( otherTriangleIndex != EMPTY_INDEX )
        {
          if ( touched[ otherTriangleIndex ] )
          {
            continue; // Triangle already processed.
          }

          touched[ otherTriangleIndex ]  = true; // Mark triangle as touched.

          triangleStack.push( otherTriangleIndex ); // Add to stack for processing.
        }
      }
    }    

    // Go back through the triangles in reverse order and delete those marked for deletion,
    // as deletion moves the last in vector to the current index.
    for ( uint32_t triangleIndex1 = static_cast< uint32_t >( output.triangles.size() ); 
          triangleIndex1 > 0; 
          --triangleIndex1 )
    {
      uint32_t triangleIndex = triangleIndex1 - 1;

      // If a triangle was marked for deletion, or wasn't touched in the walk,
      // its outside the outer bounds or in a hole.
      if ( toDelete[ triangleIndex ] || !touched[ triangleIndex ] )
      {
        output.deleteTriangle( triangleIndex );
      }
    }

    // printf( "Finished sphere, %d triangles, %d edges\n",
    //       static_cast< int >( output.triangles.size() ),
    //       static_cast< int >( output.edges.size() ) );
  }

}