#pragma once 

#include <glm/glm.hpp>
#include "triangle_contacts.h"
#include "triangle_triangle_contact_pair.h"
#include "structures/prefix_sum_map.h"
#include "structures/union_find.h"
#include "structures/winged_edge.h"
#include "structures/hash_functions.h"
#include "contact_pair.h"
#include "csg_utils.h"
#include <unordered_map>
#include <optional>
#include <execution>

#if defined (_MSC_VER)

#pragma warning( push )
#pragma warning( disable : 26812 )

#endif

#include "CDT.h"


#if defined (_MSC_VER)

#pragma warning( pop )

#endif

namespace conway::geometry {

  class CSGMesher {
  public:  

    void process(
      WingedEdgeDV3& a, // these are not const because we lazily generate dipoles.
      WingedEdgeDV3& b,
      const std::vector< TriangleTriangleContactPair >& contacts,
      const PrefixSumMap& aContactMap,
      const PrefixSumMap& bContactMap,
      const std::vector< bool > (&boundarySet)[ 2 ],
      bool aOutside,
      bool bOutside,
      bool flipBWinding,
      WingedEdgeDV3& output );

    void process(
      WingedEdgeDV3& a, // these are not const because we lazily generate dipoles.
      const std::vector< TriangleTriangleContactPair >& contacts,
      const PrefixSumMap& aContactMap,
      const std::vector< bool >& boundarySet,
      WingedEdgeDV3& output );

    void reset();

    std::string dumpNovelVertices( const std::string& preamble = "" ) const;

#if !defined( __EMSCRIPTEN__ )
    std::string dumpConstraints( const std::string& preamble = "" ) const;
#endif

  private:

    std::string dumpEdgeAndVertsToSVG() const;

    void walkAndInsertNonBoundary(
      bool outside,
      const std::vector< bool >& boundarySet,
      const WingedEdgeMesh< glm::dvec3 >& mesh,
      const WingedEdgeMesh< glm::dvec3 >& otherMesh,
      bool flippedWinding,
      uint32_t vertexOffset );

    void addEdges( std::span< const ContactPair > contactPairs, const FixedStack< uint32_t, 6 >& additionalVertices );

    template < size_t N >
    void triangulate(
      const WingedEdgeMesh< glm::dvec3 >& mesh,
      const MultiMeshVertexIndex< N >& vertices,
      uint32_t triangleInMeshIndex,
      bool flippedWinding,
      uint32_t outputStreamIndex );

    uint32_t insertLocalVertex( uint32_t inputVertex );

    uint32_t insertLocalVertexOnEdge( uint32_t inputVertex, uint32_t edgeInTriangle );

    UnionFind< uint32_t > unifiedVertices_;
    UnionFind< uint32_t > unifiedPlanes_;

    std::vector< glm::dvec3 > novelVertices_;

    std::unordered_map< std::pair< uint32_t, uint32_t >, uint32_t > edgeEdgeVertices_; 
    std::unordered_map< std::pair< uint32_t, uint32_t >, uint32_t > faceEdgeVertices_[ 2 ];
    std::unordered_map< glm::dvec3, uint32_t > duplicateVertexMap_;

    std::vector< uint32_t > onEdgeVertices_[ 3 ];

    std::vector< uint8_t > outside_[ 2 ];

    std::vector< Triangle > outputTriangleStream_;
    std::vector< std::pair< Triangle, uint32_t > > initialChartTriangles_[ 2 ];
    std::vector< CDT::Edge > edges_;
    std::unordered_map< uint32_t, uint32_t > localVertexMap_;
    std::vector< uint32_t >  localVertices_;
    std::vector< uint8_t > localVertexEdgeFlags_;
    std::vector< CDT::V2d< double > > local2DVertices_;
    std::vector< bool > walked_;
    std::vector< uint32_t > triangleStack_;
    std::vector< bool > vertexUsed_;
    std::vector< uint32_t > globalVertexMap_;

#if !defined( __EMSCRIPTEN__ )
    std::vector< std::pair< glm::vec3, glm::dvec3 > > contraintEdge_;
#endif
  };
}