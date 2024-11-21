#pragma once

#include "structures/union_find.h"
#include "csg_mesher.h"
#include "csg_utils.h"
#include "triangle_contacts.h"
#include "triangle_triangle_contact_pair.h"
#include "structures/winged_edge.h"
#include "structures/box.h"
#include "structures/aabb_tree.h"
#include "structures/prefix_sum_map.h"
#include <unordered_set>
#include <bit>

namespace conway::geometry {

  struct CSG {

    enum class Operation { 
      UNION        = 0,
      INTERSECTION = 1,
      SUBTRACTION  = 2
    };

    std::vector< TriangleTriangleContactPair > contacts;
    PrefixSumMap                               face_contact_map[ 2 ];


    std::vector< bool >                        boundary_set[ 2 ];

    /** Clean a mesh by getting its solid skin */
    void clean( WingedEdgeMesh< glm::dvec3 >& a );

    /** Index a mesh for self intersections */
    void index( WingedEdgeMesh< glm::dvec3 >& a, double tolerance = 0 );

    /**
     * Creates an index of the two meshes interactions so queries can be run.
     * 
     * This is called by run internally, and does not need to be called before calling run.
     */
    void index( WingedEdgeMesh< glm::dvec3 >& a, WingedEdgeMesh< glm::dvec3 >& b, double tolerance = 0 );

    void run( Operation operation, WingedEdgeMesh< glm::dvec3 >& a, WingedEdgeMesh< glm::dvec3 >& b, WingedEdgeMesh< glm::dvec3 >& output, double tolerance = 0 );

    void reset();

    std::string dumpNovelVertices( const std::string& preamble = "" ) const {

      return charter_.dumpNovelVertices( preamble );
    }

#if !defined( __EMSCRIPTEN__ )
    std::string dumpConstraints( const std::string& preamble = "" ) const {

      return charter_.dumpConstraints( preamble );
    }
#endif

private:

    CSGMesher charter_;

    std::vector< std::pair< uint32_t, uint32_t > > candidatePairs_;

  };

}