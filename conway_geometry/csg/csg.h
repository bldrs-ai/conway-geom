#pragma once

#include "structures/union_find.h"
#include "csg_mesher.h"
#include "csg_utils.h"
#include "triangle_contacts.h"
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
      DIFFERENCE   = 2
    };

    struct TriangleInteractions {

      int8_t vertex_face_signs[ 3 ];

      uint32_t triangle_index;
    };


    struct TriangleTriangleClassification {

      TriangleInteractions triangles[ 2 ];

      std::pair< uint32_t, uint32_t > edge_edge[ 6 ];

      FaceFace face_to_face = FaceFace::NONE;
    };


    std::vector< TriangleContacts >               contacts[ 2 ];
    
    PrefixSumMap                                  face_contact_map[ 2 ];

    std::vector< bool >                           boundary_set[ 2 ];
    std::vector< uint32_t >                       boundary[ 2 ];

    std::vector< TriangleTriangleClassification > candidate_pairs;

    void add( const TriangleContacts(&triangleContacts)[ 2 ] ) {

      for ( size_t where = 0; where < 2; ++where ) {

        const TriangleContacts& triangleContact = triangleContacts[ where ];

        if ( !triangleContact.empty() ) {

          contacts[ where ].push_back( triangleContact );

          if ( !boundary_set[ where ][ triangleContact.this_triangle_index ] ) {

            boundary_set[ where ][ triangleContact.this_triangle_index ] = true;
            boundary[ where ].push_back( triangleContact.this_triangle_index );
          }
        }
      }
    }

    FaceFace vertexFaceSigns(
        const WingedEdgeMesh< glm::dvec3 >& a,
        const WingedEdgeMesh< glm::dvec3 >& b,
        const uint32_t (&triangleIndices)[ 3 ],
        const uint32_t (&opposingTriangleIndices)[ 3 ],
        TriangleInteractions& predicates,
        double tolerance
         ) {

      int32_t signs[3] = {};

      for (uint32_t vertInTriangle = 0; vertInTriangle < 3; ++vertInTriangle) {

        // if a vertex has been unified, orient3D will *always* be 0.
        // and we should propagate that, even if orient3D should 
        // be *exact* in theory and we are paranoid about input
        // ordering to make it as deterministic as possible.
        int32_t sign =
          orient3D(
            a,
            b,
            opposingTriangleIndices[ 0 ],
            opposingTriangleIndices[ 1 ],
            opposingTriangleIndices[ 2 ],
            triangleIndices[ vertInTriangle ],
            tolerance );

        signs[ sign + 1 ] += 1;

        // Classify this vertex relative a face.
        predicates.vertex_face_signs[ vertInTriangle ] = static_cast< int8_t >( sign );
      }

      if ( signs[ 1 ] == 3 ) {

        return FaceFace::COLINEAR;
      }

      if ( signs[ 0 ] == 3 || signs[ 2 ] == 3 ) {

        return FaceFace::NONE;
      }

      return FaceFace::STRADDLES_OR_TOUCHES;
    }

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

private:

    CSGMesher charter_;

  };

}