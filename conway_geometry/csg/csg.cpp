#include "csg.h"
#include "multi_mesh_vertex_index.h"
#include "structures/thread_pool.h"


/** Clean a mesh by getting its solid skin */
void conway::geometry::CSG::clean( Geometry& a ) {

  index( a );

  Geometry output;

  charter_.process(
    a,
    contacts,
    face_contact_map[ 0 ],
    boundary_set[ 0 ],
    output
  );

  a = std::move( output );

}

/** Index a mesh for self intersections */
void conway::geometry::CSG::index( Geometry& a, double tolerance ) {

  MultiMeshVertexIndex< 2 > vertices( a.vertices );
  
  reset();

  // Note, potentially we should move this out and error on lack of bvh/dipoles,
  // and force this onto the user to keep a and b const.
  a.MakeBVH();

  AABBTree& bvhA = *a.bvh;

  if ( !bvhA.hasDipoles() ) {

    bvhA.dipoles( a );
  }

  std::vector< std::pair< uint32_t, uint32_t > >& candidatePairs = candidatePairs_;

  boundary_set[ 0 ].clear();
  boundary_set[ 0 ].resize( a.triangles.size(), false );

  candidatePairs.clear();

  bvhA.intersect( [&]( uint32_t left, uint32_t right ) {

    const Triangle& aTriangle = a.triangles[ left ];
    const Triangle& bTriangle = a.triangles[ right ];

    const TriangleEdges& aTriangleEdges = a.triangle_edges[ left ];
    const TriangleEdges& bTriangleEdges = a.triangle_edges[ right ];

    bool shareEdge = false;

    std::span< const uint32_t, 3 > edges0{ aTriangleEdges.edges };
    std::span< const uint32_t, 3 > edges1{ bTriangleEdges.edges };

    for (uint32_t edge0 : edges0) {

      for (uint32_t edge1 : edges1) {

        if ( edge0 == edge1 ) {

          shareEdge = true;
          break;
        }
      }
    }

    glm::dvec3 localVertices[ 2 ][ 3 ] = {};

    vertices.extract( 0, aTriangle.vertices, localVertices[ 0 ] );
    vertices.extract( 0, bTriangle.vertices, localVertices[ 1 ] );

    int32_t vertexFaceSigns[ 2 ][ 3 ] = {};

    int32_t(&aFaceSigns)[ 3 ] = vertexFaceSigns[ 0 ];
    int32_t(&bFaceSigns)[ 3 ] = vertexFaceSigns[ 1 ];

    FaceFace aFaceFace =
      vertex_face_signs( localVertices[ 0 ], localVertices[ 1 ], aFaceSigns, tolerance );

    if ( aFaceFace == FaceFace::NONE || shareEdge ) {

      return;
    }

    FaceFace bFaceFace =
      vertex_face_signs( localVertices[ 1 ], localVertices[ 0 ], bFaceSigns, tolerance );

    if ( bFaceFace == FaceFace::NONE ) {

      return;
    }

    candidatePairs.emplace_back(left, right);
  });

  contacts.clear();
  contacts.resize( candidatePairs.size() );

// #if defined(__EMSCRIPTEN__)
//   std::transform(
//     candidatePairs.begin(),
//     candidatePairs.end(),
//     contacts.begin(),
//     [&](const std::pair< uint32_t, uint32_t >& candidatePair) {
// #else
  // std::transform(
  //   std::execution::par,
  //   candidatePairs.begin(),
  //   candidatePairs.end(),
  //   contacts.begin(),
  //   [&](const std::pair< uint32_t, uint32_t >& candidatePair) {
//#endif

  ThreadPool::instance().parallel_for( 0, candidatePairs.size(), [&]( size_t where ) {

      const std::pair< uint32_t, uint32_t >& candidatePair = candidatePairs_[ where ];
  
      uint32_t triangleInMeshIndices[2] = { candidatePair.first, candidatePair.second };

      const Triangle& aTriangle = a.triangles[ candidatePair.first ];
      const Triangle& bTriangle = a.triangles[ candidatePair.second ];

      const TriangleEdges& aTriangleEdges = a.triangle_edges[ candidatePair.first ];
      const TriangleEdges& bTriangleEdges = a.triangle_edges[ candidatePair.second ];

      const Triangle* triangles[2] = { &aTriangle, &bTriangle };

      bool sharedEdge = false;

      std::span< const uint32_t, 3 > edges0 { aTriangleEdges.edges };
      std::span< const uint32_t, 3 > edges1 { bTriangleEdges.edges };

      for ( uint32_t edge0 : edges0 ) {

        for ( uint32_t edge1 : edges1 ) {

          if ( edge0 == edge1 ) {
            
            sharedEdge = true;
            break;
          }
        }
      }

      contacts[ where ] = find_intersections( vertices, triangleInMeshIndices, triangles, true, sharedEdge, tolerance );
    });

  candidatePairs_.clear();

  for ( const TriangleTriangleContactPair& contactPair : contacts ) {

    for (size_t where = 0; where < 2; ++where) {

      const TriangleContacts& triangleContact = contactPair.triangles[ where ];

      if (
        !triangleContact.empty()) {

        boundary_set[ 0 ][ triangleContact.this_triangle_index ] = true;
      }
    }
  }

  face_contact_map[ 0 ].construct(
    contacts,
    static_cast<uint32_t>(a.triangles.size()),
    [](const TriangleTriangleContactPair& contact) {
      return
        std::make_pair(
          contact.triangles[ 0 ].this_triangle_index,
          contact.triangles[ 1 ].this_triangle_index); });

}

/**

 * Creates an index of the two meshes interactions so queries can be run.
 * 
 * This is called by run internally, and does not need to be called before calling run.
 */
void conway::geometry::CSG::index( Geometry& a, Geometry& b, double tolerance ) {

  MultiMeshVertexIndex< 2 > vertices  = multi_mesh_vertex_index( a, b );

  reset();

  // Note, potentially we should move this out and error on lack of bvh/dipoles,
  // and force this onto the user to keep a and b const.
  AABBTree& bvhA = a.MakeBVH();
  AABBTree& bvhB = b.MakeBVH();

  if ( !bvhA.hasDipoles() ) {

    bvhA.dipoles( a );
  }

  if ( !bvhB.hasDipoles() ) {

    bvhB.dipoles( b );
  }

  std::vector< std::pair< uint32_t, uint32_t > >& candidatePairs = candidatePairs_;

  boundary_set[ 0 ].resize( a.triangles.size(), false );
  boundary_set[ 1 ].resize( b.triangles.size(), false );

  candidatePairs.clear();

  bvhA.intersect( bvhB, [&]( uint32_t left, uint32_t right ) {

    candidatePairs.emplace_back( left, right );
  } );

  contacts.clear();
  contacts.resize( candidatePairs.size() );

  // std::transform(
  //   std::execution::par,
  //   candidatePairs.begin(),
  //   candidatePairs.end(),
  //   contacts.begin(),
  //   [&]( const std::pair< uint32_t, uint32_t >& candidatePair) {
  ThreadPool::instance().parallel_for( 0, candidatePairs.size(), [&]( size_t where ) {

    const std::pair< uint32_t, uint32_t >& candidatePair = candidatePairs[ where ];

    uint32_t triangleInMeshIndices[ 2 ] = { candidatePair.first, candidatePair.second };

    const Triangle& aTriangle = a.triangles[ candidatePair.first ];
    const Triangle& bTriangle = b.triangles[ candidatePair.second ];

    const Triangle* triangles[ 2 ] = { &aTriangle, &bTriangle };

    contacts[ where ] = find_intersections( vertices, triangleInMeshIndices, triangles, false, false, tolerance );
  });

  candidatePairs_.clear();

  for ( const TriangleTriangleContactPair& contactPair : contacts ) {

    for ( size_t where = 0; where < 2; ++where ) {

      const TriangleContacts& triangleContact = contactPair.triangles[ where ];

      if ( 
        !triangleContact.empty() && 
        !boundary_set[ where ][ triangleContact.this_triangle_index ] ) {

        boundary_set[ where ][ triangleContact.this_triangle_index ] = true;
      }
    }
  }

  face_contact_map[ 0 ].construct(
    contacts,
    static_cast< uint32_t >( a.triangles.size() ),
    []( const TriangleTriangleContactPair& contact ) { return contact.triangles[ 0 ].this_triangle_index; });
  
  face_contact_map[ 1 ].construct(
    contacts,
    static_cast< uint32_t >( b.triangles.size() ),
    []( const TriangleTriangleContactPair& contact ) { return contact.triangles[ 1 ].this_triangle_index; } );
}

void conway::geometry::CSG::run( Operation operation, Geometry& a, Geometry& b, Geometry& output, double tolerance ) {

  index( a, b, tolerance );

  bool aOutside     = false;
  bool bOutside     = false;
  bool flipBWinding = false;

  switch ( operation ) {
  case Operation::UNION:

    aOutside = true;
    bOutside = true;
    break;

  case Operation::INTERSECTION:

    // default case.
    break;

  case Operation::SUBTRACTION:

    aOutside     = true;
    flipBWinding = true;
    break;
  }

  charter_.process(
    a,
    b,
    contacts,
    face_contact_map[ 0 ],
    face_contact_map[ 1 ],
    boundary_set,
    aOutside,
    bOutside,
    flipBWinding,
    output
  );

  output.MarkedCleanedup();
}


void conway::geometry::CSG::reset() {
  
  for ( uint32_t where = 0; where < 2; ++where ) {
    
    contacts.clear();
    face_contact_map[ where  ].reset();
    boundary_set[ where ].clear();
  }

  candidatePairs_.clear();
  charter_.reset();
}