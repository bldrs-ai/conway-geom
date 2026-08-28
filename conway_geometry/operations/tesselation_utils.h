#pragma once

#include <cmath>
#include <limits>
#include <glm/glm.hpp>

#include "structures/winged_edge.h"
#include "structures/scratch_arena.h"
#include "structures/alloc_telemetry.h"
#include <cstring>
#include <memory_resource>
#include <queue>
#include <unordered_set>
#include "representation/Geometry.h"
#include "representation/IfcGeometryReps.h"
#include "operations/math_utils.h"

#if defined (_MSC_VER)

#pragma warning( push )
#pragma warning( disable : 26812 )

#endif

#include "CDT.h"


#if defined (_MSC_VER)

#pragma warning( pop )

#endif

namespace conway::geometry {

  /**
   * A UV parameterized vertex on 2 parameter surface.
   */
  struct ParameterVertex {

    glm::dvec3 point;
    glm::dvec2 uv;

  };

  /** Position accessors so bound computations span both mesh vertex types. */
  inline const glm::dvec3& refinementPoint( const glm::dvec3& vertex ) {
    return vertex;
  }

  /** Position of a UV parameterized vertex. */
  inline const glm::dvec3& refinementPoint( const ParameterVertex& vertex ) {
    return vertex.point;
  }

  /**
   * How fine a deflection target the enclosing OBJECT's extent permits, as a
   * fraction of it. Read it as a display bound: at a camera that frames the
   * object on a ~1000px viewport, one pixel is 1e-3 of it, so 1e-5 of it is
   * a hundredth of a pixel — a hundred-fold zoom of headroom past the point
   * where refinement can still be seen.
   *
   * It exists because the per-face convention below, right for a part, is
   * wrong for a mosaic. `Arty_Z7.stp`'s silkscreen is 1,189 extruded-glyph
   * solids whose stroke sidewalls are 10,224 b-spline faces with a MEDIAN
   * DIAGONAL OF 0.126mm; 0.1% of such a face is a target of 0.126um, roughly
   * a thousandth of a pixel at any zoom a user reaches, and chasing it costs
   * 96% of that model's geometry payload and 89% of its geometry time
   * (bldrs-ai/conway#564). The ten thousand tiles are one visual object —
   * the printed legend — and the object, not the tile, is what sets how
   * finely it is worth resolving. That object is the defining
   * representation: `Arty_Z7_Top_Silk` is one shape representation holding
   * all 654 glyph solids, 139.03mm across against a 1.20mm median glyph.
   *
   * The floor bites only a face smaller than
   * OBJECT_DEFLECTION_FLOOR_FACTOR / RELATIVE_DEFLECTION_FACTOR = 1% of its
   * representation, and it coarsens such a face by exactly the ratio by
   * which it falls short of that 1%. Every face at or above it — which is
   * every face a mechanical part is mostly made of — keeps the target it has
   * today. That containment is the whole point: a globally coarser
   * RELATIVE_DEFLECTION_FACTOR buys the same time on Arty_Z7 but facets
   * `create-a-tube` across 29.8% of its pixels, which is precisely what
   * the 0.1% factor exists to prevent.
   */
  constexpr double OBJECT_DEFLECTION_FLOOR_FACTOR = 1e-5;

  /**
   * A mesh position, keyed for exact lookup.
   *
   * Exact, not tolerant: this is used to answer "does the mesh already have
   * this point, bitwise", which is a question about arithmetic rather than
   * about proximity. See PositionSet.
   */
  struct PositionKey {

    double x;
    double y;
    double z;

    explicit PositionKey( const glm::dvec3& point ) :
      // Normalising -0.0 keeps the two spellings of zero one key, so a
      // position cannot hide from the set behind a sign bit. Every other
      // value is stored as it arrived.
      x( point.x == 0.0 ? 0.0 : point.x ),
      y( point.y == 0.0 ? 0.0 : point.y ),
      z( point.z == 0.0 ? 0.0 : point.z ) {}

    bool operator==( const PositionKey& other ) const {

      return x == other.x && y == other.y && z == other.z;
    }
  };

  /** Hash over the bit patterns of a PositionKey's three coordinates. */
  struct PositionKeyHash {

    size_t operator()( const PositionKey& key ) const {

      uint64_t bits[ 3 ];

      std::memcpy( bits, &key, sizeof( bits ) );

      uint64_t hash = 0xcbf29ce484222325ull;

      for ( uint64_t value : bits ) {

        hash = ( hash ^ value ) * 0x100000001b3ull;
      }

      return static_cast< size_t >( hash );
    }
  };

  using PositionSet =
    std::pmr::unordered_set< PositionKey, PositionKeyHash >;

  /**
   * Every position a mesh currently holds, for `tesselate`'s duplicate guard.
   *
   * @param mesh The mesh to read.
   * @param resource Memory resource to build the set in.
   * @return The set of the mesh's vertex positions.
   */
  template< typename VertexType >
  inline PositionSet positionsOf(
    const WingedEdgeMesh< VertexType >& mesh,
    std::pmr::memory_resource*          resource ) {

    PositionSet result( mesh.vertices.size() * 2, resource );

    for ( const VertexType& vertex : mesh.vertices ) {

      result.emplace( refinementPoint( vertex ) );
    }

    return result;
  }

  /**
   * Squared deflection threshold for `tesselate`, relative to the seed
   * mesh's own extent — 0.1% of its bounding-box diagonal, squared to match
   * the squared-deflection comparison in the refinement loop — floored at
   * OBJECT_DEFLECTION_FLOOR_FACTOR of the extent of the object the face
   * belongs to.
   *
   * Unit-independence is the point. The shared absolute MAX_DEFLECTION
   * (1e-6, i.e. 1mm linear deflection) silently assumed millimetre-ish
   * numeric scale: in metre-unit STEP files (Onshape AP242 exports) every
   * curved face smaller than a metre bottomed out on the absolute constant
   * and stopped refining at ~1mm true deflection — visibly faceted
   * fingertip-sized B-spline/cylinder faces on the AmazingHand model —
   * while in millimetre-unit files the same constant was so fine it read
   * as "refine until the triangle budget runs out" (see the extrusion
   * unwrap's jet-compressor note). The relative criterion follows the
   * 0.1%-of-extent convention the revolution/extrusion/cylinder unwraps
   * already used; the tiny floor (the 2^-24 loop-point quantisation grid
   * of IfcCurve::Add3d, below which a deflection target is meaningless)
   * only guards degenerate zero-extent seeds against refine-to-zero churn.
   *
   * @param mesh The seed mesh about to be refined.
   * @param representationExtent The diagonal of the extent of the
   *                    REPRESENTATION that defines this face, in the same
   *                    units as the mesh's vertices, pinned once per
   *                    representation by the extractor that owns the parsed
   *                    file. Scoping it to the definition rather than to the
   *                    model is what makes it well defined at all — see
   *                    ParamsAddFaceToGeometry::representationExtent. Zero
   *                    means "not known" and leaves the per-face target
   *                    unfloored, which is the pre-#564 behaviour; the IFC
   *                    front end passes zero.
   * @return The squared deflection threshold to pass to `tesselate`.
   */
  template< typename VertexType >
  inline double relativeDeflectionSquared(
      const WingedEdgeMesh< VertexType >& mesh,
      double representationExtent ) {

    glm::dvec3 boxMin( std::numeric_limits< double >::max() );
    glm::dvec3 boxMax( std::numeric_limits< double >::lowest() );

    for ( const VertexType& vertex : mesh.vertices ) {

      const glm::dvec3& point = refinementPoint( vertex );

      boxMin = glm::min( boxMin, point );
      boxMax = glm::max( boxMax, point );
    }

    constexpr double RELATIVE_DEFLECTION_FACTOR = 1e-3;
    constexpr double MIN_DEFLECTION             = 0x1p-24;

    double deflection =
      mesh.vertices.empty() ?
        0.0 :
        glm::distance( boxMin, boxMax ) * RELATIVE_DEFLECTION_FACTOR;

    // A non-finite or negative extent is treated as "not known" rather than
    // propagated: std::max would carry a NaN straight into the refinement
    // comparison, where every `>` is false and the loop stops on the first
    // edge, silently under-tessellating the whole model.
    if ( std::isfinite( representationExtent ) && representationExtent > 0.0 ) {

      deflection =
        std::max(
          deflection, representationExtent * OBJECT_DEFLECTION_FLOOR_FACTOR );
    }

    return std::max( MIN_DEFLECTION * MIN_DEFLECTION, deflection * deflection );
  }

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
  template< typename VertexType > 
  struct CandidateEdge {

    double   deflection;
    uint32_t edge;

    VertexType vertex;
  };

  /**
   * Sorting operator for candidate edge priority.
   */
  template < typename VertexType >
  inline bool operator<( const CandidateEdge< VertexType >& left, const CandidateEdge< VertexType >& right ) {

    return ( left.deflection < right.deflection ) ||
      ( left.deflection == right.deflection && left.edge < right.edge );
  }

  /**
   * Append a winged edge mesh 
   */
  inline void appendMeshToGeometry( WingedEdgeMesh< ParameterVertex >& mesh, Geometry& geometry, bool sameSense ) {

    uint32_t baseVertex = geometry.vertices.size();

    for ( ConnectedTriangle& triangle : mesh.triangles ) {

      uint32_t v0 = triangle.vertices[ 0 ];
      uint32_t v1 = triangle.vertices[ 1 ];
      uint32_t v2 = triangle.vertices[ 2 ];

      if ( ( !isCCW(  
        mesh.vertices[ v0 ],
        mesh.vertices[ v1 ],
        mesh.vertices[ v2 ] ) ) != sameSense ) {

        std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
        std::swap( triangle.edges[ 0 ], triangle.edges[ 2 ] );
      }

    }

    for ( size_t vertexIndex = 0, end = mesh.vertices.size(); vertexIndex < end; ++vertexIndex ) {

      // Note, we have to have local versions of vertex and normal
      // cos addpoint isn't const correct.
      geometry.MakeVertex( mesh.vertices[ vertexIndex ].point );
    }

    for ( const ConnectedTriangle& triangle : mesh.triangles ) {   

      geometry.MakeTriangle(
        baseVertex + triangle.vertices[ 0 ],
        baseVertex + triangle.vertices[ 1 ],
        baseVertex + triangle.vertices[ 2 ] );
    }
  }

  /**
   * Append a winged edge mesh 
   */
  inline void appendMeshToGeometry( WingedEdgeMesh< glm::dvec3 >& mesh, Geometry& geometry ) {

    uint32_t baseVertex = geometry.vertices.size();

    geometry.vertices.reserve( geometry.vertices.size() + mesh.vertices.size() );
    geometry.vertices.insert(
      geometry.vertices.end(),
      mesh.vertices.begin(),
      mesh.vertices.end() );

    for ( const ConnectedTriangle& triangle : mesh.triangles ) {   

      uint32_t v0 = triangle.vertices[ 0 ];
      uint32_t v1 = triangle.vertices[ 1 ];
      uint32_t v2 = triangle.vertices[ 2 ];

      if ( conway::orient2D(  
        mesh.vertices[ v0 ],
        mesh.vertices[ v1 ],
        mesh.vertices[ v2 ]) < 0 ) {

        std::swap( v0, v2 );
      }

      geometry.MakeTriangle(
        baseVertex + v0,
        baseVertex + v1,
        baseVertex + v2 );
    }
  }

  /**
   * Append a winged edge mesh, orienting every triangle against the
   * analytic normal of the surface it lies on.
   *
   * The two-argument overload above cannot do this. It calls orient2D,
   * which projects each triangle onto ITS OWN best axis pair and forces a
   * positive sign there. On a planar face every triangle shares a plane,
   * one projection is chosen for all of them, and the result is
   * consistent. On a curved face the dominant axis changes as the surface
   * turns, so triangles end up oriented toward a fixed half-space instead
   * of consistently outward — roughly a 180 degree arc of every cylinder
   * came back wound inward, while planar faces were always correct
   * (https://github.com/bldrs-ai/conway/issues/459).
   *
   * Taking the normal from the surface removes the guess: `normalAt`
   * returns the surface's own outward normal at a point (it need not be
   * unit length, only correctly directed), and `sameSense` applies the
   * STEP advanced_face flag that says whether the face agrees with it —
   * so a boss and the bore it sits in get opposite windings from the same
   * cylinder, which is what the flag is for. The unwrap paths previously
   * computed `sameSense` and then discarded it here.
   *
   * @param mesh The tesselated mesh, in world space.
   * @param geometry Destination geometry.
   * @param sameSense The face's same_sense flag.
   * @param normalAt Outward surface normal at a point.
   */
  template< typename NormalFunction >
  inline void appendMeshToGeometry(
    WingedEdgeMesh< glm::dvec3 >& mesh,
    Geometry& geometry,
    bool sameSense,
    NormalFunction&& normalAt ) {

    uint32_t baseVertex = geometry.vertices.size();

    geometry.vertices.reserve( geometry.vertices.size() + mesh.vertices.size() );
    geometry.vertices.insert(
      geometry.vertices.end(),
      mesh.vertices.begin(),
      mesh.vertices.end() );

    for ( const ConnectedTriangle& triangle : mesh.triangles ) {

      uint32_t v0 = triangle.vertices[ 0 ];
      uint32_t v1 = triangle.vertices[ 1 ];
      uint32_t v2 = triangle.vertices[ 2 ];

      const glm::dvec3& p0 = mesh.vertices[ v0 ];
      const glm::dvec3& p1 = mesh.vertices[ v1 ];
      const glm::dvec3& p2 = mesh.vertices[ v2 ];

      glm::dvec3 winding = glm::cross( p1 - p0, p2 - p0 );
      glm::dvec3 outward = normalAt( ( p0 + p1 + p2 ) / 3.0 );

      if ( !sameSense ) {

        outward = -outward;
      }

      double agreement = glm::dot( winding, outward );

      // A zero dot means a degenerate triangle or a point where the
      // normal is undefined (on the axis). Neither has an orientation
      // worth flipping for, so it is left as the triangulator produced
      // it rather than swapped on the sign of noise.
      if ( agreement < 0.0 ) {

        std::swap( v0, v2 );
      }

      geometry.MakeTriangle(
        baseVertex + v0,
        baseVertex + v1,
        baseVertex + v2 );
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

    // AFTP: back the subdivision candidate heap with the thread scratch arena
    // too (this runs inside the per-face ScratchArenaScope of the ParameterVertex
    // tessellators). Byte-identical: heap order is set by the comparator, not the
    // allocator. Explicit std::less matches the default priority_queue ordering.
    std::priority_queue<
      CandidateEdge< ParameterVertex >,
      std::pmr::vector< CandidateEdge< ParameterVertex > >,
      std::less< CandidateEdge< ParameterVertex > > >
      candidates{
        std::less< CandidateEdge< ParameterVertex > >(),
        std::pmr::vector< CandidateEdge< ParameterVertex > >(
          conway::ThreadScratchResource() ) };

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
      conway::AllocTagScope surfaceTag( conway::AllocSite::SurfaceEval );
      glm::dvec3 newPoint     = surface( averagePoint, newUV );

      glm::dvec3 deltaNewPoint = newPoint - averagePoint;

      double deflection = glm::dot( deltaNewPoint, deltaNewPoint );

      if ( minimumDeflection > deflection ) {
        return;
      }

      candidates.push( CandidateEdge< ParameterVertex > { 
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

    // Nothing above the deflection target, so nothing below can run. Taken
    // before the position set is built, so a face that needs no refinement -
    // most of them - does not pay for one.
    if ( candidates.empty() ) {
      return;
    }

    // Every position the mesh already holds, for the duplicate guard in the
    // loop below. Arena-backed for the same reason the candidate heap is: it
    // lives exactly as long as the per-face ScratchArenaScope around this.
    PositionSet positions = positionsOf( mesh, conway::ThreadScratchResource() );

    maximumTriangles -= mesh.triangles.size();

    while ( !candidates.empty() && maximumTriangles > 0 ) {

      const CandidateEdge< ParameterVertex >& candidate = candidates.top();

      // copy edge because it mutates later
      // as may the references as the vector re-allocates.
      Edge                     edge         = mesh.edges[ candidate.edge ];

      // A queued candidate can go stale. addCandidate() rejects border edges
      // at queue time, but every subdivision below deletes two triangles and
      // deleteTriangle() clears both slots of each of their edges, so an edge
      // that was interior when queued can be a border - or, as here, fully
      // detached with EMPTY_INDEX in both slots - by the time it is popped.
      // Without this the next two lines index triangles[ 0xFFFFFFFF ] and the
      // wasm heap traps: "memory access out of bounds" on eight spherical
      // faces of Orbiter_v1.1_Gear_7.5.step, deterministically
      // (conway-geom#172). Skipping is not a loss of refinement - a border
      // edge is never subdividable in this scheme, and the triangles that
      // replaced this one were re-queued by the addCandidate() calls at the
      // end of the loop.
      if ( edge.border() ) {

        candidates.pop();
        continue;
      }

      const ConnectedTriangle& t0           = mesh.triangles[ edge.triangles[ 0 ] ];
      const ConnectedTriangle& t1           = mesh.triangles[ edge.triangles[ 1 ] ];
      uint32_t                 otherVertex0 = t0.otherVertex( edge );
      uint32_t                 otherVertex1 = t1.otherVertex( edge );

      // Drop a split that would not add a position the mesh does not already
      // have. `tesselate` splits at the midpoint of the edge's two
      // PARAMETERS, and on a trim that lands on an exactly-representable uv
      // lattice that midpoint is a FIXED POINT: newUV comes out bitwise equal
      // to the uv of a vertex already in the mesh, the surface returns that
      // vertex's point bitwise, and the "new" vertex is a copy of an old one.
      // The split then spends two triangles of budget, emits triangles that
      // are exactly degenerate for Reify's weld to discard (two of the four,
      // where the duplicate is one triangle's apex), and hands the queue a
      // sibling edge carrying the deflection the parent had - so the queue
      // never drains and the budget goes nowhere. On face `#51059` of
      // Orbiter_v1.1_Gear_7.5.step that is 78,610 of 78,833 splits; across
      // the 164 faces of that model that split more than 200 times, 658,807
      // of 718,749 (bldrs-ai/conway#625, which also records what the
      // reclaimed budget then goes on - this guard does not make a face
      // refine less, it makes what it already spends land on the surface).
      //
      // Bitwise, deliberately, and NOT "the deflection did not go down".
      // The arithmetic test sounds like the principled one and does not work:
      // measured over the same run it accepts 709,705 of those 718,749
      // splits, because each split really does hand its own children a
      // smaller deflection - it is the sibling edge it creates that comes
      // back with the parent's, which no test on the parent can see. Exact
      // equality needs no tolerance, cannot under-refine a surface that is
      // still converging, and is precisely the condition under which the
      // split provably adds nothing.
      //
      // Dropped, not re-queued. The candidate is a pure function of the
      // edge's two endpoints and the surface, and neither changes afterwards
      // (vertices are appended, never moved; an edge is destroyed by its own
      // split), so re-queuing could only pop the same rejected candidate
      // again. The vertex it duplicates also stays in the mesh, so the
      // rejection cannot go stale in the other direction either.
      if ( !positions.emplace( candidate.vertex.point ).second ) {

        candidates.pop();
        continue;
      }

      uint32_t                 newVertex    = mesh.makeVertex( candidate.vertex );

      candidates.pop();

      auto [ t0Index, t1Index ] = edge.triangles;

      if ( t0Index > t1Index ) {
        std::swap( t0Index, t1Index );
      }

      mesh.deleteTriangle( t1Index );
      mesh.deleteTriangle( t0Index );

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

  /**
   * Given a surface where a mid-point can be re-computed as point on the surface,
   * this will take a starting mesh with parameterized vertices and tesselate the internal triangles
   */
  template< typename SurfacePointFunction >
  inline void tesselate(
    WingedEdgeMesh< glm::dvec3 >& mesh,
    SurfacePointFunction surface,
    int32_t maximumTriangles,
    double minimumDeflection ) {

    std::priority_queue< CandidateEdge< glm::dvec3 > > candidates;

    auto addCandidate = [&]( uint32_t edgeIndex ) {

      if ( edgeIndex == EMPTY_INDEX  ) {
        return;
      }

      const Edge& edge = mesh.edges[ edgeIndex ];

      if ( edge.border() ) {
        return;
      }

      const glm::dvec3& v0   = mesh.vertices[ edge.vertices[ 0 ] ];
      const glm::dvec3& v1   = mesh.vertices[ edge.vertices[ 1 ] ];

      glm::dvec3 averagePoint = ( v0 + v1 ) * 0.5;
      glm::dvec3 newPoint     = surface( averagePoint );

      glm::dvec3 deltaNewPoint = newPoint - averagePoint;

      double deflection = glm::dot( deltaNewPoint, deltaNewPoint );

      if ( minimumDeflection > deflection ) {
        return;
      }

      candidates.push( ( CandidateEdge< glm::dvec3 > { 
        deflection * glm::distance( v0, v1 ),
        edgeIndex,
        newPoint 
      } ) );      
    };

    for (
      uint32_t edgeIndex = 0, end = static_cast< uint32_t >( mesh.edges.size() );
      edgeIndex < end;
      ++edgeIndex ) {

      addCandidate( edgeIndex );
    }

    // As above.
    if ( candidates.empty() ) {
      return;
    }

    // As above. This overload runs outside the parameterized tessellators'
    // scratch scope, so its set goes on the default resource, like its
    // candidate heap.
    PositionSet positions =
      positionsOf( mesh, std::pmr::get_default_resource() );

    maximumTriangles -= mesh.triangles.size();

    while ( !candidates.empty() && maximumTriangles > 0 ) {

      const CandidateEdge< glm::dvec3 >&    candidate = candidates.top();
      // copy edge because it mutates later
      // as may the references as the vector re-allocates.
      Edge                     edge         = mesh.edges[ candidate.edge ];

      // A queued candidate can go stale. addCandidate() rejects border edges
      // at queue time, but every subdivision below deletes two triangles and
      // deleteTriangle() clears both slots of each of their edges, so an edge
      // that was interior when queued can be a border - or, as here, fully
      // detached with EMPTY_INDEX in both slots - by the time it is popped.
      // Without this the next two lines index triangles[ 0xFFFFFFFF ] and the
      // wasm heap traps: "memory access out of bounds" on eight spherical
      // faces of Orbiter_v1.1_Gear_7.5.step, deterministically
      // (conway-geom#172). Skipping is not a loss of refinement - a border
      // edge is never subdividable in this scheme, and the triangles that
      // replaced this one were re-queued by the addCandidate() calls at the
      // end of the loop.
      if ( edge.border() ) {

        candidates.pop();
        continue;
      }

      const ConnectedTriangle& t0           = mesh.triangles[ edge.triangles[ 0 ] ];
      const ConnectedTriangle& t1           = mesh.triangles[ edge.triangles[ 1 ] ];
      uint32_t                 otherVertex0 = t0.otherVertex( edge );
      uint32_t                 otherVertex1 = t1.otherVertex( edge );

      // The same fixed-point guard as the ParameterVertex overload above,
      // which carries the reasoning. This overload re-projects the midpoint
      // rather than re-evaluating a parameter, but it is the same map onto
      // the same lattice and it duplicates the same way - 93 of the 164
      // faces that refine on Orbiter_v1.1_Gear_7.5.step take this path.
      if ( !positions.emplace( candidate.vertex ).second ) {

        candidates.pop();
        continue;
      }

      uint32_t                 newVertex    = mesh.makeVertex( candidate.vertex );

      candidates.pop();

      auto [ t0Index, t1Index ] = edge.triangles;

      if ( t0Index > t1Index ) {
        std::swap( t0Index, t1Index );
      }

      mesh.deleteTriangle( t1Index );
      mesh.deleteTriangle( t0Index );

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
