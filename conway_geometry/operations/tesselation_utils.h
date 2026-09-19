#pragma once

#include <cmath>
#include <limits>
#include <glm/glm.hpp>

#include "structures/winged_edge.h"
#include "structures/scratch_arena.h"
#include "structures/alloc_telemetry.h"
#include <memory_resource>
#include <queue>
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
   * Decide, once for a whole face, whether the analytic surface normal points
   * along the face's triangle winding or against it.
   *
   * Derived from the mesh rather than from `same_sense` because the two append
   * overloads reach their final winding by different routes — one negates the
   * surface normal and flips to agree, the other flips on a CCW test in the
   * parameter domain — so a shared convention read off the flag would be wrong
   * for one of them. The winding of the emitted triangles is the authoritative
   * outward direction by the time this runs, so ask it directly.
   *
   * The vote is area-weighted (the cross product is left unnormalized) and
   * taken over the whole face, which is what makes it robust: a face is one
   * surface with one orientation, and the sliver triangles along a trimmed
   * boundary — whose individual winding normals are the unreliable ones this
   * whole change exists to stop trusting — carry almost no weight.
   *
   * @param mesh The tesselated mesh, with final winding.
   * @param positionOf Maps a mesh vertex to its position.
   * @param vertexNormals Per-mesh-vertex surface normals, from
   *   evaluateVertexNormals().
   * @return +1.0 when the surface normal agrees with the winding, else -1.0.
   */
  template< typename VertexType, typename PositionFunction >
  inline double analyticNormalSign(
    const WingedEdgeMesh< VertexType >& mesh,
    PositionFunction&& positionOf,
    const std::vector< glm::dvec3 >& vertexNormals ) {

    double agreement = 0.0;

    for ( const ConnectedTriangle& triangle : mesh.triangles ) {

      uint32_t i0 = triangle.vertices[ 0 ];
      uint32_t i1 = triangle.vertices[ 1 ];
      uint32_t i2 = triangle.vertices[ 2 ];

      const glm::dvec3& p0 = positionOf( mesh.vertices[ i0 ] );
      const glm::dvec3& p1 = positionOf( mesh.vertices[ i1 ] );
      const glm::dvec3& p2 = positionOf( mesh.vertices[ i2 ] );

      glm::dvec3 winding = glm::cross( p1 - p0, p2 - p0 );

      // Summed at the corners rather than evaluated once at the centroid: a
      // centroid can fall on the axis of a cylinder or cone, where the normal
      // is undefined, even when all three corners are well away from it.
      glm::dvec3 outward =
        vertexNormals[ i0 ] + vertexNormals[ i1 ] + vertexNormals[ i2 ];

      double contribution = glm::dot( winding, outward );

      if ( std::isfinite( contribution ) ) {

        agreement += contribution;
      }
    }

    return agreement < 0.0 ? -1.0 : 1.0;
  }

  /**
   * Evaluate the analytic surface normal once per MESH VERTEX.
   *
   * Corners that share a mesh vertex share its position and parameters, so they
   * share its surface normal exactly — this is a cache, not an approximation.
   * Evaluating per corner instead costs about six times as much (the average
   * vertex valence), which on B-spline faces is the difference between a cheap
   * change and a 3.4x geometry-time regression: tinynurbs::surfaceNormal builds
   * a full derivative array and has none of the fast paths point() has.
   *
   * @param mesh The tesselated mesh.
   * @param normalOf Maps a mesh vertex to its outward surface normal.
   * @return One normal per mesh vertex, zero where undefined.
   */
  template< typename VertexType, typename NormalFunction >
  inline std::vector< glm::dvec3 > evaluateVertexNormals(
    const WingedEdgeMesh< VertexType >& mesh,
    NormalFunction&& normalOf ) {

    std::vector< glm::dvec3 > normals;

    normals.reserve( mesh.vertices.size() );

    for ( const VertexType& vertex : mesh.vertices ) {

      glm::dvec3 normal = normalOf( vertex );

      double length = glm::length( normal );

      // A pole, an apex, or a collapsed B-spline row: no usable normal. Stored
      // as zero so both consumers below read it as "fall back".
      normals.push_back(
        ( !std::isfinite( length ) || length < DBL_EPSILON ) ?
          glm::dvec3( 0.0 ) : normal / length );
    }

    return normals;
  }

  /**
   * Evaluate and store the analytic surface normal at each corner of one
   * triangle, oriented to match the face's sense.
   *
   * A point where the normal is undefined — a cone apex, a sphere pole, the
   * axis of a surface of revolution — yields a zero or non-finite vector.
   * Those corners are left at zero, which Reify() reads as "fall back to the
   * face normal", rather than being normalized into a NaN that would poison the
   * whole smoothing group.
   *
   * @param geometry Destination geometry.
   * @param triangleIndex Triangle whose corners are being recorded.
   * @param vertexNormals Per-mesh-vertex normals from evaluateVertexNormals().
   * @param i0 Mesh vertex index at corner 0.
   * @param i1 Mesh vertex index at corner 1.
   * @param i2 Mesh vertex index at corner 2.
   * @param sign +1 or -1, from analyticNormalSign() for the whole face.
   */
  inline void recordCornerNormals(
    Geometry& geometry,
    uint32_t  triangleIndex,
    const std::vector< glm::dvec3 >& vertexNormals,
    uint32_t i0,
    uint32_t i1,
    uint32_t i2,
    double sign ) {

    geometry.SetCornerNormals(
      triangleIndex,
      glm::vec3( vertexNormals[ i0 ] * sign ),
      glm::vec3( vertexNormals[ i1 ] * sign ),
      glm::vec3( vertexNormals[ i2 ] * sign ) );
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
   * UV-parameterized append that also records the analytic normal per corner.
   *
   * Same body as the three-argument overload above — the winding decision is
   * unchanged — with the shading normals captured in the second pass, after the
   * corner order is final. Split as an overload rather than a default argument
   * so the callers that have no surface normal to give keep byte-identical
   * output (bldrs-ai/conway#667).
   *
   * Takes the normal per VERTEX rather than per point, because a B-spline's
   * normal is a function of (u, v) and only the vertex carries that; the
   * quadric surfaces just read `.point` and ignore the parameters.
   *
   * @param mesh The tesselated mesh.
   * @param geometry Destination geometry.
   * @param sameSense The face's same_sense flag.
   * @param normalOf Maps a ParameterVertex to its outward surface normal.
   */
  template< typename NormalFunction >
  inline void appendMeshToGeometry(
    WingedEdgeMesh< ParameterVertex >& mesh,
    Geometry& geometry,
    bool sameSense,
    NormalFunction&& normalOf ) {

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

      geometry.MakeVertex( mesh.vertices[ vertexIndex ].point );
    }

    auto positionOf = []( const ParameterVertex& vertex ) -> const glm::dvec3& {

      return vertex.point;
    };

    std::vector< glm::dvec3 > vertexNormals =
      evaluateVertexNormals( mesh, normalOf );

    double sign = analyticNormalSign( mesh, positionOf, vertexNormals );

    for ( const ConnectedTriangle& triangle : mesh.triangles ) {

      uint32_t triangleIndex = geometry.MakeTriangle(
        baseVertex + triangle.vertices[ 0 ],
        baseVertex + triangle.vertices[ 1 ],
        baseVertex + triangle.vertices[ 2 ] );

      recordCornerNormals(
        geometry,
        triangleIndex,
        vertexNormals,
        triangle.vertices[ 0 ],
        triangle.vertices[ 1 ],
        triangle.vertices[ 2 ],
        sign );
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

    // Vertices here ARE positions, so the winding function applies unchanged.
    // Evaluated once per vertex, not once per corner — see
    // evaluateVertexNormals for why that matters.
    std::vector< glm::dvec3 > vertexNormals =
      evaluateVertexNormals( mesh, normalAt );

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

      uint32_t triangleIndex = geometry.MakeTriangle(
        baseVertex + v0,
        baseVertex + v1,
        baseVertex + v2 );

      // Record the analytic normal per corner. `outward` above is only the
      // direction at the centroid, used to pick a winding; shading wants the
      // normal AT each corner. No analyticNormalSign() vote is needed here,
      // unlike the UV overload: this loop has just forced every triangle's
      // winding to agree with `sameSense`-negated `normalAt`, so that same
      // negation is the shading sign by construction (bldrs-ai/conway#667).
      recordCornerNormals(
        geometry,
        triangleIndex,
        vertexNormals,
        v0,
        v1,
        v2,
        sameSense ? 1.0 : -1.0 );
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
   * THE BORDER EDGES `tesselate` ABOVE CANNOT TOUCH, REFINED IN PAIRS.
   *
   * `addCandidate` skips every `edge.border()`, and it has to: a border edge
   * carries one triangle, so the four-triangle replacement above has nothing
   * to build from, and a trim boundary is SHARED with the face on the other
   * side of it, which this face may not move without opening a crack.
   *
   * A periodic chart's seam is neither of those things.
   * `triangulatePeriodicUChart` cuts the chart open AFTER triangulating it,
   * by duplicating the vertices two triangles disagree about, so one interior
   * edge of the chart becomes two mesh edges with one triangle each: border
   * by the winged-edge reading, interior by the face's own geometry, and
   * shared with nobody. Left to itself it keeps the triangulation's
   * coarseness while everything beside it refines, and the single triangle on
   * each side is subdivided away from it until what is left of it is a sliver
   * lying on it - measured on `ADVANCED_FACE #19218` of Right_Hand.step as
   * 6760 of that face's 8320 triangles coming out degenerate, against 4 of
   * 1848 once the cut is refined first.
   *
   * So the two sides are split in LOCKSTEP: the same parameter on both, and
   * two new vertices carrying a BITWISE-identical position, which is what
   * keeps `Geometry::Reify`'s weld closing the seam by identity exactly as it
   * closes the duplicated corners the chart handed over. Splitting one side
   * alone would leave a T-junction against the other, which is the failure
   * this is here to avoid rather than to cause.
   *
   * THE SPLIT POINT IS EVALUATED ON THE SURFACE, unlike the layout split in
   * `triangulatePeriodicUChart`, which is placed on the segment it splits
   * because that segment is a shared trim boundary. The seam is interior:
   * there is no neighbour holding the other half of it, so there is nothing
   * to crack against, and following the surface is the whole point of
   * refining.
   *
   * Runs BEFORE `tesselate`, so the interior refinement afterwards sees the
   * seam at its final density and refines against it; both are handed the
   * same deflection target, or the seam would be refined to a different
   * fineness than its own neighbourhood.
   *
   * @param seams             Periodic copies of one chart edge, as
   *                          { a, b, a', b' } with a' the copy of a. Anything
   *                          that is no longer a pair of distinct
   *                          single-triangle border edges is skipped, so a
   *                          stale entry costs nothing.
   * @param maximumTriangles  The same budget `tesselate` takes, counted the
   *                          same way: the mesh's current size is subtracted
   *                          and each split spends two.
   */
  template< typename SurfacePointFunction >
  inline void refineSeamPairs(
    WingedEdgeMesh< ParameterVertex >& mesh,
    SurfacePointFunction surface,
    const std::vector< std::array< uint32_t, 4 > >& seams,
    int32_t maximumTriangles,
    double minimumDeflection ) {

    if ( seams.empty() ) {
      return;
    }

    maximumTriangles -= static_cast< int32_t >( mesh.triangles.size() );

    // The corners of `triangle` rotated so that ( result[ 0 ], result[ 1 ] )
    // is the edge ( a, b ) IN THE TRIANGLE'S OWN ORDER, which is what lets
    // the two halves below be emitted with the winding they replace.
    const auto orientedCorners =
      [ & ]( uint32_t triangle, uint32_t a, uint32_t b ) {

        const ConnectedTriangle& face = mesh.triangles[ triangle ];

        for ( uint32_t at = 0; at < 3; ++at ) {

          const uint32_t first  = face.vertices[ at ];
          const uint32_t second = face.vertices[ ( at + 1 ) % 3 ];

          if ( ( first == a && second == b ) ||
               ( first == b && second == a ) ) {

            return std::array< uint32_t, 3 >{
              first, second, face.vertices[ ( at + 2 ) % 3 ] };
          }
        }

        return std::array< uint32_t, 3 >{
          EMPTY_INDEX, EMPTY_INDEX, EMPTY_INDEX };
      };

    // Breadth-first bisection rather than `tesselate`'s deflection-ordered
    // heap: the seam is a path, every edge of it is measured against the same
    // floor, and there is no competition for the budget to arbitrate.
    std::vector< std::array< uint32_t, 4 > > pending( seams );

    for ( size_t head = 0;
          head < pending.size() && maximumTriangles > 0;
          ++head ) {

      const std::array< uint32_t, 4 > seam = pending[ head ];

      const std::optional< uint32_t > edge0 = mesh.getEdge( seam[ 0 ], seam[ 1 ] );
      const std::optional< uint32_t > edge1 = mesh.getEdge( seam[ 2 ], seam[ 3 ] );

      if ( !edge0.has_value() || !edge1.has_value() ||
           edge0.value() == edge1.value() ) {
        continue;
      }

      const uint32_t triangle0 = mesh.edges[ edge0.value() ].triangles[ 0 ];
      const uint32_t triangle1 = mesh.edges[ edge1.value() ].triangles[ 0 ];

      // A side that is not a live single-triangle border is not a side this
      // can split: deleteTriangle() leaves fully detached edges behind
      // (EMPTY_INDEX in both slots), and makeEdge() lets a non-manifold edge
      // carry two triangles.
      if ( !mesh.edges[ edge0.value() ].border() ||
           !mesh.edges[ edge1.value() ].border() ||
           triangle0 == EMPTY_INDEX || triangle1 == EMPTY_INDEX ||
           triangle0 == triangle1 ) {
        continue;
      }

      const ParameterVertex& from = mesh.vertices[ seam[ 0 ] ];
      const ParameterVertex& to   = mesh.vertices[ seam[ 1 ] ];

      const glm::dvec3 averagePoint = ( from.point + to.point ) * 0.5;
      const glm::dvec2 newUV        = ( from.uv + to.uv ) * 0.5;

      glm::dvec3 newPoint;

      {
        conway::AllocTagScope surfaceTag( conway::AllocSite::SurfaceEval );

        newPoint = surface( averagePoint, newUV );
      }

      const glm::dvec3 deltaNewPoint = newPoint - averagePoint;

      const double deflection = glm::dot( deltaNewPoint, deltaNewPoint );

      // The same reading, against the same floor, as addCandidate above.
      if ( minimumDeflection > deflection ) {
        continue;
      }

      // A READING THAT IS NOT ABOUT THE CHORD. Bisection drives the deflection
      // to zero only when the two ends lie on the surface the callback
      // describes; when one of them does not, the reading bottoms out at that
      // disagreement instead, and every further split halves an edge without
      // ever satisfying the floor. It happens: `triangulatePeriodicUChart`'s
      // own layout split places its point on the trim SEGMENT rather than on
      // the surface, because that segment is shared with the neighbouring
      // face - so a cut that ends on one of those points is anchored at a
      // vertex the surface does not pass through.
      //
      // Refused on the geometry rather than on a depth count. For two ends on
      // the surface the deflection of their chord is bounded by L*L*k/8 with k
      // the local curvature, so a deflection past a QUARTER of the chord needs
      // L > 2/k - a chord longer than the local diameter of curvature, which
      // two points of the surface do not span. A non-convergent one runs
      // straight into it instead: as the split points pile up at the surface
      // point the off-surface end's parameter names, the chord settles at that
      // end's own distance r from it and the deflection at r/2, which is half
      // the chord and so eight times over.
      //
      // Measured: without this, one cut of a six-sample band in
      // test/periodic_u_chart_test.cpp spent all 434 of its splits and left
      // 678 zero-length edges piled on one point. With it the same cut refines
      // and stops. Nothing in the corpus comes within two orders of magnitude
      // of the bound - `ADVANCED_FACE #19218`'s cut reads deflection 4.9e-8
      // against a chord of 1.85e-2, i.e. 1.6e-6 of the way to it.
      const glm::dvec3 chord = to.point - from.point;

      if ( ( deflection * 16.0 ) >= glm::dot( chord, chord ) ) {
        continue;
      }

      // The other side's own midpoint in uv - the same point on the surface,
      // a whole number of periods away, which is the disagreement the cut
      // exists to carry.
      const glm::dvec2 partnerUV =
        ( mesh.vertices[ seam[ 2 ] ].uv + mesh.vertices[ seam[ 3 ] ].uv ) * 0.5;

      const std::array< uint32_t, 3 > corners0 =
        orientedCorners( triangle0, seam[ 0 ], seam[ 1 ] );
      const std::array< uint32_t, 3 > corners1 =
        orientedCorners( triangle1, seam[ 2 ], seam[ 3 ] );

      if ( corners0[ 0 ] == EMPTY_INDEX || corners1[ 0 ] == EMPTY_INDEX ) {
        continue;
      }

      // Read both triangles' corners above, before either is deleted:
      // deleteTriangle() moves the back triangle into the freed slot, so the
      // second index and every reference into `triangles` goes stale the
      // moment the first one goes. Deleting the HIGHER index first is what
      // keeps the lower one addressable, exactly as tesselate() does it.
      mesh.deleteTriangle( std::max( triangle0, triangle1 ) );
      mesh.deleteTriangle( std::min( triangle0, triangle1 ) );

      const uint32_t split0 = mesh.makeVertex( { newPoint, newUV } );

      // BITWISE the same position at the other sheet's uv - the identity weld
      // the chart's own duplicates rely on, carried to the points this adds.
      const uint32_t split1 = mesh.makeVertex( { newPoint, partnerUV } );

      mesh.makeTriangle( corners0[ 0 ], split0, corners0[ 2 ] );
      mesh.makeTriangle( split0, corners0[ 1 ], corners0[ 2 ] );
      mesh.makeTriangle( corners1[ 0 ], split1, corners1[ 2 ] );
      mesh.makeTriangle( split1, corners1[ 1 ], corners1[ 2 ] );

      // Both halves, still paired: seam[ 2 ] is the copy of seam[ 0 ] and
      // split1 of split0, so the correspondence carries down the recursion.
      pending.push_back( { seam[ 0 ], split0, seam[ 2 ], split1 } );
      pending.push_back( { split0, seam[ 1 ], split1, seam[ 3 ] } );

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
