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
