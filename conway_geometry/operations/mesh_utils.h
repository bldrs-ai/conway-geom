/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/operations/mesh_utils.h
 * */

#pragma once

#include <tinynurbs/tinynurbs.h>

#include <algorithm>
#include <array>
#include <limits>
#include <glm/glm.hpp>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <ranges>

#include "geometry_utils.h"
#include "nurbs_utils.h"
#include "tesselation_utils.h"
#include "manifold_utils.h"
#include <queue>

#define CONST_PI 3.141592653589793238462643383279502884L

namespace conway::geometry {

constexpr double MAX_TRIANGLE_AMPLIFACTION = 32;


// TODO: review and simplify
// ---------------------------------------------------------------------------
// Shared helpers for unwrapping boundary loops of periodic surfaces (tori,
// surfaces of revolution) into planar CDT domains. BoundaryPoint's second
// parameter is the minor/tube angle for tori and the normalized profile
// arclength for revolutions.
// ---------------------------------------------------------------------------

namespace unwrap_detail {

constexpr double kPi    = 3.141592653589793238462643;
constexpr double kTwoPi = 2.0 * kPi;

/** Wrap an angle difference into (-pi, pi]. */
inline double wrapDeltaPi( double delta ) {

  delta = std::fmod( delta, kTwoPi );

  if ( delta > kPi ) {
    delta -= kTwoPi;
  } else if ( delta <= -kPi ) {
    delta += kTwoPi;
  }

  return delta;
}

/** Positive modulo into [0, 2pi). */
inline double positiveMod2Pi( double angle ) {

  double result = std::fmod( angle, kTwoPi );

  if ( result < 0 ) {
    result += kTwoPi;
  }

  return result;
}

/**
 * Rounding-error bound for a shoelace sum taken over `vertexCount` vertices
 * whose coordinates have been shifted to a reference point, with
 * `maxShiftedCoordinate` the largest magnitude among them.
 *
 * A shoelace result below this is indistinguishable from zero and the polygon
 * must be treated as enclosing nothing. A FIXED floor cannot serve: the charts
 * these loops live in are normalized per face, so the same absolute number
 * means "degenerate" on one face and "a real trim" on another. Measured, a
 * fixed 1e-12 classified a genuine four-corner patch of chart span 5e-7
 * (shoelace 8.4e-14) as degenerate.
 *
 * Derivation, with e = DBL_EPSILON and M = maxShiftedCoordinate:
 *
 *   - each shift fl(a - r) carries absolute error <= e*M;
 *   - a product of two shifted coordinates therefore carries
 *     M*(e*M) + M*(e*M) + e*M^2 = 3*e*M^2, and each shoelace term is a
 *     difference of two such products, so <= 8*e*M^2 including the
 *     subtraction's own rounding;
 *   - summing n terms adds running-sum rounding, bounded by
 *     n*e*max|partial sum| <= 2*n^2*e*M^2.
 *
 * 4*n^2*e*M^2 covers both with margin, and stays far below any real area: a
 * loop it can reject encloses less than the 1e-9 constraint weld can
 * represent, so it would already have welded onto a line.
 *
 * This is a "encloses literally nothing" bound, deliberately NOT a thinness
 * threshold - see bldrs-ai/conway#595 on why a shape ratio such as
 * area/perimeter^2 cannot be used, since it decays continuously as a genuine
 * lune narrows and would silently reject real trims at some width.
 */
inline double shoelaceAreaTolerance(
    size_t vertexCount, double maxShiftedCoordinate ) {

  const double n = static_cast< double >( vertexCount );

  return 4.0 * n * n * DBL_EPSILON *
         maxShiftedCoordinate * maxShiftedCoordinate;
}

struct BoundaryPoint {
  glm::dvec3 world;
  double     theta;
  double     phi;
};

/**
 * Largest empty circular gap of angles (each in [0, 2pi)).
 *
 * @return { gapSize, gapMidpoint } - for an empty input the whole circle is
 * one gap.
 */
inline std::pair< double, double > largestCircularGap( std::vector< double > &angles ) {

  if ( angles.empty() ) {
    return { kTwoPi, 0.0 };
  }

  std::sort( angles.begin(), angles.end() );

  // Wrap-around gap between the last and first sample.
  double bestGap = ( angles.front() + kTwoPi ) - angles.back();
  double bestMid = positiveMod2Pi( angles.back() + bestGap * 0.5 );

  for ( size_t where = 1; where < angles.size(); ++where ) {

    double gap = angles[ where ] - angles[ where - 1 ];

    if ( gap > bestGap ) {
      bestGap = gap;
      bestMid = angles[ where - 1 ] + gap * 0.5;
    }
  }

  return { bestGap, bestMid };
}

/**
 * Shared trimmed-unwrap triangulation used by the cylinder and cone paths
 * (the torus/revolution/extrusion unwraps predate it and keep their own
 * inlined copies of the same machinery).
 *
 * Boundary loops arrive parameterized as (theta, phi) - theta a periodic
 * angle in radians, phi non-periodic and pre-normalized to [0, 1]. They are
 * laid out in an injective 2D domain (annulus when the boundary wraps theta,
 * rectangle across the largest angular gap otherwise), welded at 1e-9
 * parameter quantization (seam-doubled boundary edges become single interior
 * constraints, which is what makes NotAllowed safe), then CDT'd with exact
 * hole nesting.
 *
 * @return true with `mesh` holding the trimmed triangulation (world
 * positions, unrefined) on success; false on any degeneracy - the caller
 * falls back to its previous behavior, so a dirty boundary can never make
 * output worse than it was.
 */
inline bool triangulateUnwrappedLoops(
    const std::vector< std::vector< BoundaryPoint > > &loops,
    WingedEdgeMesh< glm::dvec3 > &mesh,
    const char *label,
    size_t *closedLoopsReachingCdt = nullptr ) {

  if ( loops.empty() ) {
    return false;
  }

  bool windsTheta = false;

  std::vector< double > thetaSamples;

  for ( const std::vector< BoundaryPoint > &loop : loops ) {

    double sumTheta = 0.0;

    for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

      const BoundaryPoint &from = loop[ where ];
      const BoundaryPoint &to   = loop[ ( where + 1 ) % count ];

      sumTheta += wrapDeltaPi( to.theta - from.theta );

      thetaSamples.push_back( positiveMod2Pi( from.theta ) );
    }

    if ( std::llround( sumTheta / kTwoPi ) != 0 ) {
      windsTheta = true;
    }
  }

  auto [ thetaGap, thetaCut ] = largestCircularGap( thetaSamples );

  // A boundary that covers the whole circle without netting a winding
  // (seam edges walked both ways) leaves no gap for a branch cut.
  bool wrapsTheta = windsTheta || thetaGap < 1e-3;

  std::vector< std::vector< glm::dvec2 > > loops2D;

  bool layoutOk = false;

  // A rectangle edge jumping more than pi in cut-normalized theta straddles
  // the branch cut - upgrade to the wrapping (annulus) layout, which has no
  // cut to straddle, and retry once.
  for ( int attempt = 0; attempt < 2 && !layoutOk; ++attempt ) {

    loops2D.clear();
    loops2D.reserve( loops.size() );
    layoutOk = true;

    if ( wrapsTheta ) {

      // Annulus: angle = theta, radius = 1 + phi in [1, 2].
      for ( const std::vector< BoundaryPoint > &loop : loops ) {

        std::vector< glm::dvec2 > loop2D;

        loop2D.reserve( loop.size() );

        for ( const BoundaryPoint &point : loop ) {
          loop2D.emplace_back(
            ( 1.0 + point.phi ) * std::cos( point.theta ),
            ( 1.0 + point.phi ) * std::sin( point.theta ) );
        }

        loops2D.push_back( std::move( loop2D ) );
      }

    } else {

      double loTheta = std::numeric_limits< double >::max();
      double hiTheta = std::numeric_limits< double >::lowest();

      std::vector< std::vector< double > > normalized;

      normalized.reserve( loops.size() );

      for ( const std::vector< BoundaryPoint > &loop : loops ) {

        std::vector< double > loopNormalized;

        loopNormalized.reserve( loop.size() );

        double previous = 0.0;

        for ( size_t where = 0; where < loop.size(); ++where ) {

          double value = positiveMod2Pi( loop[ where ].theta - thetaCut );

          if ( where > 0 && std::abs( value - previous ) > kPi ) {
            wrapsTheta = true;
            layoutOk   = false;
            break;
          }

          previous = value;
          loTheta  = std::min( loTheta, value );
          hiTheta  = std::max( hiTheta, value );
          loopNormalized.push_back( value );
        }

        if ( !layoutOk ) {
          break;
        }

        normalized.push_back( std::move( loopNormalized ) );
      }

      if ( layoutOk ) {

        double span = std::max( hiTheta - loTheta, 1e-12 );

        for ( size_t which = 0; which < loops.size(); ++which ) {

          std::vector< glm::dvec2 > loop2D;

          loop2D.reserve( loops[ which ].size() );

          for ( size_t where = 0; where < loops[ which ].size(); ++where ) {
            loop2D.emplace_back(
              ( normalized[ which ][ where ] - loTheta ) / span,
              loops[ which ][ where ].phi );
          }

          loops2D.push_back( std::move( loop2D ) );
        }
      }
    }
  }

  if ( !layoutOk ) {
    return false;
  }

  std::vector< CDT::V2d< double > > cdtVertices;
  std::vector< glm::dvec3 >         cdtWorld;
  std::vector< CDT::Edge >          cdtEdges;

  std::map< std::pair< long long, long long >, uint32_t > weld;
  std::set< std::pair< uint32_t, uint32_t > >             edgeSet;

  bool inputFinite = true;

  auto weldVertex = [&]( const glm::dvec2 &position, const glm::dvec3 &world ) {

    std::pair< long long, long long > key(
      std::llround( position.x * 1e9 ),
      std::llround( position.y * 1e9 ) );

    auto [ found, isNew ] =
      weld.try_emplace( key, static_cast< uint32_t >( cdtVertices.size() ) );

    if ( isNew ) {
      if ( !std::isfinite( position.x ) || !std::isfinite( position.y ) ) {
        inputFinite = false;
      }
      cdtVertices.emplace_back( position.x, position.y );
      cdtWorld.push_back( world );
    }

    return found->second;
  };

  // Tallied per loop rather than globally, and counted as DISTINCT WELDED
  // VERTICES rather than as edges emitted.
  //
  // The weld above runs at 1e-9 while callers dedup at their own tolerance, so
  // a loop can arrive here intact and still lose vertices. Two ways that ends
  // badly, and only the stronger test catches both:
  //
  //   - all vertices merge to one: no edges at all, the loop vanishes from the
  //     constraints, and `cdtEdges.size() < 3` still passes on the survivors;
  //   - vertices merge to exactly TWO: the loop emits a single open edge, which
  //     is a slit rather than a boundary. The CDT triangulates straight across
  //     it, so an inner trim's hole is filled with surface.
  //
  // An edge-emitted test sees the first and misses the second - measured, that
  // second case filled a hole with 832 triangles where the caller's fallback
  // emitted none. Three distinct vertices is the least that can enclose area,
  // so that is the bar. Reported out so a caller can require every loop it
  // supplied to have cleared it; callers that do not ask (cylinder, cone) keep
  // exactly their previous behaviour.
  size_t closedLoops = 0;

  for ( size_t which = 0; which < loops.size(); ++which ) {

    const std::vector< BoundaryPoint > &loop   = loops[ which ];
    const std::vector< glm::dvec2 >    &loop2D = loops2D[ which ];

    std::set< uint32_t > distinctWelded;

    // The welded vertex sequence, so the loop's enclosed area can be taken on
    // what the CDT will actually see rather than on the points as supplied.
    std::vector< uint32_t > weldedSequence;

    for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

      size_t next = ( where + 1 ) % count;

      uint32_t v1 = weldVertex( loop2D[ where ], loop[ where ].world );
      uint32_t v2 = weldVertex( loop2D[ next ], loop[ next ].world );

      if ( v1 == v2 ) {
        continue;
      }

      // Recorded before the edgeSet dedup below: a vertex this loop shares
      // with one already inserted is still a vertex this loop reached the
      // triangulation with.
      distinctWelded.insert( v1 );
      distinctWelded.insert( v2 );

      if ( weldedSequence.empty() ) {
        weldedSequence.push_back( v1 );
      }

      weldedSequence.push_back( v2 );

      std::pair< uint32_t, uint32_t > ordered(
        std::min( v1, v2 ), std::max( v1, v2 ) );

      if ( edgeSet.insert( ordered ).second ) {
        cdtEdges.emplace_back( v1, v2 );
      }
    }

    // AREA, not cardinality, is the closing test. Three distinct vertices are
    // necessary but NOT sufficient: three points strung along one meridian are
    // distinct, survive the weld, and still enclose nothing - a slit with an
    // extra point on it, which the CDT triangulates straight across exactly as
    // it does a two-point one. Measured, a collinear inner bound filled the
    // hole at 896 triangles where the intended band is 1.2578.
    //
    // The cardinality test is kept as the cheap precondition: it is O(1)
    // against a set that has to be built anyway, and it short-circuits the
    // shoelace for the common degenerate cases.
    //
    // This is a "does it enclose literally nothing" test, deliberately NOT a
    // thinness test. #595's comment block records why: a shape ratio such as
    // area/perimeter^2 falls continuously toward zero as a genuine lune
    // narrows, so any threshold on it silently rejects real trims at some
    // width. An absolute floor at the noise level rejects only what is
    // degenerate to floating point.
    //
    // Shoelace taken relative to the loop's own first vertex, for the
    // cancellation reason recorded on conway-geom#190: on a chart offset far
    // from the origin the raw sums lose the answer to rounding.
    bool enclosesArea = false;

    if ( distinctWelded.size() >= 3 && weldedSequence.size() >= 3 ) {

      const CDT::V2d< double >& reference = cdtVertices[ weldedSequence[ 0 ] ];

      const size_t n = weldedSequence.size();

      double twiceArea  = 0.0;
      double maxShifted = 0.0;

      for ( size_t where = 0; where < n; ++where ) {

        const CDT::V2d< double >& a =
          cdtVertices[ weldedSequence[ where ] ];
        const CDT::V2d< double >& b =
          cdtVertices[ weldedSequence[ ( where + 1 ) % n ] ];

        const double ax = a.x - reference.x;
        const double ay = a.y - reference.y;
        const double bx = b.x - reference.x;
        const double by = b.y - reference.y;

        maxShifted = std::max( maxShifted,
                       std::max( std::abs( ax ), std::abs( ay ) ) );

        twiceArea += ( ax * by ) - ( bx * ay );
      }

      // The threshold is the shoelace's OWN rounding error, not a constant -
      // see shoelaceAreaTolerance for the derivation and for why a fixed floor
      // cannot work on a per-face-normalized chart.
      //
      // On the collinear case this sum is EXACTLY 0.0, because constant-theta
      // points are exactly collinear in both layouts, so any positive
      // tolerance rejects it.
      const double areaTolerance = shoelaceAreaTolerance( n, maxShifted );

      enclosesArea = std::abs( twiceArea ) > areaTolerance;
    }

    if ( enclosesArea ) {
      ++closedLoops;
    }
  }

  if ( closedLoopsReachingCdt != nullptr ) {
    *closedLoopsReachingCdt = closedLoops;
  }

  if ( !inputFinite || cdtEdges.size() < 3 ) {
    return false;
  }

  CDT::Triangulation< double > triangulation(
    CDT::VertexInsertionOrder::Auto,
    CDT::IntersectingConstraintEdges::TryResolve, 0 );

  try
  {
    conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
    triangulation.insertVertices( cdtVertices );
    triangulation.insertEdges( cdtEdges );
    triangulation.eraseOuterTrianglesAndHoles();
  }
  // std::exception, not CDT::Error (which derives from it), because the
  // failure that actually costs a model here is std::bad_alloc, not a CDT
  // diagnostic: TryResolve splits intersections into fresh intersecting
  // edges and can allocate until the wasm heap is gone (conway#473 measured
  // ~3.4GB on one face). Narrowed to CDT::Error, that escaped the fallback
  // and, before AddFaceToGeometry's per-face guard, took the whole product
  // with it. Every site below is widened for the same reason.
  catch ( const std::exception &e )
  {
    Logger::logError( "CDT Exception (%s unwrap): %s", label, e.what() );
    return false;
  }

  mesh.vertices.reserve( cdtWorld.size() );

  for ( const glm::dvec3 &world : cdtWorld ) {
    mesh.vertices.push_back( world );
  }

  // No per-triangle orientation pass: appendMeshToGeometry normalizes
  // winding by dominant-plane projection either way.
  for ( const CDT::Triangle &triangle : triangulation.triangles ) {

    auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

    if ( cdtv1 == cdtv2 || cdtv2 == cdtv3 || cdtv3 == cdtv1 ) {
      continue;
    }

    // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
    if (
      cdtv1 >= cdtWorld.size() ||
      cdtv2 >= cdtWorld.size() ||
      cdtv3 >= cdtWorld.size() ) {
      continue;
    }

    mesh.makeTriangle( cdtv1, cdtv2, cdtv3 );
  }

  return !mesh.triangles.empty();
}

}  // namespace unwrap_detail

inline void TriangulateRevolution(Geometry &geometry,
                                  std::vector<IfcBound3D> &bounds,
                                  IfcSurface &surface,
                                  double representationExtent) {

  using namespace unwrap_detail;

  glm::dvec3 cent = surface.RevolutionSurface.Direction[3];
  glm::dvec3 vecX = glm::normalize(surface.RevolutionSurface.Direction[0]);
  glm::dvec3 vecY = glm::normalize(surface.RevolutionSurface.Direction[1]);
  glm::dvec3 vecZ = glm::normalize(surface.RevolutionSurface.Direction[2]);

  // Profile in the revolution frame as (radius, height) samples plus
  // cumulative arclength - shared by the trimmed unwrap, the swept-grid
  // fallback, and the adaptive refinement's closest-point projection.

  std::vector< glm::dvec2 > profileRZ;
  std::vector< double >     profileArc;

  profileRZ.reserve( surface.RevolutionSurface.Profile.curve.points.size() );
  profileArc.reserve( surface.RevolutionSurface.Profile.curve.points.size() );

  for ( const glm::dvec3 &point : surface.RevolutionSurface.Profile.curve.points ) {

    glm::dvec3 delta = point - cent;

    double dx = glm::dot( vecX, delta );
    double dy = glm::dot( vecY, delta );
    double dz = glm::dot( vecZ, delta );

    profileRZ.emplace_back( std::sqrt( dx * dx + dy * dy ), dz );
  }

  double profileLength = 0.0;

  profileArc.push_back( 0.0 );

  for ( size_t where = 1; where < profileRZ.size(); ++where ) {
    profileLength += glm::distance( profileRZ[ where - 1 ], profileRZ[ where ] );
    profileArc.push_back( profileLength );
  }

  if ( profileRZ.empty() ) {
    return;
  }

  // Closest point on the profile polyline to a (radius, height) query:
  // the snapped point and its normalized arclength parameter.
  auto closestOnProfile = [&]( const glm::dvec2 &query ) {

    glm::dvec2 best      = profileRZ[ 0 ];
    double     bestArc   = 0.0;
    double     bestDist2 = std::numeric_limits< double >::max();

    for ( size_t where = 0; where + 1 < profileRZ.size(); ++where ) {

      const glm::dvec2 &from = profileRZ[ where ];
      const glm::dvec2 &to   = profileRZ[ where + 1 ];

      glm::dvec2 segment = to - from;
      double     length2 = glm::dot( segment, segment );

      double t =
        length2 > 0 ?
          std::clamp( glm::dot( query - from, segment ) / length2, 0.0, 1.0 ) :
          0.0;

      glm::dvec2 candidate = from + segment * t;
      double     dist2     = glm::dot( query - candidate, query - candidate );

      if ( dist2 < bestDist2 ) {
        bestDist2 = dist2;
        best      = candidate;
        bestArc   = profileArc[ where ] + t * std::sqrt( length2 );
      }
    }

    double tNormalized = profileLength > 0 ? bestArc / profileLength : 0.0;

    return std::pair< glm::dvec2, double >( best, tNormalized );
  };

  // Adaptive-refinement projection (shared by both paths): keep the query
  // point's angle around the axis, snap its (radius, height) to the profile.
  auto surfaceProjection = [&]( const glm::dvec3 &point ) {

    glm::dvec3 delta = point - cent;

    double dx = glm::dot( vecX, delta );
    double dy = glm::dot( vecY, delta );
    double dz = glm::dot( vecZ, delta );

    double dd = std::sqrt( dx * dx + dy * dy );

    if ( dd < 1e-12 ) {
      return point;  // On the axis - angle undefined, nothing to round.
    }

    double sinAngle = dx / dd;
    double cosAngle = dy / dd;

    glm::dvec2 best = closestOnProfile( glm::dvec2( dd, dz ) ).first;

    return
      vecX * ( sinAngle * best.x ) +
      vecY * ( cosAngle * best.x ) +
      vecZ * best.y + cent;
  };

  auto refineAndAppend = [&]( WingedEdgeMesh< glm::dvec3 > &mesh ) {

    tesselate(
      mesh,
      surfaceProjection,
      mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
      relativeDeflectionSquared( mesh, representationExtent ) );

    appendMeshToGeometry( mesh, geometry );
  };

  // --- Trimmed path (#149): unwrap the actual boundary loops into
  // (theta, t) - theta the periodic angle around the axis, t the normalized
  // (non-periodic) profile arclength - and CDT them, so partial sweeps and
  // interior holes are honored instead of being papered over by a swept
  // patch. Any degeneracy (no usable loops, non-finite layout, intersecting
  // constraints) falls back to the swept grid below, which reproduces the
  // pre-trim behavior - a dirty boundary can never make output worse than it
  // was. A closed profile revolved is torus-like (t periodic); that rare
  // case stays on the fallback too.

  bool profileClosed =
    profileRZ.size() > 2 &&
    glm::distance( profileRZ.front(), profileRZ.back() ) < profileLength * 1e-9;

  if ( !profileClosed && profileRZ.size() >= 2 && profileLength > 0 ) {

    std::vector< std::vector< BoundaryPoint > > loops;

    loops.reserve( bounds.size() );

    auto parameterDuplicate = []( const BoundaryPoint &a, double theta, double t ) {

      return
        std::abs( wrapDeltaPi( theta - a.theta ) ) < 1e-12 &&
        std::abs( t - a.phi ) < 1e-12;
    };

    for ( const IfcBound3D &bound : bounds ) {

      std::vector< BoundaryPoint > loop;

      loop.reserve( bound.curve.points.size() );

      for ( const glm::dvec3 &point : bound.curve.points ) {

        glm::dvec3 delta = point - cent;

        double dx = glm::dot( vecX, delta );
        double dy = glm::dot( vecY, delta );
        double dz = glm::dot( vecZ, delta );

        if ( !std::isfinite( dx ) || !std::isfinite( dy ) || !std::isfinite( dz ) ) {
          continue;
        }

        double dd = std::sqrt( dx * dx + dy * dy );

        // On-axis sample - the angle is undefined, drop it.
        if ( dd < 1e-12 ) {
          continue;
        }

        // Grid/projection convention: local x = sin(theta) * d, y = cos(theta) * d.
        double theta = std::atan2( dx, dy );
        double t     = closestOnProfile( glm::dvec2( dd, dz ) ).second;

        if ( !loop.empty() && parameterDuplicate( loop.back(), theta, t ) ) {
          continue;
        }

        loop.push_back( BoundaryPoint{ point, theta, t } );
      }

      while ( loop.size() > 1 &&
              parameterDuplicate( loop.front(), loop.back().theta, loop.back().phi ) ) {
        loop.pop_back();
      }

      if ( loop.size() >= 3 ) {
        loops.push_back( std::move( loop ) );
      }
    }

    if ( !loops.empty() ) {

      bool windsTheta = false;

      std::vector< double > thetaSamples;

      for ( const std::vector< BoundaryPoint > &loop : loops ) {

        double sumTheta = 0.0;

        for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

          const BoundaryPoint &from = loop[ where ];
          const BoundaryPoint &to   = loop[ ( where + 1 ) % count ];

          sumTheta += wrapDeltaPi( to.theta - from.theta );

          thetaSamples.push_back( positiveMod2Pi( from.theta ) );
        }

        if ( std::llround( sumTheta / kTwoPi ) != 0 ) {
          windsTheta = true;
        }
      }

      auto [ thetaGap, thetaCut ] = largestCircularGap( thetaSamples );

      // A boundary that covers the whole circle without netting a winding
      // (seam edges walked both ways) leaves no gap for a branch cut.
      bool wrapsTheta = windsTheta || thetaGap < 1e-3;

      std::vector< std::vector< glm::dvec2 > > loops2D;

      bool layoutOk = false;

      // A rectangle edge jumping more than pi in cut-normalized theta
      // straddles the branch cut - upgrade to the wrapping (annulus) layout,
      // which has no cut to straddle, and retry once.
      for ( int attempt = 0; attempt < 2 && !layoutOk; ++attempt ) {

        loops2D.clear();
        loops2D.reserve( loops.size() );
        layoutOk = true;

        if ( wrapsTheta ) {

          // Annulus: angle = theta, radius = 1 + t in [1, 2].
          for ( const std::vector< BoundaryPoint > &loop : loops ) {

            std::vector< glm::dvec2 > loop2D;

            loop2D.reserve( loop.size() );

            for ( const BoundaryPoint &point : loop ) {
              loop2D.emplace_back(
                ( 1.0 + point.phi ) * std::cos( point.theta ),
                ( 1.0 + point.phi ) * std::sin( point.theta ) );
            }

            loops2D.push_back( std::move( loop2D ) );
          }

        } else {

          double loTheta = std::numeric_limits< double >::max();
          double hiTheta = std::numeric_limits< double >::lowest();

          std::vector< std::vector< double > > normalized;

          normalized.reserve( loops.size() );

          for ( const std::vector< BoundaryPoint > &loop : loops ) {

            std::vector< double > loopNormalized;

            loopNormalized.reserve( loop.size() );

            double previous = 0.0;

            for ( size_t where = 0; where < loop.size(); ++where ) {

              double value = positiveMod2Pi( loop[ where ].theta - thetaCut );

              if ( where > 0 && std::abs( value - previous ) > kPi ) {
                wrapsTheta = true;
                layoutOk   = false;
                break;
              }

              previous = value;
              loTheta  = std::min( loTheta, value );
              hiTheta  = std::max( hiTheta, value );
              loopNormalized.push_back( value );
            }

            if ( !layoutOk ) {
              break;
            }

            normalized.push_back( std::move( loopNormalized ) );
          }

          if ( layoutOk ) {

            double span = std::max( hiTheta - loTheta, 1e-12 );

            for ( size_t which = 0; which < loops.size(); ++which ) {

              std::vector< glm::dvec2 > loop2D;

              loop2D.reserve( loops[ which ].size() );

              for ( size_t where = 0; where < loops[ which ].size(); ++where ) {
                loop2D.emplace_back(
                  ( normalized[ which ][ where ] - loTheta ) / span,
                  loops[ which ][ where ].phi );
              }

              loops2D.push_back( std::move( loop2D ) );
            }
          }
        }
      }

      if ( layoutOk ) {

        // Weld 2D-coincident vertices and dedupe edges, exactly like the
        // toroidal unwrap - seam-doubled boundary edges become single
        // interior constraints, which is what makes NotAllowed safe.
        std::vector< CDT::V2d< double > > cdtVertices;
        std::vector< glm::dvec3 >         cdtWorld;
        std::vector< CDT::Edge >          cdtEdges;

        std::map< std::pair< long long, long long >, uint32_t > weld;
        std::set< std::pair< uint32_t, uint32_t > >             edgeSet;

        bool inputFinite = true;

        auto weldVertex = [&]( const glm::dvec2 &position, const glm::dvec3 &world ) {

          std::pair< long long, long long > key(
            std::llround( position.x * 1e9 ),
            std::llround( position.y * 1e9 ) );

          auto [ found, isNew ] =
            weld.try_emplace( key, static_cast< uint32_t >( cdtVertices.size() ) );

          if ( isNew ) {
            if ( !std::isfinite( position.x ) || !std::isfinite( position.y ) ) {
              inputFinite = false;
            }
            cdtVertices.emplace_back( position.x, position.y );
            cdtWorld.push_back( world );
          }

          return found->second;
        };

        for ( size_t which = 0; which < loops.size(); ++which ) {

          const std::vector< BoundaryPoint > &loop   = loops[ which ];
          const std::vector< glm::dvec2 >    &loop2D = loops2D[ which ];

          for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

            size_t next = ( where + 1 ) % count;

            uint32_t v1 = weldVertex( loop2D[ where ], loop[ where ].world );
            uint32_t v2 = weldVertex( loop2D[ next ], loop[ next ].world );

            if ( v1 == v2 ) {
              continue;
            }

            std::pair< uint32_t, uint32_t > ordered(
              std::min( v1, v2 ), std::max( v1, v2 ) );

            if ( edgeSet.insert( ordered ).second ) {
              cdtEdges.emplace_back( v1, v2 );
            }
          }
        }

        if ( inputFinite && cdtEdges.size() >= 3 ) {

          // TryResolve does not terminate on a malformed loop set. It splits
          // each intersection into new vertices whose edges intersect again,
          // allocating until the wasm heap is exhausted - measured on ISSUE_159
          // (conway#473) at ~30s and ~3.4GB for ONE face, which then emitted
          // nothing. Entry/exit markers around insertEdges show "begin" with no
          // matching "end"; every probe that logged after the call was blind to
          // it, since the call never returns.
          //
          // The cost is real and is NOT "never worse than before": TryResolve
          // handles ordinary crossings correctly and in microseconds, and every
          // such face now falls to the swept-grid path instead, which sweeps
          // the whole profile over the union angular bbox and so fills holes
          // and squares off non-rectangular trims. Measured on
          // FM_ARC_DigitalHub: 46 faces diverted, 99,641 -> 120,947 vertices.
          // visual-diff puts that below its 0.05% whole-model threshold, which
          // is weaker evidence than it sounds - those are small fittings on a
          // large building - but it is the evidence there is.
          //
          // Gating on loop-set shape was tried and does not work. Closed loops
          // contribute one edge per vertex, so `cdtEdges.size() ==
          // cdtVertices.size()` looked like it would separate ISSUE_159's
          // runaway (92 vertices / 94 edges) from FM_ARC's resolvable
          // crossings. Measured: it does not. FM_ARC's faces are malformed by
          // that test too, take the NotAllowed branch anyway, and its vertex
          // count is unchanged at 120,947 - so the gate bought nothing and was
          // removed rather than left as dead weight in a hot path.
          //
          // That is consistent with what CDT actually branches on: a
          // transversal crossing, `fixedEdges.count( Edge( iVL, iVR ) )`, which
          // edge/vertex parity neither implies nor excludes. Separating
          // "resolves in microseconds" from "never terminates" needs a real
          // bound inside the resolver, which CDT does not expose today.
          CDT::Triangulation< double > triangulation(
            CDT::VertexInsertionOrder::Auto,
            CDT::IntersectingConstraintEdges::NotAllowed, 0 );

          bool triangulated = false;

          try
          {
            conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
            triangulation.insertVertices( cdtVertices );
            triangulation.insertEdges( cdtEdges );
            triangulation.eraseOuterTrianglesAndHoles();
            triangulated = true;
          }
          catch ( const std::exception &e )
          {
            // Fall back to the swept grid - never worse than the old output.
            Logger::logError( "CDT Exception (revolution unwrap): %s", e.what() );
          }

          if ( triangulated ) {

            WingedEdgeMesh< glm::dvec3 > mesh;

            mesh.vertices.reserve( cdtWorld.size() );

            for ( const glm::dvec3 &world : cdtWorld ) {
              mesh.vertices.push_back( world );
            }

            // No per-triangle orientation pass: appendMeshToGeometry
            // normalizes winding by dominant-plane projection either way.
            for ( const CDT::Triangle &triangle : triangulation.triangles ) {

              auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

              if ( cdtv1 == cdtv2 || cdtv2 == cdtv3 || cdtv3 == cdtv1 ) {
                continue;
              }

              // Defensive, and unreachable as written: the guard dates from
              // when this site used TryResolve, which adds split vertices
              // past the input set that have no world-space lift in
              // cdtWorld. Under NotAllowed (above) a crossing throws instead
              // of splitting, so no such vertex can reach here. Kept because
              // the cost is three compares against a value already in
              // registers and the failure it prevents is an out-of-range
              // read; the sibling copy in unwrap_detail IS live, since that
              // path is still TryResolve.
              if (
                cdtv1 >= cdtWorld.size() ||
                cdtv2 >= cdtWorld.size() ||
                cdtv3 >= cdtWorld.size() ) {
                continue;
              }

              mesh.makeTriangle( cdtv1, cdtv2, cdtv3 );
            }

            if ( !mesh.triangles.empty() ) {
              refineAndAppend( mesh );
              return;
            }
          }
        }
      }
    }
  }

  // --- Fallback: sweep the whole profile across the boundary's angular
  // bounding box (the pre-trim approximation).

  std::vector<std::vector<glm::dvec3>> newPoints;

  std::vector<glm::dvec3> bounding;
  std::vector<double> angleVec;
  std::vector<double> angleDsp;

  // Now we construct the bounding box of the boundary ...
  // ... by adding the middle point of all curves

  for (size_t i = 0; i < bounds.size(); i++) {

    const std::vector< glm::dvec3 > &curvePoints  = bounds[ i ].curve.points;
    const std::vector< uint16_t >   &curveIndices = bounds[ i ].curve.indices;

    if ( curvePoints.empty() ) {
      continue;
    }

    // indices is only populated on the edge-loop branch of GetLoop, so poly
    // loops and vertex loops arrive points-only and every read below - starting
    // with indices[0], before the loop - was out of bounds for them. Absent
    // indices are read as a single group, which is what the surrounding code
    // already does for an edge loop whose points share one edge id.
    //
    // A size match means the indices are READABLE, not that they still line up:
    // GetBound and createBound3D reverse curve.points for a .F.-oriented bound
    // and leave curve.indices alone, so those bounds group points against
    // another edge's id. Pre-existing, and untouched here.
    const bool grouped = curveIndices.size() == curvePoints.size();

    double xx = 0;
    double yy = 0;
    double zz = 0;
    double cc = 0;
    int lastTeam = grouped ? static_cast< int >( curveIndices[ 0 ] ) : 0;
    for (size_t j = 0; j < curvePoints.size(); j++) {
      const int team = grouped ? static_cast< int >( curveIndices[ j ] ) : 0;

      // If it is the first point of the group we close the previous group ...
      //  ... and create a new one. Else, the point is of the current group
      if (lastTeam != team ||
          j == (curvePoints.size() - 1)) {
        if (cc > 0) {
          xx /= cc;
          yy /= cc;
          zz /= cc;
          bounding.push_back(glm::dvec3(xx, yy, zz));
        }
        xx = curvePoints[j].x;
        yy = curvePoints[j].y;
        zz = curvePoints[j].z;
        cc = 1;

        lastTeam = team;
      } else {
        xx += curvePoints[j].x;
        yy += curvePoints[j].y;
        zz += curvePoints[j].z;
        cc++;
      }
    }
  }

  // There is a problem when points in the revolution are around 0 degrees
  // Numerical instabilities can make these points to jump from 0 to 360
  // It causes lots of trouble when drawing the boundaries in the revolution

  // The method presented here finds the angle of each point, measures the ...
  //  ... angular difference and then, if the difference is bigger than 180 ...
  //  ... corrects it to a lesser value. Finally it gets the first angle and ...
  //  ... adds the angular differences again, reconstructing a corrected
  //  boundary.

  // Now we find the angle of each point in the reference plane of the cylinder
  for (size_t j = 0; j < bounding.size(); j++) {
    double xx = bounding[j].x - cent.x;
    double yy = bounding[j].y - cent.y;
    double zz = bounding[j].z - cent.z;
    double dx = vecX.x * xx + vecX.y * yy + vecX.z * zz;
    double dy = vecY.x * xx + vecY.y * yy + vecY.z * zz;
    //				double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z *
    // zz;
    double temp = VectorToAngle(dx, dy);
    while (temp < 0) {
      temp += 360;
    }
    while (temp > 360) {
      temp -= 360;
    }
    angleVec.push_back(temp);
  }

  // angleVec.size() - 1 is unsigned, so on an empty vector it is SIZE_MAX, not
  // -1, and the loop below would push_back into angleDsp until the wasm heap
  // was exhausted. angleVec is empty when every bound carried no points; the
  // angleVec[0] reads just past the loop would be out of bounds in that case
  // too. Returning silently matches the sibling bail-outs here -
  // warnFaceAddedNothing downstream already reports a face that emits nothing.
  //
  // Scope note: this is a memory-safety fix only. The degenerate-but-bounded
  // outputs this function can still produce (a single angular sample sweeping
  // ten coincident rings, NaN angles from on-axis samples, a VERTEX_LOOP-
  // bounded full revolution emitting nothing where #461 gives the sphere a
  // full parametric grid) are pre-existing and each needs its own change and
  // fixture - attempts to fix them here kept introducing new erasures.
  // Below TWO, not merely zero. One sample leaves startDegrees == endDegrees,
  // so radSpan is 0 and the sweep lays ten coincident rings of zero-area
  // triangles - which appendMeshToGeometry appends unfiltered, suppressing
  // warnFaceAddedNothing so the face vanishes with nothing reported. Bailing
  // gives the same (empty) picture and lets that warning fire. The points-only
  // branch above produces exactly one sample, so this is its normal exit.
  constexpr size_t MINIMUM_ANGULAR_SAMPLES = 2;

  if ( angleVec.size() < MINIMUM_ANGULAR_SAMPLES ) {
    return;
  }

  for (size_t i = 0; i + 1 < angleVec.size(); i++) {
    if (angleVec[i] - angleVec[i + 1] > 180) {
      angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
    } else if (angleVec[i] - angleVec[i + 1] < -180) {
      angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
    } else {
      angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
    }
  }

  double startDegrees = angleVec[0];
  double endDegrees = angleVec[0];

  // Add angular differences starting from the first angle. We also correct the
  // start and end angles

  double temp = angleVec[0];
  for (size_t i = 0; i < angleDsp.size(); i++) {
    temp += angleDsp[i];
    if (endDegrees < temp) {
      endDegrees = temp;
    }
    if (startDegrees > temp) {
      startDegrees = temp;
    }
  }

  double startRad = startDegrees / 180 * (double)CONST_PI;
  double endRad   = endDegrees / 180 * (double)CONST_PI;
  double radSpan  = endRad - startRad;

  // The span above is measured over one CENTROID per edge, and a centroid
  // carries no angular information exactly when the edge does: a full circle
  // around the axis - a torus's equator, the natural boundary of any closed
  // profile revolved a full turn - has its centroid ON the axis, where the
  // angle is undefined. Every such face measured a span of zero here, swept a
  // grid of coincident rings, and was welded away to nothing at reification,
  // with warnFaceAddedNothing suppressed because 414 (zero-area) triangles had
  // been appended (conway#461: the torus's revolution half).
  //
  // Recover from the BOUNDARY POINTS instead, and only when the centroid span
  // collapsed - every face that produced output before produces byte-identical
  // output still, the same no-regression shape as the sphere fix (#160).
  // The covered arc is the complement of the largest angular gap between
  // boundary samples; a boundary whose points also carry no spread is the
  // degenerate-bound convention (the bound is the seam/profile itself), which
  // ISO 10303-42 reads as covering the full surface - so both cases sweep on,
  // the first over its measured arc, the second over the whole turn.
  //
  // A closed profile revolved is torus-like, and its bounds usually trim the
  // TUBE rather than the sweep: the same boundary points, mapped through
  // closestOnProfile, say which arc of the profile the face actually covers.
  // Restrict the swept rows to that arc - without it the recovery would sweep
  // the whole tube and double-cover the area the face's sibling already owns.
  const std::vector< glm::dvec2 >* sweepProfile = &profileRZ;
  std::vector< glm::dvec2 > restrictedProfile;

  constexpr double MINIMUM_RECOVERED_SPAN = 1e-6;

  if ( radSpan < MINIMUM_RECOVERED_SPAN ) {

    std::vector< double > thetaSamples;
    std::vector< double > tSamples;

    // A bound that WINDS the axis - an equator circle - covers the full turn
    // by topology, and the gap measure below cannot say so on its own: N
    // samples on a full circle always leave a largest gap of about 2*pi/N, so
    // trusting the gap alone would sweep (N-1)/N of a turn and leave an open
    // pie-slice wedge.
    //
    // The measure is the EXCURSION of the cumulative winding along each bound,
    // not its net: a torus face's single loop walks one equator forward and
    // the other back, so the net cancels to zero while the excursion reaches a
    // full turn on the way - and a partial sweep's excursion is exactly its
    // swept angle, so the same number serves both.
    bool windsFullTurn = false;

    for ( const IfcBound3D& bound : bounds ) {

      double winding       = 0.0;
      double windingLow    = 0.0;
      double windingHigh   = 0.0;
      bool   havePrevious  = false;
      double previousTheta = 0.0;

      for ( const glm::dvec3& point : bound.curve.points ) {

        glm::dvec3 delta = point - cent;

        double dx = glm::dot( vecX, delta );
        double dy = glm::dot( vecY, delta );
        double dz = glm::dot( vecZ, delta );

        if ( !std::isfinite( dx ) || !std::isfinite( dy ) || !std::isfinite( dz ) ) {
          continue;
        }

        double dd = std::sqrt( dx * dx + dy * dy );

        // On-axis sample - the angle is undefined, drop it.
        if ( dd < 1e-12 ) {
          continue;
        }

        double theta = std::atan2( dx, dy );

        thetaSamples.push_back( positiveMod2Pi( theta ) );

        if ( havePrevious ) {
          winding    += wrapDeltaPi( theta - previousTheta );
          windingLow  = std::min( windingLow, winding );
          windingHigh = std::max( windingHigh, winding );
        }

        previousTheta = theta;
        havePrevious  = true;

        if ( profileLength > 0 ) {
          tSamples.push_back( positiveMod2Pi(
            closestOnProfile( glm::dvec2( dd, dz ) ).second * kTwoPi ) );
        }
      }

      // One sample short of closing still counts: the loop's last point is
      // its first, so an excursion within one segment of a full turn is one.
      if ( windingHigh - windingLow > kTwoPi * ( 63.0 / 64.0 ) ) {
        windsFullTurn = true;
      }
    }

    if ( !thetaSamples.empty() ) {

      auto [ thetaGap, thetaMid ] = largestCircularGap( thetaSamples );

      double coveredSpan = kTwoPi - thetaGap;

      if ( windsFullTurn || coveredSpan <= MINIMUM_RECOVERED_SPAN ) {

        // Winding is full coverage by topology. No spread at all is the
        // degenerate-bound convention - the bound is the seam or the profile
        // itself, which ISO 10303-42 reads as covering the whole surface;
        // the sphere's VERTEX_LOOP (#160) is the same convention one
        // dimension down. A genuine sliver sweep whose boundary spread is
        // below MINIMUM_RECOVERED_SPAN is indistinguishable from the seam
        // case from theta alone and sweeps the full turn too - the
        // convention is real and observed, the tolerance-sliver artifact is
        // not, so the tie breaks toward the convention.
        startRad = 0.0;
        radSpan  = kTwoPi;
      } else {
        startRad = positiveMod2Pi( thetaMid + thetaGap * 0.5 );
        radSpan  = coveredSpan;
      }

      // An uncovered stretch of the tube smaller than this is sampling
      // granularity, not a trim: boundary samples land on the profile about
      // one profile-segment apart, so a genuine trim's gap (the whole missing
      // arc) dwarfs it while full coverage never reaches it.
      constexpr double MINIMUM_PROFILE_TRIM_GAP = kTwoPi / 8.0;

      if ( profileClosed && tSamples.size() >= 2 ) {

        auto [ tGap, tMid ] = largestCircularGap( tSamples );

        // Bounded on both ends. A gap below the floor is sampling
        // granularity, not a trim. A gap near the whole circle means the
        // bounds sit at ONE tube location - a single equator circle as the
        // closed surface's seam - which is "no trim information", not "the
        // face covers nothing": restricting to it would sweep a zero-width
        // ribbon and vanish the face all over again.
        const double maximumTrimGap = kTwoPi - MINIMUM_PROFILE_TRIM_GAP;

        if ( tGap > MINIMUM_PROFILE_TRIM_GAP && tGap < maximumTrimGap ) {

          // Covered arc of the profile, as arclengths along profileRZ. The
          // closing duplicate point is skipped so it cannot appear twice.
          double coveredStart =
            positiveMod2Pi( tMid + tGap * 0.5 ) / kTwoPi * profileLength;
          double coveredLength =
            ( kTwoPi - tGap ) / kTwoPi * profileLength;

          auto pointAtArc = [&]( double arc ) -> glm::dvec2 {

            arc = std::fmod( arc, profileLength );

            if ( arc < 0 ) {
              arc += profileLength;
            }

            size_t high = 1;

            while ( high + 1 < profileArc.size() && profileArc[ high ] < arc ) {
              ++high;
            }

            size_t low     = high - 1;
            double segment = profileArc[ high ] - profileArc[ low ];
            double blend   = segment > 0 ? ( arc - profileArc[ low ] ) / segment : 0;

            return glm::mix( profileRZ[ low ], profileRZ[ high ], blend );
          };

          restrictedProfile.push_back( pointAtArc( coveredStart ) );

          std::vector< std::pair< double, size_t > > interior;

          for ( size_t where = 0; where + 1 < profileRZ.size(); ++where ) {

            double lifted = profileArc[ where ] < coveredStart ?
              profileArc[ where ] + profileLength : profileArc[ where ];

            if ( lifted > coveredStart + profileLength * 1e-9 &&
                 lifted < coveredStart + coveredLength - profileLength * 1e-9 ) {
              interior.emplace_back( lifted, where );
            }
          }

          std::sort( interior.begin(), interior.end() );

          for ( const auto& [ lifted, where ] : interior ) {
            restrictedProfile.push_back( profileRZ[ where ] );
          }

          restrictedProfile.push_back( pointAtArc( coveredStart + coveredLength ) );

          if ( restrictedProfile.size() >= 2 ) {
            sweepProfile = &restrictedProfile;
          }
        }
      }
    }
  }

  // Ring count adapts to the swept span: 64 segments per full turn (sagitta
  // ~0.12% of radius) instead of the old fixed 10 rings, whose ~40-degree
  // facets read as a polygonal shaft on the Jetenginestep turbine (#149).
  // The floor keeps small spans at least as dense as before. This sets the
  // resolution of the grid BORDER (the end rings/profile rows): the interior
  // is further refined by the adaptive tesselate below, but border edges are
  // never subdivided, so they must start smooth enough on their own.
  constexpr double FULL_TURN_SEGMENTS = 64.0;

  int numRots = std::clamp(
    static_cast< int >( std::ceil(
      radSpan / ( 2.0 * (double)CONST_PI / FULL_TURN_SEGMENTS ) ) ) + 1,
    10,
    static_cast< int >( FULL_TURN_SEGMENTS ) + 1 );

  for (int r = 0; r < numRots; r++) {
    std::vector<glm::dvec3> newList;
    newPoints.push_back(newList);
  }

  double radStep  = radSpan / (numRots - 1);

  for ( const glm::dvec2 &rz : *sweepProfile ) {
    for (int r = 0; r < numRots; r++) {
      double angle = startRad + r * radStep;
      double dtempX = sin(angle) * rz.x;
      double dtempY = cos(angle) * rz.x;
      newPoints[r].push_back(
        vecX * dtempX + vecY * dtempY + vecZ * rz.y + cent );
    }
  }

  if ( newPoints[ 0 ].empty() ) {
    return;
  }

  WingedEdgeMesh< glm::dvec3 > mesh;

  size_t profileCount = newPoints[ 0 ].size();

  for ( int r = 0; r < numRots; r++ ) {
    for ( size_t s = 0; s < profileCount; s++ ) {
      mesh.vertices.push_back( newPoints[ r ][ s ] );
    }
  }

  for ( int r = 0; r < numRots - 1; r++ ) {

    uint32_t row0 = static_cast< uint32_t >( r * profileCount );
    uint32_t row1 = static_cast< uint32_t >( ( r + 1 ) * profileCount );

    for ( size_t s = 0; s + 1 < profileCount; s++ ) {

      uint32_t a = row0 + static_cast< uint32_t >( s );
      uint32_t b = row0 + static_cast< uint32_t >( s + 1 );
      uint32_t c = row1 + static_cast< uint32_t >( s );
      uint32_t d = row1 + static_cast< uint32_t >( s + 1 );

      mesh.makeTriangle( a, b, c );
      mesh.makeTriangle( c, b, d );
    }
  }

  refineAndAppend( mesh );
}

/**
 * The two halves of the sphere's dual parameterization, with the
 * projection centre filled in.
 *
 * Each half is `direction(xy) * radius`, and the radius factor vanishes at
 * the half's own centre - the point the chart projects FROM the opposite
 * pole onto the origin. That makes the centre a removable singularity: the
 * limit is the origin from every approach direction. `glm::normalize` of a
 * zero-length vector is NaN, though, so the value has to be written down
 * rather than computed.
 *
 * Without this, any spherical face with a boundary VERTEX exactly on the
 * surface's own axis produced NaN parameter coordinates and was dropped
 * whole by the non-finite CDT guard in manifold_utils.h - "conway:
 * non-finite dual-parameterization (projection pole); dropping face". Every
 * corner blend of a filleted box is such a face: an eighth-sphere whose
 * three bounding arcs meet at the three axis points, one of which is the
 * pole (conway-geom#171).
 *
 * Note the caller's `2.0 - z < DBL_EPSILON` test handles the OTHER pole -
 * the true singularity of the chart, sent to a far sentinel - and must stay
 * ahead of this, which is why this takes the already-computed radius rather
 * than recomputing it.
 *
 * @param normalFormVertex point on the unit sphere in the surface's frame.
 * @param radius the half's radius factor, 1 +/- normalFormVertex.z.
 * @return the planar parameter coordinate.
 */
inline glm::dvec2 poleSafeStereographic(
  const glm::dvec3& normalFormVertex,
  double radius ) {

  const glm::dvec2 planar( normalFormVertex );

  // Exactly zero, not a tolerance: any non-zero planar component gives
  // normalize() a finite direction, and the result is already continuous
  // into the origin because `radius` is vanishing alongside it. Widening
  // this to an epsilon would move points that are currently correct.
  if ( planar.x == 0.0 && planar.y == 0.0 ) {

    return glm::dvec2( 0.0 );
  }

  return glm::normalize( planar ) * radius;
}

inline void TriangulateSphericalSurface(Geometry &geometry,
                                        const std::vector<IfcBound3D> &bounds,
                                        IfcSurface &surface,
                                        double representationExtent) {
  if ( bounds.empty() ) {
    return;
  }

  double     radius = surface.SphericalSurface.Radius;
  glm::dvec3 cent   = surface.transformation[3];
  glm::dvec3 vecX   = glm::normalize(surface.transformation[0]);
  glm::dvec3 vecY   = glm::normalize(surface.transformation[1]);
  glm::dvec3 vecZ   = glm::normalize(surface.transformation[2]);

  // Sphere point in world space from (theta, phi): theta is longitude about
  // vecZ, phi is latitude in [-pi/2, +pi/2].
  auto spherePoint = [&]( double theta, double phi ) {

    double ring = radius * std::cos( phi );

    glm::dvec3 local(
      ring * std::cos( theta ),
      ring * std::sin( theta ),
      radius * std::sin( phi ) );

    return vecX * local.x + vecY * local.y + vecZ * local.z + cent;
  };

  // ISO 10303-42: a face bounded solely by a degenerate loop (a VERTEX_LOOP,
  // whose single vertex marks a pole rather than a trim) covers the WHOLE
  // surface. That is how OCCT writes a plain sphere: one ADVANCED_FACE on a
  // SPHERICAL_SURFACE, one FACE_BOUND, one VERTEX_LOOP at the north pole.
  // (Not how it writes fillet corner blends - those are trimmed patches with
  // real edge loops and never reach this branch.) The dual-hemisphere unwrap
  // below needs a real boundary polygon to constrain a CDT against, so one
  // point walks out as zero triangles and the sphere loads as an empty mesh
  // - silently, since nothing downstream distinguishes "trimmed away to
  // nothing" from "never triangulated".
  // See https://github.com/bldrs-ai/conway/issues/461.
  //
  // The torus already handles its own full-coverage case with a parametric
  // grid (see TriangulateToroidalSurface's wrapsTheta && wrapsPhi branch);
  // this is the same treatment for the sphere's closed domain.
  // Test the SINGLE bound, not a sum over all of them. "Bounded solely by a
  // degenerate loop" means exactly one bound, and summing loses that: a face
  // carrying a VERTEX_LOOP plus an EDGE_LOOP whose edges all failed to yield a
  // basis curve (GetLoop emits a 0-point curve in that case) sums to 1 and
  // would take this branch, which is precisely the extraction-failure case the
  // next paragraph rules out.
  //
  // One or two points cannot enclose area on any surface, so there is nothing
  // to trim with even when the loop is not literally a VERTEX_LOOP.
  //
  // ZERO points is deliberately excluded for that same reason: treating a
  // failed extraction as full coverage would turn a small trimmed patch into a
  // whole sphere at the placement centre - spurious volume, inflated bounds,
  // occlusion - and do it silently, since emitting triangles suppresses the
  // contributed-no-geometry warning. An extraction failure has to keep looking
  // like one.
  constexpr size_t MINIMUM_TRIM_POINTS = 3;

  // Multiple bounds fall through to the unwrap by reporting a count that
  // cannot be degenerate.
  const size_t boundaryPointCount =
    bounds.size() == 1 ? bounds[ 0 ].curve.points.size() : MINIMUM_TRIM_POINTS;

  // The OTHER spelling of "this face is the whole surface": a SEAM loop.
  //
  // The VERTEX_LOOP case above is how OCCT writes a plain sphere. ISO
  // 10303-42 also permits the seam form, and this file uses it - a single
  // EDGE_LOOP walking ONE great circle forward and then reversed:
  //
  //   #50626 = ADVANCED_FACE( '', ( #6222 ), #238, .T. )   <- SPHERICAL_SURFACE r=0.9
  //     EDGE_LOOP #9128:  ORIENTED_EDGE #41053 .T.  EDGE_CURVE #28750
  //                       ORIENTED_EDGE #41054 .F.  EDGE_CURVE #28750
  //
  // That loop retraces itself, so it encloses no area and cannot trim
  // anything; like the VERTEX_LOOP it means full coverage. But it arrives as
  // 47 points rather than 1, so the point-count test above does not see it,
  // and it fell through to the dual-hemisphere unwrap. There the failure is
  // total and silent: the retraced curve is a MERIDIAN, running pole to pole
  // (measured: the boundary's local z spans -1.0 to +1.0, the only spherical
  // face in the corpus that does), so the hemisphere classifier splits it
  // into two pole-to-pole arcs, neither of which encloses area in its chart,
  // and both CDTs return nothing. `_MR148ZZ Ball` on
  // step/conor/Orbiter_v1.1_Gear_7.5.step - solid 960, the bearing ball -
  // came out as 24 vertices and 0 triangles, i.e. absent
  // (bldrs-ai/conway#595, bldrs-ai/test-models-private#93).
  //
  // `bound.seam` is decided by the front end, where the ORIENTED_EDGEs are
  // still visible, and is a TOPOLOGICAL fact: every edge of the loop
  // traversed exactly twice, once in each direction.
  //
  // It deliberately replaces a geometric predicate. The first version of this
  // fix tested the loop's own vector area normalised by perimeter squared,
  // which separated the corpus by fourteen orders of magnitude - and was
  // still wrong, because it identifies LOW ENCLOSED AREA rather than
  // retracing. A thin spherical lune's area/perimeter^2 falls continuously
  // toward zero as its angular width shrinks (for width w it is about
  // w / (2*pi^2)), so below some width a GENUINE narrow trim would have been
  // silently replaced by the entire sphere - inflated bounds, spurious
  // volume, occlusion, and no warning, because emitting triangles suppresses
  // the contributed-no-geometry message. That the corpus never came within
  // five decades of the threshold was a fact about the corpus, not about the
  // code. Found by codex on bldrs-ai/conway-geom#187.
  //
  // The topological test cannot fail that way at any width: a lune is bounded
  // by two DIFFERENT edge curves, so no width makes it a retrace. It also
  // needs no tolerance, which the area test did.
  const bool retracingSeam =
    bounds.size() == 1 &&
    bounds[ 0 ].seam &&
    boundaryPointCount >= MINIMUM_TRIM_POINTS;

  if ( ( boundaryPointCount > 0 && boundaryPointCount < MINIMUM_TRIM_POINTS ) ||
       retracingSeam ) {

    // Longitude x latitude divisions. Matched to the torus grid's angular
    // density; the sphere is closed in theta but not in phi, so the poles
    // are single vertices with triangle fans rather than another ring.
    constexpr int GRID_THETA = 48;
    constexpr int GRID_PHI   = 24;

    constexpr double kHalfPi = unwrap_detail::kPi * 0.5;

    WingedEdgeMesh< glm::dvec3 > gridMesh;

    auto& gridVertices = gridMesh.vertices;

    // Layout: south pole, then GRID_PHI - 1 interior rings of GRID_THETA,
    // then north pole.
    gridVertices.push_back( spherePoint( 0.0, -kHalfPi ) );

    for ( int ring = 1; ring < GRID_PHI; ++ring ) {

      double phi = -kHalfPi + ( ring * unwrap_detail::kPi / GRID_PHI );

      for ( int step = 0; step < GRID_THETA; ++step ) {
        gridVertices.push_back( spherePoint( step * unwrap_detail::kTwoPi / GRID_THETA, phi ) );
      }
    }

    gridVertices.push_back( spherePoint( 0.0, kHalfPi ) );

    const uint32_t southPole = 0;
    const uint32_t northPole = static_cast< uint32_t >( gridVertices.size() - 1 );

    // First vertex of interior ring `ring` (1-based, matching the loop above).
    auto ringBase = []( int ring ) {

      return static_cast< uint32_t >( 1 + ( ( ring - 1 ) * GRID_THETA ) );
    };

    // Winding below is outward everywhere: cross(d/dtheta, d/dphi) points
    // away from the centre, so (lower_i, lower_i+1, upper_i+1) and
    // (lower_i, upper_i+1, upper_i) are both front-facing. The pole fans
    // follow from the same rule with one ring collapsed - which is why the
    // south cap runs i+1 before i and the north cap does not.
    for ( int step = 0; step < GRID_THETA; ++step ) {

      uint32_t next = static_cast< uint32_t >( ( step + 1 ) % GRID_THETA );

      gridMesh.makeTriangle(
        southPole,
        ringBase( 1 ) + next,
        ringBase( 1 ) + static_cast< uint32_t >( step ) );
    }

    for ( int ring = 1; ring + 1 < GRID_PHI; ++ring ) {

      uint32_t lower = ringBase( ring );
      uint32_t upper = ringBase( ring + 1 );

      for ( int step = 0; step < GRID_THETA; ++step ) {

        uint32_t here = static_cast< uint32_t >( step );
        uint32_t next = static_cast< uint32_t >( ( step + 1 ) % GRID_THETA );

        gridMesh.makeTriangle( lower + here, lower + next, upper + next );
        gridMesh.makeTriangle( lower + here, upper + next, upper + here );
      }
    }

    for ( int step = 0; step < GRID_THETA; ++step ) {

      uint32_t next = static_cast< uint32_t >( ( step + 1 ) % GRID_THETA );

      gridMesh.makeTriangle(
        northPole,
        ringBase( GRID_PHI - 1 ) + static_cast< uint32_t >( step ),
        ringBase( GRID_PHI - 1 ) + next );
    }

    // Orient against the true radial normal rather than the projection
    // heuristic, so the grid keeps the winding built above. Extractors that
    // never populate the face sense (IFC - see IfcSurface::sameSenseKnown)
    // get outward, which is the only defensible default for a closed sphere
    // and is strictly more information than the `false` the field defaults
    // to; this path emitted nothing at all before, so there is no prior
    // behaviour to preserve.
    appendMeshToGeometry(
      gridMesh,
      geometry,
      surface.sameSenseKnown ? surface.sameSense : true,
      [&]( const glm::dvec3& point ) {

        return glm::normalize( point - cent );
      } );

    return;
  }

  // --- Trimmed unwrap (bldrs-ai/conway#644): lay the boundary loops out in a
  // single injective (theta, phi) chart and CDT them once, the same treatment
  // the cylinder and cone already use via triangulateUnwrappedLoops. The
  // legacy dual-hemisphere stereographic path below stays as the fallback, so
  // any degeneracy here lands on exactly today's behaviour.
  //
  // WHY THE HEMISPHERE SPLIT FAILS, since it is not obvious from reading it.
  // That path cuts the boundary at the placement equator (normalFormVertex.z
  // <= 0), stereographically projects each half, and CDTs each separately.
  // A SPHERICAL_SURFACE's placement axis is ARBITRARY in the file, and it is
  // routinely unrelated to the part: on Orbiter's button-head dome
  // (ADVANCED_FACE #50714, express solid #964) it lands roughly perpendicular
  // to the screw, so a band 1.4 units deep on a 3.17 sphere has a boundary
  // whose normalized z spans [-0.8985, +0.8901] and is bisected LENGTHWISE by
  // an equator the geometry has no reason to respect.
  //
  // What that produces is not a clean loss of one half. Both half-CDTs run and
  // both return triangles; what differs is `edgeDiscarded`, which selects the
  // erase mode, whether the per-triangle equator filter runs, and whether
  // equator edges are collected - three things at once. On #50714 side 0
  // discards one edge of 46 and keeps 31 of 71 triangles; side 1 discards
  // NOTHING, so it takes eraseOuterTrianglesAndHoles() with no equator filter
  // and emits all 46, including triangles spanning the other hemisphere. The
  // face ships 244 triangles carrying 25.196 of a 27.95 analytic zone area
  // inside 5 of 24 azimuth sectors - about 90% of the right area in about 21%
  // of the right footprint. Past roughly 0.95 normalized z the same mechanism
  // stops emitting anything at all. Both regimes are pinned in
  // test/spherical_trim_test.cpp.
  //
  // Forcing both sides into the filtered branch was measured and makes it
  // strictly worse: the face then emits ZERO geometry. The predicate is load
  // bearing three ways, which is why this is a new chart rather than a repair
  // inside it.
  //
  // The branch cut here is placed by largestCircularGap in an empty angular
  // gap of THIS face's own boundary, or the annulus layout is used when the
  // boundary winds theta and there is no gap to place one in - which is the
  // property the fixed equator lacks. Hole nesting is exact, so an inner trim
  // (the dome's hex socket) stays a hole; #595's full-coverage grid could
  // never be widened to cover this case for exactly that reason - it would
  // pave the socket over.
  {
    using namespace unwrap_detail;

    // phi is latitude, non-periodic, normalized to [0, 1] over the face's own
    // span - the same contract the cylinder's axial coordinate satisfies.
    double loPhi = std::numeric_limits< double >::max();
    double hiPhi = std::numeric_limits< double >::lowest();

    // A boundary sample at the placement pole has no defined theta. The chart
    // cannot represent it, so the face falls through to the legacy path rather
    // than being laid out on a guess.
    bool sawPolePoint = false;

    std::vector< std::vector< BoundaryPoint > > loops;

    loops.reserve( bounds.size() );

    for ( const IfcBound3D &bound : bounds ) {

      std::vector< BoundaryPoint > loop;

      loop.reserve( bound.curve.points.size() );

      for ( const glm::dvec3 &point : bound.curve.points ) {

        glm::dvec3 delta = point - cent;

        double dx = glm::dot( vecX, delta );
        double dy = glm::dot( vecY, delta );
        double dz = glm::dot( vecZ, delta );

        if ( !std::isfinite( dx ) || !std::isfinite( dy ) ||
             !std::isfinite( dz ) ) {
          continue;
        }

        double ring = std::sqrt( dx * dx + dy * dy );

        if ( ring < radius * 1e-9 ) {
          sawPolePoint = true;
          continue;
        }

        // Latitude from the ring radius and the axial height, rather than
        // asin(dz / radius): the boundary points come off a trim and need not
        // sit exactly on the nominal radius, and atan2 stays conditioned where
        // the asin form loses precision near the poles.
        double phi = std::atan2( dz, ring );

        loPhi = std::min( loPhi, phi );
        hiPhi = std::max( hiPhi, phi );

        loop.push_back( BoundaryPoint{ point, std::atan2( dy, dx ), phi } );
      }

      if ( loop.size() >= 3 ) {
        loops.push_back( std::move( loop ) );
      }
    }

    double spanPhi = hiPhi - loPhi;

    if ( !sawPolePoint && !loops.empty() && spanPhi > 1e-12 ) {

      auto parameterDuplicate = []( const BoundaryPoint &a,
                                    const BoundaryPoint &b ) {

        return
          std::abs( wrapDeltaPi( b.theta - a.theta ) ) < 1e-12 &&
          std::abs( b.phi - a.phi ) < 1e-12;
      };

      for ( std::vector< BoundaryPoint > &loop : loops ) {

        for ( BoundaryPoint &point : loop ) {
          point.phi = ( point.phi - loPhi ) / spanPhi;
        }

        std::vector< BoundaryPoint > cleaned;

        cleaned.reserve( loop.size() );

        for ( const BoundaryPoint &point : loop ) {
          if ( cleaned.empty() || !parameterDuplicate( cleaned.back(), point ) ) {
            cleaned.push_back( point );
          }
        }

        while ( cleaned.size() > 1 &&
                parameterDuplicate( cleaned.front(), cleaned.back() ) ) {
          cleaned.pop_back();
        }

        loop = std::move( cleaned );
      }

      // ONE invariant, checked twice because a loop can be lost in two
      // different places: EVERY BOUND SUPPLIED MUST REACH THE TRIANGULATION.
      // Losing one and triangulating the survivors paves an inner trim's hole
      // over and adds volume silently - the failure that ruled #595's
      // full-coverage grid out of this job in the first place, since it would
      // have replaced the dome's hex socket with surface.
      //
      // Where a loop can vanish:
      //
      //   1. AT COLLECTION. A non-finite sample is skipped, and a loop left
      //      under three points is never pushed into `loops` at all. Counting
      //      collapses among the loops that survived cannot see this, because
      //      the casualty is already gone - so compare against bounds.size().
      //
      //   2. AT THE CONSTRAINT WELD, inside triangulateUnwrappedLoops. It
      //      welds at 1e-9 while the dedup above works at 1e-12, so points
      //      that legitimately survive here can merge there. A loop merged
      //      onto ONE vertex emits no edges and disappears from the
      //      constraints; a loop merged onto TWO emits a single open edge,
      //      which is a slit rather than a boundary and gets triangulated
      //      straight across. Either way `cdtEdges.size() < 3` still passes on
      //      the outer loop alone. Only the helper can see this, so it reports
      //      how many loops survived as CLOSED constraints - three distinct
      //      welded vertices, the least that can enclose area - and the
      //      equality is checked after it returns.
      //
      // Failing any of these sends the face down the legacy path exactly as if
      // this code were not here, which is what keeps the "can never make
      // output worse" contract literally true rather than true-except-here.
      //
      // Worth being precise about what that buys, because it differs by case:
      // for a loop collapsing under this function's own dedup the legacy path
      // paves too, so the check only stops this chart OWNING the failure; but
      // for the two-vertex case the legacy path DECLINES, emitting nothing
      // where this chart filled a hole with 832 triangles. That one was a
      // genuine regression, not a shared defect.
      //
      // All found by review on bldrs-ai/conway-geom#191 and #192.
      const bool everyBoundSurvivedCollection = loops.size() == bounds.size();

      WingedEdgeMesh< glm::dvec3 > unwrapMesh;

      size_t closedLoopsReachingCdt = 0;

      if ( everyBoundSurvivedCollection && !loops.empty() &&
           triangulateUnwrappedLoops(
             loops, unwrapMesh, "sphere", &closedLoopsReachingCdt ) &&
           closedLoopsReachingCdt == loops.size() ) {

        tesselate(
          unwrapMesh,
          [&]( const glm::dvec3 &point ) {

            return glm::normalize( point - cent ) * radius + cent;
          },
          unwrapMesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
          relativeDeflectionSquared( unwrapMesh, representationExtent ) );

        // A sphere's outward normal is the offset from the centre, so unlike
        // the cylinder there is no degenerate direction to guard: a triangle
        // centroid can only land on the centre if the triangle spans the whole
        // sphere, which a trimmed loop cannot produce. Sense handling matches
        // the cylinder path - only orient against the surface normal when the
        // extractor actually populated the flag (IfcSurface::sameSenseKnown),
        // and use surface.sameSense rather than the local copy, which carries
        // a uv-space handedness correction that would be applied twice here.
        if ( surface.sameSenseKnown ) {

          appendMeshToGeometry(
            unwrapMesh,
            geometry,
            surface.sameSense,
            [&]( const glm::dvec3& point ) {

              return point - cent;
            } );

        } else {

          appendMeshToGeometry( unwrapMesh, geometry );
        }

        return;
      }
    }
  }

  WingedEdgeMesh< glm::dvec3 > mesh;

  tesselateDualParametrization(
    mesh,
    bounds,
    [&]( const glm::dvec3& point ) {
      // Produce a normalized vector from the centroid to the point.
      glm::dvec3 deltaCentroid = glm::normalize( point - cent );

      // we can normalize first because rotation is invariant
      // relative the centroid
      double dx = glm::dot( vecX, deltaCentroid );
      double dy = glm::dot( vecY, deltaCentroid );
      double dz = glm::dot (vecZ, deltaCentroid );

      // Project the point onto the unit sphere surface
      return glm::normalize( glm::dvec3( dx, dy, dz ) );
    },
    [&]( const glm::dvec3& normalFormVertex ) {

      const double z = ( 1 + normalFormVertex.z );

      if ( 2.0 - z < DBL_EPSILON ) {
        return glm::dvec2( 4, 4 );
      }

      return poleSafeStereographic( normalFormVertex, z );
    },
    [&]( const glm::dvec3& normalFormVertex ) {

      const double z = ( 1 - normalFormVertex.z );

      if ( 2.0 - z < DBL_EPSILON ) {
        return glm::dvec2( 4, 4 );
      }

      return poleSafeStereographic( normalFormVertex, z );
    },
    []( const glm::dvec3& normalFormVertex ) {

      return ( normalFormVertex.z <= 0.0 ) ? 0 : 1;
    },
    []( const glm::dvec3& normalFormVertex1,
        const glm::dvec3& normalFormVertex2,
        const glm::dvec2& paramVertex1,
        const glm::dvec2& paramVertex2 ) {

      if ( normalFormVertex1.z <= 0.0 && normalFormVertex2.z <= 0.0 ) {
        return false;
      }

      if ( normalFormVertex1.z > ( 1.0 - DBL_EPSILON ) || normalFormVertex2.z > ( 1.0 - DBL_EPSILON ) ) {
        return true;
      }

      return glm::distance( paramVertex1, paramVertex2 ) > 1.0;
    },
    []( const glm::dvec3& normalFormVertex1,
        const glm::dvec3& normalFormVertex2,
        const glm::dvec2& paramVertex1,
        const glm::dvec2& paramVertex2 ) {

        if ( normalFormVertex1.z > 0.0 && normalFormVertex2.z > 0.0 ) {
          return false;
        }

        if ( normalFormVertex1.z < ( -1.0 + DBL_EPSILON ) || normalFormVertex2.z < ( -1.0 + DBL_EPSILON ) ) {
          return true;
        }

        return glm::distance( paramVertex1, paramVertex2 ) > 1.0;
    } );

  tesselate(
    mesh,
    [&]( const glm::dvec3& point ) { 
      
      return glm::normalize( point - cent ) * radius + cent;
    },
    mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
    relativeDeflectionSquared( mesh, representationExtent ) );

  appendMeshToGeometry( mesh, geometry );
}


// ---------------------------------------------------------------------------
// Toroidal faces: topology-aware (theta, phi) unwrap.
//
// A torus's parameter domain is doubly periodic (major angle theta around the
// axis, minor angle phi around the tube), so no single planar unwrap is
// injective for every face. The previous dual-hemisphere machinery unwrapped
// each tube half onto an annulus keyed on sin(phi) (2-to-1 per half, folding
// the inner/outer equators together) and stitched the halves with a third CDT
// pass fed by distance-thresholded "leak" edges. That fold made equator-
// touching boundaries coincide in 2D (CDT "intersecting constraint edges" /
// "duplicate vertex" - dropped jet-engine shaft fillets, test-models#47), and
// the discard thresholds were calibrated to the exact annulus scale, so any
// reparameterization broke hole nesting elsewhere (the nist_ftc_11 washer).
//
// Instead, classify each face by its boundary winding/coverage per axis and
// pick a 2D domain that is injective for that class:
//
//   - wraps neither axis  -> rectangle in (theta', phi'), branch cuts placed
//                            in empty angular gaps of the boundary;
//   - wraps theta only    -> annulus: angle = theta, radius = phi' interval
//                            normalized to [1, 2] (washer-style rims);
//   - wraps phi only      -> annulus: angle = phi, radius = theta' interval
//                            normalized to [1, 2] (tube segments / fillets);
//   - wraps both          -> the boundary covers the whole torus (seam-edge
//                            faces): emit a parametric grid directly.
//
// Every boundary loop stays closed in 2D, so hole nesting is handled exactly
// by CDT's eraseOuterTrianglesAndHoles - no edge discarding, no hull fill, no
// equator stitching. Seam-doubled boundary edges dedupe to a single interior
// constraint. The coarse CDT triangles are then refined by the shared
// adaptive on-surface subdivision (tesselate), which restores metric accuracy
// regardless of the unwrap's distortion.
// ---------------------------------------------------------------------------

inline void TriangulateToroidalSurface(
    Geometry &geometry,
    const std::vector<IfcBound3D> &bounds,
    IfcSurface &surface,
    double representationExtent) {

  using namespace unwrap_detail;

  if ( bounds.empty() ) {
    return;
  }

  double     majorRadius = surface.ToroidalSurface.MajorRadius;
  double     minorRadius = surface.ToroidalSurface.MinorRadius;
  glm::dvec3 cent        = surface.transformation[3];
  glm::dvec3 vecX        = glm::normalize(surface.transformation[0]);
  glm::dvec3 vecY        = glm::normalize(surface.transformation[1]);
  glm::dvec3 vecZ        = glm::normalize(surface.transformation[2]);

  // Torus point in world space from (theta, phi).
  auto torusPoint = [&]( double theta, double phi ) {

    double planar = majorRadius + minorRadius * std::cos( phi );

    glm::dvec3 local(
      planar * std::cos( theta ),
      planar * std::sin( theta ),
      minorRadius * std::sin( phi ) );

    return vecX * local.x + vecY * local.y + vecZ * local.z + cent;
  };

  // Outward tube normal in world space at an arbitrary world-space point,
  // used only to orient output triangles consistently.
  auto tubeNormal = [&]( const glm::dvec3 &point ) {

    glm::dvec3 delta = point - cent;

    double dx = glm::dot( vecX, delta );
    double dy = glm::dot( vecY, delta );
    double dz = glm::dot( vecZ, delta );

    double planar = std::sqrt( dx * dx + dy * dy );

    if ( planar < 1e-12 ) {
      return glm::dvec3( 0.0 );
    }

    glm::dvec3 ringCenter = glm::dvec3( dx / planar, dy / planar, 0.0 ) * majorRadius;
    glm::dvec3 local      = glm::dvec3( dx, dy, dz ) - ringCenter;

    double length = glm::length( local );

    if ( length < 1e-12 ) {
      return glm::dvec3( 0.0 );
    }

    local /= length;

    return vecX * local.x + vecY * local.y + vecZ * local.z;
  };

  // --- 1. Boundary loops -> (theta, phi) samples ---------------------------

  std::vector< std::vector< BoundaryPoint > > loops;

  loops.reserve( bounds.size() );

  for ( const IfcBound3D &bound : bounds ) {

    const std::vector< glm::dvec3 > &points = bound.curve.points;

    std::vector< BoundaryPoint > loop;

    loop.reserve( points.size() );

    auto parameterDuplicate = []( const BoundaryPoint &a, double theta, double phi ) {

      return
        std::abs( wrapDeltaPi( theta - a.theta ) ) < 1e-12 &&
        std::abs( wrapDeltaPi( phi - a.phi ) ) < 1e-12;
    };

    for ( const glm::dvec3 &point : points ) {

      glm::dvec3 delta = point - cent;

      double dx = glm::dot( vecX, delta );
      double dy = glm::dot( vecY, delta );
      double dz = glm::dot( vecZ, delta );

      if ( !std::isfinite( dx ) || !std::isfinite( dy ) || !std::isfinite( dz ) ) {
        continue;
      }

      double planar = std::sqrt( dx * dx + dy * dy );

      // On-axis sample (degenerate spindle torus) - theta is undefined, drop.
      if ( planar < 1e-12 ) {
        continue;
      }

      double theta = std::atan2( dy, dx );
      double phi   = std::atan2( dz, planar - majorRadius );

      if ( !loop.empty() && parameterDuplicate( loop.back(), theta, phi ) ) {
        continue;
      }

      loop.push_back( BoundaryPoint{ point, theta, phi } );
    }

    // Drop the closing duplicate if the sampler emitted first == last.
    while ( loop.size() > 1 &&
            parameterDuplicate( loop.front(), loop.back().theta, loop.back().phi ) ) {
      loop.pop_back();
    }

    if ( loop.size() >= 3 ) {
      loops.push_back( std::move( loop ) );
    }
  }

  if ( loops.empty() ) {
    return;
  }

  // --- 2. Per-axis winding + coverage classification -----------------------

  bool windsTheta = false;
  bool windsPhi   = false;

  std::vector< double > thetaSamples;
  std::vector< double > phiSamples;

  for ( const std::vector< BoundaryPoint > &loop : loops ) {

    double sumTheta = 0.0;
    double sumPhi   = 0.0;

    for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

      const BoundaryPoint &from = loop[ where ];
      const BoundaryPoint &to   = loop[ ( where + 1 ) % count ];

      sumTheta += wrapDeltaPi( to.theta - from.theta );
      sumPhi   += wrapDeltaPi( to.phi - from.phi );

      thetaSamples.push_back( positiveMod2Pi( from.theta ) );
      phiSamples.push_back( positiveMod2Pi( from.phi ) );
    }

    if ( std::llround( sumTheta / kTwoPi ) != 0 ) {
      windsTheta = true;
    }

    if ( std::llround( sumPhi / kTwoPi ) != 0 ) {
      windsPhi = true;
    }
  }

  auto [ thetaGap, thetaCut ] = largestCircularGap( thetaSamples );
  auto [ phiGap, phiCut ]     = largestCircularGap( phiSamples );

  // Seam-doubled boundaries (a closed torus written with its seam edges
  // walked in both directions) have zero net winding but leave no empty
  // angular gap for a branch cut - treat a covered axis as wrapping.
  constexpr double MINIMUM_CUT_GAP = 1e-3;

  bool wrapsTheta = windsTheta || thetaGap < MINIMUM_CUT_GAP;
  bool wrapsPhi   = windsPhi || phiGap < MINIMUM_CUT_GAP;

  WingedEdgeMesh< glm::dvec3 > mesh;

  // pmr::vector (AFTP arena-backing), mirrors tesselateDualParametrization.
  auto &meshVertices = mesh.vertices;

  if ( wrapsTheta && wrapsPhi ) {

    // --- Full-coverage face: parametric grid -------------------------------
    // The boundary covers both periods (closed torus with seam edges, or an
    // exotic doubly-winding trim, which this approximates as the full torus).

    constexpr int GRID_THETA = 48;
    constexpr int GRID_PHI   = 24;

    for ( int j = 0; j < GRID_PHI; ++j ) {
      for ( int i = 0; i < GRID_THETA; ++i ) {
        meshVertices.push_back(
          torusPoint( i * kTwoPi / GRID_THETA, j * kTwoPi / GRID_PHI ) );
      }
    }

    for ( int j = 0; j < GRID_PHI; ++j ) {
      for ( int i = 0; i < GRID_THETA; ++i ) {

        uint32_t i00 = static_cast< uint32_t >( j * GRID_THETA + i );
        uint32_t i10 = static_cast< uint32_t >( j * GRID_THETA + ( i + 1 ) % GRID_THETA );
        uint32_t i01 = static_cast< uint32_t >( ( ( j + 1 ) % GRID_PHI ) * GRID_THETA + i );
        uint32_t i11 = static_cast< uint32_t >( ( ( j + 1 ) % GRID_PHI ) * GRID_THETA + ( i + 1 ) % GRID_THETA );

        // (+theta, +phi) quad order gives outward winding; honor reversed
        // face sense like the CDT path's senseSign does.
        if ( surface.sameSense ) {
          mesh.makeTriangle( i00, i10, i11 );
          mesh.makeTriangle( i00, i11, i01 );
        } else {
          mesh.makeTriangle( i00, i11, i10 );
          mesh.makeTriangle( i00, i01, i11 );
        }
      }
    }

  } else {

    // --- 3. Injective 2D layout for the chosen topology --------------------

    std::vector< std::vector< glm::dvec2 > > loops2D;

    // Interior Steiner points on the analytic (theta, phi) grid, in layout
    // coordinates with exact world positions. CDT sees only boundary
    // vertices otherwise, and for a face whose rails sit diametrically
    // opposite on the tube (a half-tube patch has both long edges at z = 0)
    // its triangles chord straight through the tube - and the closest-point
    // refinement can never recover the arc, because every chord midpoint
    // projects back onto the same z = 0 equators. Seeding the interior with
    // on-surface points makes the initial triangulation follow the tube;
    // refinement then only polishes. (Seen on the jet engine's fuel
    // manifold: a thin R=538/r=17.5 torus split into half-tube segments
    // rendered as flat creased ribbons.)
    std::vector< std::pair< glm::dvec2, glm::dvec3 > > steinerPoints;

    constexpr double THETA_STEP = kTwoPi / 48.0;
    constexpr double PHI_STEP   = kTwoPi / 24.0;
    constexpr int    MAX_GRID   = 96;

    auto gridSteps = []( double span, double step, int maxSteps ) {
      return std::clamp(
        static_cast< int >( std::ceil( span / step ) ), 1, maxSteps );
    };

    // A boundary edge whose cut-normalized coordinate jumps by more than pi
    // straddles the branch cut, i.e. the face actually passes through the
    // "empty" gap - upgrade that axis to wrapping and re-layout (each retry
    // strictly increases the wrap flags, so this terminates).
    for ( int attempt = 0; attempt < 3; ++attempt ) {

      if ( wrapsTheta && wrapsPhi ) {
        // Both upgraded by straddle detection - triangulated below as a grid
        // is not possible anymore at this point; drop the face like the CDT
        // failure path would. In practice straddles this deep mean a
        // malformed boundary.
        Logger::logWarning( "Toroidal unwrap: boundary straddles both branch cuts, dropping face." );
        return;
      }

      loops2D.assign( loops.size(), {} );
      steinerPoints.clear();

      double normMin  = 0.0;
      double normSpan = 1.0;

      // Cut-normalized coordinate of the non-wrapping axis (or both, for the
      // rectangle case), computed per point: cut + ((raw - cut) mod 2pi).
      auto normalizeTheta = [&]( double raw ) {
        return positiveMod2Pi( raw - thetaCut );
      };
      auto normalizePhi = [&]( double raw ) {
        return positiveMod2Pi( raw - phiCut );
      };

      if ( wrapsTheta ) {

        // Annulus: angle = theta, radius = normalized phi in [1, 2].
        double lo = std::numeric_limits< double >::max();
        double hi = std::numeric_limits< double >::lowest();

        for ( const std::vector< BoundaryPoint > &loop : loops ) {
          for ( const BoundaryPoint &bp : loop ) {
            double value = normalizePhi( bp.phi );
            lo = std::min( lo, value );
            hi = std::max( hi, value );
          }
        }

        normMin  = lo;
        normSpan = hi - lo;

        if ( normSpan < 1e-9 ) {
          return;  // Degenerate band - no interior to triangulate.
        }

        for ( size_t which = 0; which < loops.size(); ++which ) {
          for ( const BoundaryPoint &bp : loops[ which ] ) {

            double radius = 1.0 + ( normalizePhi( bp.phi ) - normMin ) / normSpan;

            loops2D[ which ].push_back(
              glm::dvec2( radius * std::cos( bp.theta ), radius * std::sin( bp.theta ) ) );
          }
        }

        int gridTheta = gridSteps( kTwoPi, THETA_STEP, MAX_GRID );
        int gridPhi   = gridSteps( normSpan, PHI_STEP, MAX_GRID );

        for ( int i = 0; i < gridTheta; ++i ) {
          for ( int j = 1; j < gridPhi; ++j ) {

            double theta  = i * kTwoPi / gridTheta;
            double phiCutNormalized =
              normMin + normSpan * j / static_cast< double >( gridPhi );
            double radius = 1.0 + ( phiCutNormalized - normMin ) / normSpan;

            steinerPoints.emplace_back(
              glm::dvec2( radius * std::cos( theta ), radius * std::sin( theta ) ),
              torusPoint( theta, phiCut + phiCutNormalized ) );
          }
        }

      } else if ( wrapsPhi ) {

        // Annulus: angle = phi, radius = normalized theta in [1, 2].
        double lo = std::numeric_limits< double >::max();
        double hi = std::numeric_limits< double >::lowest();

        for ( const std::vector< BoundaryPoint > &loop : loops ) {
          for ( const BoundaryPoint &bp : loop ) {
            double value = normalizeTheta( bp.theta );
            lo = std::min( lo, value );
            hi = std::max( hi, value );
          }
        }

        normMin  = lo;
        normSpan = hi - lo;

        if ( normSpan < 1e-9 ) {
          return;
        }

        for ( size_t which = 0; which < loops.size(); ++which ) {
          for ( const BoundaryPoint &bp : loops[ which ] ) {

            double radius = 1.0 + ( normalizeTheta( bp.theta ) - normMin ) / normSpan;

            loops2D[ which ].push_back(
              glm::dvec2( radius * std::cos( bp.phi ), radius * std::sin( bp.phi ) ) );
          }
        }

        int gridPhi   = gridSteps( kTwoPi, PHI_STEP, MAX_GRID );
        int gridTheta = gridSteps( normSpan, THETA_STEP, MAX_GRID );

        for ( int j = 0; j < gridPhi; ++j ) {
          for ( int i = 1; i < gridTheta; ++i ) {

            double phi = j * kTwoPi / gridPhi;
            double thetaCutNormalized =
              normMin + normSpan * i / static_cast< double >( gridTheta );
            double radius = 1.0 + ( thetaCutNormalized - normMin ) / normSpan;

            steinerPoints.emplace_back(
              glm::dvec2( radius * std::cos( phi ), radius * std::sin( phi ) ),
              torusPoint( thetaCut + thetaCutNormalized, phi ) );
          }
        }

      } else {

        // Rectangle in (theta', phi'), normalized per axis to [0, 1].
        double loT = std::numeric_limits< double >::max();
        double hiT = std::numeric_limits< double >::lowest();
        double loP = std::numeric_limits< double >::max();
        double hiP = std::numeric_limits< double >::lowest();

        for ( const std::vector< BoundaryPoint > &loop : loops ) {
          for ( const BoundaryPoint &bp : loop ) {
            double tU = normalizeTheta( bp.theta );
            double pU = normalizePhi( bp.phi );
            loT = std::min( loT, tU );
            hiT = std::max( hiT, tU );
            loP = std::min( loP, pU );
            hiP = std::max( hiP, pU );
          }
        }

        if ( hiT - loT < 1e-9 || hiP - loP < 1e-9 ) {
          return;
        }

        for ( size_t which = 0; which < loops.size(); ++which ) {
          for ( const BoundaryPoint &bp : loops[ which ] ) {
            loops2D[ which ].push_back(
              glm::dvec2(
                ( normalizeTheta( bp.theta ) - loT ) / ( hiT - loT ),
                ( normalizePhi( bp.phi ) - loP ) / ( hiP - loP ) ) );
          }
        }

        int gridTheta = gridSteps( hiT - loT, THETA_STEP, MAX_GRID );
        int gridPhi   = gridSteps( hiP - loP, PHI_STEP, MAX_GRID );

        for ( int i = 1; i < gridTheta; ++i ) {
          for ( int j = 1; j < gridPhi; ++j ) {

            double thetaCutNormalized =
              loT + ( hiT - loT ) * i / static_cast< double >( gridTheta );
            double phiCutNormalized =
              loP + ( hiP - loP ) * j / static_cast< double >( gridPhi );

            steinerPoints.emplace_back(
              glm::dvec2(
                i / static_cast< double >( gridTheta ),
                j / static_cast< double >( gridPhi ) ),
              torusPoint(
                thetaCut + thetaCutNormalized, phiCut + phiCutNormalized ) );
          }
        }
      }

      // Straddle check: consecutive samples on a non-wrapping axis must stay
      // within half a period of each other after cut-normalization.
      bool straddleTheta = false;
      bool straddlePhi   = false;

      for ( size_t which = 0; which < loops.size(); ++which ) {

        const std::vector< BoundaryPoint > &loop = loops[ which ];

        for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

          const BoundaryPoint &from = loop[ where ];
          const BoundaryPoint &to   = loop[ ( where + 1 ) % count ];

          if ( !wrapsTheta &&
               std::abs( normalizeTheta( to.theta ) - normalizeTheta( from.theta ) ) > kPi ) {
            straddleTheta = true;
          }

          if ( !wrapsPhi &&
               std::abs( normalizePhi( to.phi ) - normalizePhi( from.phi ) ) > kPi ) {
            straddlePhi = true;
          }
        }
      }

      if ( straddleTheta || straddlePhi ) {
        wrapsTheta = wrapsTheta || straddleTheta;
        wrapsPhi   = wrapsPhi || straddlePhi;
        continue;
      }

      break;
    }

    if ( wrapsTheta && wrapsPhi ) {
      return;  // Straddle upgrades exhausted the planar cases (logged above).
    }

    // --- 4. CDT input: weld 2D-coincident vertices, dedupe edges ------------
    // Seam-doubled boundary edges (the same physical edge walked twice) and
    // shared loop corners map to identical 2D points; welding them up front
    // is what makes CDT's NotAllowed constraint mode safe here.

    std::vector< CDT::V2d< double > > cdtVertices;
    std::vector< glm::dvec3 >         cdtWorld;
    std::vector< CDT::Edge >          cdtEdges;

    std::map< std::pair< long long, long long >, uint32_t > weld;
    std::set< std::pair< uint32_t, uint32_t > >             edgeSet;

    auto weldVertex = [&]( const glm::dvec2 &position, const glm::dvec3 &world ) {

      // 2D layouts are normalized to O(1) extents, so a fixed 1e-9 quantum
      // welds only numerically-identical parameter points.
      std::pair< long long, long long > key(
        std::llround( position.x * 1e9 ),
        std::llround( position.y * 1e9 ) );

      auto [ found, isNew ] = weld.try_emplace( key, static_cast< uint32_t >( cdtVertices.size() ) );

      if ( isNew ) {
        cdtVertices.emplace_back( position.x, position.y );
        cdtWorld.push_back( world );
      }

      return found->second;
    };

    for ( size_t which = 0; which < loops.size(); ++which ) {

      const std::vector< BoundaryPoint > &loop   = loops[ which ];
      const std::vector< glm::dvec2 >    &loop2D = loops2D[ which ];

      for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

        size_t next = ( where + 1 ) % count;

        uint32_t v1 = weldVertex( loop2D[ where ], loop[ where ].world );
        uint32_t v2 = weldVertex( loop2D[ next ], loop[ next ].world );

        if ( v1 == v2 ) {
          continue;
        }

        std::pair< uint32_t, uint32_t > ordered( std::min( v1, v2 ), std::max( v1, v2 ) );

        if ( edgeSet.insert( ordered ).second ) {
          cdtEdges.emplace_back( v1, v2 );
        }
      }
    }

    if ( cdtEdges.size() < 3 ) {
      return;
    }

    // Interior grid points join as plain (unconstrained) vertices after the
    // boundary, so constraint edge indices stay valid. Skip any that weld
    // onto an existing boundary vertex.
    size_t boundaryVertexCount = cdtVertices.size();

    for ( const std::pair< glm::dvec2, glm::dvec3 > &steiner : steinerPoints ) {

      std::pair< long long, long long > key(
        std::llround( steiner.first.x * 1e9 ),
        std::llround( steiner.first.y * 1e9 ) );

      if ( weld.find( key ) != weld.end() ) {
        continue;
      }

      cdtVertices.emplace_back( steiner.first.x, steiner.first.y );
      cdtWorld.push_back( steiner.second );
    }

    for ( const CDT::V2d< double > &v : cdtVertices ) {
      if ( !std::isfinite( v.x ) || !std::isfinite( v.y ) ) {
        throw std::runtime_error(
          "conway: non-finite toroidal unwrap; dropping face" );
      }
    }

    // Try the seeded triangulation first; if CDT rejects the input (e.g. a
    // grid point landing exactly on a constraint edge), retry with the
    // boundary alone before dropping the face - only that final failure is
    // worth logging.
    auto runCdt = [&]( size_t vertexCount )
        -> std::optional< std::vector< CDT::Triangle > > {

      CDT::Triangulation< double > triangulation(
        CDT::VertexInsertionOrder::Auto,
        CDT::IntersectingConstraintEdges::TryResolve, 0);

      bool logFailure = vertexCount == boundaryVertexCount;

      try
      {
        conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
        triangulation.insertVertices(
          std::vector< CDT::V2d< double > >(
            cdtVertices.begin(),
            cdtVertices.begin() + static_cast< ptrdiff_t >( vertexCount ) ) );
        triangulation.insertEdges( cdtEdges );
        triangulation.eraseOuterTrianglesAndHoles();
      }
      catch ( const CDT::IntersectingConstraintsError &e )
      {
        if ( logFailure ) {
          const CDT::V2d< double > &ev1 = cdtVertices[ e.e1().v1() ];
          const CDT::V2d< double > &ev2 = cdtVertices[ e.e1().v2() ];
          const CDT::V2d< double > &ev3 = cdtVertices[ e.e2().v1() ];
          const CDT::V2d< double > &ev4 = cdtVertices[ e.e2().v2() ];

          Logger::logError( "CDT Exception (torus unwrap) ((%f,%f),(%f,%f)) -> ((%f,%f),(%f,%f)): %s",
            ev1.x, ev1.y, ev2.x, ev2.y, ev3.x, ev3.y, ev4.x, ev4.y, e.what() );
        }
        return std::nullopt;
      }
      catch ( const std::exception &e )
      {
        if ( logFailure ) {
          Logger::logError( "CDT Exception (torus unwrap): %s", e.what() );
        }
        return std::nullopt;
      }

      return std::vector< CDT::Triangle >(
        triangulation.triangles.begin(), triangulation.triangles.end() );
    };

    size_t usedVertexCount = cdtVertices.size();

    std::optional< std::vector< CDT::Triangle > > cdtTriangles =
      runCdt( usedVertexCount );

    if ( !cdtTriangles.has_value() && usedVertexCount > boundaryVertexCount ) {
      usedVertexCount = boundaryVertexCount;
      cdtTriangles   = runCdt( usedVertexCount );
    }

    if ( !cdtTriangles.has_value() ) {
      return;
    }

    meshVertices.reserve( usedVertexCount );

    for ( size_t where = 0; where < usedVertexCount; ++where ) {
      meshVertices.push_back( cdtWorld[ where ] );
    }

    // --- 5. Emit triangles, oriented by the analytic tube normal -----------

    double senseSign = surface.sameSense ? 1.0 : -1.0;

    for ( const CDT::Triangle &triangle : *cdtTriangles ) {

      auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

      if ( cdtv1 == cdtv2 || cdtv2 == cdtv3 || cdtv3 == cdtv1 ) {
        continue;
      }

      // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
      if (
        cdtv1 >= usedVertexCount ||
        cdtv2 >= usedVertexCount ||
        cdtv3 >= usedVertexCount ) {
        continue;
      }

      const glm::dvec3 &w1 = meshVertices[ cdtv1 ];
      const glm::dvec3 &w2 = meshVertices[ cdtv2 ];
      const glm::dvec3 &w3 = meshVertices[ cdtv3 ];

      glm::dvec3 reference =
        tubeNormal( ( w1 + w2 + w3 ) / 3.0 ) * senseSign;

      if ( glm::dot( glm::cross( w2 - w1, w3 - w1 ), reference ) < 0.0 ) {
        mesh.makeTriangle( cdtv1, cdtv3, cdtv2 );
      } else {
        mesh.makeTriangle( cdtv1, cdtv2, cdtv3 );
      }
    }
  }

  // --- 6. Shared adaptive on-surface refinement -----------------------------

  // Scale-aware refinement threshold (relativeDeflectionSquared): the old
  // absolute MAX_DEFLECTION was unreachable at large model scale and read as
  // "refine until the budget runs out" - tolerable when the budget was tiny
  // (boundary-only CDT), runaway now that interior seeding raises it.
  tesselate(
    mesh,
    [&]( const glm::dvec3& point ) {

      // Produce a normalized vector from the centroid to the point.
      glm::dvec3 deltaCentroid = point - cent;

      // we can normalize first because rotation is invariant
      // relative the centroid
      double dx = glm::dot( vecX, deltaCentroid );
      double dy = glm::dot( vecY, deltaCentroid );
      double dz = glm::dot (vecZ, deltaCentroid );

      glm::dvec2 planar = glm::normalize( glm::dvec2( dx, dy ) );

      // Centroid on the ring.
      glm::dvec3 ringCenter =
        glm::dvec3( glm::normalize( planar ) * majorRadius, 0.0 );

      glm::dvec3 normalPointOnRing = glm::normalize( glm::dvec3( dx, dy, dz ) - ringCenter );

      glm::dvec3 pointOnIdentityRing = ringCenter + ( normalPointOnRing * minorRadius );

      // Move back to the original coordinate frame.
      return
        vecX * pointOnIdentityRing.x +
        vecY * pointOnIdentityRing.y +
        vecZ * pointOnIdentityRing.z + cent;
    },
    mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
    relativeDeflectionSquared( mesh, representationExtent ) );

  appendMeshToGeometry( mesh, geometry );
}


// TODO: review and simplify
inline void TriangulateConicalSurface(
  Geometry &geometry,
  const std::vector<IfcBound3D> &bounds,
  IfcSurface &surface,
  double representationExtent) {
  // First we get the cylinder data

  if ( bounds.empty() ) {
    return;
  }

  double radius    = surface.ConicalSurface.Radius;
  double semiAngle = surface.ConicalSurface.SemiAngle;

  double sinSemiAngle = tan( fabs( semiAngle ) );
  
  glm::dvec3 cent = surface.transformation[3];
  glm::dvec3 vecX = glm::normalize(surface.transformation[0]);
  glm::dvec3 vecY = glm::normalize(surface.transformation[1]);
  glm::dvec3 vecZ = glm::normalize(surface.transformation[2]);
  
  bool sameSense = surface.sameSense;

  // A mirrored (left-handed) placement flips the surface normal, so the
  // face sense has to flip with it. The previous test, dot(vecZ, vecX),
  // compared two ORTHONORMAL columns of that placement — always ~1e-17,
  // so its sign was rounding noise and an entire face's winding turned on
  // a coin toss. See https://github.com/bldrs-ai/conway/issues/463.
  if ( glm::dot( glm::cross( vecX, vecY ), vecZ ) < 0 ) {

    sameSense = !sameSense;
  }

  // --- Trimmed unwrap (jet compressor rings, follow-up to #149): lay the
  // boundary loops out as (theta, t) - theta the angle around the cone
  // axis, t the normalized axial coordinate - and CDT them. The legacy path
  // below projects boundaries axially ((dx, dy) / maxR, dropping z
  // entirely), which collapses steep, near-cylindrical cone bands to
  // slivers (inner and outer rims land on nearly the same circle); earcut
  // then emits chord triangles that the budgeted refinement inflates into
  // lumpy off-surface shelves. Any degeneracy here falls back to that
  // legacy path unchanged.
  {
    using namespace unwrap_detail;

    // Generator line of the cone in (z, r), least-squares fitted from the
    // boundary samples themselves (they lie exactly on the surface). The
    // fit sidesteps the SemiAngle sign convention - the legacy projection
    // uses tan(fabs(semiAngle)), which points the wrong way for narrowing
    // cones - and tolerates a Radius given at a different reference plane.
    // Fitted in centered form: raw moments (n*sumZZ - sumZ^2) cancel
    // catastrophically for a face at a large z-offset with a small span.
    std::vector< glm::dvec2 > generatorSamples;

    double loT = std::numeric_limits< double >::max();
    double hiT = std::numeric_limits< double >::lowest();

    // A boundary sample on the axis is the cone's apex: its theta is
    // undefined, so the unwrap cannot represent the fan around it - the
    // legacy path (which keeps the apex at the projection center and fans
    // to the rim) handles that case better. Bail out to it.
    bool sawAxisPoint = false;

    std::vector< std::vector< BoundaryPoint > > loops;

    loops.reserve( bounds.size() );

    for ( const IfcBound3D &bound : bounds ) {

      std::vector< BoundaryPoint > loop;

      loop.reserve( bound.curve.points.size() );

      for ( const glm::dvec3 &point : bound.curve.points ) {

        glm::dvec3 delta = point - cent;

        double dx = glm::dot( vecX, delta );
        double dy = glm::dot( vecY, delta );
        double dz = glm::dot( vecZ, delta );

        if ( !std::isfinite( dx ) || !std::isfinite( dy ) ||
             !std::isfinite( dz ) ) {
          continue;
        }

        double dd = std::sqrt( dx * dx + dy * dy );

        if ( dd < 1e-12 ) {
          sawAxisPoint = true;
          continue;
        }

        generatorSamples.emplace_back( dz, dd );

        loT = std::min( loT, dz );
        hiT = std::max( hiT, dz );

        loop.push_back( BoundaryPoint{ point, std::atan2( dy, dx ), dz } );
      }

      if ( loop.size() >= 3 ) {
        loops.push_back( std::move( loop ) );
      }
    }

    double spanT = hiT - loT;

    double meanZ = 0.0;
    double meanR = 0.0;

    for ( const glm::dvec2 &sample : generatorSamples ) {
      meanZ += sample.x;
      meanR += sample.y;
    }

    if ( !generatorSamples.empty() ) {
      meanZ /= static_cast< double >( generatorSamples.size() );
      meanR /= static_cast< double >( generatorSamples.size() );
    }

    double covZZ = 0.0;
    double covZR = 0.0;

    for ( const glm::dvec2 &sample : generatorSamples ) {
      covZZ += ( sample.x - meanZ ) * ( sample.x - meanZ );
      covZR += ( sample.x - meanZ ) * ( sample.y - meanR );
    }

    if ( !sawAxisPoint && !loops.empty() && spanT > 1e-12 && covZZ > 0.0 ) {

      double slope = covZR / covZZ;

      // Normalize t into [0, 1] and drop parameter-duplicate neighbors
      // (the weld would collapse them anyway; empty edges upset CDT less
      // when they never exist).
      auto parameterDuplicate = []( const BoundaryPoint &a,
                                    const BoundaryPoint &b ) {

        return
          std::abs( wrapDeltaPi( b.theta - a.theta ) ) < 1e-12 &&
          std::abs( b.phi - a.phi ) < 1e-12;
      };

      for ( std::vector< BoundaryPoint > &loop : loops ) {

        for ( BoundaryPoint &point : loop ) {
          point.phi = ( point.phi - loT ) / spanT;
        }

        std::vector< BoundaryPoint > cleaned;

        cleaned.reserve( loop.size() );

        for ( const BoundaryPoint &point : loop ) {
          if ( cleaned.empty() || !parameterDuplicate( cleaned.back(), point ) ) {
            cleaned.push_back( point );
          }
        }

        while ( cleaned.size() > 1 &&
                parameterDuplicate( cleaned.front(), cleaned.back() ) ) {
          cleaned.pop_back();
        }

        loop = std::move( cleaned );
      }

      std::erase_if( loops, []( const std::vector< BoundaryPoint > &loop ) {
        return loop.size() < 3;
      } );

      WingedEdgeMesh< glm::dvec3 > mesh;

      if ( triangulateUnwrappedLoops( loops, mesh, "cone" ) ) {

        tesselate(
          mesh,
          [&]( const glm::dvec3 &point ) {

            glm::dvec3 delta = point - cent;

            double dx = glm::dot( vecX, delta );
            double dy = glm::dot( vecY, delta );
            double dz = glm::dot( vecZ, delta );
            double dd = std::sqrt( dx * dx + dy * dy );

            if ( dd < 1e-12 ) {
              return point;  // On the axis - angle undefined.
            }

            double targetRadius =
              std::max( meanR + slope * ( dz - meanZ ), 0.0 );

            return cent +
              vecX * ( dx / dd * targetRadius ) +
              vecY * ( dy / dd * targetRadius ) +
              vecZ * dz;
          },
          mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
          relativeDeflectionSquared( mesh, representationExtent ) );

        // Same idea as the cylinder, tilted by the cone's taper: the
        // radius model this unwrap tesselates against is
        // `meanR + slope * (dz - meanZ)`, so the generator runs along
        // (radial, slope) and the outward normal along (radial, -slope).
        // Gated on sameSenseKnown, and reading surface.sameSense rather
        // than the handedness-flipped local, for the same reasons as the
        // cylinder path above.
        if ( surface.sameSenseKnown ) {

          appendMeshToGeometry(
            mesh,
            geometry,
            surface.sameSense,
            [&]( const glm::dvec3& point ) {

              glm::dvec3 delta  = point - cent;
              glm::dvec3 radial = delta - vecZ * glm::dot( delta, vecZ );
              double     length = glm::length( radial );

              // Radius-relative for the same reason as the cylinder
              // path; an absolute epsilon misjudges very small and very
              // large cones alike.
              if ( length < fabs( meanR ) * 1e-9 ) {

                return glm::dvec3( 0.0, 0.0, 0.0 );  // On the axis (apex).
              }

              return ( radial / length ) - vecZ * slope;
            } );

        } else {

          appendMeshToGeometry( mesh, geometry );
        }

        return;
      }
    }
  }

  std::vector<std::vector<glm::dvec3>> newPoints;

  double minR = DBL_MAX;
  double maxR = -DBL_MAX;

  std::priority_queue< std::pair< double, size_t > > outsideMostBoundaries;

  // Find the relative coordinates of each curve point in the cylinder reference
  // plane Only retain the max and min relative Z
  for (size_t i = 0; i < bounds.size(); i++) {

    double localMaxR = -DBL_MAX;
    double localMinR = DBL_MAX;

    for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {

      glm::dvec3 pt = bounds[ i ].curve.points[ j ];
      glm::dvec3 vv = pt - cent;

      double dx = glm::dot(vecX, vv);
      double dy = glm::dot(vecY, vv);

      double dr = glm::length( glm::dvec2( dx, dy ) );
      
      localMaxR = std::max( localMaxR, dr );
      localMinR = std::min( localMinR, dr );
    }

    outsideMostBoundaries.push( std::make_pair( localMaxR, i ) );

    maxR = std::max( maxR, localMaxR );
    minR = std::min( minR, localMinR );
  }

  using Point = std::array<double, 2>;
  std::vector<std::vector<Point>> uvBoundaryValues;
  std::vector<ParameterVertex> vertices;

  // AFTP: back this per-face tessellation mesh (esp. edge_map, which allocates
  // a node per edge) with the thread scratch arena, rewound at function exit.
  // Byte-identical: the arena changes only where the mesh's nodes live.
  conway::ScratchArenaScope arenaScope;
  WingedEdgeMesh< ParameterVertex > mesh{ conway::ThreadScratchResource() };

  while ( !outsideMostBoundaries.empty() ) {

    std::vector<Point> points;

    size_t boundsIndex = outsideMostBoundaries.top().second;

    outsideMostBoundaries.pop();
   const IfcBound3D& bound = bounds[ boundsIndex ];

    if ( bound.curve.points.empty() ) {
      continue;
    }

    for ( const glm::dvec3& pt : bound.curve.points ) {

      glm::dvec3 vv = pt - cent;

      double dx = glm::dot(vecX, vv);
      double dy = glm::dot(vecY, vv);
      double dz = glm::dot(vecZ, vv);

      glm::dvec2 pInv = glm::dvec2( dx, dy ) / maxR;

      points.push_back({pInv.x, pInv.y});
      mesh.makeVertex( { pt, pInv } );
    }

    uvBoundaryValues.push_back( points );
  }

#if (OUTPUT_SVG_DEBUG == 1) 

    // atomic: may run concurrently on the thread pool during staged face
    // finalization.
    static std::atomic< size_t > svgIndex = 0;

    size_t outputIndex = svgIndex++;

    std::ofstream svgFile( "cone_" + std::to_string( outputIndex ) + ".svg" );

    auto svgScale = []( double value ) {

      return 50 + ( 1024.0 * ( value + 1.0 ) / 2.0 ); 

    };

    svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
              << "width=\"1124\" height=\"1124\">\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"512\" style=\"stroke:rgb(255, 132, 0);stroke-width:2\" fill=\"none\"/>\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"256\" style=\"stroke:rgb(0, 0, 255);stroke-width:2\" fill=\"none\"/>\n";

    for ( const std::vector< Point >& loop  : uvBoundaryValues )
    {
      bool firstInLoop = true;
      
      glm::dvec2 lastPoint;

      for ( const Point &point : loop )
      {
        glm::dvec2 svgPoint = glm::dvec2( point[ 0 ], point[ 1 ] );

        if ( firstInLoop )
        {
          svgFile << "<polyline points=\"";
          firstInLoop = false;
        }
        else
        {
          svgFile << " ";
        }

        svgFile << svgScale( svgPoint.x ) << "," << svgScale( svgPoint.y );

        lastPoint = svgPoint;
      }

      svgFile << "\" style=\"fill:none;stroke:rgb(0,0,0);stroke-width:2\" />\n";
      
      for ( const Point &point : loop )
      {
        glm::dvec2 svgPoint = glm::dvec2( point[ 0 ], point[ 1 ] );
        
        svgFile << "<circle cx=\"" << svgScale( svgPoint.x ) << "\" cy=\"" << svgScale( svgPoint.y )
            << "\" r=\"3\" fill=\"red\"/>\n";
      }
    }

    svgFile << "</svg>\n";

    svgFile.close();

#endif

  // Triangulate projected boundary
  // Subdivide resulting triangles to increase definition
  // r indicates the level of subdivision, currently 3 you can increase it to
  // 5

  std::vector<uint32_t> indices;
  {
    conway::AllocTagScope earcutTag( conway::AllocSite::Earcut );
    indices = mapbox::earcut<uint32_t>(uvBoundaryValues);
  }

  for (size_t i = 0; i < indices.size(); i += 3) {

    mesh.makeTriangle( 
      indices[ i  + 0 ], 
      indices[ i  + 1 ], 
      indices[ i  + 2 ] );
  }

  tesselate(
    mesh,
    [&]( const glm::dvec3& point, [[maybe_unused]]const glm::dvec2& from ) { 
      
      glm::dvec3 vv = point - cent;
      double     dx = glm::dot( vecX, vv );
      double     dy = glm::dot( vecY, vv );    
      double     dz = glm::dot( vecZ, vv );

      glm::dvec3 coneSpacePoint = glm::dvec3( ( radius + dz * sinSemiAngle ) * glm::normalize( glm::dvec2( dx, dy ) ), dz );

      return cent + coneSpacePoint.x * vecX + coneSpacePoint.y * vecY + coneSpacePoint.z * vecZ;
    },
    mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
    relativeDeflectionSquared( mesh, representationExtent ) );

  appendMeshToGeometry( mesh, geometry, sameSense );
}

// TODO: review and simplify
inline void TriangulateCylindricalSurface(Geometry &geometry,
                                          const std::vector<IfcBound3D> &bounds,
                                          IfcSurface &surface,
                                          double representationExtent) {
  // First we get the cylinder data

  if ( bounds.empty() ) {
    return;
  }

  double radius = surface.CylinderSurface.Radius;
  glm::dvec3 cent = surface.transformation[3];
  glm::dvec3 vecX = glm::normalize(surface.transformation[0]);
  glm::dvec3 vecY = glm::normalize(surface.transformation[1]);
  glm::dvec3 vecZ = glm::normalize(surface.transformation[2]);

  bool sameSense = surface.sameSense;

  // A mirrored (left-handed) placement flips the surface normal, so the
  // face sense has to flip with it. The previous test, dot(vecZ, vecX),
  // compared two ORTHONORMAL columns of that placement — always ~1e-17,
  // so its sign was rounding noise and an entire face's winding turned on
  // a coin toss. See https://github.com/bldrs-ai/conway/issues/463.
  if ( glm::dot( glm::cross( vecX, vecY ), vecZ ) < 0 ) {

    sameSense = !sameSense;
  }

  // --- Trimmed unwrap (jet compressor rings, follow-up to #149): lay the
  // boundary loops out as (theta, t) - theta the angle around the cylinder
  // axis, t the normalized axial coordinate - and CDT them with exact hole
  // nesting. The legacy path below earcuts an annulus projection whose
  // ring ordering assumes strict radial nesting by max-z; its chord
  // triangles then get inflated by the budgeted refinement into lumpy
  // off-surface shelves. Any degeneracy here falls back to that legacy
  // path unchanged.
  {
    using namespace unwrap_detail;

    double loT = std::numeric_limits< double >::max();
    double hiT = std::numeric_limits< double >::lowest();

    // An on-axis boundary sample can't lie on the cylinder - it flags
    // degenerate input the legacy path already has behavior for.
    bool sawAxisPoint = false;

    std::vector< std::vector< BoundaryPoint > > loops;

    loops.reserve( bounds.size() );

    for ( const IfcBound3D &bound : bounds ) {

      std::vector< BoundaryPoint > loop;

      loop.reserve( bound.curve.points.size() );

      for ( const glm::dvec3 &point : bound.curve.points ) {

        glm::dvec3 delta = point - cent;

        double dx = glm::dot( vecX, delta );
        double dy = glm::dot( vecY, delta );
        double dz = glm::dot( vecZ, delta );

        if ( !std::isfinite( dx ) || !std::isfinite( dy ) ||
             !std::isfinite( dz ) ) {
          continue;
        }

        double dd = std::sqrt( dx * dx + dy * dy );

        if ( dd < 1e-12 ) {
          sawAxisPoint = true;
          continue;
        }

        loT = std::min( loT, dz );
        hiT = std::max( hiT, dz );

        loop.push_back( BoundaryPoint{ point, std::atan2( dy, dx ), dz } );
      }

      if ( loop.size() >= 3 ) {
        loops.push_back( std::move( loop ) );
      }
    }

    double spanT = hiT - loT;

    if ( !sawAxisPoint && !loops.empty() && spanT > 1e-12 ) {

      auto parameterDuplicate = []( const BoundaryPoint &a,
                                    const BoundaryPoint &b ) {

        return
          std::abs( wrapDeltaPi( b.theta - a.theta ) ) < 1e-12 &&
          std::abs( b.phi - a.phi ) < 1e-12;
      };

      for ( std::vector< BoundaryPoint > &loop : loops ) {

        for ( BoundaryPoint &point : loop ) {
          point.phi = ( point.phi - loT ) / spanT;
        }

        std::vector< BoundaryPoint > cleaned;

        cleaned.reserve( loop.size() );

        for ( const BoundaryPoint &point : loop ) {
          if ( cleaned.empty() || !parameterDuplicate( cleaned.back(), point ) ) {
            cleaned.push_back( point );
          }
        }

        while ( cleaned.size() > 1 &&
                parameterDuplicate( cleaned.front(), cleaned.back() ) ) {
          cleaned.pop_back();
        }

        loop = std::move( cleaned );
      }

      std::erase_if( loops, []( const std::vector< BoundaryPoint > &loop ) {
        return loop.size() < 3;
      } );

      WingedEdgeMesh< glm::dvec3 > unwrapMesh;

      if ( triangulateUnwrappedLoops( loops, unwrapMesh, "cylinder" ) ) {

        tesselate(
          unwrapMesh,
          [&]( const glm::dvec3 &point ) {

            glm::dvec3 delta = point - cent;

            double dx = glm::dot( vecX, delta );
            double dy = glm::dot( vecY, delta );
            double dz = glm::dot( vecZ, delta );
            double dd = std::sqrt( dx * dx + dy * dy );

            if ( dd < 1e-12 ) {
              return point;  // On the axis - angle undefined.
            }

            return cent +
              vecX * ( dx / dd * radius ) +
              vecY * ( dy / dd * radius ) +
              vecZ * dz;
          },
          unwrapMesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
          relativeDeflectionSquared( unwrapMesh, representationExtent ) );

        // A cylinder's outward normal is the radial component of the
        // offset from the axis; the axial part carries no orientation.
        // Only orient against the surface normal when the face sense is
        // actually known — see IfcSurface::sameSenseKnown. Extractors that
        // never populated it (IFC) keep the previous projection-based
        // behaviour rather than being handed a default as if it were an
        // answer.
        //
        // surface.sameSense, NOT the local `sameSense`: that one carries
        // the legacy path's handedness flip (see the det(vecX,vecY,vecZ)
        // test above), which corrects a uv-space winding decision. The
        // normal handed to appendMeshToGeometry below is computed in world
        // space, so it already points outward whichever way the placement
        // is handed. Feeding the flipped local in would apply the
        // correction twice and invert every face on a mirrored placement.
        if ( surface.sameSenseKnown ) {

          appendMeshToGeometry(
            unwrapMesh,
            geometry,
            surface.sameSense,
            [&]( const glm::dvec3& point ) {

              glm::dvec3 delta  = point - cent;
              glm::dvec3 radial = delta - vecZ * glm::dot( delta, vecZ );

              // A centroid can land on the axis even when all three
              // vertices are on the surface — the midpoint of two
              // diametrically opposed points is exactly there. The radial
              // direction is then pure noise and the triangle would be
              // swapped on its sign, which is the failure this overload
              // exists to remove. Scaled by the radius so a millimetre
              // bore and a metre drum behave alike.
              if ( glm::length( radial ) < radius * 1e-9 ) {

                return glm::dvec3( 0.0, 0.0, 0.0 );
              }

              return radial;
            } );

        } else {

          appendMeshToGeometry( unwrapMesh, geometry );
        }

        return;
      }
    }
  }

  std::vector<std::vector<glm::dvec3>> newPoints;

  double minZ = DBL_MAX;
  double maxZ = -DBL_MAX;

  std::priority_queue< std::pair< double, size_t > > outsideMostBoundaries;

  if ( bounds.size() == 1 && bounds[0].curve.points.size() < 3 ) {
    // If there is no curve, we can not triangulate
    return;
  }

  // Find the relative coordinates of each curve point in the cylinder reference
  // plane Only retain the max and min relative Z
  for (size_t i = 0; i < bounds.size(); i++) {

    double localMaxZ = -DBL_MAX;
    double localMinZ = DBL_MAX;

    for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {
      glm::dvec3 vv = bounds[i].curve.points[j] - cent;
      //					double dx = glm::dot(vecX, vv);
      //					double dy = glm::dot(vecY, vv);
      double dz = glm::dot(vecZ, vv);
      
      localMaxZ = std::max( localMaxZ, dz );
      localMinZ = std::min( localMinZ, dz );
    }

    outsideMostBoundaries.push( std::make_pair( localMaxZ, i ) );

    maxZ = std::max( maxZ, localMaxZ );
    minZ = std::min( minZ, localMinZ );
  }

   using Point = std::array<double, 2>;
  std::vector<std::vector<Point>> uvBoundaryValues;
  std::vector<ParameterVertex> vertices;

  // AFTP: back this per-face tessellation mesh (esp. edge_map, which allocates
  // a node per edge) with the thread scratch arena, rewound at function exit.
  // Byte-identical: the arena changes only where the mesh's nodes live.
  conway::ScratchArenaScope arenaScope;
  WingedEdgeMesh< ParameterVertex > mesh{ conway::ThreadScratchResource() };

  double zScale = 0.5 / ( maxZ - minZ );

  while ( !outsideMostBoundaries.empty() ) {

    std::vector<Point> points;

    size_t boundsIndex = outsideMostBoundaries.top().second;

    outsideMostBoundaries.pop();

    const IfcBound3D& bound = bounds[ boundsIndex ];

    if ( bound.curve.points.empty() ) {
      continue;
    }

    for ( const glm::dvec3& pt : bound.curve.points ) {

      glm::dvec3 vv = pt - cent;

      double dx = glm::dot(vecX, vv);
      double dy = glm::dot(vecY, vv);
      double dz = glm::dot(vecZ, vv);

      glm::dvec2 pInv =
        glm::normalize( glm::dvec2( dx, dy ) ) *
        ( 0.5 + ( dz - minZ ) * zScale );

      points.push_back({pInv.x, pInv.y});
      mesh.makeVertex( { pt, pInv } );
    } 

    uvBoundaryValues.push_back( points );
  }

#if (OUTPUT_SVG_DEBUG == 1) 

    // atomic: may run concurrently on the thread pool during staged face
    // finalization.
    static std::atomic< size_t > svgIndex = 0;

    size_t outputIndex = svgIndex++;

    std::ofstream svgFile( "cylinder_" + std::to_string( outputIndex ) + ".svg" );

    auto svgScale = []( double value ) {

      return 50 + ( 1024.0 * ( value + 1.0 ) / 2.0 ); 

    };

    svgFile << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
              << "width=\"1124\" height=\"1124\">\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"512\" style=\"stroke:rgb(255, 132, 0);stroke-width:2\" fill=\"none\"/>\n";

    svgFile << "<circle cx=\"" << svgScale( 0 ) << "\" cy=\"" << svgScale( 0 )
            << "\" r=\"256\" style=\"stroke:rgb(0, 0, 255);stroke-width:2\" fill=\"none\"/>\n";

    for ( const std::vector< Point >& loop  : uvBoundaryValues )
    {
      bool firstInLoop = true;
      
      glm::dvec2 lastPoint;

      for ( const Point &point : loop )
      {
        glm::dvec2 svgPoint = glm::dvec2( point[ 0 ], point[ 1 ] );

        if ( firstInLoop )
        {
          svgFile << "<polyline points=\"";
          firstInLoop = false;
        }
        else
        {
          svgFile << " ";
        }

        svgFile << svgScale( svgPoint.x ) << "," << svgScale( svgPoint.y );

        lastPoint = svgPoint;
      }

      svgFile << "\" style=\"fill:none;stroke:rgb(0,0,0);stroke-width:2\" />\n";
      
      for ( const Point &point : loop )
      {
        glm::dvec2 svgPoint = glm::dvec2( point[ 0 ], point[ 1 ] );
        
        svgFile << "<circle cx=\"" << svgScale( svgPoint.x ) << "\" cy=\"" << svgScale( svgPoint.y )
            << "\" r=\"3\" fill=\"red\"/>\n";
      }
    }

    svgFile << "</svg>\n";

    svgFile.close();

#endif

  std::vector<uint32_t> indices;
  {
    conway::AllocTagScope earcutTag( conway::AllocSite::Earcut );
    indices = mapbox::earcut<uint32_t>(uvBoundaryValues);
  }

  for (size_t i = 0; i < indices.size(); i += 3) {

    mesh.makeTriangle( 
      indices[ i  + 0 ], 
      indices[ i  + 1 ], 
      indices[ i  + 2 ] );
  }

  tesselate(
    mesh,
    [&]( const glm::dvec3& point, [[maybe_unused]]const glm::dvec2& from ) { 
      
      glm::dvec3 vv                = point - cent;
      double     dx                = glm::dot(vecX, vv);
      double     dy                = glm::dot(vecY, vv);
      double     dz                = glm::dot(vecZ, vv);
      glm::dvec2 inPlane           = glm::dvec2( dx, dy );
      glm::dvec2 normalizedInPlane = glm::normalize( from );

      glm::dvec2 newInPlane = normalizedInPlane * radius;

      return cent + newInPlane.x * vecX + newInPlane.y * vecY + vecZ * dz;
    },
    mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
    relativeDeflectionSquared( mesh, representationExtent ) );

  appendMeshToGeometry( mesh, geometry, sameSense );
}

// TODO: review and simplify
// ---------------------------------------------------------------------------
// Extrusion surfaces: exact developable unwrap with CDT trimming.
//
// The previous implementation ignored the trim loops, dropped the profile's
// z-coordinate, and swept in the wrong frame ("SIMPLE EXTRUSION, NOT
// TRIMMED") - in practice worse than the TriangulateBounds fallback these
// faces used to hit when the extractor failed to supply a profile. An
// extrusion surface is developable, so the unwrap is exact: s = arclength
// along the profile projected perpendicular to the axis, t = distance along
// the axis. Trim loops CDT in (s, t) with exact hole nesting; a closed
// profile makes s periodic and gets the same winding/gap/annulus treatment
// as the revolution unwrap. Any degeneracy falls back to TriangulateBounds,
// which reproduces the old fail-soft appearance.
// ---------------------------------------------------------------------------

inline void TriangulateExtrusion(Geometry &geometry,
                                 std::vector<IfcBound3D> &bounds,
                                 IfcSurface &surface,
                                 double representationExtent) {

  using namespace unwrap_detail;

  const std::vector< glm::dvec3 > &profilePoints =
    surface.ExtrusionSurface.Profile.curve.points;

  glm::dvec3 axis = surface.ExtrusionSurface.Direction;

  double axisLength = glm::length( axis );

  if ( profilePoints.size() < 2 || axisLength < 1e-12 ) {
    TriangulateBounds( geometry, bounds );
    return;
  }

  axis /= axisLength;

  // Profile projected onto the plane through the origin perpendicular to the
  // axis - the developable surface is fully described by (base curve, axis).
  std::vector< glm::dvec3 > basePoints;
  std::vector< double >     baseArc;

  basePoints.reserve( profilePoints.size() );
  baseArc.reserve( profilePoints.size() );

  double baseLength = 0.0;

  for ( const glm::dvec3 &point : profilePoints ) {

    glm::dvec3 projected = point - axis * glm::dot( point, axis );

    if ( !basePoints.empty() ) {
      baseLength += glm::distance( basePoints.back(), projected );
    }

    basePoints.push_back( projected );
    baseArc.push_back( baseLength );
  }

  if ( baseLength < 1e-12 ) {
    TriangulateBounds( geometry, bounds );
    return;
  }

  bool profileClosed =
    glm::distance( basePoints.front(), basePoints.back() ) <
      baseLength * 1e-9;

  // Closest point on the projected profile polyline: the snapped base point
  // and its normalized arclength parameter.
  auto closestOnBase = [&]( const glm::dvec3 &query ) {

    glm::dvec3 best      = basePoints[ 0 ];
    double     bestArc   = 0.0;
    double     bestDist2 = std::numeric_limits< double >::max();

    for ( size_t where = 0; where + 1 < basePoints.size(); ++where ) {

      const glm::dvec3 &from = basePoints[ where ];
      const glm::dvec3 &to   = basePoints[ where + 1 ];

      glm::dvec3 segment = to - from;
      double     length2 = glm::dot( segment, segment );

      double t =
        length2 > 0 ?
          std::clamp( glm::dot( query - from, segment ) / length2, 0.0, 1.0 ) :
          0.0;

      glm::dvec3 candidate = from + segment * t;
      double     dist2     = glm::dot( query - candidate, query - candidate );

      if ( dist2 < bestDist2 ) {
        bestDist2 = dist2;
        best      = candidate;
        bestArc   = baseArc[ where ] + t * std::sqrt( length2 );
      }
    }

    return std::pair< glm::dvec3, double >( best, bestArc / baseLength );
  };

  // Exact on-surface projection for the adaptive refinement: keep the axial
  // coordinate, snap the cross-axis position to the profile.
  auto surfaceProjection = [&]( const glm::dvec3 &point ) {

    double t = glm::dot( point, axis );

    glm::dvec3 foot = closestOnBase( point - axis * t ).first;

    return foot + axis * t;
  };

  // Scale-aware refinement threshold (relativeDeflectionSquared): stop
  // subdividing once midpoint deviation drops below 0.1% of the face's
  // extent (the same visual criterion as the revolution ring count). The
  // old shared absolute MAX_DEFLECTION read as "refine until the budget
  // runs out" at large model scales, and with 1,370 extrusion faces on the
  // jet compressor alone that doubled the model's triangle count for
  // invisible gains.
  auto refineAndAppend = [&]( WingedEdgeMesh< glm::dvec3 > &mesh ) {

    tesselate(
      mesh,
      surfaceProjection,
      mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
      relativeDeflectionSquared( mesh, representationExtent ) );

    appendMeshToGeometry( mesh, geometry );
  };

  // --- Boundary loops -> (s, t): s the normalized profile arclength
  // (periodic for a closed profile), t the axial distance.

  std::vector< std::vector< BoundaryPoint > > loops;

  loops.reserve( bounds.size() );

  auto parameterDuplicate = [&]( const BoundaryPoint &a, double s, double t ) {

    double deltaS =
      profileClosed ? wrapDeltaPi( s - a.theta ) : ( s - a.theta );

    return std::abs( deltaS ) < 1e-12 && std::abs( t - a.phi ) < 1e-12;
  };

  for ( const IfcBound3D &bound : bounds ) {

    std::vector< BoundaryPoint > loop;

    loop.reserve( bound.curve.points.size() );

    for ( const glm::dvec3 &point : bound.curve.points ) {

      if ( !std::isfinite( point.x ) || !std::isfinite( point.y ) ||
           !std::isfinite( point.z ) ) {
        continue;
      }

      double t = glm::dot( point, axis );
      double s = closestOnBase( point - axis * t ).second;

      // A closed profile's s is periodic - store it as an angle so the
      // shared winding/gap helpers apply verbatim.
      if ( profileClosed ) {
        s = positiveMod2Pi( s * kTwoPi );
      }

      if ( !loop.empty() && parameterDuplicate( loop.back(), s, t ) ) {
        continue;
      }

      loop.push_back( BoundaryPoint{ point, s, t } );
    }

    while ( loop.size() > 1 &&
            parameterDuplicate( loop.front(), loop.back().theta, loop.back().phi ) ) {
      loop.pop_back();
    }

    if ( loop.size() >= 3 ) {
      loops.push_back( std::move( loop ) );
    }
  }

  if ( loops.empty() ) {
    TriangulateBounds( geometry, bounds );
    return;
  }

  // Axial range for normalizing t in the 2D layouts.
  double loT = std::numeric_limits< double >::max();
  double hiT = std::numeric_limits< double >::lowest();

  for ( const std::vector< BoundaryPoint > &loop : loops ) {
    for ( const BoundaryPoint &point : loop ) {
      loT = std::min( loT, point.phi );
      hiT = std::max( hiT, point.phi );
    }
  }

  double spanT = std::max( hiT - loT, 1e-12 );

  bool wrapsS = false;

  if ( profileClosed ) {

    bool windsS = false;

    std::vector< double > sSamples;

    for ( const std::vector< BoundaryPoint > &loop : loops ) {

      double sumS = 0.0;

      for ( size_t where = 0, count = loop.size(); where < count; ++where ) {
        sumS += wrapDeltaPi( loop[ ( where + 1 ) % count ].theta - loop[ where ].theta );
        sSamples.push_back( loop[ where ].theta );
      }

      if ( std::llround( sumS / kTwoPi ) != 0 ) {
        windsS = true;
      }
    }

    auto [ sGap, sCut ] = largestCircularGap( sSamples );

    wrapsS = windsS || sGap < 1e-3;

    if ( !wrapsS ) {
      // Open up the covered arc: rebase every s off the branch cut so the
      // rectangle layout below sees a contiguous interval.
      for ( std::vector< BoundaryPoint > &loop : loops ) {
        for ( BoundaryPoint &point : loop ) {
          point.theta = positiveMod2Pi( point.theta - sCut );
        }
      }
    }
  }

  std::vector< std::vector< glm::dvec2 > > loops2D;

  loops2D.reserve( loops.size() );

  if ( wrapsS ) {

    // Annulus: angle = s, radius = 1 + normalized t in [1, 2].
    for ( const std::vector< BoundaryPoint > &loop : loops ) {

      std::vector< glm::dvec2 > loop2D;

      loop2D.reserve( loop.size() );

      for ( const BoundaryPoint &point : loop ) {

        double radius = 1.0 + ( point.phi - loT ) / spanT;

        loop2D.emplace_back(
          radius * std::cos( point.theta ),
          radius * std::sin( point.theta ) );
      }

      loops2D.push_back( std::move( loop2D ) );
    }

  } else {

    double loS = std::numeric_limits< double >::max();
    double hiS = std::numeric_limits< double >::lowest();

    for ( const std::vector< BoundaryPoint > &loop : loops ) {
      for ( const BoundaryPoint &point : loop ) {
        loS = std::min( loS, point.theta );
        hiS = std::max( hiS, point.theta );
      }
    }

    double spanS = std::max( hiS - loS, 1e-12 );

    for ( const std::vector< BoundaryPoint > &loop : loops ) {

      std::vector< glm::dvec2 > loop2D;

      loop2D.reserve( loop.size() );

      for ( const BoundaryPoint &point : loop ) {
        loop2D.emplace_back(
          ( point.theta - loS ) / spanS,
          ( point.phi - loT ) / spanT );
      }

      loops2D.push_back( std::move( loop2D ) );
    }
  }

  // --- Weld, dedupe, CDT - identical machinery to the torus/revolution
  // unwraps; failures fall back to the pre-existing bounds triangulation.

  std::vector< CDT::V2d< double > > cdtVertices;
  std::vector< glm::dvec3 >         cdtWorld;
  std::vector< CDT::Edge >          cdtEdges;

  std::map< std::pair< long long, long long >, uint32_t > weld;
  std::set< std::pair< uint32_t, uint32_t > >             edgeSet;

  bool inputFinite = true;

  auto weldVertex = [&]( const glm::dvec2 &position, const glm::dvec3 &world ) {

    std::pair< long long, long long > key(
      std::llround( position.x * 1e9 ),
      std::llround( position.y * 1e9 ) );

    auto [ found, isNew ] =
      weld.try_emplace( key, static_cast< uint32_t >( cdtVertices.size() ) );

    if ( isNew ) {
      if ( !std::isfinite( position.x ) || !std::isfinite( position.y ) ) {
        inputFinite = false;
      }
      cdtVertices.emplace_back( position.x, position.y );
      cdtWorld.push_back( world );
    }

    return found->second;
  };

  for ( size_t which = 0; which < loops.size(); ++which ) {

    const std::vector< BoundaryPoint > &loop   = loops[ which ];
    const std::vector< glm::dvec2 >    &loop2D = loops2D[ which ];

    for ( size_t where = 0, count = loop.size(); where < count; ++where ) {

      size_t next = ( where + 1 ) % count;

      uint32_t v1 = weldVertex( loop2D[ where ], loop[ where ].world );
      uint32_t v2 = weldVertex( loop2D[ next ], loop[ next ].world );

      if ( v1 == v2 ) {
        continue;
      }

      std::pair< uint32_t, uint32_t > ordered(
        std::min( v1, v2 ), std::max( v1, v2 ) );

      if ( edgeSet.insert( ordered ).second ) {
        cdtEdges.emplace_back( v1, v2 );
      }
    }
  }

  if ( !inputFinite || cdtEdges.size() < 3 ) {
    TriangulateBounds( geometry, bounds );
    return;
  }

  CDT::Triangulation< double > triangulation(
    CDT::VertexInsertionOrder::Auto,
    CDT::IntersectingConstraintEdges::TryResolve, 0 );

  try
  {
    conway::AllocTagScope cdtTag( conway::AllocSite::Cdt );
    triangulation.insertVertices( cdtVertices );
    triangulation.insertEdges( cdtEdges );
    triangulation.eraseOuterTrianglesAndHoles();
  }
  catch ( const std::exception &e )
  {
    Logger::logError( "CDT Exception (extrusion unwrap): %s", e.what() );
    TriangulateBounds( geometry, bounds );
    return;
  }

  WingedEdgeMesh< glm::dvec3 > mesh;

  mesh.vertices.reserve( cdtWorld.size() );

  for ( const glm::dvec3 &world : cdtWorld ) {
    mesh.vertices.push_back( world );
  }

  // No per-triangle orientation pass: appendMeshToGeometry normalizes
  // winding by dominant-plane projection either way.
  for ( const CDT::Triangle &triangle : triangulation.triangles ) {

    auto [ cdtv1, cdtv2, cdtv3 ] = triangle.vertices;

    if ( cdtv1 == cdtv2 || cdtv2 == cdtv3 || cdtv3 == cdtv1 ) {
      continue;
    }

    // TryResolve can add split vertices past the input set; they have no world-space lift here, so skip the sliver triangles that touch them.
    if (
      cdtv1 >= cdtWorld.size() ||
      cdtv2 >= cdtWorld.size() ||
      cdtv3 >= cdtWorld.size() ) {
      continue;
    }

    mesh.makeTriangle( cdtv1, cdtv2, cdtv3 );
  }

  if ( mesh.triangles.empty() ) {
    TriangulateBounds( geometry, bounds );
    return;
  }

  refineAndAppend( mesh );
}

constexpr size_t INVERSE_GRID_SIDE   = 8.0; 
constexpr double INVERSE_GRID_SIZE_D = static_cast< double >( INVERSE_GRID_SIDE );
constexpr double INVERSE_GRID_FACTOR = 1.0 / ( INVERSE_GRID_SIZE_D - 1.0 );
/**
 * How close, relative to the surface's own extent, the NURBS inverse solve
 * has to land before it calls a uv converged.
 *
 * This used to be an absolute 0.001 model units, the last of the absolute
 * tolerances the relative-deflection work replaced elsewhere (see
 * relativeCurveDeflectionSquared and relativeDeflectionSquared, both of
 * which are 0.1% of the entity's own extent). Absolute is wrong here for
 * the same reason it was wrong there, but the failure is louder, because
 * the number this one feeds is a NOISE FLOOR under a target that IS
 * relative: `tesselate` refines an edge while
 *
 *     | surface( uvMid ) - chordMid |  >  0.001 * faceDiagonal
 *
 * and the seed vertices carry uv from this solve, so that difference can
 * never fall below the solve's own residual. On a face smaller than
 * MAX_ERROR / 0.001 - i.e. anything under a metre in millimetre units -
 * the residual is the larger of the two and the criterion is unsatisfiable:
 * the loop refines until the 32x MAX_TRIANGLE_AMPLIFACTION budget runs out,
 * then stops, having spent all of it chasing noise.
 *
 * Arty_Z7.stp is the case that surfaced it (bldrs-ai/conway#564). Its
 * silkscreen is modelled as 1,189 extruded-glyph solids whose stroke
 * sidewalls are 10,224 degree-(3,1) b-spline faces with a MEDIAN DIAGONAL
 * OF 0.126 mm. Target 0.126 um, residual ~1 um: the residual runs 7.7x the
 * target at p50 and 62x at worst, every one of those faces amplifies
 * exactly 32x, and 345,645 seed triangles become 10,902,597 - 96% of the
 * model's geometry payload and 89% of its geometry time. The bilinear
 * (degree 1,1) faces in the same file are the control: same code path,
 * residual 4e-5 of the target, amplification 1.0.
 *
 * 1e-6 puts the residual three decades under the 1e-3 deflection target, so
 * the target is reachable at every scale and the refinement stops on real
 * curvature. The floor is the 2^-24 quantisation grid of IfcCurve::Add3d,
 * matching the deflection helpers, and also guards a degenerate
 * zero-extent surface.
 */
constexpr double RELATIVE_INVERSE_ERROR = 1e-6;
constexpr double MIN_INVERSE_ERROR      = 0x1p-24;

/** Armijo sufficient-decrease coefficient for the solve's line search. */
constexpr double ARMIJO_COEFFICIENT  = 0.001;

constexpr double ALPHA_ERROR         = 1e-6;
constexpr double MIN_STEP            = 1e-9;

/**
 * Iteration ceiling for the solve.
 *
 * Worth knowing what it now means: against the old absolute target every call
 * converged (382,831 of 382,919 on Arty_Z7, mean 8.8 iterations). Against a
 * target scaled to the surface, 72.8% of calls reach this ceiling still
 * improving - they land around 5% of the tessellation deflection target,
 * which is all the refinement needs, but three decades under the surface is
 * more than 50 Gauss-Newton steps will reach. That is 15.3M iterations where
 * there were 3.4M, and it is why the net win on Arty_Z7 is 16% rather than
 * more. Raising this ceiling would cost time and buy nothing;
 * RELATIVE_INVERSE_ERROR is the dial worth measuring if the solve ever
 * needs to be cheaper.
 */
constexpr size_t MAX_ITERATION       = 50;

struct RationalNurbsInverseMethod {

  glm::dvec3 grid[ INVERSE_GRID_SIDE ][ INVERSE_GRID_SIDE ];
  const tinynurbs::RationalSurface3d& surface;

  /** Allocation-free sampler shared by the inverse solve and tessellation. */
  RationalSurfaceEvaluator evaluator;

  glm::dvec2 min_extent;
  glm::dvec2 max_extent;

  /** Convergence target for the solve, scaled to this surface's extent. */
  double convergence_error;

  /**
   * The uv returned by the last operator() call, offered as a seed candidate
   * to the next one - see the continuity-seed comment in operator().
   *
   * Per-face state on a per-face object: TriangulateBspline constructs one
   * of these per face and inverts that face's boundary on one thread, so
   * there is no sharing to guard.
   */
  glm::dvec2 previousSolution_    = glm::dvec2( 0.0 );
  bool       hasPreviousSolution_ = false;

  /**
   * Is the support surface CLOSED in u / in v - does the first row (column)
   * of control points coincide with the last?
   *
   * When it is, the two ends of that parameter's domain name the SAME 3D
   * curve, so a point on the seam has two exact preimages and no residual can
   * tell them apart. See seamContinuousSolution().
   */
  bool closedU_ = false;
  bool closedV_ = false;

  /**
   * Drop the continuity memo, so the next query is seeded from the grid
   * alone.
   *
   * Call this at every boundary-loop transition. The continuity seed's whole
   * justification is that consecutive queries are a CHORD apart on one
   * ordered polyline; across two different trim loops of the same face there
   * is no such relationship, and the last point of one loop can still be
   * near enough in 3D to beat every grid sample and win the seed for the
   * first point of the next.
   *
   * That is not merely a worse seed, and it is the one failure the retry in
   * operator() cannot catch (codex, second pass on conway-geom#186). On a
   * folded, self-near or multiply-covered support surface the cross-loop
   * seed can converge cleanly onto a DIFFERENT sheet or preimage - and since
   * the retry is gated on non-convergence, a successful-but-wrong result
   * never triggers it. The bad uv then seeds the rest of that loop, and
   * earcut receives a misplaced hole, which on an inner trim loop is exactly
   * the defect conway#594 is about, relocated rather than removed.
   *
   * Measured on the locally materialised corpus, the cross-loop seed never
   * actually won: see the count in conway-geom#186. So this is a latent bug
   * being closed, not an observed one - but nothing about an unseen model
   * bounds it, and the invariant is cheap to enforce and was previously only
   * asserted.
   */
  void resetContinuity() {
    hasPreviousSolution_ = false;
  }

  RationalNurbsInverseMethod( const tinynurbs::RationalSurface3d& srf )
    : surface( srf ), evaluator( srf ) {

    size_t degreeU = static_cast< size_t >( srf.degree_u );
    size_t degreeV = static_cast< size_t >( srf.degree_v );
   
    if ( srf.knots_u.size() == 2 ) {
   
      degreeU = 0;

    }

    if ( srf.knots_v.size() == 2 ) {

      degreeV = 0;

    }

    // The valid parameter domain of a clamped NURBS is
    // [ knots[ degree ], knots[ knotCount - degree - 1 ] ] per axis. These
    // comparisons were previously inverted (size < degree, which is never
    // true for a valid knot vector), pinning the domain to [0,1]^2. IFC
    // exporters normalise knots to 0..1 so it went unnoticed, but STEP
    // (OCCT) writes real-valued knot ranges (e.g. u in [0, 200] down the
    // length of a cylinder), and the [0,1] pin collapsed every rational
    // b-spline face to a sliver - see conway#350 (AS1 rod invisible).
    min_extent = glm::dvec2(
      srf.knots_u.size() > degreeU ? srf.knots_u[ degreeU ] : 0.0,
      srf.knots_v.size() > degreeV ? srf.knots_v[ degreeV ] : 0.0 );

    // Unsigned underflow when knots.size() < degree + 1 makes uM/vM huge,
    // so the bounds check below also rejects that malformed case.
    size_t uM = srf.knots_u.size() - ( degreeU + 1 );
    size_t vM = srf.knots_v.size() - ( degreeV + 1 );

    max_extent = glm::dvec2(
      uM < srf.knots_u.size() ? srf.knots_u[ uM ] : 1.0,
      vM < srf.knots_v.size() ? srf.knots_v[ vM ] : 1.0 );


    for ( size_t i = 0; i < INVERSE_GRID_SIDE; ++i ) {
      for ( size_t j = 0; j < INVERSE_GRID_SIDE; ++j ) {

        glm::dvec2 uv = 
          min_extent + 
          ( max_extent - min_extent ) * 
          glm::dvec2( 
            static_cast< double >( i ) * INVERSE_GRID_FACTOR,
            static_cast< double >( j ) * INVERSE_GRID_FACTOR );

        grid[ i ][ j ] = evaluator.point( uv.x, uv.y );
      }
    }

    // Seeded off the grid, which spans the WHOLE surface, rather than off the
    // trimmed face: the solve is constructed before any bound is projected,
    // so the face's extent is not known yet. That over-estimates a face that
    // is an island on a larger surface, and over-estimating is NOT
    // automatically safe - too loose is its own failure mode. Two thresholds
    // bound it, both in terms of surfaceDiagonal / trimDiagonal:
    //
    //   >= 1e3  the residual reaches the tessellation deflection target
    //           (RELATIVE_INVERSE_ERROR / RELATIVE_DEFLECTION_FACTOR), i.e.
    //           conway-geom#178's own bug returning at a bigger ratio;
    //   >= 1e6  the tolerance reaches the trim's own diameter, so every
    //           boundary point starts within one tolerance of the same grid
    //           sample, operator() iterates zero times, returns one uv for the
    //           entire loop, and the uv earcut degenerates.
    //
    // Measured over every b-spline face in the locally materialised corpus -
    // 10,422 faces across Arty_Z7, Right_Hand, nist_ctc_02 and
    // FM_ARC_DigitalHub - the maximum surfaceDiagonal / trimDiagonal is
    // 2.273: three to six decades of margin, so no corpus model reaches
    // either threshold. That is corpus absence, not safety, and a valid
    // unseen file can sit anywhere on that ratio - which is why
    // TriangulateBspline calls boundConvergenceToTrim() below before it
    // projects anything, capping this seed on the face actually being
    // tessellated.
    glm::dvec3 gridMin( std::numeric_limits< double >::max() );
    glm::dvec3 gridMax( std::numeric_limits< double >::lowest() );

    for ( size_t i = 0; i < INVERSE_GRID_SIDE; ++i ) {
      for ( size_t j = 0; j < INVERSE_GRID_SIDE; ++j ) {

        gridMin = glm::min( gridMin, grid[ i ][ j ] );
        gridMax = glm::max( gridMax, grid[ i ][ j ] );
      }
    }

    convergence_error =
      std::max(
        MIN_INVERSE_ERROR,
        glm::distance( gridMin, gridMax ) * RELATIVE_INVERSE_ERROR );

    // Closure, from the control grid rather than by sampling: a b-spline
    // surface is closed in a parameter exactly when the first and last rows
    // (columns) of control points define the SAME curve, and that is a
    // property of the definition rather than a measurement of it.
    //
    // On a RATIONAL surface the Cartesian control points alone do not settle
    // that. The endpoint curves share this surface's knot vector and degree,
    // so with control points P and weights w they are
    //
    //   C( t ) = sum_j N_j( t ) w_j P_j / sum_j N_j( t ) w_j
    //
    // and two rows of coincident P with DIFFERENT w are different rational
    // curves - they agree only where the basis is interpolating, i.e. at the
    // ends. Comparing P alone therefore asserts more than it checks, and
    // getting it wrong is not cosmetic: closedU_/closedV_ now gate two
    // separate paths - the seam crossing in solveFromSeed and the
    // full-coverage grid in tryFullCoverageSeamGrid - so a false positive
    // makes the descent wrap across a genuine geometric discontinuity and
    // return the wrong branch. Weights are compared as well.
    //
    // The weight test is RELATIVE, against the same RELATIVE_INVERSE_ERROR
    // the distance tolerance below is built from, because a weight is
    // dimensionless and cannot be judged against a length. No new constant
    // enters.
    //
    // Equal ( P, w ) is sufficient for the curves to coincide; it is not
    // necessary, because scaling a whole row of weights by a common factor
    // leaves the rational curve unchanged. That direction is deliberately not
    // chased: the test errs toward "not closed", which costs a face nothing
    // more than the behaviour it had before either path existed.
    //
    // "Coincide" is judged at the same tolerance the solve uses to call two
    // points the same point, so no new constant is introduced there either.
    // On the surface this was written for the rows are bitwise equal - the
    // seam of `B_SPLINE_SURFACE_WITH_KNOTS #1553` matches to 0.0 across all
    // 84 columns, and its weights match to 0.0 as well - but an exporter that
    // round-trips through floats needs the slack, and a surface that is
    // merely NEAR closed is one where the wrap below would be rejected by its
    // own residual test anyway.
    const double closureTolerance =
      std::max( MIN_INVERSE_ERROR,
                glm::distance( gridMin, gridMax ) * RELATIVE_INVERSE_ERROR );

    const size_t rows = surface.control_points.rows();
    const size_t cols = surface.control_points.cols();

    // Absent or mis-sized weights mean the surface is evaluated as
    // polynomial, so every weight is 1 and the weight test is vacuous.
    const bool haveWeights =
      surface.weights.rows() == rows && surface.weights.cols() == cols;

    const auto sameWeight =
      [ & ]( double left, double right ) {

        if ( !haveWeights ) {
          return true;
        }

        return std::abs( left - right ) <=
               std::max( std::abs( left ), std::abs( right ) ) *
                 RELATIVE_INVERSE_ERROR;
      };

    if ( rows > 1 && cols > 0 ) {

      closedU_ = true;

      for ( size_t col = 0; col < cols; ++col ) {

        if ( glm::distance( surface.control_points( 0, col ),
                            surface.control_points( rows - 1, col ) ) >
               closureTolerance ||
             !sameWeight( haveWeights ? surface.weights( 0, col ) : 1.0,
                          haveWeights ? surface.weights( rows - 1, col ) : 1.0 ) ) {

          closedU_ = false;
          break;
        }
      }
    }

    if ( cols > 1 && rows > 0 ) {

      closedV_ = true;

      for ( size_t row = 0; row < rows; ++row ) {

        if ( glm::distance( surface.control_points( row, 0 ),
                            surface.control_points( row, cols - 1 ) ) >
               closureTolerance ||
             !sameWeight( haveWeights ? surface.weights( row, 0 ) : 1.0,
                          haveWeights ? surface.weights( row, cols - 1 ) : 1.0 ) ) {

          closedV_ = false;
          break;
        }
      }
    }

  }

  /**
   * Cap the convergence target on the trimmed face rather than the whole
   * support surface.
   *
   * Called by TriangulateBspline once the face's 3D bounds are known and
   * before the first inversion, so every call to operator() on this face
   * sees the capped target. The cap only ever tightens: a face larger than
   * its own support surface is not a thing, and loosening would reintroduce
   * exactly what conway-geom#178 fixed.
   *
   * Capping rather than raising MIN_INVERSE_ERROR is deliberate. The floor
   * is a noise floor - the 2^-24 quantisation grid of IfcCurve::Add3d - and
   * raising it would loosen the solve on every face at every scale, whereas
   * this narrows it only on the faces whose support surface over-states
   * them. `trimDiagonal` costs one bbox pass over the same bound points the
   * projection loop walks next, and no extra inversion.
   *
   * Worth keeping from #178's measurement, because it says this cap is not
   * the whole story: the nearest any corpus face came to the >= 1e3
   * threshold was 0.121 of it, a 0.49mm face in metre-unit Right_Hand - and
   * that came from MIN_INVERSE_ERROR, not from a surface/trim mismatch (that
   * face's ratio is 0.9996). This cap does nothing for that case. If the
   * floor ever needs revisiting, that is the face to revisit it against.
   *
   * @param trimDiagonal 3D diagonal of the bounding box of the face's bound
   *                     curve points, in the same (post-`scaling`) units the
   *                     points are inverted in. Zero or non-finite leaves the
   *                     target at the MIN_INVERSE_ERROR floor, which is what
   *                     a zero-extent face gets from the grid path too.
   */
  void boundConvergenceToTrim( double trimDiagonal ) {

    convergence_error =
      std::max(
        MIN_INVERSE_ERROR,
        std::min( convergence_error, trimDiagonal * RELATIVE_INVERSE_ERROR ) );
  }

  glm::dvec2 operator()( const glm::dvec3& point ) {

    glm::dvec2 gridGuess;
    glm::dvec3 gridPoint;
    double     gridDistance2 = DBL_MAX;

    // Take initial guess from the grid.
    for ( size_t i = 0; i < INVERSE_GRID_SIDE; ++i ) {
      for ( size_t j = 0; j < INVERSE_GRID_SIDE; ++j ) {

        glm::dvec3 deltaP = grid[ i ][ j ] - point;

        double distance2 = glm::dot( deltaP, deltaP );

        if ( distance2 < gridDistance2 ) {

          gridGuess =
            min_extent + 
            ( max_extent - min_extent ) * 
            glm::dvec2( 
              static_cast< double >( i ) * INVERSE_GRID_FACTOR,
              static_cast< double >( j ) * INVERSE_GRID_FACTOR );

          gridPoint     = grid[ i ][ j ];
          gridDistance2 = distance2;
        }
      }
    }

    glm::dvec2 bestGuess    = gridGuess;
    glm::dvec3 bestPoint    = gridPoint;
    double     minDistance2 = gridDistance2;

    bool usedContinuitySeed = false;

    // Continuity seed: the uv this method returned for the PREVIOUS query,
    // offered as one more seed candidate and kept only when it is nearer in
    // 3D than every grid sample.
    //
    // The callers invert an ORDERED boundary polyline, so consecutive
    // queries are a chord apart - 0.116 model units on the Orbiter thread
    // ribbons - while the 8x8 grid spans the whole support surface, whose
    // own diagonal is 20. A seed three orders of magnitude closer is not a
    // speed optimisation; on a swept surface it is the difference between
    // converging and not.
    //
    // Why the grid alone fails there (bldrs-ai/conway#594). A thread
    // modelled as a swept NURBS ribbon runs ~7 turns down v, so an 8x8 grid
    // places barely one sample per turn and the nearest one routinely
    // belongs to a NEIGHBOURING turn - which passes within 0.118 of the
    // query. Gauss-Newton then descends in that turn's basin, and 50
    // iterations of a basin that does not contain the answer end at
    // MAX_ITERATION with a residual of ~0.9 rather than the 2e-5 target.
    // Measured on face #51059 of solid 971: median residual 0.880, only
    // 18.1% of 5,089 boundary points reaching the target, 238 single-turn
    // jumps and 137 v-direction reversals in what should be a ribbon
    // boundary with two. Model-wide, 71% of 53,981 solves failed to
    // converge. With this seed the same face reaches median 1.3e-5, 100%
    // converged, and ONE reversal.
    //
    // Kept as a candidate rather than used unconditionally: a bare
    // previous-point chain has no way back once one query lands badly, and
    // measured worse tails than the grid on two of twelve faces.
    //
    // Note what this comparison does and does NOT establish. It picks the
    // nearer SEED, which is not the same as the nearer answer - see the
    // retry below, which is where the grid actually becomes an escape hatch
    // rather than a discarded alternative.
    //
    // Scoped to ONE ordered polyline: callers must call resetContinuity()
    // at every boundary-loop transition, and TriangulateBspline does. The
    // first commit left the memo running across loops, reasoning that an
    // unrelated candidate is simply further away and loses the comparison.
    // That reasoning is wrong, and it is wrong in the one way the retry
    // below cannot cover - see resetContinuity().
    if ( hasPreviousSolution_ ) {

      glm::dvec3 previousPoint =
        evaluator.point( previousSolution_.x, previousSolution_.y );

      glm::dvec3 deltaP = previousPoint - point;

      double distance2 = glm::dot( deltaP, deltaP );

      if ( distance2 < minDistance2 ) {

        bestGuess          = previousSolution_;
        bestPoint          = previousPoint;
        minDistance2       = distance2;
        usedContinuitySeed = true;
      }
    }

    double     solvedDistance2 = 0.0;
    glm::dvec2 solution        =
      solveFromSeed( point, bestGuess, bestPoint, minDistance2, solvedDistance2 );

    // Retry from the grid seed when the continuity seed did not get there,
    // and keep whichever solve ended NEARER - by final residual, not by seed
    // distance.
    //
    // This is what makes "the grid is still the escape hatch" a mechanism
    // rather than a claim, and it is a correctness fix rather than a
    // refinement (codex on conway-geom#186). A nearer seed does not imply a
    // nearer answer: Armijo guarantees decrease only from the seed actually
    // chosen, so a previous solution that is closer in 3D but sits in a
    // poorer basin displaces a grid seed whose outcome is never computed and
    // therefore never compared. Choosing on seed distance alone made that
    // trade blind, and on an unseen model nothing bounds how bad the
    // discarded alternative was.
    //
    // Measured before this retry existed: 178 of 721 b-spline faces on
    // Orbiter ended at a HIGHER residual than the grid-only solve. That
    // population is exactly this failure, and comparing final residuals
    // retires it.
    //
    // Cost is bounded to the calls that need it - only a non-converged solve
    // that took the continuity seed retries, and it retries once.
    if ( usedContinuitySeed &&
         solvedDistance2 > convergence_error * convergence_error ) {

      double     gridSolvedDistance2 = 0.0;
      glm::dvec2 gridSolution        =
        solveFromSeed(
          point, gridGuess, gridPoint, gridDistance2, gridSolvedDistance2 );

      if ( gridSolvedDistance2 < solvedDistance2 ) {

        solution        = gridSolution;
        solvedDistance2 = gridSolvedDistance2;
      }
    }

    solution = seamContinuousSolution( point, solution );

    previousSolution_    = solution;
    hasPreviousSolution_ = true;

    return solution;
  }

 private:

  /**
   * On a surface closed in u or v, re-express a solved uv on the branch that
   * is CONTINUOUS with the previous boundary point.
   *
   * Why a converged solve can still be the wrong answer. Where a surface is
   * closed, the two ends of that parameter's domain are the same 3D curve:
   * S( uMin, v ) is identically S( uMax, v ). Every point on that seam
   * therefore has TWO exact preimages, and both have the same residual - zero
   * - so the retry in operator() cannot see the difference. Its gate asks
   * "did this converge?", and the wrong answer converges perfectly. That is
   * the failure resetContinuity()'s comment names as the one the retry cannot
   * catch, met here in its own right rather than across loops.
   *
   * It matters because the consumer is EARCUT, which reads the uv polyline as
   * a polygon. A trimmed face on a closed surface is bounded by the seam
   * TWICE - once at each end of the domain - and in this file that is not an
   * inference: `ADVANCED_FACE #50977` of `MANIFOLD_SOLID_BREP #970` ("Spring
   * v1") on `step/conor/Orbiter_v1.1_Gear_7.5.step` has an EDGE_LOOP of four
   * edges of which two are the SAME `EDGE_CURVE #29690`, once forward and
   * once reversed. Those two legs are identical in 3D and can only be told
   * apart by which side of the seam they are on. Solved independently they
   * both land on whichever side the seed happened to favour, the rectangle
   * earcut should have received collapses onto its own edge, and the
   * triangulation that comes back joins opposite sides of the coil: 404 edges
   * longer than 3.0 on a spring whose outer diameter is 6.9, and 54 zero-area
   * triangles out of 287 (bldrs-ai/conway#611).
   *
   * The correction is a branch choice, not a tolerance. A jump of more than
   * HALF the parameter's own period between consecutive points on one ordered
   * polyline cannot be a chord; the shorter way round is the continuous one,
   * so the solution is re-expressed there. Half a period is the standard
   * branch-continuity criterion and is read from the surface's own knot
   * domain, so nothing here is scaled to a model.
   *
   * Two things keep it from inventing geometry:
   *
   * - It runs only where the control grid says the surface is CLOSED, so the
   *   two representatives really are the same point.
   * - It is accepted only if the re-expressed uv is itself a CONVERGED fit -
   *   within the same convergence_error the solve uses to call a uv solved,
   *   or no worse than the solved uv, whichever is the weaker demand. So this
   *   can move a uv onto a different branch, but never onto one that does not
   *   answer the query.
   *
   * That last test is deliberately not "at least as good as the solved uv".
   * The descent CLAMPS into the domain and stops at a local optimum, so a
   * seam point typically converges a whisker INSIDE the bound - u =
   * -3.141592633 against a domain end of -3.141592654 - and the wrapped
   * representative, which lands exactly ON the far bound, is then worse by
   * about 1e-8. Demanding "no worse" rejected the wrap on precisely the
   * points it exists for: measured on `#50977`, the strict form left one
   * uncorrected flip at boundary point 32 of 175, and one is enough, because
   * every later point on that leg takes its branch from it.
   *
   * Scoped to one ordered polyline by the same hasPreviousSolution_ memo the
   * continuity seed uses, and cleared by the same resetContinuity().
   *
   * @param point    The 3D point being inverted.
   * @param solution The uv the descent returned.
   * @return The same point's uv, on the branch continuous with the previous
   *   boundary point, or `solution` unchanged when there is nothing to fix.
   */
  glm::dvec2 seamContinuousSolution(
      const glm::dvec3& point,
      glm::dvec2        solution ) const {

    if ( !hasPreviousSolution_ || ( !closedU_ && !closedV_ ) ) {

      return solution;
    }

    const glm::dvec2 period = max_extent - min_extent;

    glm::dvec2 candidate = solution;

    for ( int axis = 0; axis < 2; ++axis ) {

      if ( ( axis == 0 && !closedU_ ) || ( axis == 1 && !closedV_ ) ) {
        continue;
      }

      if ( !( period[ axis ] > 0.0 ) ) {
        continue;
      }

      const double delta = solution[ axis ] - previousSolution_[ axis ];

      if ( std::abs( delta ) <= period[ axis ] * 0.5 ) {
        continue;
      }

      // Toward the previous point by exactly one period, then back into the
      // domain. On a closed surface the clamp lands on the opposite bound,
      // which names the identical 3D point - that is what makes this exact
      // rather than approximate.
      candidate[ axis ] =
        glm::clamp(
          solution[ axis ] - std::copysign( period[ axis ], delta ),
          min_extent[ axis ],
          max_extent[ axis ] );
    }

    if ( candidate == solution ) {

      return solution;
    }

    const glm::dvec3 solvedPoint    = evaluator.point( solution.x, solution.y );
    const glm::dvec3 candidatePoint = evaluator.point( candidate.x, candidate.y );

    return glm::distance( candidatePoint, point ) <=
           std::max( glm::distance( solvedPoint, point ), convergence_error ) ?
             candidate : solution;
  }

 public:

  /**
   * The displacement a trial step actually made, per axis.
   *
   * Armijo compares a MEASURED decrease against a PREDICTED one, and the two
   * have to describe the same step. The two axis kinds differ:
   *
   * - On an OPEN axis the domain edge truncates the step, so the movement is
   *   whatever the clamp left. Crediting the full request there is what made
   *   the descent stall at a domain corner and report its own seed as the
   *   answer.
   * - On a CLOSED axis nothing is truncated. S( u + period, v ) IS S( u, v ),
   *   so evaluating at the wrapped representative evaluates the requested
   *   point exactly, and the movement is the whole request.
   *
   * Reducing a closed axis to its shortest representative instead - taking
   * `from - to` and subtracting whole periods - is wrong for a step wider than
   * half a period, and `deltaUV` is uncapped so those steps are reachable. The
   * reduction hands back a component pointing the OTHER WAY from the one
   * requested, `dot( displacement, jte )` goes negative, the Armijo threshold
   * rises ABOVE the current residual, and the trial that gets accepted is one
   * that INCREASES it. Returning the request removes that case rather than
   * guarding it, because on a closed axis the request is what happened.
   *
   * Public, and split out of the line search, so it can be tested directly.
   * Good seeding keeps `deltaUV` small enough that the wide-step case is hard
   * to reach through operator() - a residual bounded by the grid seed gives a
   * step of |e| / |dS/du|, which on a uniformly parameterised closed axis
   * cannot exceed 2r / r - so an end-to-end test does not reliably discriminate
   * it. See test/nurbs_seam_test.cpp.
   *
   * @param from      The uv the trial started from.
   * @param to        The uv it landed on, already brought into the domain.
   * @param requested The unclamped, unwrapped step that was asked for.
   * @return The displacement to measure the predicted decrease over.
   */
  glm::dvec2 trialDisplacement(
      const glm::dvec2& from,
      const glm::dvec2& to,
      const glm::dvec2& requested ) const {

    glm::dvec2 displacement = from - to;

    for ( int axis = 0; axis < 2; ++axis ) {

      if ( ( axis == 0 ? closedU_ : closedV_ ) &&
           ( max_extent[ axis ] - min_extent[ axis ] ) > 0.0 ) {

        displacement[ axis ] = requested[ axis ];
      }
    }

    return displacement;
  }

 private:

  /**
   * Bring a uv back into the parameter domain: by WRAPPING on an axis the
   * surface is closed in, by clamping on one it is not.
   *
   * The descent below has to leave the domain to reach some answers. Where a
   * surface is closed in u, uMin and uMax name the same 3D curve, so stepping
   * off one end is stepping onto the other and the step is not leaving the
   * surface at all - it is crossing the seam. Clamping treats that crossing
   * as a wall: the step is truncated to zero length, the Armijo trial
   * measures no movement, every alpha is rejected, damping grows, and the
   * solve returns its own seed.
   *
   * Measured on `ADVANCED_FACE #50977` of `MANIFOLD_SOLID_BREP #970` ("Spring
   * v1", step/conor/Orbiter_v1.1_Gear_7.5.step): boundary point 88 of 174 is
   * seeded at ( -pi, 43.982297 ) - the domain corner, which is where the
   * preceding seam point sits - and its Gauss-Newton direction points to
   * u < -pi. All 21 Armijo trials evaluate the seed unchanged, and the solve
   * reports ( -pi, 43.982297 ) at a residual of 0.1226 against a target of
   * 1.19e-05. The answer is at u ~ 2.8687, one seam crossing away, and the
   * surface is smooth all the way there: sampling u down from pi walks the
   * residual 0.1226 -> 0.0104 with the minimum in between. Nothing about the
   * geometry stops the descent; only the clamp does.
   *
   * Closure is the same fact seamContinuousSolution uses, read from the
   * control grid rather than measured, so this introduces no tolerance and no
   * constant. On a surface closed in neither parameter every call is the
   * previous glm::clamp, unchanged.
   *
   * @param uv A uv that may lie outside the domain.
   * @return The equivalent uv inside it.
   */
  glm::dvec2 domainRepresentative( glm::dvec2 uv ) const {

    for ( int axis = 0; axis < 2; ++axis ) {

      const double period = max_extent[ axis ] - min_extent[ axis ];

      if ( ( axis == 0 ? closedU_ : closedV_ ) && period > 0.0 ) {

        const double offset = uv[ axis ] - min_extent[ axis ];

        uv[ axis ] =
          min_extent[ axis ] + ( offset - ( period * std::floor( offset / period ) ) );

      } else {

        uv[ axis ] = glm::clamp( uv[ axis ], min_extent[ axis ], max_extent[ axis ] );
      }
    }

    return uv;
  }

  /**
   * Damped Gauss-Newton descent from a caller-supplied seed.
   *
   * Split out of operator() so the same descent can be run twice from
   * different seeds and the results compared on their FINAL residual - see
   * the retry in operator().
   *
   * @param point          The 3D point being inverted.
   * @param bestGuess      Seed uv, taken by value and refined in place.
   * @param bestPoint      surface( bestGuess ), so the caller's existing
   *                       evaluation is not repeated.
   * @param minDistance2   | bestPoint - point |^2, likewise.
   * @param solvedDistance2 Out: the squared residual the descent ended at.
   * @return The uv it ended at.
   */
  glm::dvec2 solveFromSeed(
      const glm::dvec3& point,
      glm::dvec2        bestGuess,
      glm::dvec3        bestPoint,
      double            minDistance2,
      double&           solvedDistance2 ) const {

    size_t iteration = 0;

    double damping = 1e-6;

    // glm::dvec2 alphaUV    = max_extent - min_extent;
    // double     startAlpha = 1.0 / std::max( alphaUV.x, alphaUV.y );

    const double convergenceError2 = convergence_error * convergence_error;

    while ( minDistance2 > convergenceError2 && iteration++ < MAX_ITERATION ) {

      glm::dvec3 deltaP = bestPoint - point;

      if ( minDistance2 <= convergenceError2 ) {
        break;
      }

      auto [dU, dV] = evaluator.tangent( bestGuess.x, bestGuess.y );

      glm::dmat2x3 jacobianT( dU, dV ); // Jacobian
      glm::dmat3x2 jacobian = glm::transpose( jacobianT );   // Transposed Jacobian

      glm::dmat2x2 jtj = jacobian * jacobianT; // J^T * J
      glm::dvec2   jte = jacobian * deltaP;    // J^T * e

      jtj[ 0 ][ 0 ] += damping;
      jtj[ 1 ][ 1 ] += damping;

      glm::dvec2 deltaUV = robust_2x2_solve( jtj, jte );

      double alpha = 1.0;//startAlpha;
      double phi   = 0.5 * minDistance2;

      bool success = false;

      while ( alpha > ALPHA_ERROR ) {

        const glm::dvec2 requested = deltaUV * alpha;

        glm::dvec2 newGuessUV = domainRepresentative( bestGuess - requested );

        // The displacement actually made, which is what Armijo's predicted
        // decrease has to be measured against.
        //
        // The two axis kinds differ, and the difference is not cosmetic. On an
        // OPEN axis the domain edge truncates the step, so the movement is
        // what the clamp left: crediting the full request there is what made
        // the descent stall at a corner. On a CLOSED axis nothing is
        // truncated - S( u + period, v ) IS S( u, v ), so evaluating at the
        // wrapped representative evaluates the requested point exactly, and
        // the movement is the full request.
        //
        // Reducing a closed axis to its shortest representative instead - as
        // the first version of this did - is wrong for a step wider than half
        // a period. `deltaUV` is uncapped, so those steps are reachable, and
        // the reduction can hand back a component pointing the OTHER WAY from
        // the one requested. `glm::dot( ..., jte )` then goes negative, the
        // threshold below rises ABOVE minDistance2, and the trial that gets
        // accepted is one that INCREASES the residual. Using the request
        // removes that case rather than guarding it, because on a closed axis
        // the request is what happened.
        const glm::dvec2 displacement =
          trialDisplacement( bestGuess, newGuessUV, requested );

        glm::dvec3 newPoint =
          evaluator.point( newGuessUV.x, newGuessUV.y );

        glm::dvec3 newDeltaP = newPoint - point;

        double newDistance2 = glm::dot( newDeltaP, newDeltaP );

        // Floored at zero so the threshold can never exceed minDistance2, and
        // an accepted trial therefore always strictly decreases the residual.
        //
        // Before this file wrapped or truncated anything, that was free:
        // `deltaUV` solves an SPD system, so dot( deltaUV, jte ) =
        // jte^T (J^T J + damping I)^-1 jte is positive by construction and the
        // predicted decrease could not come out negative. Truncating a
        // component on an open axis breaks that guarantee - the surviving
        // components are no longer the SPD solution - so the property is
        // asserted here rather than assumed. Where the predicted decrease is
        // positive this changes nothing at all.
        const double predictedDecrease =
          ARMIJO_COEFFICIENT * std::max( 0.0, glm::dot( displacement, jte ) );

        if ( newDistance2 < minDistance2 - predictedDecrease ) {

          bestPoint = newPoint;
          bestGuess = newGuessUV;
          damping  *= 0.1;
          
          minDistance2 = std::min( newDistance2, newDistance2 );
          success = true;
          break;
        }

        alpha *= 0.5;
      }

      // Tested BEFORE the failure path, not after it. `damping *= 10.0`
      // shrinks the next deltaUV monotonically, so a step already under
      // MIN_STEP can only get smaller: no later iteration can move the guess,
      // and Armijo's strict decrease has nothing left to accept. Downstream
      // of that `continue` the test was unreachable on exactly the calls that
      // need it - a stationary point, or a guess whose step clamps back to
      // the same uv - and each one ran out the whole MAX_ITERATION budget
      // failing the same way. Rare, and a termination guard rather than a
      // speed win: 293 of 382,919 solves on Arty_Z7 (0.08%) end this way, and
      // fixing it moves total solve iterations by under 0.1%. The population
      // exists at all only because the target above is now tight enough that
      // some points cannot reach it.
      //
      // Bounded, and measured as such: the steps it skips are under MIN_STEP
      // (1e-9 in uv), so the uv it returns differs by at most ~1e-9 - three
      // decades under the convergence target itself. Corpus effect is 6 rows
      // of 178,273, all of them rows the tolerance change had already moved,
      // and 9 pixels of 1,440,000 on Right_Hand with Arty's silkscreen and
      // nist_ctc_02 byte-identical.
      if ( glm::dot( deltaUV, deltaUV ) < MIN_STEP * MIN_STEP ) {
        break;
      }

      if ( !success ) {
        damping *= 10.0;
        continue;
      }
    }

    solvedDistance2 = minDistance2;

    return bestGuess;
  }

 public:

};


/*inline double InverseMethod(glm::dvec3 pt, const tinynurbs::RationalSurface3d& srf,
                            double pr, double rotations, double minError,
                            double maxError, double &fU, double &fV,
                            double &divisor, double maxDistance) {
  while (maxDistance > maxError && divisor < 10000) {
    for (double r = 1; r < 5; r++) {
      int round = 0;
      while (maxDistance > minError && round < 3) {
      //  printf("maxError: %.3f\n", maxError);
       // printf("minError: %.3f\n", minError);
      //  printf("round: %i\n", round);
        for (double i = 0; i < rotations; i++) {
          double rads = (i / rotations) * (double)CONST_PI * 2;
          double incU = glm::sin(rads) / (r * r * divisor);
          double incV = glm::cos(rads) / (r * r * divisor);
          if (pr > 1) {
            incV *= pr;
          } else {
            incU /= pr;
          }
          bool repeat = true;
          while (repeat) {
            double ffU = fU + incU;
            double ffV = fV + incV;
            glm::highp_dvec3 pt00 = tinynurbs::surfacePoint(srf, ffU, ffV);
            double di = glm::distance(pt00, pt);
            if (di < maxDistance) {
              maxDistance = di;
              fU = ffU;
              fV = ffV;
            } else {
              repeat = false;
            }
          }
        }
        round++;
      }
    }
    divisor *= 3;
    // printf("divisor: %.3f\n", divisor);
    // printf("maxError: %.3f\n", maxError);
    // printf("minError: %.3f\n", minError);
  }
  return maxDistance;
}*/

// inline glm::dvec2 BSplineInverseEvaluation(glm::dvec3 pt,
//                                            const tinynurbs::RationalSurface3d& srf,
//                                            double scaling) {
//   glm::highp_dvec3 ptc = tinynurbs::surfacePoint(srf, 0.0, 0.0);
//   glm::highp_dvec3 pth = tinynurbs::surfacePoint(srf, 1.0, 0.0);
//   glm::highp_dvec3 ptv = tinynurbs::surfacePoint(srf, 0.0, 1.0);

//   double dh = glm::distance(ptc, pth);
//   double dv = glm::distance(ptc, ptv);
//   double pr = (dh + 1) / (dv + 1);

//   double minError = 0.00001;
//   double maxError = 0.001;
//   double rotations = 6;

//   double fU = 0.5;
//   double fV = 0.5;
//   double divisor = 100.0;
//   double maxDistance = 1e+100;

//   //printf("scaling: %.3f\n", scaling);
//   maxDistance =
//     InverseMethod(
//         pt,
//         srf,
//         pr,
//         rotations,
//         minError / scaling,
//         maxError / scaling,
//         fU,
//         fV,
//         divisor,
//         maxDistance );

//   return glm::dvec2(fU, fV);
// }


/**
 * Triangulate a b-spline face that wraps its surface's seam as a structured
 * grid over the whole parametric chart, instead of ear-clipping the chart.
 *
 * WHY EARCUT CANNOT DO THIS FACE. When a face wraps a closed parameter, its
 * uv chart is a rectangle whose two long edges are THE SAME 3D CURVE - the
 * seam, walked down one side and back up the other. earcut is an ear-clipping
 * triangulator: it has only the boundary vertices to work with and cannot
 * introduce interior ones, so it fills that rectangle with chords that join
 * far-apart boundary points. Measured on `ADVANCED_FACE #50977`, the coil
 * spring of solid 970 in `step/conor/Orbiter_v1.1_Gear_7.5.step`, against a
 * chart that was by then provably correct (enclosed area 276.3487 of a
 * 276.3489 chart, coverage 1.0000):
 *
 *   earcut        171 tri, max 3D edge 9.7199, 128 of 171 with an edge > 3.0
 *   tesselated   5473 tri, max 3D edge 9.7199, 5232 near-zero-area
 *   shipped       299 tri  (the weld in Geometry::Reify discards the rest)
 *
 * 9.72 on a spring whose outer diameter is 6.9. Every chord joining u = uMin
 * to u = uMax is degenerate in 3D, because those are the same point, and
 * `tesselate` cannot repair it: it only SPLITS edges, and splitting a
 * degenerate triangle leaves it degenerate. This is the seed-triangulation
 * half of bldrs-ai/conway#608, which lists this face at 103x its deflection
 * target.
 *
 * WHAT REPLACES IT. Every other surface kind in this file already has a
 * full-coverage path - TriangulateSphericalSurface's parametric grid, reached
 * for a degenerate loop, is the direct precedent (conway-geom#187). This is
 * that path for b-splines, and it is a grid over the chart with the trim
 * curve's own points on its borders.
 *
 * HOW COVERAGE IS DECIDED. Not here, and not from the geometry: the caller
 * passes `seamPair`, which the AP214 front end read off the ORIENTED_EDGEs
 * (see hasSeamEdgePair). The file STATES that two of loop `#9507`'s four
 * edges are `EDGE_CURVE #29690` with opposite senses; a face only needs that
 * spelling when it wraps the closed parameter completely. The alternative -
 * testing whether the projected uv boundary encloses the whole chart - was
 * rejected: it infers topology from a measurement and needs a tolerance, and
 * a tolerance tuned until one model passes is how this defect family got
 * here.
 *
 * HOW THE GRID IS BUILT, also without a tolerance. The boundary polyline is
 * read structurally, by EXACT point identity. The two seam legs are the same
 * curve evaluated twice, so their points are bitwise equal - verified on
 * `#50977`: 63 of 63 with a maximum coordinate difference of exactly 0. So:
 *
 *   - a point that occurs more than once in the loop is ON the seam;
 *   - circularly, those flags form exactly four alternating runs,
 *     seam / ring / seam / ring;
 *   - the two seam runs must be equal in length and exact reverses;
 *   - the two ring runs must be equal in length.
 *
 * Anything else returns false and the caller ear-clips exactly as before, so
 * an unfamiliar spelling degrades to today's behaviour rather than to a
 * guess. On `#50977` this yields 2 seam runs of 65 and 2 ring runs of 22,
 * 174 points in all.
 *
 * The grid's rows are the seam run's own v values and its columns the end
 * ring's own u values, so the RESOLUTION IS THE MODEL'S OWN: whatever the
 * adaptive curve tessellator decided those trim curves needed. Nothing is
 * chosen here. The two border rows are the ring points THEMSELVES rather
 * than re-evaluations, which is what keeps this watertight against the two
 * PLANE cap faces that share those circles, and the seam column is likewise
 * the seam curve's own points. Interior points are evaluated on the surface.
 * The grid wraps in u, so the seam is interior to the mesh and only the two
 * end rings are borders - which is correct, since those are exactly the
 * edges that must not move.
 *
 * @param mesh          The face mesh, already holding the projected boundary
 *                      vertices as its first `boundaryCount` entries.
 * @param solve         The face's inverse solve, for surface evaluation and
 *                      for its closure flags.
 * @param boundaryCount Number of boundary vertices in `mesh`.
 * @return True when the grid was built and the caller must not ear-clip.
 */
inline bool tryFullCoverageSeamGrid(
    WingedEdgeMesh< ParameterVertex >& mesh,
    const RationalNurbsInverseMethod&  solve,
    size_t                             boundaryCount,
    double                             deflectionSquared ) {

  // `seamPair` says the TOPOLOGY wraps a seam, but not WHICH parameter it
  // wraps - the predicate is true for either - and everything below reads the
  // repeated legs as the u seam: columns come from the ring and rows from the
  // leg. So this has to establish that the only seam available is the u one.
  //
  // Closed in u and NOT in v does that. A surface closed in neither cannot
  // have produced a retraced pair, and is rejected. A surface closed in v only
  // has a v seam, which this grid would read as a u seam - pulling nearly
  // constant u from the ring and nearly constant v from the leg, and paving
  // the chart with collapsed geometry. A surface closed in BOTH is ambiguous
  // from here and equally unsafe, so it is rejected too rather than guessed
  // at: handling a v seam properly means carrying which parameter wrapped
  // through IfcBound3D and generalising this grid to either axis, which is
  // real scope and not something to infer at this line.
  //
  // Rejecting is cheap and cannot regress anything. The caller ear-clips,
  // which is exactly what every one of these faces got before this function
  // existed. Measured across the STEP corpus reachable here - Orbiter,
  // Right_Hand, nist_ctc_01/02, create-a-tube, supercap, driver board - there
  // is not one surface closed in v at all, with or without u, so this costs
  // nothing today and is a guard against a file that has one.
  if ( !solve.closedU_ || solve.closedV_ ) {
    return false;
  }

  size_t count = boundaryCount;

  // A loop that repeats its first point at the end would make that point look
  // seam-like for the wrong reason. Drop it before counting occurrences.
  if ( count > 1 && mesh.vertices[ count - 1 ].point == mesh.vertices[ 0 ].point ) {
    --count;
  }

  // Four runs of at least two points each is the smallest thing this can be.
  if ( count < 8 ) {
    return false;
  }

  std::unordered_map< glm::dvec3, uint32_t > occurrences;

  for ( size_t at = 0; at < count; ++at ) {
    ++occurrences[ mesh.vertices[ at ].point ];
  }

  std::vector< bool > onSeam( count );

  for ( size_t at = 0; at < count; ++at ) {
    onSeam[ at ] = occurrences[ mesh.vertices[ at ].point ] > 1;
  }

  // Circular runs: rotate to a point where the flag changes, so a run that
  // straddles index 0 - which the seam leg does on this face - is one run
  // rather than two.
  size_t rotation = 0;

  while ( rotation < count &&
          onSeam[ rotation ] == onSeam[ ( rotation + count - 1 ) % count ] ) {
    ++rotation;
  }

  if ( rotation == count ) {
    // Every point has the same flag: no alternation, so no seam structure.
    return false;
  }

  std::vector< std::vector< uint32_t > > runs;

  for ( size_t step = 0; step < count; ++step ) {

    const size_t at = ( rotation + step ) % count;

    if ( step == 0 ||
         onSeam[ at ] != onSeam[ ( rotation + step - 1 ) % count ] ) {

      runs.emplace_back();
    }

    runs.back().push_back( static_cast< uint32_t >( at ) );
  }

  if ( runs.size() != 4 ) {
    return false;
  }

  // Runs alternate by construction, so the seam pair is either {0,2} or
  // {1,3}; `rotation` lands on a run start, whose flag says which.
  const bool seamFirst = onSeam[ runs[ 0 ][ 0 ] ];

  const std::vector< uint32_t >& seamA = seamFirst ? runs[ 0 ] : runs[ 1 ];
  const std::vector< uint32_t >& ringA = seamFirst ? runs[ 1 ] : runs[ 2 ];
  const std::vector< uint32_t >& seamB = seamFirst ? runs[ 2 ] : runs[ 3 ];
  const std::vector< uint32_t >& ringB = seamFirst ? runs[ 3 ] : runs[ 0 ];

  const size_t rows = seamA.size();

  if ( rows < 2 || seamB.size() != rows ||
       ringA.size() != ringB.size() || ringA.empty() ) {
    return false;
  }

  // The two legs must be the same curve walked both ways. Exact, because they
  // ARE the same evaluation; a near-match would mean this is not a seam.
  for ( size_t at = 0; at < rows; ++at ) {

    if ( !( mesh.vertices[ seamA[ at ] ].point ==
            mesh.vertices[ seamB[ rows - 1 - at ] ].point ) ) {

      return false;
    }
  }

  // A ring row is the boundary from one seam leg's end to the next leg's
  // start, INCLUSIVE of those two seam points - they are the row's ends, at
  // u = uMin and u = uMax. Walking the loop in order, ringA follows seamA, so
  // seamA's last point opens it and seamB's first point closes it.
  const auto ringRow =
    [ & ]( const std::vector< uint32_t >& ring,
           uint32_t                       opening,
           uint32_t                       closing ) {

      std::vector< uint32_t > row;

      row.reserve( ring.size() + 2 );
      row.push_back( opening );
      row.insert( row.end(), ring.begin(), ring.end() );
      row.push_back( closing );

      return row;
    };

  // seamA runs one way in v and seamB the other; orient both so index 0 is
  // the end that ringA sits at.
  const std::vector< uint32_t > rowFirst =
    ringRow( ringA, seamA.back(), seamB.front() );
  const std::vector< uint32_t > rowLast =
    ringRow( ringB, seamB.back(), seamA.front() );

  if ( rowFirst.size() != rowLast.size() ) {
    return false;
  }

  // Both ends of a ring row are the same 3D point - the seam - and BOTH are
  // kept as columns rather than wrapping the grid around.
  //
  // Wrapping is the obvious thing and it is wrong, because `tesselate` refines
  // by averaging the uv of an edge's two endpoints. Across a wrap that average
  // is meaningless: on this face the last column sits at u = -2.868745 and the
  // first at u = +3.141593, which are neighbours in 3D but 0.136 - the far side
  // of the wire - when averaged in parameter. Refinement then inserts vertices
  // diametrically opposite where they belong and folds the mesh over itself.
  // Measured: a wrapped grid gave 453.668 of surface area against a true
  // 373.434, with single edges shared by up to 66 triangles.
  //
  // Keeping both seam columns makes u monotone across the whole chart, so every
  // midpoint is meaningful and the fold cannot arise. The two columns carry the
  // same 3D points, so the mesh is seam-duplicated exactly as a cylinder chart
  // normally is, and Geometry::Reify's weld merges them.
  const size_t columns = rowFirst.size();

  // Both ends of each ring row must be the SAME 3D point: that is what makes
  // them the seam, and it is the structural check that this really is a chart
  // that closes rather than a loop that happens to have a repeated edge.
  if ( columns < 4 ||
       !( mesh.vertices[ rowFirst.front() ].point ==
          mesh.vertices[ rowFirst.back() ].point ) ||
       !( mesh.vertices[ rowLast.front() ].point ==
          mesh.vertices[ rowLast.back() ].point ) ) {

    return false;
  }

  // rowFirst is ordered along u from the seam; rowLast is the far ring and
  // runs the other way round the tube, so reverse it to share the ordering.
  std::vector< uint32_t > rowLastAligned( rowLast.rbegin(), rowLast.rend() );

  // seamA carries v from rowFirst's end to rowLast's end. Its last point is
  // rowFirst's opening column, so index it from the back.
  const auto seamAt =
    [ & ]( size_t row ) {
      return seamA[ rows - 1 - row ];
    };

  // Column parameters from the near ring - the model's own adaptive sampling
  // of the curve that BOUNDS this face, which is also what the two PLANE cap
  // faces are built from, so keeping it is what keeps the solid closed.
  std::vector< double > columnU( columns );

  for ( size_t column = 0; column < columns; ++column ) {
    columnU[ column ] = mesh.vertices[ rowFirst[ column ] ].uv.x;
  }

  // -------------------------------------------------------------------------
  // THE CHART CONTRACT.
  //
  // Everything below builds the tensor product columnU x rowV and evaluates
  // interior vertices at ( columnU[ column ], rowV[ row ] ). That is correct
  // on a RECTANGULAR ISOPARAMETRIC PATCH and on nothing else. The precondition
  // is therefore stated here in full and asserted, rather than inferred from
  // proxies that happen to correlate with it.
  //
  // It is written this way because it had to be. Three review rounds each
  // found a DIFFERENT proxy standing in for the same missing assertion - the
  // seam axis, then constant-v rings, then two rings with equal point counts
  // but unrelated u partitions - and each time the run analysis above was
  // satisfied, because it only counts runs, checks they alternate, checks they
  // are pairwise equal in length, and compares seam points. None of that says
  // anything about where the trim actually lies in parameter space.
  // (bldrs-ai/conway#611.)
  //
  // The face must satisfy all of:
  //
  //   1. the support surface is closed in u and NOT in v, so that "the seam"
  //      names one axis unambiguously        [asserted at the top of this
  //                                           function]
  //   2. the loop is four alternating runs - seam, ring, seam, ring - whose
  //      seam legs are exact reverses of one another and whose rings are equal
  //      in length                           [asserted above, run analysis]
  //   3. the columns are a STRICTLY MONOTONE u partition
  //   4. the rows are a STRICTLY MONOTONE v partition
  //   5. the columns span the surface's ENTIRE u period - what makes this full
  //      coverage rather than a strip that merely touches the seam twice
  //   6. both ring runs lie on their constant-v lines, sampled AT THE COLUMNS
  //      the grid will use - which also forces the two rings' u partitions to
  //      correspond, rather than merely to have equal counts
  //   7. both seam legs lie on their constant-u lines, at the first and last
  //      columns
  //
  // What that guarantees to the code below, so a reader does not have to
  // reconstruct it: columnU is an ordered partition of the full
  // [ min_extent.x, max_extent.x ]; rowV is an ordered partition of the face's
  // v range; and every boundary vertex the grid reuses verbatim already sits
  // within `deflectionSquared` of where the tensor product would put it. So
  // substituting S( columnU[ c ], rowV[ r ] ) for the real trim is a
  // substitution the tessellation tolerance already covers.
  //
  // Failing any of it returns false and the caller ear-clips, which is exactly
  // what these faces got before this function existed - so the gate is free to
  // be strict, and is. Narrow and provably correct beats broad and hopeful.
  //
  // TOLERANCE. Only 6 and 7 need one, and it is `deflectionSquared`, the
  // tessellation target already passed in, in the units it already has.
  // Measured on `#50977`, the alternatives are wrong rather than merely
  // stricter:
  //
  //   ring deviation from its isocurve   8.29e-04
  //   seam deviation from its isocurve   1.41e-05
  //   solve->convergence_error           1.19e-05
  //   sqrt( deflectionSquared )          1.54e-02
  //
  // The rings miss constant v by 8.29e-04 while their v values span only
  // 2.2e-05 - a thousandth of that - so the deviation is not the face being
  // non-rectangular. It is the inverse solve's own residual: those ring points
  // converge to ~8e-04, seventy times its target, because they hit
  // MAX_ITERATION. Testing against convergence_error would reject this face
  // for the SOLVER's accuracy rather than for its own shape, which measures
  // the wrong thing - and would take the coil spring with it. Deflection is
  // the honest bound: a boundary lying on its isocurve to within the tolerance
  // the mesh is being built to is one the tessellation already cannot
  // distinguish from the rectangle. It is also what the refinement below aims
  // at, so both halves of this function are judged against a single target.
  //
  // Checks 3, 4 and 5 need no tolerance at all - they are exact comparisons on
  // ordering and on the domain bounds, which the descent assigns exactly
  // because it clamps to them.
  // -------------------------------------------------------------------------
  {
    const double vFirst = mesh.vertices[ seamAt( 0 ) ].uv.y;
    const double vLast  = mesh.vertices[ seamAt( rows - 1 ) ].uv.y;

    // 3. Strictly monotone columns, in either direction. Without this the
    //    chart is not a bijection and cells fold back over one another.
    bool ascending = true, descending = true;

    for ( size_t column = 0; column + 1 < columns; ++column ) {
      ascending  = ascending  && columnU[ column ] < columnU[ column + 1 ];
      descending = descending && columnU[ column ] > columnU[ column + 1 ];
    }

    if ( !ascending && !descending ) {
      return false;
    }

    // 4. Likewise the rows, which the refinement below bisects: a non-monotone
    //    v partition makes those midpoints meaningless. This also makes
    //    vFirst != vLast, so the v range cannot be degenerate.
    bool risingV = true, fallingV = true;

    for ( size_t row = 0; row + 1 < rows; ++row ) {

      const double lower = mesh.vertices[ seamAt( row ) ].uv.y;
      const double upper = mesh.vertices[ seamAt( row + 1 ) ].uv.y;

      risingV  = risingV  && lower < upper;
      fallingV = fallingV && lower > upper;
    }

    if ( !risingV && !fallingV ) {
      return false;
    }

    // 5. The columns must run bound to bound. Monotonicity already rules out
    //    the two ends coinciding, but not a partition that covers only PART of
    //    the period - a strip whose ends happen to be the same 3D point
    //    because the surface passes through it twice. Full coverage is the
    //    entire claim this function makes, so it is asserted rather than
    //    assumed. Exact: the descent clamps onto the bounds, so a seam leg
    //    that genuinely sits on one lands on it exactly.
    const double lowColumn  = ascending ? columnU[ 0 ] : columnU[ columns - 1 ];
    const double highColumn = ascending ? columnU[ columns - 1 ] : columnU[ 0 ];

    if ( lowColumn != solve.min_extent.x || highColumn != solve.max_extent.x ) {
      return false;
    }

    const auto onIsoLine =
      [ & ]( uint32_t vertex, double u, double v ) {

        const glm::dvec3 delta =
          solve.evaluator.point( u, v ) - mesh.vertices[ vertex ].point;

        return glm::dot( delta, delta ) <= deflectionSquared;
      };

    // 6. Both ring rows lie on their constant-v lines, sampled at the columns
    //    the grid will actually use. Testing the FAR ring at columnU - rather
    //    than at its own u values - is what makes this catch two rings that
    //    share a point count but not a partition.
    for ( size_t column = 0; column < columns; ++column ) {

      if ( !onIsoLine( rowFirst[ column ], columnU[ column ], vFirst ) ||
           !onIsoLine( rowLastAligned[ column ], columnU[ column ], vLast ) ) {

        return false;
      }
    }

    // 7. And both seam legs on their constant-u lines, which are the first and
    //    last columns.
    for ( size_t row = 0; row < rows; ++row ) {

      if ( !onIsoLine( seamAt( row ), columnU[ 0 ],
                       mesh.vertices[ seamAt( row ) ].uv.y ) ||
           !onIsoLine( seamB[ row ], columnU[ columns - 1 ],
                       mesh.vertices[ seamB[ row ] ].uv.y ) ) {

        return false;
      }
    }
  }

  // Rows start as the seam curve's own v samples, then are refined IN
  // PARAMETER SPACE until they meet the same deflection target `tesselate`
  // would have used.
  //
  // Refining here rather than handing the grid to `tesselate` is the point.
  // `tesselate` refines topologically, by splitting the edge with the largest
  // deflection and re-projecting its midpoint, and on this mesh that does not
  // converge: measured, it ran to the full 32x amplification cap - 2,944
  // triangles to 94,208 - and took the surface area to 396.89 against a true
  // 373.6, i.e. it folded the mesh over itself, after which Geometry::Reify's
  // weld discarded 90% of what it built. That is the same non-convergence
  // conway#608 reports on this face.
  //
  // Subdividing the v partition cannot fold, because the parameterisation
  // stays monotone: every row is a distinct v, every column a distinct u, and
  // the grid remains a bijection onto the chart no matter how fine it gets.
  // The criterion is the deflection target itself, so the resolution is still
  // derived rather than chosen.
  //
  // Only v is refined. The columns are the trim curve's own points and must
  // not move - see above - and on this face they already meet the target,
  // which is what one would expect of a curve the exporter tessellated to its
  // own deflection tolerance.
  std::vector< double >   rowV;
  std::vector< uint32_t > rowSeam;

  rowV.reserve( rows );
  rowSeam.reserve( rows );

  for ( size_t row = 0; row < rows; ++row ) {

    rowV.push_back( mesh.vertices[ seamAt( row ) ].uv.y );
    rowSeam.push_back( static_cast< uint32_t >( row ) );
  }

  // Bounded by the same budget the topological refinement had, so this cannot
  // cost more than the path it replaces.
  const size_t maximumRows =
    rows * static_cast< size_t >( MAX_TRIANGLE_AMPLIFACTION );

  for ( bool refined = true; refined && rowV.size() < maximumRows; ) {

    refined = false;

    std::vector< double >   nextV;
    std::vector< uint32_t > nextSeam;

    nextV.reserve( rowV.size() * 2 );
    nextSeam.reserve( rowV.size() * 2 );

    for ( size_t row = 0; row + 1 < rowV.size(); ++row ) {

      nextV.push_back( rowV[ row ] );
      nextSeam.push_back( rowSeam[ row ] );

      if ( nextV.size() + ( rowV.size() - row ) >= maximumRows ) {
        continue;
      }

      const double middle = ( rowV[ row ] + rowV[ row + 1 ] ) * 0.5;

      // Worst column decides, so a row is only kept coarse where the whole
      // row is within target.
      //
      // Both the v EDGE and the cell DIAGONAL are tested, because the diagonal
      // is a real edge of the mesh - makeTriangle below splits every cell
      // along it - and it is the longer chord, so testing only the v edge
      // understates the cell. Measured on `#50977` before the diagonal was
      // included: the v edges converged to 0.89x the target while the
      // diagonals sat at 1.16x, i.e. the grid missed the tolerance it claims
      // to aim at on the very edges it emits.
      //
      // Refining v answers both, and terminates: as the v spacing shrinks,
      // the diagonal's deflection falls toward the u-direction deflection at
      // that column, which on this face is 0.29x the target. That floor is
      // also the honest limit of this loop - it can only subdivide v. The
      // columns are the trim curve's own points and cannot be subdivided
      // here, because the two PLANE cap faces are built from those same
      // points and a column inserted between them would put a vertex on the
      // end rings that the caps do not have. So a surface whose u direction
      // ALONE misses the target needs the caps refined with it, which is
      // cross-face coordination this architecture does not have. Nothing in
      // the corpus is in that state, and the u term is measured here so the
      // criterion is at least aware of it rather than silently ignoring it.
      double worst = 0.0;

      for ( size_t column = 0; column < columns; ++column ) {

        const glm::dvec3 lower = solve.evaluator.point( columnU[ column ], rowV[ row ] );
        const glm::dvec3 upper =
          solve.evaluator.point( columnU[ column ], rowV[ row + 1 ] );

        const glm::dvec3 edgeDelta =
          solve.evaluator.point( columnU[ column ], middle ) - ( ( lower + upper ) * 0.5 );

        worst = std::max( worst, glm::dot( edgeDelta, edgeDelta ) );

        if ( column + 1 >= columns ) {
          continue;
        }

        // The diagonal makeTriangle uses runs ( row, column ) to
        // ( row + 1, column + 1 ).
        const glm::dvec3 far =
          solve.evaluator.point( columnU[ column + 1 ], rowV[ row + 1 ] );

        const glm::dvec3 diagonalDelta =
          solve.evaluator.point( ( columnU[ column ] + columnU[ column + 1 ] ) * 0.5,
                                 middle ) - ( ( lower + far ) * 0.5 );

        worst = std::max( worst, glm::dot( diagonalDelta, diagonalDelta ) );
      }

      if ( worst > deflectionSquared ) {

        nextV.push_back( middle );
        nextSeam.push_back( EMPTY_INDEX );
        refined = true;
      }
    }

    nextV.push_back( rowV.back() );
    nextSeam.push_back( rowSeam.back() );

    rowV.swap( nextV );
    rowSeam.swap( nextSeam );
  }

  const size_t gridRows = rowV.size();


  std::vector< uint32_t > grid( gridRows * columns );

  for ( size_t row = 0; row < gridRows; ++row ) {

    const uint32_t seamRow = rowSeam[ row ];
    const bool     onSeam  = seamRow != EMPTY_INDEX;

    for ( size_t column = 0; column < columns; ++column ) {

      uint32_t vertex;

      if ( onSeam && seamRow == 0 ) {

        vertex = rowFirst[ column ];

      } else if ( onSeam && seamRow == rows - 1 ) {

        vertex = rowLastAligned[ column ];

      } else if ( onSeam && column == 0 ) {

        // Both seam columns take the seam curve's own points rather than
        // re-evaluations. seamB is seamA reversed, so the point matching
        // seam row `seamRow` on the far side is seamB[ seamRow ].
        vertex = seamAt( seamRow );

      } else if ( onSeam && column == columns - 1 ) {

        vertex = seamB[ seamRow ];

      } else {

        const glm::dvec2 uv( columnU[ column ], rowV[ row ] );

        vertex = mesh.makeVertex(
          ParameterVertex{ solve.evaluator.point( uv.x, uv.y ), uv } );
      }

      grid[ ( row * columns ) + column ] = vertex;
    }
  }

  for ( size_t row = 0; row + 1 < gridRows; ++row ) {
    for ( size_t column = 0; column + 1 < columns; ++column ) {

      const size_t nextColumn = column + 1;

      const uint32_t a = grid[ ( row * columns ) + column ];
      const uint32_t b = grid[ ( row * columns ) + nextColumn ];
      const uint32_t c = grid[ ( ( row + 1 ) * columns ) + nextColumn ];
      const uint32_t d = grid[ ( ( row + 1 ) * columns ) + column ];

      // Guard against a cell that collapses onto a shared vertex. The grid's
      // construction should preclude it - every cell's four corners are
      // either four distinct boundary vertices or freshly made ones - so this
      // is defensive, and on `#50977` it never fires.
      //
      // Note what index identity CANNOT see: a cell whose corners are
      // distinct VERTICES sharing a PARAMETER. Two grid columns at the same u
      // evaluate to the same point on every interior row, Geometry::Reify's
      // weld merges them, and the cells between them vanish - leaving
      // duplicated half-edges along a whole generatrix that no test here
      // would catch. That is not hypothetical: it is what this face did until
      // domainRepresentative() taught the descent to cross a closed surface's
      // seam, and it cost 249 duplicated half-edges of 34,047 plus two
      // zero-area triangles. A parametric duplicate is a defect in the
      // inverse solve, and it has to be fixed there - synthesising a
      // replacement parameter here would invent the geometry this path exists
      // to read. See bldrs-ai/conway#611.
      if ( a != b && b != c && c != a ) {
        mesh.makeTriangle( a, b, c );
      }

      if ( a != c && c != d && d != a ) {
        mesh.makeTriangle( a, c, d );
      }
    }
  }

  return true;
}


/**
 * Sign changes of dv around a closed uv boundary, with `dv == 0` treated as
 * NO INFORMATION rather than as a reversal.
 *
 * The obvious test - `previous * next < 0` on the two adjacent differences -
 * silently misses a turning point that sits on a PLATEAU, and on these
 * boundaries the loop's start vertex is exactly such a plateau. Measured on
 * `step/conor/Orbiter_v1.1_Gear_7.5.step`: the 5089-, 5057-, 1017- and
 * 1025-point b-spline faces each report one plateau, at the wrap, so the
 * strict test finds ONE reversal where there are two and rejects a boundary
 * that is a textbook ribbon. That cost a full implementation cycle, so the
 * plateau handling is the reason this helper exists rather than being written
 * inline (bldrs-ai/conway#608).
 *
 * @param mesh  The mesh holding the boundary as its first `count` vertices.
 * @param count Number of boundary vertices.
 * @return Indices at which the v direction reverses, ascending.
 */
inline std::vector< size_t > uvMonotoneBreaks(
    const WingedEdgeMesh< ParameterVertex >& mesh,
    size_t                                   count ) {

  std::vector< size_t > breaks;

  if ( count < 3 ) {
    return breaks;
  }

  const auto deltaAt =
    [ & ]( size_t at ) {
      return mesh.vertices[ ( at + 1 ) % count ].uv.y - mesh.vertices[ at ].uv.y;
    };

  // Start the scan ON a non-zero delta and seed the running sign from it.
  //
  // Seeding with "no sign yet" instead loses a reversal: the first non-zero
  // delta would only establish the sign, so a boundary whose index 0 falls in
  // the middle of a run reports one break rather than two - which is the whole
  // failure this helper exists to avoid, arrived at from the other direction.
  size_t start = count;

  for ( size_t at = 0; at < count; ++at ) {

    if ( deltaAt( at ) != 0.0 ) {
      start = at;
      break;
    }
  }

  if ( start == count ) {
    // Constant v all the way round: no runs, so no ribbon.
    return breaks;
  }

  int lastSign = deltaAt( start ) > 0.0 ? 1 : -1;

  // Exactly one full turn, so every delta including the wrap is seen once.
  for ( size_t step = 1; step <= count; ++step ) {

    const size_t at = ( start + step ) % count;

    const double delta = deltaAt( at );

    if ( delta == 0.0 ) {
      continue;
    }

    const int sign = delta > 0.0 ? 1 : -1;

    if ( sign != lastSign ) {
      breaks.push_back( at );
    }

    lastSign = sign;
  }

  std::sort( breaks.begin(), breaks.end() );
  breaks.erase( std::unique( breaks.begin(), breaks.end() ), breaks.end() );

  return breaks;
}

/**
 * Triangulate a trimmed b-spline face that is a RIBBON - a boundary of two
 * v-monotone legs meeting at two turning points - with a monotone sweep over
 * its own boundary vertices, instead of ear-clipping its parametric chart.
 *
 * WHY EARCUT CANNOT DO THIS FACE. The thread flanks of
 * `step/conor/Orbiter_v1.1_Gear_7.5.step` are trimmed strips whose uv chart
 * has an aspect ratio of 7489:1. Ear clipping introduces no interior vertices
 * and has no quality criterion, so on a strip that narrow it fills the chart
 * with slivers running its length - and a sliver long in v is a CHORD ACROSS
 * THE COIL. Measured on `ADVANCED_FACE #51059`, whose chart is provably sound
 * (inverse solve converged on 5,087 of 5,089 boundary points, earcut coverage
 * 1.0000, zero inverted seed triangles):
 *
 *   seed triangles spanning > 98% of the v range     1
 *   seed slivers (> 1% of v, < 50% of u)         2,456 of 5,086
 *   longest seed edge                          19.4998  (face diagonal 19.94)
 *   shipped p90 deviation                       2.3519  (118x its target)
 *
 * The seed is a VALID tiling of the chart and a useless triangulation of the
 * surface at the same time, which is why per-axis uv normalisation before
 * earcut measured as a byte-identical no-op (conway#601). Refinement cannot
 * repair it either: `tesselate` only splits edges, so it inherits the chords;
 * with the livelock guard removed and the budget raised sixteenfold the p90
 * still sits at 11.8x target while the mesh area inflates 45x, i.e. it folds.
 *
 * WHAT THIS DOES INSTEAD. A v-monotone polygon has a textbook linear-time
 * triangulation (de Berg et al.), and its decisive property here is that it
 * emits triangles between EXISTING boundary vertices only. So:
 *
 *   - it cannot create a T-vertex, and the trim polyline reaches the mesh
 *     exactly as the neighbouring faces build it - watertight topologically,
 *     not merely geometrically;
 *   - it is correct for ANY monotone polygon, concave included, with no
 *     convexity test and no tolerance;
 *   - a horizontal edge is just a tie in the sweep order, so it needs no
 *     restriction on flat stretches of a leg.
 *
 * Two earlier designs are recorded here because each was refuted by
 * measurement rather than by argument, and both look reasonable on paper.
 *
 * A TENSOR GRID on shared rows, mirroring tryFullCoverageSeamGrid, needs the
 * two legs to CORRESPOND. They do not: leg lengths 2544/2545, 2528/2529,
 * 512/513, 498/519, and across all 21 lofted Orbiter faces the two legs share
 * exactly TWO v values - their turning points - out of a 5,088-value union.
 *
 * A LOFT joining the legs at one v fixes the geometry but not the topology.
 * Wherever the legs' samples differ, one rung endpoint must be interpolated
 * onto a leg; the point is collinear so nothing leaks, but Reify's weld does
 * not subdivide the neighbour's edge, leaving two half-edges facing one.
 * Measured on the shipped GLB with half-edges matched on exact world position,
 * Orbiter's unpaired half-edges went 1,448 -> 70,698 - a 48x increase against
 * the metric conway-geom#188 used to certify the spring fix by driving it to
 * zero (codex on conway-geom#190).
 *
 * WHERE QUALITY COMES FROM. Not the seed. `tesselate` runs over this mesh as
 * it does over the ear-clipped one, and it splits INTERIOR edges only - it
 * returns early on edge.border() - so every vertex it adds is strictly
 * interior and the boundary stays verbatim however finely the face resolves.
 * The amplification budget is therefore the same one every other b-spline face
 * gets, and MAX_TRIANGLE_AMPLIFACTION stays the real limit.
 *
 * SCOPE. Genuine two-leg ribbons only. A boundary with more than two
 * v-monotone runs returns false and ear-clips exactly as it does today,
 * including the four faces whose uv chart is self-overlapping because their
 * inverse solve did not converge (conway#642). Narrow and provably correct
 * beats broad and hopeful; this is the posture tryFullCoverageSeamGrid takes.
 *
 * @param mesh              The face mesh, holding the projected boundary
 *                          vertices as its first `boundaryCount` entries.
 * @param solve             The face's inverse solve, unused by the sweep and
 *                          kept so the signature matches the seam grid's.
 * @param boundaryCount     Number of boundary vertices in `mesh`.
 * @return True when the sweep triangulated the face.
 */
inline bool tryRibbonLoft(
    WingedEdgeMesh< ParameterVertex >& mesh,
    const RationalNurbsInverseMethod&  solve,
    size_t                             boundaryCount ) {

  (void)solve;

  size_t count = boundaryCount;

  // A loop that repeats its first point at the end would read as a plateau for
  // the wrong reason. Drop it before looking at the v profile.
  if ( count > 1 &&
       mesh.vertices[ count - 1 ].point == mesh.vertices[ 0 ].point ) {
    --count;
  }

  if ( count < 8 ) {
    return false;
  }

  const std::vector< size_t > breaks = uvMonotoneBreaks( mesh, count );

  // Exactly two turning points is what "ribbon" means. Everything else - a
  // notched or oscillating boundary - is out of scope and ear-clips.
  if ( breaks.size() != 2 ) {
    return false;
  }

  const size_t first  = breaks[ 0 ];
  const size_t second = breaks[ 1 ];

  // Leg A is the arc first -> second; leg B is the rest of the loop. Both
  // INCLUDE the two turning points, which they share.
  // Carried as INDICES into the mesh's own boundary vertices, not as copies:
  // the emitted mesh has to share those vertices, which is what makes the trim
  // polylines bit-identical rather than merely equal.
  std::vector< uint32_t > legA;
  std::vector< uint32_t > legB;

  for ( size_t at = first; at <= second; ++at ) {
    legA.push_back( static_cast< uint32_t >( at ) );
  }

  for ( size_t at = second; at <= second + ( count - ( second - first ) ); ++at ) {
    legB.push_back( static_cast< uint32_t >( at % count ) );
  }

  if ( legA.size() < 2 || legB.size() < 2 ) {
    return false;
  }

  // Structural coverage: every boundary vertex belongs to exactly one leg, and
  // the two turning points to both. Nothing may be dropped - a mesh built over
  // part of a face still sits close to the surface, so the deflection metric
  // cannot see a truncated boundary and this count is the only thing that can.
  if ( legA.size() + legB.size() != count + 2 ) {
    return false;
  }

  // Orient both legs ascending in v so the staircase below is monotone.
  const auto vertexOf =
    [ & ]( uint32_t index ) -> const ParameterVertex& {
      return mesh.vertices[ index ];
    };

  if ( vertexOf( legA.front() ).uv.y > vertexOf( legA.back() ).uv.y ) {
    std::reverse( legA.begin(), legA.end() );
  }

  if ( vertexOf( legB.front() ).uv.y > vertexOf( legB.back() ).uv.y ) {
    std::reverse( legB.begin(), legB.end() );
  }

  // The legs must meet at both ends, exactly. That is what makes this one
  // closed ribbon rather than two unrelated curves, and it is an identity
  // between vertices the loop already shares, so it needs no tolerance.
  if ( !( vertexOf( legA.front() ).point == vertexOf( legB.front() ).point ) ||
       !( vertexOf( legA.back() ).point == vertexOf( legB.back() ).point ) ) {
    return false;
  }

  // The polygon's own winding, so emitted triangles can be oriented to match
  // what earcut would have produced on this same boundary.
  //
  // Relative to a boundary point for the same reason the validity gate's sums
  // are - see there. This one only needs a SIGN, but a sign is exactly what
  // cancellation destroys first, and getting it wrong inverts every triangle
  // the face emits.
  const glm::dvec2 windingReference = mesh.vertices[ 0 ].uv;

  double shoelace = 0.0;

  for ( size_t at = 0; at < count; ++at ) {

    const glm::dvec2 here = mesh.vertices[ at ].uv - windingReference;
    const glm::dvec2 next =
      mesh.vertices[ ( at + 1 ) % count ].uv - windingReference;

    shoelace += ( here.x * next.y ) - ( next.x * here.y );
  }

  const double winding = shoelace >= 0.0 ? 1.0 : -1.0;

  // MONOTONE TWO-CHAIN SWEEP.
  //
  // The textbook triangulation of a v-monotone polygon (de Berg et al.), and
  // the reason it is used here rather than anything cleverer is that it emits
  // triangles between EXISTING boundary vertices only. It never introduces a
  // vertex on the boundary, so it cannot create a T-vertex, and the face's
  // trim polyline reaches the mesh exactly as the neighbouring faces build it.
  //
  // That is not a refinement of the previous design, it is the correction of
  // it. The loft this replaced joined the two legs at one v, which meant that
  // wherever the legs' v samples differed - which is everywhere, they are
  // independently tessellated trim curves, and measured they share only their
  // two turning points out of 5,088 - one rung endpoint had to be interpolated
  // ONTO a leg. That point is collinear, so the solid stayed geometrically
  // gap-free, but Reify's weld does not subdivide the neighbour's edge, so the
  // result was two half-edges facing one. Measured on the shipped GLB with
  // half-edges matched on exact world position: Orbiter's unpaired half-edges
  // went 1,448 -> 70,698, a 48x increase, against a metric this epic used to
  // certify the spring fix by driving it to zero (codex on conway-geom#190).
  //
  // The sweep also subsumes the concavity fix that preceded it: it is correct
  // for ANY v-monotone polygon, concave or not, and needs no plateau
  // restriction, because a horizontal edge is just a tie in the sweep order.
  //
  // Quality does not come from the seed here. It comes from `tesselate`
  // afterwards, which splits INTERIOR edges only - it returns early on
  // edge.border() - so every vertex refinement adds is strictly interior and
  // the boundary stays verbatim no matter how finely the face is resolved.
  enum class Chain : uint8_t { A, B };

  struct SweepVertex {

    uint32_t index;
    Chain    chain;
  };

  std::vector< SweepVertex > sweep;

  sweep.reserve( count );

  // Descending in v: the shared top, then both chains' interiors merged, then
  // the shared bottom. Both legs are ascending, so they are walked backwards.
  sweep.push_back( { legA.back(), Chain::A } );

  {
    size_t onA = legA.size() >= 2 ? legA.size() - 2 : 0;
    size_t onB = legB.size() >= 2 ? legB.size() - 2 : 0;

    const bool hasA = legA.size() >= 2;
    const bool hasB = legB.size() >= 2;

    size_t leftA = hasA ? legA.size() - 2 : 0;
    size_t leftB = hasB ? legB.size() - 2 : 0;

    while ( leftA > 0 || leftB > 0 ) {

      // Descending in v, and AT EQUAL v ordered by u - a geometric ordering,
      // not a chain preference.
      //
      // `>=` alone silently means "chain A first", and that is wrong wherever
      // the two chains meet at one v. A horizontal edge on chain A puts two A
      // vertices at the same v; if a B vertex shares that v, chain preference
      // emits both A vertices before it, the sweep then reaches three
      // consecutive vertices at one v, the triangle it forms has zero area and
      // is dropped, and the face finishes one triangle short of `count - 2`.
      // The validity gate then rejects a perfectly good v-monotone ribbon back
      // to the ear-clipper - restoring exactly the long chords this path
      // exists to remove (codex on conway-geom#190).
      //
      // Which u comes first depends on which side is which, and that is what
      // `winding` carries: the two chains are traversed in opposite senses
      // around the loop, so the boundary's own orientation fixes the order
      // without anyone having to decide which leg is "left".
      bool takeA = leftB == 0;

      if ( !takeA && leftA > 0 ) {

        const glm::dvec2& onChainA = vertexOf( legA[ onA ] ).uv;
        const glm::dvec2& onChainB = vertexOf( legB[ onB ] ).uv;

        takeA =
          onChainA.y != onChainB.y ?
            ( onChainA.y > onChainB.y ) :
            ( ( onChainA.x - onChainB.x ) * winding > 0.0 );
      }

      if ( takeA ) {

        sweep.push_back( { legA[ onA ], Chain::A } );

        --leftA;

        if ( onA > 0 ) { --onA; }

      } else {

        sweep.push_back( { legB[ onB ], Chain::B } );

        --leftB;

        if ( onB > 0 ) { --onB; }
      }
    }
  }

  sweep.push_back( { legA.front(), Chain::B } );

  if ( sweep.size() < 3 ) {
    return false;
  }

  // Emit with the boundary's own winding, so the face's normals match what the
  // ear-clipped path produced and `appendMeshToGeometry` flips them the same
  // way.
  std::vector< std::array< uint32_t, 3 > > swept;

  swept.reserve( count );

  const auto emit =
    [ & ]( uint32_t a, uint32_t b, uint32_t c ) {

      if ( a == b || b == c || c == a ) {
        return;
      }

      const glm::dvec2& p = mesh.vertices[ a ].uv;
      const glm::dvec2& q = mesh.vertices[ b ].uv;
      const glm::dvec2& r = mesh.vertices[ c ].uv;

      const double area =
        ( ( q.x - p.x ) * ( r.y - p.y ) ) - ( ( r.x - p.x ) * ( q.y - p.y ) );

      if ( area == 0.0 ) {
        return;
      }

      if ( ( area * winding ) > 0.0 ) {
        swept.push_back( { a, b, c } );
      } else {
        swept.push_back( { a, c, b } );
      }
    };

  // Is the diagonal from `current` to `next` interior? On one chain that asks
  // for a left turn and on the other a right turn; `winding` folds in which
  // way round the boundary runs.
  const auto diagonalInside =
    [ & ]( const SweepVertex& current,
           const SweepVertex& last,
           const SweepVertex& next ) {

      const glm::dvec2& p = mesh.vertices[ next.index ].uv;
      const glm::dvec2& q = mesh.vertices[ last.index ].uv;
      const glm::dvec2& r = mesh.vertices[ current.index ].uv;

      const double turn =
        ( ( q.x - p.x ) * ( r.y - p.y ) ) - ( ( r.x - p.x ) * ( q.y - p.y ) );

      const double want = current.chain == Chain::A ? -1.0 : 1.0;

      return ( turn * winding * want ) > 0.0;
    };

  std::vector< SweepVertex > stack;

  stack.reserve( sweep.size() );

  stack.push_back( sweep[ 0 ] );
  stack.push_back( sweep[ 1 ] );

  for ( size_t at = 2; ( at + 1 ) < sweep.size(); ++at ) {

    const SweepVertex current = sweep[ at ];

    if ( current.chain != stack.back().chain ) {

      while ( stack.size() > 1 ) {

        const SweepVertex top = stack.back();

        stack.pop_back();

        emit( current.index, top.index, stack.back().index );
      }

      stack.pop_back();
      stack.push_back( sweep[ at - 1 ] );
      stack.push_back( current );

    } else {

      SweepVertex last = stack.back();

      stack.pop_back();

      while ( !stack.empty() && diagonalInside( current, last, stack.back() ) ) {

        const SweepVertex next = stack.back();

        stack.pop_back();

        emit( current.index, last.index, next.index );

        last = next;
      }

      stack.push_back( last );
      stack.push_back( current );
    }
  }

  {
    const SweepVertex bottom = sweep.back();

    while ( stack.size() > 1 ) {

      const SweepVertex top = stack.back();

      stack.pop_back();

      emit( bottom.index, top.index, stack.back().index );
    }
  }

  // VALIDITY GATE. The sweep is only committed if it actually triangulated the
  // polygon, and that is checked exactly rather than trusted.
  //
  // Any triangulation of a simple n-gon has exactly n - 2 triangles, so a
  // different count means the sweep terminated early or fanned - it dropped
  // part of the face. And the triangles' summed |uv area| must equal the
  // boundary's shoelace area, or they overlap. Both are properties of the
  // OUTPUT, so neither depends on trusting the sweep's own reasoning.
  //
  // This is not belt-and-braces: measured across 284 lofted faces on the
  // corpus, 3 faces violate the count (one 45-point boundary emitted 21
  // triangles where 43 are required) and 9 overlap by more than 1e-6, up to
  // 4.5%. Every one of them is a SMALL face - 15 to 67 boundary points - and
  // every one of the large thread faces this path exists for passes both
  // checks exactly. Rather than ship a construction that is right on the
  // faces that motivated it and quietly wrong on a handful of others, a face
  // that fails either check is handed back to the ear-clipper, which is
  // exactly what it got before this function existed.
  //
  // The tolerance on the area check is relative and loose (1e-6) because the
  // sum runs over thousands of triangles: at 1e-9 the honest faces fail on
  // double-precision accumulation alone. The count check needs no tolerance.
  {
    if ( swept.size() != count - 2 ) {
      return false;
    }

    // Both sums are taken RELATIVE TO A BOUNDARY POINT, never on the raw
    // parameters.
    //
    // A shoelace term is a difference of products, and a b-spline's uv comes
    // straight from the file's knot domain, which carries whatever offset the
    // exporter used. On a domain like [1e4, 1e4 + 0.03] every product is ~1e8
    // while their difference is ~1e-3, so the leading digits cancel and the
    // area arrives as noise - it can come out zero, or wrong by more than the
    // 1e-6 the check below allows. Either way a perfectly good sweep is
    // rejected to the ear-clipper and the long chords come back, which is a
    // failure that would have been invisible: the gate is designed to reject,
    // so it rejecting looks like it working (codex on conway-geom#190).
    //
    // Subtracting a vertex of the polygon itself makes every coordinate small
    // relative to the shape, and area is translation invariant so the value is
    // unchanged. The triangle sum is already written as differences and is
    // sound as it stands, but it is shifted by the same reference so the two
    // quantities the ratio compares are computed the same way rather than
    // merely being equal in exact arithmetic.
    const glm::dvec2 reference = mesh.vertices[ 0 ].uv;

    double shoelaceArea = 0.0;

    for ( size_t at = 0; at < count; ++at ) {

      const glm::dvec2 here = mesh.vertices[ at ].uv - reference;
      const glm::dvec2 next =
        mesh.vertices[ ( at + 1 ) % count ].uv - reference;

      shoelaceArea += ( here.x * next.y ) - ( next.x * here.y );
    }

    shoelaceArea = fabs( shoelaceArea ) * 0.5;

    double covered = 0.0;

    for ( const std::array< uint32_t, 3 >& triangle : swept ) {

      const glm::dvec2 a = mesh.vertices[ triangle[ 0 ] ].uv - reference;
      const glm::dvec2 b = mesh.vertices[ triangle[ 1 ] ].uv - reference;
      const glm::dvec2 c = mesh.vertices[ triangle[ 2 ] ].uv - reference;

      covered +=
        fabs( ( ( b.x - a.x ) * ( c.y - a.y ) ) -
              ( ( c.x - a.x ) * ( b.y - a.y ) ) ) * 0.5;
    }

    if ( shoelaceArea <= 0.0 ||
         fabs( ( covered / shoelaceArea ) - 1.0 ) > 1e-6 ) {
      return false;
    }
  }

  for ( const std::array< uint32_t, 3 >& triangle : swept ) {
    mesh.makeTriangle( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );
  }

  return true;
}


// TODO: review and simplify
/**
 * Re-solve the head of a CLOSED trim loop, seeded by its own last point.
 *
 * RationalNurbsInverseMethod::resetContinuity() scopes the continuity seed to
 * one bound, which is correct - a cross-loop seed can converge silently onto a
 * different sheet. But it leaves the FIRST point of every bound with no seed at
 * all, cold-starting from the grid. On a surface that passes near itself (a
 * multi-turn thread flank, the conway#594 configuration) that grid can land in
 * the wrong basin, and the bad uv then seeds its successors until the descent
 * claws its way back. Measured across three models, 10 of 848 b-spline faces
 * cold-start into a residual more than 100x their own loop median and above
 * 0.1% of the face extent; the worst lands 26% of the face away from the
 * surface and one returns the domain corner (0,0) outright.
 *
 * A closed loop's first point is adjacent to its last, so after the caller's
 * first pass the solver's continuity seed is exactly the right one. Re-solving
 * the head from there replaces a grid guess with a genuine neighbour. Stop at
 * the first point whose answer does not change: from there the caller's pass
 * was already running on a seed of its own.
 *
 * Measured, the case is even plainer than "adjacent": across 859 trim bounds in
 * the three b-spline-carrying regression models, the closing gap is EXACTLY
 * zero in every one - the last point of a trim polyline is bit-identical to the
 * first. So the head and the tail are not merely neighbours, they are the SAME
 * QUERY, and a cold start that disagrees with the continuity-seeded solve of
 * that identical point is self-evidently the wrong one of the two. This pass
 * does not choose between two defensible answers; it removes an inconsistency
 * where one point had two.
 *
 * THE CLOSURE GATE IS NOT OPTIONAL. Everything above depends on the last point
 * genuinely neighbouring the first. An OPEN bound - a trim that failed to
 * close, an extraction that dropped its tail - has a tail somewhere unrelated,
 * and seeding the head from it is precisely the cross-sheet failure
 * resetContinuity exists to prevent, arriving by a different door. So the gate
 * measures the closing gap against the loop's own median segment and declines
 * anything that does not look like a closure, leaving those bounds on the cold
 * start they have today.
 *
 * @return the number of head points rewritten; 0 if the gate declined.
 */
inline size_t reSolveClosedTrimHead(
    RationalNurbsInverseMethod      &solve,
    const std::vector< glm::dvec3 > &scaledPoints,
    std::vector< glm::dvec2 >       &solved,
    bool                             closedByConstruction ) {

  const size_t count = scaledPoints.size();

  if ( count < 3 || solved.size() != count ) {
    return 0;
  }

  // EXPLICIT CLOSURE ONLY: the last point must repeat the first.
  //
  // This gate was first written as a proximity heuristic - closing gap against
  // the loop's own median segment - and that approach is not merely
  // mis-tuned, it is unusable, because the gap is fixed by geometry while the
  // median scales with sampling density. The same open arc spanning one coil
  // turn, whose tail sits one pitch from its head on a different sheet,
  // measures:
  //
  //     samples   120    60    40    30    20
  //     ratio    4.05  2.02  1.35  1.01  0.68
  //
  // against an adjacency-closed loop's ~1.0. At 30 samples the open arc is
  // INDISTINGUISHABLE from a closed loop, and at 20 it looks more closed than
  // one. No threshold separates the two classes; a coarser mesh walks an open
  // bound through any gate you pick. Found on bldrs-ai/conway#655.
  //
  // Explicit closure has no such failure mode. Measured, it also costs nothing:
  // all 859 trim bounds in the three b-spline-carrying regression models close
  // by EXACT duplication, the last point bit-identical to the first. So the
  // head and tail are not merely neighbours but the SAME QUERY, and a cold
  // start that disagrees with the continuity-seeded solve of that identical
  // point is self-evidently the wrong one of the two - this pass removes an
  // inconsistency rather than choosing between two defensible answers.
  //
  // DELIBERATELY DECLINED: loops that close by ADJACENCY - last point one
  // sample step from the first, no duplicate. An earlier revision accepted
  // those, on a proximity test that has since been shown unusable. No corpus
  // bound exhibits adjacency closure; the case was synthetic, built for a
  // fixture. If a real producer of it appears, it needs a closure signal from
  // the front end (where the ORIENTED_EDGEs are still visible and closure is a
  // topological fact) rather than a geometric guess here.
  //
  // The tolerance is float noise, not a proximity budget: head and tail are
  // multiplied by the same `scaling` and so stay bit-identical through it, but
  // a few ULP of the coordinate magnitude costs nothing to allow and keeps the
  // test from depending on that staying true.
  // Two ways a bound can be known to close, and BOTH are facts rather than
  // measurements:
  //
  //   - the producer says so. A point-list loop is a closed polygon that does
  //     not repeat its head (IfcCurve::closedByConstruction, set in GetLoop
  //     where the entity type makes it certain). Without this the whole
  //     poly-loop producer class - every IFCPOLYLOOP face bound - would be
  //     declined and silently keep the cold-start error. Found by review on
  //     bldrs-ai/conway-geom#195; it is the negative answer to "is there a
  //     front-end closure signal", namely that the producer knows and simply
  //     was not encoding it.
  //
  //   - the points say so, by repeating the head as the tail. This is how every
  //     one of the 859 edge-loop trim bounds in the regression corpus arrives.
  //
  // The tolerance below is float noise, not a proximity budget: head and tail
  // take the same `scaling` multiplication and stay bit-identical through it,
  // but a few ULP costs nothing to allow and keeps this from depending on that.
  if ( !closedByConstruction ) {

    double maxMagnitude = 0.0;

    for ( const glm::dvec3& point : scaledPoints ) {
      maxMagnitude = std::max( maxMagnitude,
        std::max( std::abs( point.x ),
                  std::max( std::abs( point.y ), std::abs( point.z ) ) ) );
    }

    const double closingGap =
      glm::distance( scaledPoints[ count - 1 ], scaledPoints[ 0 ] );

    if ( closingGap > 8.0 * DBL_EPSILON * std::max( maxMagnitude, 1.0 ) ) {
      return 0;
    }
  }

  size_t rewritten = 0;

  for ( size_t i = 0; i < count; ++i ) {

    glm::dvec2 again;
    {
      conway::AllocTagScope inverseTag( conway::AllocSite::NurbsInverse );
      again = solve( scaledPoints[ i ] );
    }

    if ( again == solved[ i ] ) {
      break;
    }

    solved[ i ] = again;
    ++rewritten;
  }

  return rewritten;
}

inline void TriangulateBspline(Geometry &geometry,
                               const std::vector<IfcBound3D> &bounds,
                               IfcSurface &surface, double scaling,
                               double representationExtent) {

//  printf( "Triangulating BSpline Surface\n" );

  tinynurbs::RationalSurface3d srf;
  
  srf.degree_u = surface.BSplineSurface.UDegree;
  srf.degree_v = surface.BSplineSurface.VDegree;
  size_t num_u = surface.BSplineSurface.ControlPoints.size();
  size_t num_v = surface.BSplineSurface.ControlPoints[0].size();

  std::vector<glm::dvec3> controlPoints;

  controlPoints.reserve( num_u * num_v );

  for ( const std::vector<glm::dvec3>& row : surface.BSplineSurface.ControlPoints ) {
    for (const glm::dvec3& point : row) {
      controlPoints.push_back({point.x, point.y, point.z});
    }
  }

  srf.control_points = tinynurbs::array2(num_u, num_v, controlPoints);

  std::vector<double> weights;
  weights.reserve( num_u * num_v );
  // Read WeightPoints, not Weights: WeightPoints is the field the embind
  // surface parameters (and the IFC GetSurface path) actually populate.
  // Weights has no writer anywhere, so reading it always fell through to
  // the all-1.0 default below and rational surfaces (e.g. STEP cylinders
  // written as weighted Bezier arcs) evaluated as plain polynomials,
  // bulging the profile - part of conway#350.
  for (const std::vector<double>& row : surface.BSplineSurface.WeightPoints) {
    for (double weight : row) {
      weights.push_back(weight);
    }
  }

  if (weights.size() != num_u * num_v) {
    for (size_t i = 0; i < num_u * num_v; i++) {
      weights.push_back(1.0);
    }
  }
  
  srf.weights = tinynurbs::array2(num_u, num_v, weights);

  for (size_t i = 0; i < surface.BSplineSurface.UMultiplicity.size(); i++) {
    for (size_t r = 0; r < surface.BSplineSurface.UMultiplicity[i]; r++) {
      srf.knots_u.push_back(surface.BSplineSurface.UKnots[i]);
    }
  }

  for (size_t i = 0; i < surface.BSplineSurface.VMultiplicity.size(); i++) {
    for (size_t r = 0; r < surface.BSplineSurface.VMultiplicity[i]; r++) {
      srf.knots_v.push_back(surface.BSplineSurface.VKnots[i]);
    }
  }

  // If the NURBS surface is valid we continue
  

//  printf( "Evaluating inverse parameter space\n" );

  if (tinynurbs::surfaceIsValid(srf)) {

    // Constructed only for valid surfaces: this builds the homogeneous
    // control grid and samples the inverse-evaluation seed grid, all of
    // which would be wasted on the invalid-surface path below.
    RationalNurbsInverseMethod bSplineInverseEvaluation( srf );

    // Cap the solve's convergence target on this face before inverting
    // anything - see boundConvergenceToTrim. The bounds are already 3D and in
    // hand, so this is one pass over the same points the projection loop
    // below walks, and it must run first: operator() reads convergence_error
    // on every call.
    {
      glm::dvec3 trimMin( std::numeric_limits< double >::max() );
      glm::dvec3 trimMax( std::numeric_limits< double >::lowest() );

      for ( const IfcBound3D& bound : bounds ) {
        for ( const glm::dvec3& point : bound.curve.points ) {

          // Same `scaling` the projection loop applies; the target has to be
          // in the units the points are inverted in.
          const glm::dvec3 scaled = point * scaling;

          trimMin = glm::min( trimMin, scaled );
          trimMax = glm::max( trimMax, scaled );
        }
      }

      // trimMin > trimMax when no bound carried a point, which leaves the
      // distance non-finite; boundConvergenceToTrim floors that case.
      bSplineInverseEvaluation.boundConvergenceToTrim(
        glm::all( glm::lessThanEqual( trimMin, trimMax ) )
          ? glm::distance( trimMin, trimMax )
          : 0.0 );
    }

    // Find projected boundary using NURBS inverse evaluation

    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> uvBoundaryValues;

    std::vector<ParameterVertex> vertices;

    // AFTP: arena-back this per-face mesh (esp. edge_map); rewound at scope exit.
    conway::ScratchArenaScope arenaScope;
    WingedEdgeMesh< ParameterVertex > mesh{ conway::ThreadScratchResource() };

    for ( size_t i = 0; i < bounds.size(); ++i ) {

      std::vector<Point> points;

      // Each bound is its own ordered polyline; the continuity seed is only
      // meaningful within one. See RationalNurbsInverseMethod::resetContinuity.
      bSplineInverseEvaluation.resetContinuity();

      const uint32_t firstVertex =
        static_cast< uint32_t >( mesh.vertices.size() );

      const size_t boundPointCount = bounds[ i ].curve.points.size();

      std::vector< glm::dvec3 > scaledPoints;
      std::vector< glm::dvec2 > solvedUv;

      scaledPoints.reserve( boundPointCount );
      solvedUv.reserve( boundPointCount );

      for (size_t j = 0; j < boundPointCount; j++) {
        glm::dvec3 pt = bounds[i].curve.points[j];

        //hack 
        pt.x *= scaling;
        pt.y *= scaling;
        pt.z *= scaling;

        glm::dvec2 pInv;
        {
          conway::AllocTagScope inverseTag( conway::AllocSite::NurbsInverse );
          pInv = bSplineInverseEvaluation( pt );
        }

        scaledPoints.push_back( pt );
        solvedUv.push_back( pInv );

        points.push_back({pInv.x, pInv.y});
        mesh.makeVertex( { pt, pInv } );
      }

      // The head of a closed loop has no continuity seed - see
      // reSolveClosedTrimHead, which also carries the closure gate that keeps
      // an OPEN bound on its cold start rather than seeding it cross-sheet.
      const size_t rewritten =
        reSolveClosedTrimHead( bSplineInverseEvaluation, scaledPoints, solvedUv,
                               bounds[ i ].curve.closedByConstruction );

      for ( size_t j = 0; j < rewritten; ++j ) {
        mesh.vertices[ firstVertex + j ].uv = solvedUv[ j ];
        points[ j ] = { solvedUv[ j ].x, solvedUv[ j ].y };
      }

      uvBoundaryValues.push_back(points);
    }

  //  printf( "Earcutting parameter space %zu\n", mesh.vertices.size() );

    // Triangulate projected boundary
    // Subdivide resulting triangles to increase definition
    // r indicates the level of subdivision, currently 3 you can increase it to
    // 5

    // A face that wraps the surface's seam covers the whole chart, and
    // ear-clipping that chart produces chords across the solid rather than a
    // tessellation of it - see tryFullCoverageSeamGrid, which builds the grid
    // instead. `seamPair` is a topological fact from the front end, so this
    // costs one bool test on every other face and changes none of them.
    //
    // Single-bound only: an inner trim loop is a hole, and a grid over the
    // whole chart would pave straight over it.
    //
    // The target is computed here, from the boundary vertices alone, because
    // the grid builder refines against it and `tesselate` below is handed the
    // same expression - the two paths must aim at the same tolerance.
    const double seamGridDeflection2 =
      relativeDeflectionSquared( mesh, representationExtent * scaling );

    const bool builtSeamGrid =
      bounds.size() == 1 &&
      bounds[ 0 ].seamPair &&
      tryFullCoverageSeamGrid(
        mesh, bSplineInverseEvaluation, mesh.vertices.size(), seamGridDeflection2 );

    // A trimmed RIBBON - two v-monotone boundary legs - is lofted rather than
    // ear-clipped, for the reasons in tryRibbonLoft. Single-bound only, like
    // the seam grid: an inner trim loop is a hole, and a loft spanning the
    // whole chart would pave over it. Anything that is not a ribbon returns
    // false here and takes exactly the path it takes today.
    const bool builtRibbonLoft =
      !builtSeamGrid &&
      bounds.size() == 1 &&
      tryRibbonLoft( mesh, bSplineInverseEvaluation, mesh.vertices.size() );

    const bool builtStructured = builtSeamGrid || builtRibbonLoft;

    std::vector<uint32_t> indices;

    if ( !builtStructured ) {
  {
    conway::AllocTagScope earcutTag( conway::AllocSite::Earcut );
    indices = mapbox::earcut<uint32_t>(uvBoundaryValues);
  }

    for ( size_t i = 0; i < indices.size(); i += 3 ) {

      mesh.makeTriangle( 
        indices[ i  + 0 ], 
        indices[ i  + 1 ], 
        indices[ i  + 2 ] );
    }
    }
    
  //  printf( "Tesselating BSpline Surface\n" );

    // Skipped for the seam grid, which arrives already refined to this same
    // target in parameter space. Running the topological refiner over it does
    // not converge - it spends the whole 32x budget and folds the mesh - see
    // tryFullCoverageSeamGrid.
    // Runs over the ribbon sweep too: it splits interior edges only, so it is
    // where the ribbon's quality comes from, and it cannot touch the boundary.
    if ( !builtSeamGrid )
    tesselate(
      mesh,
      [&bSplineInverseEvaluation]( [[maybe_unused]]const glm::dvec3&, const glm::dvec2& from ) {
        return bSplineInverseEvaluation.evaluator.point( from.x, from.y );
      },
      mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
      // `representationExtent * scaling`, not the raw extent. This is the ONE
      // triangulator that does not tessellate in the units its bound points
      // arrive in: the projection loop above multiplies every point by
      // `scaling` before seeding the mesh, so the mesh's own bounding box -
      // the other half of the comparison relativeDeflectionSquared makes -
      // is in post-`scaling` units too. Handing it a raw extent would put a
      // millimetre-to-metre load's floor 1000x too coarse and a
      // metre-to-millimetre one's 1000x too fine, i.e. off entirely. Same
      // reason boundConvergenceToTrim above is fed a scaled trim diagonal.
      //
      // Both front ends pass scaling = 1 today, so this multiply is a no-op
      // in every shipping path; it is here so the invariant on
      // ParamsAddFaceToGeometry::representationExtent ("in the units the
      // bound points arrive in") stays true if that ever changes.
      relativeDeflectionSquared( mesh, representationExtent * scaling ) );

    appendMeshToGeometry( mesh, geometry, !surface.sameSense );

  //  printf( "Tesselated BSpline Surface with %zu triangles\n", mesh.triangles.size() );


  } else {
    Logger::logError( "Surface was not valid!\n");
  }
}
}  // namespace conway::geometry
