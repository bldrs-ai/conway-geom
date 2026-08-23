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
    const char *label ) {

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
                                  IfcSurface &surface) {

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
      relativeDeflectionSquared( mesh ) );

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
                                        IfcSurface &surface) {
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

  if ( boundaryPointCount > 0 && boundaryPointCount < MINIMUM_TRIM_POINTS ) {

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
    relativeDeflectionSquared( mesh ) );

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
    IfcSurface &surface) {

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
    relativeDeflectionSquared( mesh ) );

  appendMeshToGeometry( mesh, geometry );
}


// TODO: review and simplify
inline void TriangulateConicalSurface(
  Geometry &geometry,
  const std::vector<IfcBound3D> &bounds,
  IfcSurface &surface) {
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
          relativeDeflectionSquared( mesh ) );

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
    relativeDeflectionSquared( mesh ) );

  appendMeshToGeometry( mesh, geometry, sameSense );
}

// TODO: review and simplify
inline void TriangulateCylindricalSurface(Geometry &geometry,
                                          const std::vector<IfcBound3D> &bounds,
                                          IfcSurface &surface) {
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
          relativeDeflectionSquared( unwrapMesh ) );

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
    relativeDeflectionSquared( mesh ) );

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
                                 IfcSurface &surface) {

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
      relativeDeflectionSquared( mesh ) );

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

    // Scaled off the grid, which spans the WHOLE surface, rather than off the
    // trimmed face: the solve is constructed before any bound is projected.
    // That over-estimates a face that is an island on a larger surface, and
    // over-estimating is NOT automatically safe - too loose is its own failure
    // mode. Two thresholds bound it, both in terms of surfaceDiagonal /
    // trimDiagonal:
    //
    //   >= 1e3  the residual reaches the tessellation deflection target
    //           (RELATIVE_INVERSE_ERROR / RELATIVE_DEFLECTION_FACTOR), i.e.
    //           the bug this whole change fixes, returning at a bigger ratio;
    //   >= 1e6  the tolerance reaches the trim's own diameter, so every
    //           boundary point starts within one tolerance of the same grid
    //           sample, operator() iterates zero times, returns one uv for the
    //           entire loop, and the uv earcut degenerates.
    //
    // Measured over every b-spline face in the locally materialised corpus -
    // 10,422 faces across Arty_Z7, Right_Hand, nist_ctc_02 and
    // FM_ARC_DigitalHub - the maximum surfaceDiagonal / trimDiagonal is
    // 2.273. Zero faces reach either threshold, zero collapse their uv
    // boundary, zero produce an empty seed mesh: these exporters fit the
    // surface to the face.
    //
    // The nearest approach to the first threshold is 0.121 of it, on a 0.49mm
    // face in metre-unit Right_Hand - and it comes from MIN_INVERSE_ERROR
    // below, not from any surface/trim mismatch (that face's ratio is 0.9996).
    // If a file ever does need this capped, cap it on the trim rather than
    // loosening the floor: bounds[].curve.points are 3D and already in hand
    // before the projection loop, so their bounding box costs one pass over
    // points the next loop walks anyway, and no extra inversion.
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
  }

  glm::dvec2 operator()( const glm::dvec3& point ) const {

    glm::dvec2 bestGuess;
    glm::dvec3 bestPoint;
    double     minDistance2 = DBL_MAX;
    
    // Take initial guess from the grid.
    for ( size_t i = 0; i < INVERSE_GRID_SIDE; ++i ) {
      for ( size_t j = 0; j < INVERSE_GRID_SIDE; ++j ) {

        glm::dvec3 deltaP = grid[ i ][ j ] - point;

        double distance2 = glm::dot( deltaP, deltaP );

        if ( distance2 < minDistance2 ) {

          bestGuess =
            min_extent + 
            ( max_extent - min_extent ) * 
            glm::dvec2( 
              static_cast< double >( i ) * INVERSE_GRID_FACTOR,
              static_cast< double >( j ) * INVERSE_GRID_FACTOR );

          bestPoint    = grid[ i ][ j ];
          minDistance2 = distance2;
        }
      }
    }

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

        glm::dvec2 newGuessUV = bestGuess - deltaUV * alpha;

        newGuessUV = glm::clamp( newGuessUV, min_extent, max_extent );

        glm::dvec3 newPoint =
          evaluator.point( newGuessUV.x, newGuessUV.y );

        glm::dvec3 newDeltaP = newPoint - point;

        double newDistance2 = glm::dot( newDeltaP, newDeltaP );

        if ( newDistance2 <
             minDistance2 - ARMIJO_COEFFICIENT * alpha * glm::dot( deltaUV, jte ) ) {

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

    return bestGuess;
  }

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

// TODO: review and simplify
inline void TriangulateBspline(Geometry &geometry,
                               const std::vector<IfcBound3D> &bounds,
                               IfcSurface &surface, double scaling) {

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

    // Find projected boundary using NURBS inverse evaluation

    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> uvBoundaryValues;

    std::vector<ParameterVertex> vertices;

    // AFTP: arena-back this per-face mesh (esp. edge_map); rewound at scope exit.
    conway::ScratchArenaScope arenaScope;
    WingedEdgeMesh< ParameterVertex > mesh{ conway::ThreadScratchResource() };

    for ( size_t i = 0; i < bounds.size(); ++i ) {

      std::vector<Point> points;

      for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {
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

        points.push_back({pInv.x, pInv.y});
        mesh.makeVertex( { pt, pInv } );
      }

      uvBoundaryValues.push_back(points);
    }

  //  printf( "Earcutting parameter space %zu\n", mesh.vertices.size() );

    // Triangulate projected boundary
    // Subdivide resulting triangles to increase definition
    // r indicates the level of subdivision, currently 3 you can increase it to
    // 5

    std::vector<uint32_t> indices;
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
    
  //  printf( "Tesselating BSpline Surface\n" );

    tesselate(
      mesh,
      [&bSplineInverseEvaluation]( [[maybe_unused]]const glm::dvec3&, const glm::dvec2& from ) {
        return bSplineInverseEvaluation.evaluator.point( from.x, from.y );
      },
      mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
      relativeDeflectionSquared( mesh ) );

    appendMeshToGeometry( mesh, geometry, !surface.sameSense );

  //  printf( "Tesselated BSpline Surface with %zu triangles\n", mesh.triangles.size() );


  } else {
    Logger::logError( "Surface was not valid!\n");
  }
}
}  // namespace conway::geometry
