/*
 * Tests for triangulatePeriodicUChart in mesh_utils.h.
 *
 * THIS FILE REPLACES A FILE OF REFUSAL TESTS, and most of what it stopped
 * testing stopped existing. Its predecessor pinned the gates of
 * `tryPeriodicUStrip`, which cut the periodic annulus open into a simple
 * polygon because `mapbox::earcut` needs one: the rim classification, the
 * per-ring nearest-image unwrap, the cut-crossing tests, the per-hole
 * containment test, and the late-refusal leak each had a test. The CDT path
 * has none of those constructions, so it has none of those gates, and a test
 * that pins a gate that cannot exist is not coverage - it is a fossil.
 *
 * Every one of those cases is still here, as an input. What changed is the
 * assertion, and the file is organised by what happened to it:
 *
 *   PORTED UNCHANGED - the two gates that are facts about the SURFACE rather
 *   than about a construction, and which the CDT path keeps verbatim: the
 *   declared `u_closed`, and the per-knot-span closure check with enough
 *   samples per span to settle each span's polynomial.
 *
 *   NOW BUILD, where the old code refused a spelling it could not read. A rim
 *   closed by ADJACENCY rather than by repeating its head was refused because
 *   `netDelta` measured P - P/n; here every ring's closing segment is a
 *   constraint edge like any other, so the two spellings are the same input
 *   and both build. A rim that WANDERS back across its own shortest cut
 *   needed a 60-candidate cut search; here there is no cut to search for. A
 *   ring whose two ends are one 3D point the solve answered twice, 0.003 of a
 *   period apart, needed the winding to be rounded rather than measured; here
 *   nothing reads a winding per ring at all.
 *
 *   REFUSED BY THE ONE POST-HOC CHECK, where the old code had a gate per
 *   failure mode. A hole reaching out through a rim, a hole whose unlisted
 *   closing segment crossed the cut, and a ring that crosses itself are all
 *   the same defect - rings that intersect - and all three are caught by the
 *   same reading: `CDT::IntersectingConstraintEdges::TryResolve` RESOLVES
 *   intersections by adding vertices, so a vertex set that came back larger
 *   than it went in is the signal, and it refuses the whole face. The two
 *   spellings of the crossing hole (head repeated, head not repeated) are
 *   pinned as behaving IDENTICALLY, which is the defect class disappearing
 *   rather than one more gate covering it.
 *
 *   CHARACTERISED, where the removal of a gate legitimately changed what is
 *   accepted: a bound that lies clear of every other bound used to be refused
 *   by hole containment and is now triangulated as its own island, because
 *   `eraseOuterTrianglesAndHoles` decides interior by parity and a closed loop
 *   in the chart has an interior. That is the same semantics
 *   `triangulateUnwrappedLoops` has shipped with on the cylinder and cone
 *   paths; it is recorded here rather than gated.
 *
 *   NEW, for what the CDT path builds that the cut did not. The lift is the
 *   whole reason a cut was ever needed - `tesselate`'s ParameterVertex
 *   overload midpoints uv, so an edge spanning u ~ uMax and u ~ uMin lands on
 *   the far side of the surface and folds the mesh - so the operative
 *   assertion is that NO emitted triangle spans half a period in the lifted
 *   chart. The emergent seam's duplicates are asserted to be BITWISE copies of
 *   their originals' positions, which is what welds them in Geometry::Reify.
 *   Chord subdivision is asserted to fire on a coarse rim and to keep the
 *   triangulated area right, because the polar layout's chords dip inward.
 *   And every refusal is asserted to leave the mesh exactly as it was handed
 *   over, which is the one property of the old file that carries over
 *   untouched.
 *
 * The surface is synthetic and minimal on purpose: a degree-1 tube whose last
 * control row IS its first, which is closed in u by evaluation - the only
 * property of the surface this reads, besides its knot domain.
 * `makeSplitSeamTube` is the same tube with more v knot spans and a seam that
 * only closes on the spans the old five probes landed on.
 *
 * Standalone by design: it includes mesh_utils.h directly and links nothing
 * but the Logger stubs below, matching outer_bound_order_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <array>
#include <limits>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// mesh_utils.h's error paths call these; the rest of the header is header-only.
void Logger::logError( const char*, ... ) {}
void Logger::logWarning( const char*, ... ) {}

namespace {

int failures = 0;

void check( bool condition, const std::string& what ) {

  if ( condition ) {
    printf( "  ok    %s\n", what.c_str() );
    return;
  }

  printf( "  FAIL  %s\n", what.c_str() );
  ++failures;
}

using conway::geometry::ParameterVertex;
using conway::geometry::WingedEdgeMesh;

constexpr double TUBE_RADIUS = 10.0;
constexpr double TUBE_HEIGHT = 20.0;

/** Segments round the tube. The last control row repeats the first. */
constexpr size_t TUBE_SEGMENTS = 12;

constexpr double TWO_PI = 6.283185307179586;

/**
 * A degree-1 tube, closed in u because its last control row IS its first, with
 * the u knot domain [0, 1] so a period reads as 1 everywhere below.
 */
tinynurbs::RationalSurface3d makeClosedTube() {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  const size_t rows = TUBE_SEGMENTS + 1;

  std::vector< glm::dvec3 > control;

  control.reserve( rows * 2 );

  for ( size_t row = 0; row < rows; ++row ) {

    // The wrap is written as a modulus rather than as a second cos/sin call so
    // the closing row is BIT-identical to the first; the closure gate compares
    // evaluated points, and a 1-ulp difference there is a different test.
    const double angle =
      TWO_PI * ( static_cast< double >( row % TUBE_SEGMENTS ) /
                 static_cast< double >( TUBE_SEGMENTS ) );

    control.push_back( { TUBE_RADIUS * std::cos( angle ),
                         TUBE_RADIUS * std::sin( angle ), 0.0 } );
    control.push_back( { TUBE_RADIUS * std::cos( angle ),
                         TUBE_RADIUS * std::sin( angle ), TUBE_HEIGHT } );
  }

  surface.control_points = tinynurbs::array2( rows, 2, control );
  surface.weights = tinynurbs::array2( rows, 2, std::vector< double >( rows * 2, 1.0 ) );

  surface.knots_u.push_back( 0.0 );

  for ( size_t at = 0; at <= TUBE_SEGMENTS; ++at ) {
    surface.knots_u.push_back( static_cast< double >( at ) /
                               static_cast< double >( TUBE_SEGMENTS ) );
  }

  surface.knots_u.push_back( 1.0 );

  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  return surface;
}

/** v knot spans in `makeSplitSeamTube`: eight, at 0, 1/8 ... 1. */
constexpr size_t SPLIT_SEAM_SPANS = 8;

/** How far the split seam opens, in model units. Far above closureTolerance. */
constexpr double SPLIT_SEAM_GAP = 1.0;

/**
 * The same tube, with eight v knot spans instead of one and a seam that closes
 * only at EVEN v control columns - which is to say, only at v = 0, 0.25, 0.5,
 * 0.75 and 1.
 *
 * That is exactly the set of parameters the old five-probe closure gate read,
 * so it certified this surface closed. It is not: at every odd column the last
 * u control row stands SPLIT_SEAM_GAP further out than the first, and the two
 * u-end isocurves - degree-1 polylines through those columns - are a full
 * SPLIT_SEAM_GAP apart there. A strip built on it has cut edges that land on
 * different 3D curves and weld to nothing.
 */
tinynurbs::RationalSurface3d makeSplitSeamTube() {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  const size_t rows    = TUBE_SEGMENTS + 1;
  const size_t columns = SPLIT_SEAM_SPANS + 1;

  std::vector< glm::dvec3 > control;

  control.reserve( rows * columns );

  for ( size_t row = 0; row < rows; ++row ) {

    const double angle =
      TWO_PI * ( static_cast< double >( row % TUBE_SEGMENTS ) /
                 static_cast< double >( TUBE_SEGMENTS ) );

    for ( size_t column = 0; column < columns; ++column ) {

      const double v =
        static_cast< double >( column ) /
        static_cast< double >( SPLIT_SEAM_SPANS );

      // Only the CLOSING row moves, and only where the five probes did not
      // look. Row 0 is untouched, so the seam is shut at every even column.
      const double radius =
        ( row == TUBE_SEGMENTS && ( column % 2 ) == 1 ) ?
          ( TUBE_RADIUS + SPLIT_SEAM_GAP ) : TUBE_RADIUS;

      control.push_back( { radius * std::cos( angle ),
                           radius * std::sin( angle ),
                           TUBE_HEIGHT * v } );
    }
  }

  surface.control_points = tinynurbs::array2( rows, columns, control );
  surface.weights =
    tinynurbs::array2( rows, columns,
                       std::vector< double >( rows * columns, 1.0 ) );

  surface.knots_u.push_back( 0.0 );

  for ( size_t at = 0; at <= TUBE_SEGMENTS; ++at ) {
    surface.knots_u.push_back( static_cast< double >( at ) /
                               static_cast< double >( TUBE_SEGMENTS ) );
  }

  surface.knots_u.push_back( 1.0 );

  surface.knots_v.push_back( 0.0 );

  for ( size_t at = 0; at <= SPLIT_SEAM_SPANS; ++at ) {
    surface.knots_v.push_back( static_cast< double >( at ) /
                               static_cast< double >( SPLIT_SEAM_SPANS ) );
  }

  surface.knots_v.push_back( 1.0 );

  return surface;
}

/** The tube point a (u, v) names, with u taken as a fraction of the period. */
glm::dvec3 tubePoint( double u, double v ) {

  const double angle = TWO_PI * u;

  return { TUBE_RADIUS * std::cos( angle ),
           TUBE_RADIUS * std::sin( angle ),
           TUBE_HEIGHT * v };
}

/**
 * One trim ring, as the u samples it visits at a fixed v. `repeatHead` appends
 * the ring's first sample again, which is how every edge-loop bound in the
 * corpus arrives and what makes the head-to-tail winding a full traversal.
 */
struct RingSpec {

  std::vector< double > us;
  double                v;
  bool                  repeatHead;

  // Non-zero makes the ring straddle two v values instead of sitting on one,
  // so a four-sample ring encloses area. Rims leave it at zero: their shape in
  // v is not what any gate under test reads.
  double                halfHeight = 0.0;

  // When non-empty this IS the ring, verbatim, and everything above is
  // ignored. A ring whose shape in BOTH parameters is the point of the test
  // cannot be spelled as `us` at one v.
  std::vector< std::array< double, 2 > > explicitPoints;

  // When non-empty, one more entry is appended carrying THIS uv and the HEAD'S
  // 3D POINT. That is a closed loop whose two ends are one point the inverse
  // solve answered twice and disagreed about - which is what
  // `ADVANCED_FACE #19215`'s seam-straddling hole is: head ( 0.933400,
  // 0.124924 ) against tail ( 0.936026, 0.125770 ) for one point, 0.0026 of a
  // period apart. Spelled as an override because nothing derived from a uv can
  // express it: the two ends have to disagree in uv while agreeing exactly in
  // 3D.
  std::vector< std::array< double, 2 > > tailUv;
};

/** A rim going once round the tube, evenly, in the given direction. */
RingSpec evenRim( size_t samples, double v, bool forward, bool repeatHead ) {

  RingSpec spec;

  spec.v          = v;
  spec.repeatHead = repeatHead;

  for ( size_t at = 0; at < samples; ++at ) {

    const double fraction =
      static_cast< double >( at ) / static_cast< double >( samples );

    // Wrapped into [0, 1), which is where the inverse solve returns u.
    spec.us.push_back( forward ? fraction : std::fmod( 1.0 - fraction, 1.0 ) );
  }

  return spec;
}

/** An ordinary non-wrapping ring: a small box in the chart. */
RingSpec smallRing( double centreU, double centreV, double halfSize ) {

  RingSpec spec;

  spec.v          = centreV;
  spec.repeatHead = true;
  spec.halfHeight = halfSize;
  spec.us         = { centreU - halfSize, centreU + halfSize,
                      centreU + halfSize, centreU - halfSize };

  return spec;
}

/** Built rings, and the mesh whose vertices they index in ring order. */
struct Built {

  WingedEdgeMesh< ParameterVertex >                     mesh;
  std::vector< std::vector< std::array< double, 2 > > > rings;
};

/**
 * Lays the rings out the way TriangulateBspline does: every boundary point
 * becomes a mesh vertex, in ring order, starting at vertex 0.
 *
 * `smallRing` carries a v half-size as well, so its four samples make a box
 * rather than a degenerate segment; rims are at one v, which is all the gates
 * under test read.
 */
Built build( const std::vector< RingSpec >& specs ) {

  Built built;

  for ( const RingSpec& spec : specs ) {

    std::vector< std::array< double, 2 > > ring;

    if ( !spec.explicitPoints.empty() ) {

      for ( const std::array< double, 2 >& point : spec.explicitPoints ) {

        ring.push_back( point );
        built.mesh.makeVertex( { tubePoint( point[ 0 ], point[ 1 ] ),
                                 glm::dvec2( point[ 0 ], point[ 1 ] ) } );
      }

      if ( !spec.tailUv.empty() ) {

        ring.push_back( spec.tailUv.front() );
        built.mesh.makeVertex(
          { tubePoint( spec.explicitPoints.front()[ 0 ],
                       spec.explicitPoints.front()[ 1 ] ),
            glm::dvec2( spec.tailUv.front()[ 0 ], spec.tailUv.front()[ 1 ] ) } );
      }

      built.rings.push_back( std::move( ring ) );
      continue;
    }

    const size_t count = spec.us.size();

    for ( size_t at = 0; at <= count; ++at ) {

      if ( at == count ) {

        if ( !spec.repeatHead ) {
          break;
        }

        // The SAME sample again: same uv, same point, bit-identical.
        ring.push_back( ring.front() );
        built.mesh.makeVertex(
          { tubePoint( ring.front()[ 0 ], ring.front()[ 1 ] ),
            glm::dvec2( ring.front()[ 0 ], ring.front()[ 1 ] ) } );
        break;
      }

      const double v =
        spec.halfHeight == 0.0 ?
          spec.v :
          ( ( at == 0 || at == 1 ) ? spec.v - spec.halfHeight
                                   : spec.v + spec.halfHeight );

      ring.push_back( { spec.us[ at ], v } );
      built.mesh.makeVertex( { tubePoint( spec.us[ at ], v ),
                               glm::dvec2( spec.us[ at ], v ) } );
    }

    built.rings.push_back( std::move( ring ) );
  }

  return built;
}

/**
 * Runs the function under test and reports everything the assertions below
 * need: whether it built, what it did to the mesh, and the emitted triangles.
 */
struct Outcome {

  bool                                     built;
  size_t                                   vertexGrowth;
  std::vector< std::array< uint32_t, 3 > > triangles;
  double                                   period;
  double                                   uMin;
};

Outcome run( Built& state, const tinynurbs::RationalSurface3d& surface,
             bool declaredClosedU = true ) {

  std::vector< std::array< uint32_t, 3 > > triangles;

  double period = 0.0;
  double uMin   = 0.0;

  const size_t before = state.mesh.vertices.size();

  const bool built =
    conway::geometry::triangulatePeriodicUChart(
      state.mesh, surface, declaredClosedU, state.rings,
      [ &surface ]( double u, double v ) {

        return tinynurbs::surfacePoint( surface, u, v );
      },
      triangles, period, uMin );

  return { built, state.mesh.vertices.size() - before, std::move( triangles ),
           period, uMin };
}

/**
 * THE OPERATIVE PROPERTY. `tesselate`'s ParameterVertex overload refines by
 * `newUV = ( v0.uv + v1.uv ) * 0.5`, so a triangle whose corners sit on
 * different periodic images of the chart has edge midpoints on the far side
 * of the surface, and the mesh folds when it is refined. The whole reason a
 * cut was ever needed here is to make sure no such triangle exists; the lift
 * replaces it, and this is the assertion that says the lift worked.
 *
 * @return The widest u span of any emitted triangle, as a fraction of the
 *         period. Anything at or above 0.5 is a folded triangle.
 */
double worstLiftedSpan( const Built& state, const Outcome& outcome ) {

  double worst = 0.0;

  for ( const std::array< uint32_t, 3 >& triangle : outcome.triangles ) {

    double low  = std::numeric_limits< double >::max();
    double high = std::numeric_limits< double >::lowest();

    for ( size_t corner = 0; corner < 3; ++corner ) {

      const double u = state.mesh.vertices[ triangle[ corner ] ].uv.x;

      low  = std::min( low, u );
      high = std::max( high, u );
    }

    worst = std::max( worst, high - low );
  }

  return worst / outcome.period;
}

/** The emitted mesh's area in 3D, which is what the face actually ships. */
double emittedArea( const Built& state, const Outcome& outcome ) {

  double total = 0.0;

  for ( const std::array< uint32_t, 3 >& triangle : outcome.triangles ) {

    const glm::dvec3& a = state.mesh.vertices[ triangle[ 0 ] ].point;
    const glm::dvec3& b = state.mesh.vertices[ triangle[ 1 ] ].point;
    const glm::dvec3& c = state.mesh.vertices[ triangle[ 2 ] ].point;

    total += glm::length( glm::cross( b - a, c - a ) ) * 0.5;
  }

  return total;
}

/**
 * Every vertex the run added beyond `boundaryCount` must carry a position
 * that is BITWISE some existing vertex's position, or it is a subdivision
 * point rather than a seam duplicate. Counted separately so a test can say
 * which kind it expected.
 */
struct Growth {

  size_t duplicates;
  size_t introduced;
};

Growth classifyGrowth( const Built& state, size_t boundaryCount ) {

  Growth growth { 0, 0 };

  for ( size_t at = boundaryCount; at < state.mesh.vertices.size(); ++at ) {

    bool isCopy = false;

    for ( size_t other = 0; other < boundaryCount && !isCopy; ++other ) {

      isCopy = state.mesh.vertices[ at ].point ==
               state.mesh.vertices[ other ].point;
    }

    if ( isCopy ) {
      ++growth.duplicates;
    } else {
      ++growth.introduced;
    }
  }

  return growth;
}

}  // namespace


int main() {

  const tinynurbs::RationalSurface3d surface = makeClosedTube();

  // The band between two rims at v = 0.25 and v = 0.75, which most of the
  // cases below triangulate: half the tube's height, all the way round. The
  // boundary is a polygon, so the reference is the inscribed n-gon's perimeter
  // times the height and not the circle's - an EXACT quantity, which is what
  // lets the area assertions run at 1e-9 rather than at a few percent.
  const auto bandArea =
    []( size_t sides ) {

      const double perimeter =
        static_cast< double >( sides ) * 2.0 * TUBE_RADIUS *
        std::sin( TWO_PI / ( 2.0 * static_cast< double >( sides ) ) );

      return perimeter * TUBE_HEIGHT * 0.5;
    };

  printf( "=== the tube is closed in u by evaluation ===\n" );

  check( glm::distance( tinynurbs::surfacePoint( surface, 0.0, 0.5 ),
                        tinynurbs::surfacePoint( surface, 1.0, 0.5 ) ) == 0.0,
         "the two ends of the u domain evaluate to the same point" );

  printf( "=== two rims make a chart, with no cut to place ===\n" );

  {
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( outcome.built, "the chart is triangulated" );

    // The band is a ribbon between two 24-gons inscribed in the tube, so its
    // area is that 24-gon's perimeter times the height - exactly, not to a
    // tolerance. Anything that covers it twice, or drops a rim, misses this by
    // a factor rather than by an epsilon.
    const double polygonal = bandArea( 24 );

    check( std::abs( emittedArea( state, outcome ) - polygonal ) <
             ( polygonal * 1e-9 ),
           "and covers the whole band between the rims, exactly once" );

    check( worstLiftedSpan( state, outcome ) < 0.5,
           "no emitted triangle spans half a period in the lifted chart" );

    // The seam is emergent: the triangulation chose where it runs, and the
    // only cost is one duplicated vertex per chart vertex it passes through.
    const Growth growth = classifyGrowth( state, boundary );

    check( growth.introduced == 0,
           "nothing is introduced - the mesh grows only by seam duplicates" );
    check( growth.duplicates > 0 && growth.duplicates <= 4,
           "and the seam is a short chain of them" );
  }

  printf( "=== the seam's duplicates weld by identity, not by periodicity ===\n" );

  {
    // The old cut's two edges welded in Geometry::Reify only because the
    // surface really was closed in u and they were evaluated a period apart. A
    // duplicate made here is a BITWISE copy of its original's world position,
    // so welder.weld( *this, DBL_EPSILON ) closes the seam by identity. That is
    // strictly stronger, and this is what says so.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( outcome.built, "the chart is triangulated" );

    size_t exact = 0;

    for ( size_t at = boundary; at < state.mesh.vertices.size(); ++at ) {
      for ( size_t other = 0; other < boundary; ++other ) {

        if ( state.mesh.vertices[ at ].point ==
               state.mesh.vertices[ other ].point ) {

          ++exact;
          break;
        }
      }
    }

    check( exact == ( state.mesh.vertices.size() - boundary ),
           "every added vertex is a bit-identical copy of a boundary vertex" );

    // And the uv it carries is a whole number of periods from its original's,
    // which is what makes it the same point on the surface.
    size_t lifted = 0;

    for ( size_t at = boundary; at < state.mesh.vertices.size(); ++at ) {
      for ( size_t other = 0; other < boundary; ++other ) {

        if ( state.mesh.vertices[ at ].point !=
               state.mesh.vertices[ other ].point ) {
          continue;
        }

        const double gap = std::abs( state.mesh.vertices[ at ].uv.x -
                                     state.mesh.vertices[ other ].uv.x );

        if ( std::abs( gap - std::round( gap / outcome.period ) *
                               outcome.period ) < 1e-12 ) {
          ++lifted;
        }

        break;
      }
    }

    check( lifted == ( state.mesh.vertices.size() - boundary ),
           "and sits a whole number of periods away in u" );
  }

  printf( "=== a rim closed by adjacency builds, at any sampling ===\n" );

  {
    // WAS A REFUSAL. `netDelta` was measured head-to-tail over the sampled
    // points, so a rim that does not repeat its head measured P - P/n: refused
    // at ordinary sampling, and past ~1000 samples P/n fell inside periodSlack
    // and it was ACCEPTED with cut edges that were not periodic copies and
    // welded to nothing. Neither reading exists now - a ring's closing segment
    // is a constraint edge like every other segment - so both spellings are
    // the same input and both build.
    Built coarse = build( { evenRim( 24, 0.25, true, false ),
                            evenRim( 24, 0.75, false, false ) } );

    const Outcome coarseOutcome = run( coarse, surface );

    check( coarseOutcome.built,
           "a coarsely sampled rim that does not repeat its head builds" );

    Built fine = build( { evenRim( 1200, 0.25, true, false ),
                          evenRim( 1200, 0.75, false, false ) } );

    check( glm::distance( fine.mesh.vertices.front().point,
                          fine.mesh.vertices[ 1199 ].point ) > 0.01,
           "the finely sampled rim's first and last points really are apart" );

    const Outcome fineOutcome = run( fine, surface );

    check( fineOutcome.built,
           "and so does the 1200-sample spelling that used to be accepted "
           "wrongly" );
    check( worstLiftedSpan( fine, fineOutcome ) < 0.5,
           "with no triangle spanning half a period" );

    // The one that mattered: those cut edges welded nothing. There is no cut,
    // and every vertex this added is a bit-identical copy.
    const Growth growth = classifyGrowth( fine, 2400 );

    check( growth.introduced == 0,
           "and every vertex it adds is a copy that welds in Reify" );
  }

  printf( "=== a rim that wanders back over itself builds, with no cut to "
          "search for ===\n" );

  {
    // WAS A 60-CANDIDATE SEARCH. A rim that is not u-monotone can wander back
    // across its own shortest cut, so the old code sorted every rotation of
    // the opposite rim by cut length and took the first that crossed nothing.
    // `ADVANCED_FACE #19215` of Right_Hand.step has exactly one clean candidate
    // out of 60. There is no cut here, so there is nothing to be crossed and
    // nothing to search.
    RingSpec wandering;

    wandering.explicitPoints = { { 0.00, 0.90 }, { 0.20, 0.90 }, { 0.40, 0.90 },
                                 { 0.60, 0.90 }, { 0.75, 0.90 }, { 0.85, 0.70 },
                                 { 0.95, 0.50 }, { 0.04, 0.48 }, { 0.94, 0.52 },
                                 { 0.96, 0.70 }, { 0.00, 0.90 } };

    Built state = build( { wandering, evenRim( 24, 0.10, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "a rim that is not u-monotone builds" );
    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
  }

  printf( "=== a ring whose two ends disagree in uv is an ordinary ring ===\n" );

  {
    // WAS A ROUNDED WINDING READING. The two ends are ONE 3D point the solve
    // answered twice, 0.003 of a period apart in u - three times the 1e-3 slack
    // the old code measured against, which is why that reading had to be
    // rounded rather than measured. `ADVANCED_FACE #19215`'s seam-straddling
    // hole is this shape. Nothing reads a winding per ring now: the ring is a
    // closed loop with one short closing segment, like any other.
    RingSpec hole;

    hole.explicitPoints = { { 0.45, 0.45 }, { 0.55, 0.45 },
                            { 0.55, 0.55 }, { 0.45, 0.55 } };
    hole.tailUv         = { { 0.453, 0.4508 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           hole } );

    const size_t head = 50;
    const size_t tail = state.mesh.vertices.size() - 1;

    check( state.mesh.vertices[ head ].point ==
             state.mesh.vertices[ tail ].point,
           "the hole's two ends are the same 3D point" );

    check( std::abs( state.rings[ 2 ].back()[ 0 ] -
                     state.rings[ 2 ].front()[ 0 ] ) > 0.001,
           "and disagree in u by more than the old fixed slack" );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "the chart is triangulated anyway" );

    // And the hole is a hole: the band minus the box, not the band.
    check( emittedArea( state, outcome ) < ( bandArea( 24 ) * 0.995 ),
           "with the hole carved out of it" );
  }

  printf( "=== an ordinary hole is carved, and the two spellings agree ===\n" );

  {
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           smallRing( 0.5, 0.5, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "a hole inside the band still builds the chart" );
    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
  }

  printf( "=== rings that intersect are refused, by one reading ===\n" );

  {
    // THREE OLD GATES, ONE CHECK. `IntersectingConstraintEdges::TryResolve`
    // resolves intersections by ADDING vertices rather than by rejecting the
    // input, so a vertex set that came back larger than it went in is the
    // signal that the rings this was handed cross - and it refuses the whole
    // face. What used to need a cut-crossing test on listed segments, a second
    // one on the implicit closing segment, and a per-point hole-containment
    // test is this one reading of the triangulation's own output.
    //
    // A hole reaching out through the top rim. Old code: refused by hole
    // containment, at every candidate cut.
    RingSpec reachingOut;

    reachingOut.explicitPoints = { { 0.45, 0.45 }, { 0.55, 0.45 },
                                   { 0.55, 0.95 }, { 0.45, 0.55 },
                                   { 0.45, 0.45 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           reachingOut } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built, "a hole that reaches out through a rim is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  {
    // The crossing hole in BOTH spellings. Old code refused the head-repeated
    // one by the listed segment and BUILT the other, because the gate walked
    // only listed segments and a point-list bound does not repeat its head.
    // Here the closing segment is a constraint edge either way, so the two
    // spellings are one input - which is the defect class disappearing rather
    // than a wider gate covering it.
    RingSpec unlisted;

    unlisted.explicitPoints = { { 0.9996, 0.40 },
                                { 0.90,   0.90 },
                                { 1.05,   0.90 },
                                { 1.0004, 0.45 } };

    RingSpec listed = unlisted;

    listed.explicitPoints.push_back( { 0.9996, 0.40 } );

    Built unlistedState = build( { evenRim( 24, 0.25, true, true ),
                                   evenRim( 24, 0.75, false, true ),
                                   unlisted } );

    Built listedState = build( { evenRim( 24, 0.25, true, true ),
                                 evenRim( 24, 0.75, false, true ),
                                 listed } );

    const Outcome unlistedOutcome = run( unlistedState, surface );
    const Outcome listedOutcome   = run( listedState, surface );

    check( !unlistedOutcome.built,
           "a hole crossing a rim with its closing segment unlisted is "
           "refused" );
    check( !listedOutcome.built,
           "and so is the same hole with its head repeated" );
    check( unlistedOutcome.built == listedOutcome.built,
           "the two spellings are one input now, not two cases" );
    check( unlistedOutcome.vertexGrowth == 0 &&
             listedOutcome.vertexGrowth == 0,
           "and neither touches the mesh" );
  }

  {
    // A ring that crosses ITSELF. Old code refused it for disagreeing with its
    // own head by 0.3 of a period, which was a reading of the winding; here it
    // is refused for what is actually wrong with it.
    RingSpec selfCrossing;

    selfCrossing.explicitPoints = { { 0.45, 0.45 }, { 0.55, 0.45 },
                                    { 0.55, 0.55 }, { 0.45, 0.55 } };
    selfCrossing.tailUv         = { { 0.75, 0.4508 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           selfCrossing } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built, "a ring that crosses itself is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== the declared u_closed is the first gate ===\n" );

  {
    // PORTED UNCHANGED. A periodic chart is a consequence of the author's
    // topology, and no amount of sampling makes a surface closed that its
    // author did not close.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, surface, false );

    check( !outcome.built,
           "a surface the file does not declare closed in u is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== five probes were not a closure test ===\n" );

  {
    // PORTED UNCHANGED. The two u-end isocurves are splines; they can agree at
    // any five chosen v and part company between them. A degree-1 surface with
    // eight v knot spans does it exactly. A chart built on such a surface has
    // an emergent seam whose two sides are NOT one point on the surface, so it
    // welds nothing and the face ships open.
    const tinynurbs::RationalSurface3d split = makeSplitSeamTube();

    double worstProbe = 0.0;

    for ( size_t at = 0; at < 5; ++at ) {

      const double v = static_cast< double >( at ) / 4.0;

      worstProbe =
        std::max( worstProbe,
                  glm::distance( tinynurbs::surfacePoint( split, 0.0, v ),
                                 tinynurbs::surfacePoint( split, 1.0, v ) ) );
    }

    check( worstProbe == 0.0,
           "the split-seam tube closes EXACTLY at the five parameters the old "
           "gate read" );

    check( glm::distance( tinynurbs::surfacePoint( split, 0.0, 0.125 ),
                          tinynurbs::surfacePoint( split, 1.0, 0.125 ) ) >
             ( SPLIT_SEAM_GAP * 0.5 ),
           "and stands a seam gap apart between them, so it is not closed" );

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, split );

    check( !outcome.built,
           "a surface that closes only where the old probes landed is "
           "refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a boundary that does not wrap is left to earcut ===\n" );

  {
    // The reading that chooses between two triangulators, and the only thing
    // it can cost: a boundary that nets no winding is an ordinary outer ring
    // with holes, which earcut reads correctly and more cheaply.
    Built state = build( { smallRing( 0.3, 0.5, 0.10 ),
                           smallRing( 0.7, 0.5, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built, "a boundary that does not wrap the closure is "
                           "refused, for earcut to take" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== chord subdivision keeps a coarse rim off the inner rim ===\n" );

  {
    // BUILT IN, NOT GATED. The polar layout distorts long chords: a segment at
    // radius r spanning dTheta dips to r * cos( dTheta / 2 ). A four-sample rim
    // steps a quarter period - 90 degrees - so without subdivision its chords
    // dip to 0.707 of their radius, well inside the inner rim, and the CDT
    // triangulates a shape that is not the band. With MAX_CHART_CHORD_FRACTION
    // at 1/12 each of those steps becomes three, and the area comes out right.
    Built state = build( { evenRim( 4, 0.25, true, true ),
                           evenRim( 4, 0.75, false, true ) } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( outcome.built, "a rim sampled every quarter period builds" );

    const Growth growth = classifyGrowth( state, boundary );

    check( growth.introduced > 0,
           "and subdivision introduced points, because 90 degrees is past the "
           "threshold" );

    // The band the subdivided boundary bounds is the 12-gon's, not the
    // 4-gon's: a quarter-period step splits into three 30-degree ones and the
    // added points are evaluated ON the tube, at the twelfths that are its
    // control rows. Exact, because those land on the prism's own corners.

    check( std::abs( emittedArea( state, outcome ) - bandArea( 12 ) ) <
             ( bandArea( 12 ) * 1e-9 ),
           "and the emitted area is the subdivided band's, exactly" );

    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
  }

  printf( "=== a boundary that touches itself without crossing ===\n" );

  {
    // Codex 4049291629 on bldrs-ai/conway-geom#207 attacked the cut search's
    // crossing test for detecting only PROPER crossings: a candidate cut that
    // touched a boundary vertex, or ran collinear along a boundary segment,
    // made the rewritten polygon self-touching, and earcut has no answer for
    // one. There is no candidate cut any more and no polygon is built, so the
    // question that remains is what a self-touching BOUNDARY does here.
    //
    // It is a legal constraint configuration. A hole whose corner sits exactly
    // on the lower rim shares that vertex with it after the 1e-9 weld, and the
    // parity peel gives the region those edges bound - pinched at the contact,
    // which is what the input says. Nothing is dropped and nothing is paved:
    // the band minus the hole, to nine digits. Measured rather than argued,
    // because "CDT does not need a simple polygon" is the claim the whole
    // change rests on and a touching contact is where it is least obvious.
    RingSpec touching;

    // The rims sit at v = 0.25 and v = 0.75; this hole's bottom edge lies ON
    // the lower rim, sharing two of its sample points exactly.
    touching.explicitPoints = { { 5.0 / 24.0, 0.25 }, { 7.0 / 24.0, 0.25 },
                                { 7.0 / 24.0, 0.45 }, { 5.0 / 24.0, 0.45 },
                                { 5.0 / 24.0, 0.25 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           touching } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "a hole touching a rim without crossing it builds" );

    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );

    // The hole is carved, not ignored: the band less the 2/24-wide, 0.2-tall
    // patch it takes out of it. Both are exact on the 24-gon.
    const double removed =
      2.0 * 2.0 * TUBE_RADIUS * std::sin( TWO_PI / 48.0 ) * TUBE_HEIGHT * 0.20;

    check( std::abs( emittedArea( state, outcome ) -
                     ( bandArea( 24 ) - removed ) ) < ( bandArea( 24 ) * 1e-9 ),
           "and exactly that hole is missing from the band" );
  }

  printf( "=== a boundary that winds but does not bound a strip is refused ===\n" );

  {
    // THE SECOND POST-HOC READING, and the state that reaches it. Two rims
    // bound a strip and the disc inside the inner one is erased, which is what
    // keeps the origin out of every kept triangle and therefore what makes
    // each triangle's nearest-image lift determined. ONE wrapping rim plus an
    // ordinary hole is not a strip: nothing erases the centre, the disc inside
    // the rim is kept by parity, and its triangles span the origin - for which
    // there is no nearest image, because both are equally far.
    //
    // Read off the triangulation rather than argued: a triangle spans less
    // than half a period exactly when its three pairwise nearest-image offsets
    // sum to zero, and a triangle containing the origin cannot. Confirmed to
    // be THIS check that refuses, by instrumenting it - it is a state the
    // construction produces, not a defense against one it cannot.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           smallRing( 0.5, 0.75, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "one wrapping rim with no second rim to close the strip is "
           "refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a bound clear of every other bound is its own island ===\n" );

  {
    // CHARACTERISATION, and a genuine change. This used to be refused by hole
    // containment: a ring at v = 4.0 is nowhere near the band between the rims.
    // `eraseOuterTrianglesAndHoles` decides interior by parity, and a closed
    // loop in the chart has an interior, so it is now triangulated as its own
    // patch. That is the same semantics triangulateUnwrappedLoops has shipped
    // with on the cylinder and cone paths. Recorded, not gated: no reading of
    // the input can tell a stray bound from a legitimate disjoint one, and the
    // gate that used to try cost more than it caught.
    Built banded = build( { evenRim( 24, 0.25, true, true ),
                            evenRim( 24, 0.75, false, true ) } );

    const Outcome bandedOutcome = run( banded, surface );

    Built strayed = build( { evenRim( 24, 0.25, true, true ),
                             evenRim( 24, 0.75, false, true ),
                             smallRing( 0.5, 4.0, 0.05 ) } );

    const Outcome strayedOutcome = run( strayed, surface );

    check( strayedOutcome.built,
           "a bound clear of the band is triangulated rather than refused" );

    check( strayedOutcome.triangles.size() > bandedOutcome.triangles.size(),
           "as extra triangles, not as a replacement for the band" );

    check( emittedArea( strayed, strayedOutcome ) >
             emittedArea( banded, bandedOutcome ),
           "and it adds area rather than inverting what was kept" );
  }

  printf( failures == 0 ? "PASS\n" : "FAIL (%d)\n", failures );

  return failures == 0 ? 0 : 1;
}
