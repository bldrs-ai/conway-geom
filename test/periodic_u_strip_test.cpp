/*
 * Refusal tests for tryPeriodicUStrip in mesh_utils.h.
 *
 * The function's safety claim is that every one of its gates refuses by
 * returning false and leaving the caller exactly the state it was handed, so an
 * unfamiliar face degrades to today's ear-clipping rather than to a guess. Three
 * ways that claim was not true, all found by review on bldrs-ai/conway-geom#207
 * and all pinned here:
 *
 *   1. A RIM CLOSED BY ADJACENCY. `netDelta` is measured head-to-tail over the
 *      sampled points, so it is the loop's full winding only when the closing
 *      edge is degenerate - i.e. when the rim repeats its first point as its
 *      last. A point-list bound is a closed polygon that does NOT repeat its
 *      head (IfcCurve::closedByConstruction, GetLoop), and such a rim measures
 *      P - P/n: refused at ordinary sampling, but past ~1000 samples P/n falls
 *      inside periodSlack and the rim is ACCEPTED with distinct endpoints. The
 *      two synthesized cut edges are then not periodic copies of one segment,
 *      weld nothing in Reify, and leave the strip open along the cut.
 *
 *   2. AN AMBIGUOUS UNWRAP STEP. Nearest-image unwrapping is a reading only
 *      while the runner-up image is decisively further. Two steps whose genuine
 *      motion exceeds half a period are each read a full period wrong, and if
 *      they are opposite in sign the errors CANCEL - `netDelta` still measures
 *      exactly one period, and the strip is built round a boundary routed
 *      through chart the face never visits. The cut-crossing check cannot see
 *      it: that tests only the two synthesized cuts against the rings.
 *
 *   3. A HOLE'S IMPLICIT CLOSING SEGMENT. The cut-crossing gate walked each
 *      ring's LISTED segments, which for a point-list bound - one that does
 *      not repeat its head - leaves `back() -> front()` untested. The
 *      net-delta gate does not cover it either: an ordinary ring is admitted
 *      whenever its head-to-tail u change is within periodSlack, which is
 *      exactly what a finely sampled ring has, so its closing edge is short in
 *      u, near-parallel to the cut, and in the best orientation to cross it.
 *      A hole crossing the cut is handed to earcut as a hole crossing the
 *      outer polygon.
 *
 *   4. FIVE PROBES ARE NOT A CLOSURE TEST. The two u-end isocurves are
 *      splines; they can agree at any five chosen v and part company between
 *      them. A degree-1 surface with eight v knot spans does it exactly. The
 *      gate now takes the file's own `u_closed` declaration first, and checks
 *      it per knot span with enough samples to settle each span's polynomial.
 *
 *   5. A LATE REFUSAL THAT MUTATES. The chain's closing duplicate was added
 *      with mesh.makeVertex before the hole-containment and cut-crossing gates
 *      had run, so refusing at either left an extra vertex referenced by no
 *      triangle. This is not latent: `ADVANCED_FACE #19215` of
 *      `Right_Hand.step` reaches the chain build and then refuses, leaking one
 *      vertex per call on the corpus as it stands.
 *
 * Red-proven by reverting rather than by reading: with mesh_utils.h alone rolled
 * back to the reviewed revision, the three tests named above go red (the
 * adjacency-closed rim is built, the ambiguous rim is built, and the refused
 * face grows the mesh by one vertex) and the run exits 1. The remaining tests
 * are characterisations that hold on both revisions and say so where they stand.
 *
 * The surface is synthetic and minimal on purpose: a degree-1 tube whose last
 * control row IS its first, which is closed in u by evaluation - the only
 * property of the surface tryPeriodicUStrip reads, besides its knot domain.
 * `makeSplitSeamTube` is the same tube with more v knot spans and a seam that
 * only closes on the spans the old five probes landed on.
 *
 * Standalone by design: it includes mesh_utils.h directly and links nothing but
 * the Logger stubs below, matching outer_bound_order_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <array>
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

/** Runs the function under test and reports what it did to the mesh. */
struct Outcome {

  bool   built;
  size_t vertexGrowth;
};

Outcome run( Built& state, const tinynurbs::RationalSurface3d& surface,
             bool declaredClosedU = true ) {

  std::vector< uint32_t > flatToVertex;

  double period = 0.0;
  double uMin   = 0.0;

  const size_t before = state.mesh.vertices.size();

  const bool built =
    conway::geometry::tryPeriodicUStrip(
      state.mesh, surface, declaredClosedU, state.rings, flatToVertex, period,
      uMin );

  return { built, state.mesh.vertices.size() - before };
}

}  // namespace


int main() {

  const tinynurbs::RationalSurface3d surface = makeClosedTube();

  printf( "=== the tube is closed in u by evaluation ===\n" );

  check( glm::distance( tinynurbs::surfacePoint( surface, 0.0, 0.5 ),
                        tinynurbs::surfacePoint( surface, 1.0, 0.5 ) ) == 0.0,
         "the two ends of the u domain evaluate to the same point" );

  printf( "=== two rims that repeat their head make a strip ===\n" );

  {
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "the strip is built" );

    // The cut closes chain B on a duplicate of its opening vertex, and that is
    // the ONLY vertex this function may ever add.
    check( outcome.vertexGrowth == 1,
           "success adds exactly one vertex - the chain's closing duplicate" );
  }

  printf( "=== a rim closed by adjacency is refused ===\n" );

  {
    // 1200 samples: the missing closing interval is P/1200 = 8.3e-4, INSIDE the
    // 1e-3 periodSlack, so the head-to-tail winding reads as a full period
    // while the rim's first and last points are 0.05 radians apart in 3D. This
    // is the case that used to be built, and the one whose cut edges would not
    // have welded. RED on the reviewed revision.
    Built state = build( { evenRim( 1200, 0.25, true, false ),
                           evenRim( 1200, 0.75, false, false ) } );

    check( glm::distance( state.mesh.vertices.front().point,
                          state.mesh.vertices[ 1199 ].point ) > 0.01,
           "the adjacency-closed rim's first and last points really are apart" );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a finely sampled rim that does not repeat its head is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  {
    // The same spelling at ordinary sampling, where the missing interval is far
    // outside periodSlack. Characterisation: refused on both revisions, by the
    // net-delta gate rather than by the closure test.
    Built state = build( { evenRim( 24, 0.25, true, false ),
                           evenRim( 24, 0.75, false, false ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a coarsely sampled rim that does not repeat its head is refused" );
  }

  printf( "=== ambiguous unwrap steps are refused ===\n" );

  {
    // A rim whose genuine motion advances 0.6 of a period and then retreats by
    // the same amount before completing the period in small steps. Each of
    // those two steps is read a full period wrong and the errors CANCEL, so the
    // head-to-tail winding still measures exactly +1 period and the net-delta
    // gate sees nothing. RED on the reviewed revision.
    RingSpec ambiguous;

    ambiguous.v          = 0.25;
    ambiguous.repeatHead = true;
    ambiguous.us         = { 0.0, 0.6, 0.0 };

    for ( size_t at = 1; at < 20; ++at ) {
      ambiguous.us.push_back( 0.05 * static_cast< double >( at ) );
    }

    Built state = build( { ambiguous, evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a rim with two cancelling half-period steps is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a late refusal leaves the mesh untouched ===\n" );

  {
    // Two ordinary rims and a ring that lands OUTSIDE the strip in v. Every
    // gate up to and including the chain build passes; hole containment is what
    // refuses. The chain's closing duplicate is built by then, so on the
    // reviewed revision the mesh comes back one vertex longer than it went in.
    // RED there on the growth assertion.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           smallRing( 0.5, 4.0, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built, "a hole outside the strip refuses the face" );
    check( outcome.vertexGrowth == 0,
           "a refusal after the chain is built adds no vertex" );
  }

  printf( "=== the declared u_closed is the first gate ===\n" );

  {
    // The same two rims on the same closed tube, with the file declaring the
    // surface OPEN in u. Nothing this builds - the cut, the two edges that have
    // to be periodic copies of one another - means anything on a surface whose
    // author did not close it, so the declaration is read before any of it.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, surface, false );

    check( !outcome.built,
           "a surface the file does not declare closed in u is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== five probes were not a closure test ===\n" );

  {
    const tinynurbs::RationalSurface3d split = makeSplitSeamTube();

    // The five parameters the old gate read, spread over the whole v domain.
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

    // One span in, where no probe landed.
    check( glm::distance( tinynurbs::surfacePoint( split, 0.0, 0.125 ),
                          tinynurbs::surfacePoint( split, 1.0, 0.125 ) ) >
             ( SPLIT_SEAM_GAP * 0.5 ),
           "and stands a seam gap apart between them, so it is not closed" );

    // Declared closed - a file that says .T. about a surface that is not. This
    // is what the numerical half of the gate is for. RED on the reviewed
    // revision, which built the strip.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ) } );

    const Outcome outcome = run( state, split );

    check( !outcome.built,
           "a surface that closes only where the old probes landed is "
           "refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a hole's implicit closing segment crosses the cut ===\n" );

  {
    // A four-point hole that does NOT repeat its head, so its closing segment
    // (1.0004, 0.45) -> (0.9996, 0.40) is not listed. That segment crosses the
    // cut at u = 1 between the two rims; none of the listed three does, because
    // the one that reaches over the cut does so at v = 0.9, clear of the cut's
    // own span. Its head-to-tail u change is 8e-4 of a period, inside
    // periodSlack, so the net-delta gate reads it as an ordinary ring and the
    // centroid lands inside the strip. RED on the reviewed revision, which
    // built the strip and handed earcut a hole crossing its outer polygon.
    RingSpec hole;

    hole.explicitPoints = { { 0.9996, 0.40 },
                            { 0.90,   0.90 },
                            { 1.05,   0.90 },
                            { 1.0004, 0.45 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           hole } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a hole whose unlisted closing segment crosses the cut is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  {
    // The SAME ring, with its head repeated so the crossing segment is listed.
    // Refused on both revisions - which is what isolates the mechanism above:
    // the geometry did not change, only whether the gate could see it.
    RingSpec hole;

    hole.explicitPoints = { { 0.9996, 0.40 },
                            { 0.90,   0.90 },
                            { 1.05,   0.90 },
                            { 1.0004, 0.45 },
                            { 0.9996, 0.40 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           hole } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "the same hole with its head repeated is refused by the listed "
           "segment" );
  }

  {
    // And an ordinary hole that crosses nothing still builds, so the wider walk
    // did not simply refuse everything. A ring that repeats its head adds a
    // zero-length closing segment; a degenerate segment cannot report a
    // crossing, and this is what says so.
    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           smallRing( 0.5, 0.5, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "a hole inside the strip still builds the strip" );
    check( outcome.vertexGrowth == 1,
           "and still adds exactly one vertex" );
  }


  printf( "=== a ring whose two ends disagree in uv is still an ordinary ring ===\n" );

  {
    // The two ends are ONE 3D point the solve answered twice, 0.003 of a period
    // apart in u - three times the 1e-3 periodSlack this used to measure
    // against, and the reading that refused `ADVANCED_FACE #19215` outright
    // once its seam-straddling hole stopped being clamped flat. A ring that
    // closes on its own head winds a whole number of periods, so the winding is
    // ROUNDED rather than measured. RED before that change.
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

    check( outcome.built,
           "the strip is built anyway - the winding rounds to zero" );
    check( outcome.vertexGrowth == 1, "and still adds exactly one vertex" );
  }

  {
    // The rounding is bounded, not unbounded: a ring whose ends disagree by
    // 0.3 of a period is not decisively any winding - the same quarter-period
    // margin the unwrap step demands - and is still refused. Characterisation
    // on the old revision, which refused it for being outside 1e-3.
    RingSpec hole;

    hole.explicitPoints = { { 0.45, 0.45 }, { 0.55, 0.45 },
                            { 0.55, 0.55 }, { 0.45, 0.55 } };
    hole.tailUv         = { { 0.75, 0.4508 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           hole } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a ring whose ends disagree by 0.3 of a period is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== the shortest cut is not the only cut tried ===\n" );

  {
    // A rim that is not u-monotone: it runs forward along v = 0.90, dips to
    // v ~ 0.50 while overshooting the seam to u = 1.04, comes back to u = 0.94,
    // and finishes at u = 1.00. Net winding is exactly one period and its two
    // ends are one point, so it is a rim - but its excursion straddles u = 1.00
    // at v ~ 0.49, which is precisely where the SHORTEST cut runs.
    //
    // That is `ADVANCED_FACE #19215`'s shape: a 213-point rim advancing 141
    // steps and retreating 71 over 1.066 periods, whose shortest cut is crossed
    // by the rim's own segment and exactly one of whose 60 candidate cuts
    // crosses nothing. RED before the candidate search.
    RingSpec wandering;

    wandering.explicitPoints = { { 0.00, 0.90 }, { 0.20, 0.90 }, { 0.40, 0.90 },
                                 { 0.60, 0.90 }, { 0.75, 0.90 }, { 0.85, 0.70 },
                                 { 0.95, 0.50 }, { 0.04, 0.48 }, { 0.94, 0.52 },
                                 { 0.96, 0.70 }, { 0.00, 0.90 } };

    Built state = build( { wandering, evenRim( 24, 0.10, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built,
           "a strip whose shortest cut is crossed is built on a longer one" );
    check( outcome.vertexGrowth == 1, "and still adds exactly one vertex" );
  }

  {
    // The search does not bypass the gates it is searching against. The same
    // wandering rim, with a hole whose centroid is inside the strip but one of
    // whose corners reaches out through the top rim: refused at every
    // candidate, and the mesh left exactly as it was handed over.
    RingSpec wandering;

    wandering.explicitPoints = { { 0.00, 0.90 }, { 0.20, 0.90 }, { 0.40, 0.90 },
                                 { 0.60, 0.90 }, { 0.75, 0.90 }, { 0.85, 0.70 },
                                 { 0.95, 0.50 }, { 0.04, 0.48 }, { 0.94, 0.52 },
                                 { 0.96, 0.70 }, { 0.00, 0.90 } };

    RingSpec reachingOut;

    // Centroid v = 0.60, inside; the third corner stands at v = 0.95, above the
    // top rim at v = 0.90.
    reachingOut.explicitPoints = { { 0.45, 0.45 }, { 0.55, 0.45 },
                                   { 0.55, 0.95 }, { 0.45, 0.55 },
                                   { 0.45, 0.45 } };

    Built state = build( { wandering, evenRim( 24, 0.10, false, true ),
                           reachingOut } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a hole that reaches out through a rim is refused at every cut" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( failures == 0 ? "PASS\n" : "FAIL (%d)\n", failures );

  return failures == 0 ? 0 : 1;
}
