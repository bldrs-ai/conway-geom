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
#include <map>
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

constexpr double CONST_PI_TEST = 3.141592653589793;

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

/**
 * THE SURFACE A STEP CAN REDUCE THE WRONG WAY ON, and the reason the
 * nearest-image reduction needs confirming against something.
 *
 * A closed degree-1 prism whose cross-section is the rectangle
 * [-10, 10] x [0, 6], with its four corners at u = 0, 0.8, 0.9 and 0.95 of
 * the period. Nothing is odd about it except the knot spacing, which no
 * exporter is obliged to make uniform: the long bottom edge occupies 0.8 of
 * the u period and the other three share the remaining 0.2.
 *
 * WHY THAT MATTERS. The bottom edge is STRAIGHT in 3D, so a trim curve
 * running along it is sampled EXACTLY by its two endpoints - a tessellator
 * has no reason to put a point in the middle of a straight line. The
 * resulting polyline is a perfect representation of the boundary, and its one
 * segment steps +0.80 of a period. Reduced to its nearest image that reads
 * -0.20, which is a route round the OTHER three sides, and every gate that
 * reads only u is content: -0.20 is well inside
 * MAX_DECISIVE_IMAGE_STEP_FRACTION.
 *
 * On the tube above the same spelling would be a coarse sampling of a curved
 * rim, and the short route would be the honest reading of it. Here the chord
 * lies ON the surface, the long route is the one it draws, and the two
 * readings are a 0.80-period band against a 0.20-period sliver of one.
 */
constexpr double PRISM_HALF_WIDTH = 10.0;
constexpr double PRISM_DEPTH      = 6.0;

/** The prism's profile corners, as fractions of the u period. */
constexpr double PRISM_CORNER_U[ 5 ] = { 0.0, 0.8, 0.9, 0.95, 1.0 };

tinynurbs::RationalSurface3d makeFlatSidedPrism() {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  const glm::dvec2 corners[ 5 ] = {
    {  PRISM_HALF_WIDTH, 0.0 },
    { -PRISM_HALF_WIDTH, 0.0 },
    { -PRISM_HALF_WIDTH, PRISM_DEPTH },
    {  PRISM_HALF_WIDTH, PRISM_DEPTH },
    {  PRISM_HALF_WIDTH, 0.0 } };

  std::vector< glm::dvec3 > control;

  control.reserve( 10 );

  for ( size_t row = 0; row < 5; ++row ) {

    control.push_back( { corners[ row ].x, corners[ row ].y, 0.0 } );
    control.push_back( { corners[ row ].x, corners[ row ].y, TUBE_HEIGHT } );
  }

  surface.control_points = tinynurbs::array2( 5, 2, control );
  surface.weights =
    tinynurbs::array2( 5, 2, std::vector< double >( 10, 1.0 ) );

  // Degree 1, five control rows: the clamped knot vector repeats the two ends
  // and carries the corner parameters between them.
  surface.knots_u = { PRISM_CORNER_U[ 0 ], PRISM_CORNER_U[ 0 ],
                      PRISM_CORNER_U[ 1 ], PRISM_CORNER_U[ 2 ],
                      PRISM_CORNER_U[ 3 ],
                      PRISM_CORNER_U[ 4 ], PRISM_CORNER_U[ 4 ] };

  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  return surface;
}

/** The prism point a (u, v) names, by the same linear interpolation. */
glm::dvec3 prismPoint( double u, double v ) {

  const double wrapped = u - std::floor( u );

  size_t corner = 0;

  while ( corner < 3 && wrapped >= PRISM_CORNER_U[ corner + 1 ] ) {
    ++corner;
  }

  const double fraction =
    ( wrapped - PRISM_CORNER_U[ corner ] ) /
    ( PRISM_CORNER_U[ corner + 1 ] - PRISM_CORNER_U[ corner ] );

  const glm::dvec2 profile[ 5 ] = {
    {  PRISM_HALF_WIDTH, 0.0 },
    { -PRISM_HALF_WIDTH, 0.0 },
    { -PRISM_HALF_WIDTH, PRISM_DEPTH },
    {  PRISM_HALF_WIDTH, PRISM_DEPTH },
    {  PRISM_HALF_WIDTH, 0.0 } };

  const glm::dvec2 at =
    profile[ corner ] +
    ( ( profile[ corner + 1 ] - profile[ corner ] ) * fraction );

  return { at.x, at.y, TUBE_HEIGHT * v };
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
using PointOfUv = glm::dvec3 ( * )( double, double );

Built build( const std::vector< RingSpec >& specs,
             PointOfUv                      point = tubePoint ) {

  Built built;

  for ( const RingSpec& spec : specs ) {

    std::vector< std::array< double, 2 > > ring;

    if ( !spec.explicitPoints.empty() ) {

      for ( const std::array< double, 2 >& explicitPoint :
              spec.explicitPoints ) {

        ring.push_back( explicitPoint );
        built.mesh.makeVertex(
          { point( explicitPoint[ 0 ], explicitPoint[ 1 ] ),
            glm::dvec2( explicitPoint[ 0 ], explicitPoint[ 1 ] ) } );
      }

      if ( !spec.tailUv.empty() ) {

        ring.push_back( spec.tailUv.front() );
        built.mesh.makeVertex(
          { point( spec.explicitPoints.front()[ 0 ],
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
          { point( ring.front()[ 0 ], ring.front()[ 1 ] ),
            glm::dvec2( ring.front()[ 0 ], ring.front()[ 1 ] ) } );
        break;
      }

      const double v =
        spec.halfHeight == 0.0 ?
          spec.v :
          ( ( at == 0 || at == 1 ) ? spec.v - spec.halfHeight
                                   : spec.v + spec.halfHeight );

      ring.push_back( { spec.us[ at ], v } );
      built.mesh.makeVertex( { point( spec.us[ at ], v ),
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
  std::vector< std::array< uint32_t, 4 > > seamEdges;
  double                                   period;
  double                                   uMin;
};

Outcome run( Built& state, const tinynurbs::RationalSurface3d& surface,
             bool declaredClosedU = true ) {

  std::vector< std::array< uint32_t, 3 > > triangles;
  std::vector< std::array< uint32_t, 4 > > seamEdges;

  double period = 0.0;
  double uMin   = 0.0;

  const size_t before = state.mesh.vertices.size();

  const bool built =
    conway::geometry::triangulatePeriodicUChart(
      state.mesh, surface, declaredClosedU, state.rings,
      triangles, seamEdges, period, uMin );

  return { built, state.mesh.vertices.size() - before, std::move( triangles ),
           std::move( seamEdges ), period, uMin };
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


/**
 * The tube with a radius that bulges in v.
 *
 * `refineSeamPairs` is handed THIS and not `tubePoint`, because the test
 * tube is ruled in v: a cut runs from one rim to the other at nearly constant
 * u, so its chord lies exactly ON the flat tube and the deflection reading is
 * zero however coarse the cut is. Nothing would split, and a refinement test
 * that cannot make anything split pins nothing. The bulge is the smallest
 * change that gives a radial chord something to deviate from; the chart gates
 * above still read the flat tube, which is the surface they are about.
 *
 * THE BULGE VANISHES AT THE TWO RIMS the cases below use, so every boundary
 * vertex the chart hands over lies exactly on this surface too. A surface
 * that disagreed with its own boundary would make the deflection reading
 * bottom out at that disagreement instead of at zero, and the refinement
 * would subdivide until it ran out of budget - which is a property of the
 * mismatch and not of the pass under test.
 */
constexpr double BULGE_RIM_LOW  = 0.2;
constexpr double BULGE_RIM_HIGH = 0.8;

glm::dvec3 bulgedTubePoint( double u, double v ) {

  const double angle = TWO_PI * u;

  const double alongBand =
    ( v - BULGE_RIM_LOW ) / ( BULGE_RIM_HIGH - BULGE_RIM_LOW );

  const double radius =
    TUBE_RADIUS * ( 1.0 + ( 0.1 * std::sin( CONST_PI_TEST * alongBand ) ) );

  return { radius * std::cos( angle ), radius * std::sin( angle ),
           TUBE_HEIGHT * v };
}

/** Border edges of `mesh`, keyed by their two endpoint POSITIONS. */
std::map< std::array< double, 6 >, size_t > borderEdgesByPosition(
    const WingedEdgeMesh< ParameterVertex >& mesh ) {

  std::map< std::array< double, 6 >, size_t > keys;

  for ( const conway::geometry::Edge& edge : mesh.edges ) {

    // A fully detached edge is neither a border nor an interior edge: it is
    // what deleteTriangle() leaves behind when both its triangles go.
    if ( edge.triangles[ 0 ] == conway::geometry::EMPTY_INDEX ||
         !edge.border() ) {
      continue;
    }

    const glm::dvec3& first  = mesh.vertices[ edge.vertices[ 0 ] ].point;
    const glm::dvec3& second = mesh.vertices[ edge.vertices[ 1 ] ].point;

    std::array< double, 3 > low  = { first.x, first.y, first.z };
    std::array< double, 3 > high = { second.x, second.y, second.z };

    if ( high < low ) {
      std::swap( low, high );
    }

    ++keys[ { low[ 0 ], low[ 1 ], low[ 2 ],
              high[ 0 ], high[ 1 ], high[ 2 ] } ];
  }

  return keys;
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

  printf( "=== two bounds clear of each other are two faces, not one ===\n" );

  {
    // THE REASON THIS IS REFUSED CHANGED. It used to be the winding reading:
    // neither ring wraps, so the face was handed to earcut. That reading is
    // gone (see the non-wrapping outer-and-hole case below, which now
    // builds), and what refuses this is the component count - two rings that
    // enclose nothing of each other are two regions, and an ADVANCED_FACE is
    // one. Same answer, from a reading that cannot be wrong about it.
    Built state = build( { smallRing( 0.3, 0.5, 0.10 ),
                           smallRing( 0.7, 0.5, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built, "two bounds enclosing nothing of each other are "
                           "refused as two regions" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a coarse rim keeps the boundary its neighbour has ===\n" );

  {
    // THE LAYOUT IS REFINED; THE BOUNDARY IS NOT. A four-sample rim steps a
    // quarter period, past the 30-degree threshold, so the layout splits each
    // step into three. Those split points used to be EVALUATED ON THE
    // SURFACE, and this test used to pin the consequence as a feature: "the
    // emitted area is the subdivided band's, exactly" - the 12-gon's, where
    // the boundary handed over was a 4-gon. That is a face-local rewrite of a
    // SHARED trim segment, so the face on the other side of that edge keeps
    // the 4-gon and the two disagree by the sagitta. Codex 4051389895 on
    // bldrs-ai/conway-geom#207 read the test as demonstrating the mismatch,
    // which is exactly what it was doing.
    //
    // The split point is now the linear interpolation of the segment's own
    // two endpoints, so it lies ON the edge the neighbour kept: the split is
    // invisible in 3D, the emitted area is the 4-gon band's to the last
    // digit, and no crack is geometrically possible. What remains is a
    // T-junction, which cannot be closed from inside this function.
    Built state = build( { evenRim( 4, 0.25, true, true ),
                           evenRim( 4, 0.75, false, true ) } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( outcome.built, "a rim sampled every quarter period builds" );

    const Growth growth = classifyGrowth( state, boundary );

    check( growth.introduced > 0,
           "and the layout still splits it, because 90 degrees is past the "
           "threshold" );

    check( std::abs( emittedArea( state, outcome ) - bandArea( 4 ) ) <
             ( bandArea( 4 ) * 1e-9 ),
           "but the emitted area is the boundary polyline's own, exactly - "
           "the 4-gon band, not the 12-gon's" );

    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
  }

  printf( "=== every split point lies on the segment it splits ===\n" );

  {
    // THE PROPERTY THE AREA ASSERTION IS A CONSEQUENCE OF, stated directly
    // and over a boundary whose segments are not all the same length: every
    // vertex this function invents must be collinear with - and between - the
    // two boundary vertices of the segment it came from, because that segment
    // is the neighbouring face's edge. Read off the mesh rather than the
    // chart: any introduced vertex must lie on SOME boundary segment.
    Built state = build( { evenRim( 5, 0.25, true, true ),
                           evenRim( 3, 0.75, false, true ) } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( outcome.built, "a band with two differently sampled rims builds" );

    size_t introduced = 0;
    size_t onSegment  = 0;

    for ( size_t at = boundary; at < state.mesh.vertices.size(); ++at ) {

      const glm::dvec3& point = state.mesh.vertices[ at ].point;

      bool duplicate = false;

      for ( size_t other = 0; other < boundary && !duplicate; ++other ) {
        duplicate = state.mesh.vertices[ other ].point == point;
      }

      if ( duplicate ) {
        continue;
      }

      ++introduced;

      // Against every segment of every ring: collinear, and strictly between.
      for ( const std::vector< std::array< double, 2 > >& ring : state.rings ) {

        const size_t count = ring.size();

        for ( size_t step = 0; step < count; ++step ) {

          const glm::dvec3 a =
            tubePoint( ring[ step ][ 0 ], ring[ step ][ 1 ] );
          const glm::dvec3 b =
            tubePoint( ring[ ( step + 1 ) % count ][ 0 ],
                       ring[ ( step + 1 ) % count ][ 1 ] );

          const glm::dvec3 along = b - a;
          const double     length2 = glm::dot( along, along );

          if ( length2 == 0.0 ) {
            continue;
          }

          const double t = glm::dot( point - a, along ) / length2;

          if ( t <= 0.0 || t >= 1.0 ) {
            continue;
          }

          if ( glm::distance( point, a + ( along * t ) ) <
                 ( TUBE_RADIUS * 1e-12 ) ) {

            ++onSegment;
            step  = count;
            break;
          }
        }
      }
    }

    check( introduced > 0, "and the layout splits at least one of its steps" );
    check( introduced == onSegment,
           "with every invented vertex lying on the boundary segment it came "
           "from, which is the edge the neighbouring face keeps" );
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

    // The rims sit at v = 0.25 and v = 0.75; this hole's bottom edge IS one
    // segment of the lower rim, sharing both of its endpoints exactly. One
    // segment rather than two, so that no step here is wide enough for the
    // layout split and the contact is the only thing under test.
    touching.explicitPoints = { { 5.0 / 24.0, 0.25 }, { 6.0 / 24.0, 0.25 },
                                { 6.0 / 24.0, 0.45 }, { 5.0 / 24.0, 0.45 },
                                { 5.0 / 24.0, 0.25 } };

    Built state = build( { evenRim( 24, 0.25, true, true ),
                           evenRim( 24, 0.75, false, true ),
                           touching } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "a hole touching a rim without crossing it builds" );

    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );

    // The hole is carved, not ignored: the band less the 1/24-wide, 0.2-tall
    // patch it takes out of it. Both are exact on the 24-gon.
    const double removed =
      2.0 * TUBE_RADIUS * std::sin( TWO_PI / 48.0 ) * TUBE_HEIGHT * 0.20;

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
    // sum to zero, and a triangle containing the origin cannot.
    //
    // THE HOLE IS INSIDE THE RIM, not outside it, and that is load-bearing
    // rather than incidental. Spelled the other way round - the wrapping rim
    // at the LOW v and the hole above it - the hole lands outside the rim's
    // circle in the annulus layout, the kept set is two parity islands, and
    // the CONNECTIVITY refusal fires first; the case would still be refused
    // and would pin nothing about this reading. Here the hole sits inside the
    // rim, the kept region is one component, and the only thing left to refuse
    // it is the well-posedness sum. Confirmed by red proof: disabling the
    // well-posedness reading alone makes this case build.
    Built state = build( { evenRim( 24, 0.75, true, true ),
                           smallRing( 0.5, 0.25, 0.05 ) } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "one wrapping rim with no second rim to close the strip is "
           "refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a bound clear of every other bound is refused ===\n" );

  {
    // THE THIRD POST-HOC READING, and the hole the other two leave. A ring at
    // v = 4.0 is nowhere near the band between the rims, so it intersects
    // nothing and `TryResolve` adds no vertex to resolve anything - the
    // vertex-growth check sees a clean input. `eraseOuterTrianglesAndHoles`
    // decides interior by parity, a closed loop in the chart has an interior,
    // and the face would ship with that loop's area welded onto it: a
    // WRONG-BUT-DEFINED triangulation reached without growing the vertex set.
    //
    // This file used to CHARACTERISE that outcome - "its own island ... adds
    // area rather than inverting what was kept" - which is a defect written
    // down as a feature. Codex 4051389899 on bldrs-ai/conway-geom#207 read it
    // that way and was right to. The reading that refuses it is the lift's
    // own: the breadth-first walk is defined over one connected component, a
    // stray bound is a second, and a face is one connected region by
    // definition.
    Built banded = build( { evenRim( 24, 0.25, true, true ),
                            evenRim( 24, 0.75, false, true ) } );

    const Outcome bandedOutcome = run( banded, surface );

    check( bandedOutcome.built, "the band alone still builds" );

    Built strayed = build( { evenRim( 24, 0.25, true, true ),
                             evenRim( 24, 0.75, false, true ),
                             smallRing( 0.5, 4.0, 0.05 ) } );

    const size_t boundary = strayed.mesh.vertices.size();
    const Outcome strayedOutcome = run( strayed, surface );

    check( !strayedOutcome.built,
           "a bound clear of the band is refused, not welded on as extra "
           "area" );
    check( strayedOutcome.vertexGrowth == 0, "and the mesh is untouched" );
    check( boundary == strayed.mesh.vertices.size(),
           "leaving the caller exactly the boundary it handed over, to "
           "ear-clip" );
  }

  printf( "=== a stray bound is not caught by the vertex-growth reading ===\n" );

  {
    // WHY THE COMPONENT COUNT IS NOT REDUNDANT. The same stray-bound input,
    // asserted against the OTHER reading: it must be the connectivity refusal
    // that fires and not the vertex-growth one, or the new check is pinning
    // nothing. A ring at v = 4.0 crosses no other constraint, so CDT resolves
    // nothing and the vertex set comes back exactly as it went in.
    //
    // Read by asking the same question of a case that DOES intersect: a hole
    // reaching out through a rim grows the vertex set and is refused by the
    // older reading. Both refuse; only one of them can refuse the stray
    // bound, and this pins which.
    RingSpec reaching;

    reaching.explicitPoints = { { 0.40, 0.50 }, { 0.60, 0.50 },
                                { 0.60, 0.10 }, { 0.40, 0.10 },
                                { 0.40, 0.50 } };

    Built crossing = build( { evenRim( 24, 0.25, true, true ),
                              evenRim( 24, 0.75, false, true ),
                              reaching } );

    const Outcome crossingOutcome = run( crossing, surface );

    check( !crossingOutcome.built,
           "a hole reaching out through a rim is refused, as it was before" );
    check( crossingOutcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== the cut is reported, and it is a pair of border edges ===\n" );

  {
    // A coarse band: six samples a rim, so the cut across it is long enough
    // for a refinement to have something to do.
    Built state = build( { evenRim( 6, BULGE_RIM_LOW, true, true ),
                           evenRim( 6, BULGE_RIM_HIGH, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "the coarse band is triangulated" );
    check( !outcome.seamEdges.empty(),
           "and the chart says where the cut runs" );

    for ( const std::array< uint32_t, 3 >& triangle : outcome.triangles ) {
      state.mesh.makeTriangle( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );
    }

    // THE PROPERTY THE FINDING IS ABOUT. Each side of the cut carries one
    // triangle, so `edge.border()` is true of both - and `tesselate` skips
    // every border edge, which is why the cut cannot refine by itself.
    size_t borderBoth  = 0;
    size_t bitwiseCopy = 0;

    for ( const std::array< uint32_t, 4 >& seam : outcome.seamEdges ) {

      const std::optional< uint32_t > first =
        state.mesh.getEdge( seam[ 0 ], seam[ 1 ] );
      const std::optional< uint32_t > second =
        state.mesh.getEdge( seam[ 2 ], seam[ 3 ] );

      if ( first.has_value() && second.has_value() &&
           state.mesh.edges[ first.value() ].border() &&
           state.mesh.edges[ second.value() ].border() ) {
        ++borderBoth;
      }

      // The partner carries the SAME 3D point - bitwise - a whole number of
      // periods away in u, which is what welds the cut in Geometry::Reify.
      const bool positionsMatch =
        state.mesh.vertices[ seam[ 0 ] ].point ==
          state.mesh.vertices[ seam[ 2 ] ].point &&
        state.mesh.vertices[ seam[ 1 ] ].point ==
          state.mesh.vertices[ seam[ 3 ] ].point;

      const double firstShift =
        state.mesh.vertices[ seam[ 2 ] ].uv.x -
        state.mesh.vertices[ seam[ 0 ] ].uv.x;

      if ( positionsMatch &&
           std::abs( firstShift -
                     ( outcome.period *
                       std::round( firstShift / outcome.period ) ) ) < 1e-12 &&
           firstShift != 0.0 ) {
        ++bitwiseCopy;
      }
    }

    check( borderBoth == outcome.seamEdges.size(),
           "every edge of the cut is a border edge on both sides, which is "
           "what tesselate skips" );
    check( bitwiseCopy == outcome.seamEdges.size(),
           "and the two sides are bitwise the same point, a whole period "
           "apart in u" );
  }

  printf( "=== the cut refines on both sides at once ===\n" );

  {
    Built state = build( { evenRim( 6, BULGE_RIM_LOW, true, true ),
                           evenRim( 6, BULGE_RIM_HIGH, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built, "the coarse band is triangulated" );

    for ( const std::array< uint32_t, 3 >& triangle : outcome.triangles ) {
      state.mesh.makeTriangle( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );
    }

    const std::map< std::array< double, 6 >, size_t > before =
      borderEdgesByPosition( state.mesh );

    size_t trimBefore = 0;
    size_t cutBefore  = 0;

    for ( const auto& entry : before ) {
      ( entry.second == 1 ? trimBefore : cutBefore ) += 1;
    }

    const size_t trianglesBefore = state.mesh.triangles.size();
    const size_t verticesBefore  = state.mesh.vertices.size();

    // A floor far below the bulge, so the cut is refined to it rather than
    // stopping on the first reading.
    constexpr double FLOOR = 1e-6;

    conway::geometry::refineSeamPairs(
      state.mesh,
      []( const glm::dvec3&, const glm::dvec2& uv ) {
        return bulgedTubePoint( uv.x, uv.y );
      },
      outcome.seamEdges,
      static_cast< int32_t >( state.mesh.triangles.size() * 32 ),
      FLOOR );

    const size_t splits =
      ( state.mesh.triangles.size() - trianglesBefore ) / 2;

    check( splits > 0, "the cut is split" );

    // Each split is TWO triangles and TWO vertices - one of each per side -
    // which is the lockstep stated as arithmetic.
    check( state.mesh.triangles.size() - trianglesBefore == splits * 2 &&
           state.mesh.vertices.size() - verticesBefore == splits * 2,
           "adding one triangle and one vertex to each side, never to one" );

    size_t pairedNewVertices = 0;

    for ( size_t at = verticesBefore; at < state.mesh.vertices.size(); ++at ) {
      for ( size_t other = verticesBefore;
            other < state.mesh.vertices.size();
            ++other ) {

        if ( other != at &&
             state.mesh.vertices[ at ].point ==
               state.mesh.vertices[ other ].point ) {

          ++pairedNewVertices;
          break;
        }
      }
    }

    check( pairedNewVertices == splits * 2,
           "and every point it adds is bitwise the same on both sides, so "
           "Reify welds the refined cut as it welded the coarse one" );

    // THE CUT IS STILL CLOSED. Every border edge of the cut still has exactly
    // one partner carrying the same two positions; the trim boundary, which
    // has no partner and must not be touched, is unchanged in count.
    const std::map< std::array< double, 6 >, size_t > after =
      borderEdgesByPosition( state.mesh );

    size_t trimAfter   = 0;
    size_t cutAfter    = 0;
    size_t unpairedCut = 0;

    for ( const auto& entry : after ) {

      if ( entry.second == 1 ) {
        ++trimAfter;
      } else if ( entry.second == 2 ) {
        ++cutAfter;
      } else {
        ++unpairedCut;
      }
    }

    check( unpairedCut == 0,
           "no border edge of the refined mesh is anything but a trim "
           "segment or one half of a cut pair" );
    check( cutAfter == cutBefore + splits,
           "the cut gains exactly one paired edge per split" );
    check( trimAfter == trimBefore,
           "and the trim boundary, which is shared with the neighbouring "
           "face, is not touched" );

    // THE SPLIT POINT IS ON THE SURFACE, not on the chord it replaces. That
    // is the opposite of the layout split above, and for the opposite
    // reason: the cut is interior to this face, so no neighbour holds the
    // other half of it and there is nothing to crack against.
    size_t onSurface = 0;
    size_t offChord  = 0;

    for ( size_t at = verticesBefore; at < state.mesh.vertices.size(); ++at ) {

      const ParameterVertex& added = state.mesh.vertices[ at ];

      if ( glm::distance( added.point,
                          bulgedTubePoint( added.uv.x, added.uv.y ) ) < 1e-12 ) {
        ++onSurface;
      }
    }

    for ( const std::array< uint32_t, 4 >& seam : outcome.seamEdges ) {

      const glm::dvec3 chordMid =
        ( state.mesh.vertices[ seam[ 0 ] ].point +
          state.mesh.vertices[ seam[ 1 ] ].point ) * 0.5;

      for ( size_t at = verticesBefore;
            at < state.mesh.vertices.size();
            ++at ) {

        if ( glm::distance( state.mesh.vertices[ at ].point, chordMid ) >
               1e-9 ) {
          continue;
        }

        ++offChord;
      }
    }

    check( onSurface == splits * 2,
           "every point the refinement adds lies on the surface" );
    check( offChord == 0,
           "and none of them is the chord midpoint the coarse cut ran "
           "through" );
  }

  printf( "=== a mesh with no cut is not touched ===\n" );

  {
    Built state = build( { evenRim( 6, BULGE_RIM_LOW, true, true ),
                           evenRim( 6, BULGE_RIM_HIGH, false, true ) } );

    const Outcome outcome = run( state, surface );

    for ( const std::array< uint32_t, 3 >& triangle : outcome.triangles ) {
      state.mesh.makeTriangle( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );
    }

    const size_t triangles = state.mesh.triangles.size();
    const size_t vertices  = state.mesh.vertices.size();

    conway::geometry::refineSeamPairs(
      state.mesh,
      []( const glm::dvec3&, const glm::dvec2& uv ) {
        return bulgedTubePoint( uv.x, uv.y );
      },
      {},
      static_cast< int32_t >( state.mesh.triangles.size() * 32 ),
      1e-6 );

    check( state.mesh.triangles.size() == triangles &&
           state.mesh.vertices.size() == vertices,
           "an empty seam list leaves the mesh exactly as it was" );
  }

  printf( "=== a boundary that does not wrap the closure is triangulated too ===\n" );

  {
    // THE READING THAT USED TO STAND HERE, AND WHY IT IS GONE. This face's
    // rings net no winding, and the old code sent it to earcut on that
    // reading alone. codex 4051389902 on bldrs-ai/conway-geom#207 showed the
    // reading cannot be made sound: the lifted steps +0.60P, -0.10P, +0.49P,
    // +0.01P sum to one turn, and nearest-image reduction turns the first
    // into -0.40P and the sum into zero, so a rim that DOES wrap reads as one
    // that does not and is ear-clipped - the sliver this path exists to
    // remove. No threshold closes it either: +0.80P, +0.10P, +0.05P, +0.05P
    // also sum to one turn and every one of them reduces below 0.20P.
    //
    // So nothing chooses a triangulator from the winding any more, and this
    // case - an ordinary outer ring with a hole in it, on a surface that
    // happens to be closed in u - is triangulated here rather than refused.
    // The property that says the winding is not being read is the CUT: a
    // boundary that does not wrap has no monodromy, so there is nothing to
    // duplicate and `seamEdges` comes back empty.
    RingSpec outer;

    outer.explicitPoints = { { 0.20, 0.30 }, { 0.45, 0.30 },
                             { 0.45, 0.70 }, { 0.20, 0.70 },
                             { 0.20, 0.30 } };

    RingSpec hole;

    hole.explicitPoints = { { 0.28, 0.45 }, { 0.28, 0.55 },
                            { 0.37, 0.55 }, { 0.37, 0.45 },
                            { 0.28, 0.45 } };

    Built state = build( { outer, hole } );

    const Outcome outcome = run( state, surface );

    check( outcome.built,
           "an outer ring with a hole, netting no winding, is triangulated "
           "rather than handed to earcut on a reading of that winding" );
    check( outcome.seamEdges.empty(),
           "and it needs no cut, because a boundary that does not wrap has "
           "no monodromy to cut" );
    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
    check( emittedArea( state, outcome ) > 0.0,
           "and it emits area" );
  }

  printf( "=== a step past the decisive bound is refused, not guessed at ===\n" );

  {
    // codex's own counterexample, spelled as the solve would return it - every
    // u wrapped into [0, 1). The true steps are +0.60, -0.10, +0.49, +0.01 of
    // the period and the rim wraps once; reduced to nearest images they read
    // -0.40, -0.10, +0.49, +0.01 and net zero.
    //
    // The 0.49 and 0.40 are past MAX_DECISIVE_IMAGE_STEP_FRACTION, so this
    // face is refused HERE - at the reading every later one depends on -
    // rather than laid out from a step whose direction the reduction has
    // reversed. The caller ear-clips, which is what it does today; what is
    // new is that the construction says it cannot read this boundary instead
    // of proceeding as if it could.
    RingSpec ambiguous;

    ambiguous.explicitPoints = { { 0.00, 0.70 }, { 0.60, 0.80 },
                                 { 0.50, 0.70 }, { 0.99, 0.80 } };

    Built state = build( { evenRim( 24, 0.25, true, true ), ambiguous } );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a rim whose steps reduce past the decisive bound is refused" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
    check( boundary == state.mesh.vertices.size(),
           "leaving the caller exactly the boundary it handed over" );
  }

  printf( "=== a wide step is refused even when it reads decisively ===\n" );

  {
    // THE CASE THE BOUND IS FOR, and the one that shows it is load-bearing.
    // This rim steps 0.40 of a period three times over - a shape the layout
    // would lay out and the CDT would accept - and 0.40 is exactly what
    // codex's +0.60 reduces to. Nothing recoverable from u separates the two
    // spellings, so both are refused, and the caller ear-clips as it does
    // today. The bound is conservative on purpose: it costs a coarse rim that
    // happened to be spelled correctly, and it buys never laying out a
    // reversed one.
    RingSpec wide;

    wide.explicitPoints = { { 0.00, 0.75 }, { 0.40, 0.75 }, { 0.80, 0.75 } };

    Built state = build( { evenRim( 24, 0.25, true, true ), wide } );

    const Outcome outcome = run( state, surface );

    check( !outcome.built,
           "a rim stepping 0.40 of a period is refused, because that is what "
           "a reversed 0.60 reduces to" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
  }

  printf( "=== a step that reduced the wrong way is refused ===\n" );

  {
    // CODEX 4051991542 ON bldrs-ai/conway-geom#207 AND bldrs-ai/conway#711:
    // the decisive-step bound is necessary and not sufficient, because a LARGE
    // true step can reduce INTO the decisive band. A rim stepping +0.80 of a
    // period reads as -0.20, every later reading takes that -0.20 as fact, and
    // the three post-hoc refusals on the triangulation's own output cannot
    // recover a route that was discarded before the triangulation was built.
    //
    // THE FACE THIS IS SPELLED ON is the prism's side wall - two rims, each
    // sampled twelve times, which read correctly - with a SLOT cut along its
    // straight bottom face. The slot spans u from 0.05 to 0.75, and because
    // the face it lies on is flat its two long edges are represented exactly
    // by their endpoints: a tessellator has no reason to put a point in the
    // middle of a straight line. So the polyline is not coarse, it is the
    // slot, and its long edge steps +0.70 of a period.
    //
    // MEASURED AGAINST THE HEADER BEFORE THE ROUTE CONFIRMATION, this face
    // BUILDS: 49 triangles and 476.29 of area, with CDT adding no vertex
    // (nothing intersects), the kept triangles one connected component (the
    // misread slot is a box inside the band, not a second island) and every
    // triangle's offsets summing to zero. All three post-hoc refusals pass on
    // a face whose slot covers 0.30 of the period where the boundary it was
    // handed describes one covering 0.70 - 0.40 of a period of slot paved
    // over. That is codex's claim, and it is right.
    //
    // WHAT THE ROUTE READING SEES. The +0.70 image runs along the flat face,
    // so its surface image IS the chord: length 17.50, which is the chord's
    // own length and the floor no route can go under. The -0.30 image leaves
    // the bottom face for the right, top and left ones and measures 29.93. The
    // reduction is not the shortest route, so it is not the route the chord
    // stands for, and the face is refused. The caller ear-clips, as it does
    // today.
    const tinynurbs::RationalSurface3d prism = makeFlatSidedPrism();

    RingSpec slot;

    slot.explicitPoints = { { 0.05, 0.45 }, { 0.75, 0.45 },
                            { 0.75, 0.55 }, { 0.05, 0.55 } };

    Built state = build( { evenRim( 12, 0.20, true, true ),
                           evenRim( 12, 0.80, false, true ), slot },
                         prismPoint );

    const size_t boundary = state.mesh.vertices.size();
    const Outcome outcome = run( state, prism );

    check( !outcome.built,
           "a segment whose nearest image is not the route its own chord "
           "draws is refused, though the image reads decisively" );
    check( outcome.vertexGrowth == 0, "and the mesh is untouched" );
    check( boundary == state.mesh.vertices.size(),
           "leaving the caller exactly the boundary it handed over" );
  }

  printf( "=== the same face, sampled so the reduction IS the route ===\n" );

  {
    // THE CONTROL THE REFUSAL ABOVE NEEDS, or it would be indistinguishable
    // from refusing the prism outright. Same surface, same band, same slot -
    // only the slot's two long edges carry one more sample each, so no step
    // exceeds 0.35 of a period and every one of them IS the shortest route.
    // Measured: each half of the slot edge reads 8.75, its own chord exactly,
    // against 34.77 for the nearest rival image. Nothing about the surface
    // changed; what changed is that the polyline no longer needs a route the
    // reduction cannot reach.
    const tinynurbs::RationalSurface3d prism = makeFlatSidedPrism();

    RingSpec slot;

    slot.explicitPoints = { { 0.05, 0.45 }, { 0.40, 0.45 }, { 0.75, 0.45 },
                            { 0.75, 0.55 }, { 0.40, 0.55 }, { 0.05, 0.55 } };

    Built state = build( { evenRim( 12, 0.20, true, true ),
                           evenRim( 12, 0.80, false, true ), slot },
                         prismPoint );

    const Outcome outcome = run( state, prism );

    check( outcome.built,
           "the same band, sampled so every step's nearest image is the route "
           "its chord draws, is triangulated" );
    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
    check( emittedArea( state, outcome ) > 0.0, "and it emits area" );
  }

  printf( "=== the coarsest legitimate rim still builds ===\n" );

  {
    // THE FLOOR THE BOUND HAS TO CLEAR. Three samples is the fewest a closed
    // rim can have and its steps are exactly a third of a period, an ulp
    // either side. A bound at P/3 refuses this; 3/8 does not.
    Built state = build( { evenRim( 3, 0.25, true, true ),
                           evenRim( 3, 0.75, false, true ) } );

    const Outcome outcome = run( state, surface );

    check( outcome.built,
           "a three-sample rim, whose every step is a third of a period, is "
           "read rather than refused" );
    check( worstLiftedSpan( state, outcome ) < 0.5,
           "with no triangle spanning half a period" );
  }

  printf( failures == 0 ? "PASS\n" : "FAIL (%d)\n", failures );

  return failures == 0 ? 0 : 1;
}
