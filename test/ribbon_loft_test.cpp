/*
 * Soundness tests for tryRibbonLoft's boundary decomposition in mesh_utils.h.
 *
 * The three properties pinned here are the three that a deviation metric
 * CANNOT see, which is why they need a test rather than a corpus digest:
 *
 *   1. Plateau-tolerant reversal detection. A turning point sitting on a
 *      plateau (dv == 0) is invisible to the obvious `previous * next < 0`
 *      test, and on every Orbiter thread flank the loop's start vertex is
 *      exactly such a plateau. With the strict test the gate sees ONE reversal
 *      where there are two and rejects a textbook ribbon - the faces this
 *      whole path exists for silently keep their old behaviour, and nothing
 *      downstream reports anything wrong.
 *
 *   2. Full boundary coverage. A mesh built over PART of a face still sits
 *      close to the surface, so deflection cannot distinguish it from a mesh
 *      built over the whole face. An earlier prototype selected the two
 *      longest monotone runs and silently discarded up to 8% of the boundary
 *      while reporting deviation inside target. Coverage has to be a count.
 *
 *   3. The gate is narrow. A boundary with more than two v-monotone runs must
 *      be rejected outright rather than approximated, so it ear-clips exactly
 *      as it does today.
 *
 * Each test fails if the corresponding guard is reverted; that was verified by
 * reverting, not by reading. See bldrs-ai/conway#608.
 *
 * Standalone by design: it includes mesh_utils.h directly and links nothing,
 * matching nurbs_seam_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

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

constexpr double PI = 3.14159265358979323846;

using conway::geometry::ParameterVertex;
using conway::geometry::WingedEdgeMesh;
using conway::geometry::uvMonotoneBreaks;

/**
 * A helical ribbon: v runs along the coil, u across its width. The same shape
 * as an Orbiter thread flank, at a size a test can hold.
 */
glm::dvec3 ribbonPoint( double u, double v ) {

  constexpr double RADIUS = 10.0;
  constexpr double PITCH  = 3.0;

  return glm::dvec3(
    ( RADIUS + u ) * std::cos( 2.0 * PI * v ),
    ( RADIUS + u ) * std::sin( 2.0 * PI * v ),
    PITCH * v );
}

ParameterVertex ribbonVertex( double u, double v ) {

  return ParameterVertex{ ribbonPoint( u, v ), glm::dvec2( u, v ) };
}

/** The degree-1 surface the ribbon lies on, for the inverse solve. */
tinynurbs::RationalSurface3d flatSurface() {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  std::vector< glm::dvec3 > controlPoints = {
    ribbonPoint( 0.0, 0.0 ), ribbonPoint( 0.0, 1.0 ),
    ribbonPoint( 1.0, 0.0 ), ribbonPoint( 1.0, 1.0 ) };

  surface.control_points = tinynurbs::array2( 2, 2, controlPoints );
  surface.weights        = tinynurbs::array2( 2, 2,
    std::vector< double >{ 1.0, 1.0, 1.0, 1.0 } );

  return surface;
}

/**
 * A closed ribbon boundary: up one long edge at u = 0, across, back down the
 * other at u = uWidth, and across again.
 *
 * `plateauAtWrap` repeats the turning point's v on the vertex either side of
 * index 0, which is what the real trim curves do and what defeats the strict
 * reversal test.
 */
WingedEdgeMesh< ParameterVertex > buildRibbon(
    size_t steps,
    double uWidth,
    bool   plateauAtWrap ) {

  WingedEdgeMesh< ParameterVertex > mesh;

  // Index 0 is the TOP turning point, so the wrap is where a plateau would
  // sit - the arrangement measured on the Orbiter faces. `plateauAtWrap`
  // duplicates that apex, which is enough to defeat the strict reversal test.
  if ( plateauAtWrap ) {
    mesh.makeVertex( ribbonVertex( 0.0, 1.0 ) );
  }

  // Down one side, v strictly decreasing, ending ON the bottom turning point.
  for ( size_t at = 0; at <= steps; ++at ) {

    const double v =
      1.0 - ( static_cast< double >( at ) / static_cast< double >( steps ) );

    mesh.makeVertex( ribbonVertex( 0.0, v ) );
  }

  // Back up the other side, STRICTLY BETWEEN the two turning points: the
  // apexes are single vertices shared by both legs. Repeating them here would
  // put a plateau INSIDE a leg, where it means the leg runs flat and the
  // region is no longer an interval - which the loft rejects, by design.
  // Offset by half a step so the two legs' v samples INTERLEAVE rather than
  // coincide. That is the realistic case - the two trim curves are tessellated
  // independently - and it is what makes the rung count approach the boundary
  // count, which is the regime where a constant column cap overruns the
  // amplification budget.
  for ( size_t at = 1; at < steps; ++at ) {

    const double v =
      ( static_cast< double >( at ) - 0.5 ) / static_cast< double >( steps );

    mesh.makeVertex( ribbonVertex( uWidth, v ) );
  }

  return mesh;
}

/**
 * The plateau at the wrap must not hide a turning point.
 *
 * With the strict `previous * next < 0` test this reports one break and the
 * ribbon is rejected.
 */
void testPlateauAtWrapIsNotAReversal() {

  printf( "\n== plateau at the wrap ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 24, 1.0, true );

  const std::vector< size_t > breaks =
    uvMonotoneBreaks( mesh, mesh.vertices.size() );

  printf( "  breaks found: %zu\n", breaks.size() );

  check( breaks.size() == 2,
         "a ribbon whose turning point sits on a plateau reports two breaks" );

  // And the strict test really would have missed it - stated as arithmetic so
  // the test does not merely assert its own implementation.
  const size_t count = mesh.vertices.size();

  size_t strict = 0;

  for ( size_t at = 0; at < count; ++at ) {

    const double previous =
      mesh.vertices[ at ].uv.y - mesh.vertices[ ( at + count - 1 ) % count ].uv.y;
    const double next =
      mesh.vertices[ ( at + 1 ) % count ].uv.y - mesh.vertices[ at ].uv.y;

    if ( previous * next < 0.0 ) {
      ++strict;
    }
  }

  printf( "  the strict previous*next<0 test would find: %zu\n", strict );

  check( strict < 2,
         "the strict test is genuinely blind here, so this test can fail" );
}

/** Without the plateau the decomposition still finds exactly two breaks. */
void testPlainRibbonDecomposes() {

  printf( "\n== plain ribbon ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 24, 1.0, false );

  const std::vector< size_t > breaks =
    uvMonotoneBreaks( mesh, mesh.vertices.size() );

  printf( "  breaks found: %zu\n", breaks.size() );

  check( breaks.size() == 2, "a plain ribbon reports two breaks" );
}

/**
 * The loft consumes every boundary vertex, and reuses both trim polylines
 * verbatim.
 *
 * Coverage is asserted as a COUNT because deflection cannot see a truncated
 * boundary - the failure an earlier prototype shipped undetected.
 */
void testLoftCoversAndReusesTheBoundary() {

  printf( "\n== coverage and verbatim boundary ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 24, 1.0, false );

  const size_t boundaryCount = mesh.vertices.size();

  // Keep a copy of the boundary to compare against after the loft appends.
  const std::vector< ParameterVertex > before(
    mesh.vertices.begin(), mesh.vertices.end() );

  conway::geometry::RationalNurbsInverseMethod* noSolve = nullptr;

  (void)noSolve;

  // The loft needs a solve only to evaluate interior columns. Build the
  // surface the ribbon lies on so the evaluator is the real one.
  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  std::vector< glm::dvec3 > controlPoints = {
    ribbonPoint( 0.0, 0.0 ), ribbonPoint( 0.0, 1.0 ),
    ribbonPoint( 1.0, 0.0 ), ribbonPoint( 1.0, 1.0 ) };

  surface.control_points = tinynurbs::array2( 2, 2, controlPoints );
  surface.weights        = tinynurbs::array2( 2, 2,
    std::vector< double >{ 1.0, 1.0, 1.0, 1.0 } );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  check( built, "the ribbon lofts" );

  if ( !built ) {
    return;
  }

  printf( "  triangles: %zu, vertices %zu -> %zu\n",
          mesh.triangles.size(), boundaryCount, mesh.vertices.size() );

  // Every original boundary vertex is still there, unmoved: the loft appends
  // interior vertices and never rewrites the ones it was handed.
  size_t moved = 0;

  for ( size_t at = 0; at < boundaryCount; ++at ) {

    if ( !( mesh.vertices[ at ].point == before[ at ].point ) ) {
      ++moved;
    }
  }

  check( moved == 0, "no boundary vertex was moved by the loft" );

  // Every boundary vertex is referenced by the emitted mesh. This is the count
  // that a deflection metric cannot replace.
  std::vector< bool > used( boundaryCount, false );

  for ( const conway::geometry::ConnectedTriangle& triangle : mesh.triangles ) {
    for ( uint32_t corner : triangle.vertices ) {
      if ( corner < boundaryCount ) {
        used[ corner ] = true;
      }
    }
  }

  const size_t consumed =
    static_cast< size_t >( std::count( used.begin(), used.end(), true ) );

  printf( "  boundary vertices consumed: %zu of %zu\n",
          consumed, boundaryCount );

  check( consumed == boundaryCount,
         "every boundary vertex is used by the emitted mesh" );
}

/**
 * A boundary with more than two v-monotone runs is rejected, not approximated.
 *
 * This is the notched case, and the whole point of the gate: it must fall
 * through to ear-clipping rather than silently lofting over the notch.
 */
void testNotchedBoundaryIsRejected() {

  printf( "\n== notched boundary is rejected ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh;

  // Up, back down a little, up again, then all the way down: four runs.
  const double profile[] = {
    0.0, 0.2, 0.4, 0.3, 0.2, 0.5, 0.8, 1.0, 0.7, 0.4, 0.1 };

  for ( double v : profile ) {
    mesh.makeVertex( ribbonVertex( 0.0, v ) );
  }

  for ( double v : profile ) {
    mesh.makeVertex( ribbonVertex( 1.0, 1.0 - v ) );
  }

  const std::vector< size_t > breaks =
    uvMonotoneBreaks( mesh, mesh.vertices.size() );

  printf( "  breaks found: %zu\n", breaks.size() );

  check( breaks.size() > 2, "the notched boundary really has more than two runs" );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  std::vector< glm::dvec3 > controlPoints = {
    ribbonPoint( 0.0, 0.0 ), ribbonPoint( 0.0, 1.0 ),
    ribbonPoint( 1.0, 0.0 ), ribbonPoint( 1.0, 1.0 ) };

  surface.control_points = tinynurbs::array2( 2, 2, controlPoints );
  surface.weights        = tinynurbs::array2( 2, 2,
    std::vector< double >{ 1.0, 1.0, 1.0, 1.0 } );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  const size_t before = mesh.triangles.size();

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, mesh.vertices.size() );

  check( !built, "a notched boundary is rejected" );
  check( mesh.triangles.size() == before,
         "and the mesh is left untouched for the ear-clipper" );
}

/**
 * A CONCAVE two-leg monotone ribbon still tiles its chart exactly once.
 *
 * Two v-direction sign changes prove each leg is a graph over v and nothing
 * more - the legs may wander arbitrarily in u between their samples. The first
 * implementation paired legA[i] with legB[j] at DIFFERENT v, and such a
 * slanted rung can leave the trim while consecutive slanted rungs cross, so
 * the quads between them overlap. Measured on the corpus at the time: 7 of 284
 * lofted faces covered part of their chart up to 1.045x over, and a planar
 * face reports zero deflection so nothing downstream rejected them.
 *
 * The oracle is uv area: for a tiling of a simple polygon the sum of the
 * triangles' |area| equals the boundary's shoelace area. Anything above it is
 * overlap. Deflection cannot see this - every vertex is on the surface.
 *
 * Rejecting concave legs was measured to be the wrong fix: all seven of the
 * large thread faces the path exists for are concave in u (up to 702 direction
 * changes on one leg), while 5 of 84 u-monotone faces overlapped anyway.
 */
void testConcaveRibbonTilesExactlyOnce() {

  printf( "\n== concave ribbon tiles exactly once ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh;

  // Leg A is just the two turning points, so the staircase has nothing to
  // advance on that side and fans from one vertex across the whole face. Leg B
  // wanders hard in u - out to 3.0, back to 0.1, repeatedly - while staying
  // strictly monotone in v. A fan tiles a polygon only if it is star-shaped
  // from the fan's apex, and this one is not: the triangle spanning from the
  // bottom apex to a far-out leg-B vertex crosses the narrow waist where the
  // face is only 0.1 wide, so it lies far outside the trim.
  mesh.makeVertex( ribbonVertex( 0.0, 0.0 ) );
  mesh.makeVertex( ribbonVertex( 0.0, 1.0 ) );

  const double legBU[] = { 3.0, 0.1, 3.0, 0.1, 3.0, 0.1 };
  const double legBV[] = { 0.85, 0.7, 0.55, 0.4, 0.25, 0.1 };

  for ( size_t at = 0; at < sizeof( legBU ) / sizeof( legBU[ 0 ] ); ++at ) {
    mesh.makeVertex( ribbonVertex( legBU[ at ], legBV[ at ] ) );
  }

  const size_t boundaryCount = mesh.vertices.size();

  const std::vector< size_t > breaks = uvMonotoneBreaks( mesh, boundaryCount );

  printf( "  boundary %zu, breaks %zu\n", boundaryCount, breaks.size() );

  check( breaks.size() == 2,
         "the concave ribbon is still v-monotone with two breaks" );

  // Shoelace of the uv boundary, before the loft appends anything.
  double shoelace = 0.0;

  for ( size_t at = 0; at < boundaryCount; ++at ) {

    const glm::dvec2& p = mesh.vertices[ at ].uv;
    const glm::dvec2& q = mesh.vertices[ ( at + 1 ) % boundaryCount ].uv;

    shoelace += ( p.x * q.y ) - ( q.x * p.y );
  }

  shoelace = std::fabs( shoelace ) * 0.5;

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  check( built, "the concave ribbon lofts" );

  if ( !built ) {
    return;
  }

  double absolute = 0.0;

  for ( const conway::geometry::ConnectedTriangle& triangle : mesh.triangles ) {

    const glm::dvec2& a = mesh.vertices[ triangle.vertices[ 0 ] ].uv;
    const glm::dvec2& b = mesh.vertices[ triangle.vertices[ 1 ] ].uv;
    const glm::dvec2& c = mesh.vertices[ triangle.vertices[ 2 ] ].uv;

    absolute += std::fabs(
      0.5 * ( ( ( b.x - a.x ) * ( c.y - a.y ) ) -
              ( ( c.x - a.x ) * ( b.y - a.y ) ) ) );
  }

  printf( "  shoelace %.9f, emitted |area| %.9f, coverage %.9f\n",
          shoelace, absolute, absolute / shoelace );

  check( std::fabs( ( absolute / shoelace ) - 1.0 ) < 1e-9,
         "the emitted triangles cover the chart exactly once - no overlap" );
}

/**
 * The sweep emits EXACTLY `boundary - 2` triangles, which is what keeps the
 * amplification budget honest.
 *
 * Any triangulation of a simple n-gon has exactly n - 2 triangles, so a sweep
 * that emits that many has triangulated the polygon and nothing else - it has
 * not fanned, doubled back, or dropped an ear. That is also precisely earcut's
 * seed size on the same boundary, so the face reaches `tesselate` with the
 * budget it would have had anyway and MAX_TRIANGLE_AMPLIFACTION stays the real
 * limit.
 *
 * Asserted as an equality rather than as "within budget": under the sweep the
 * inequality cannot fail - n - 2 is always under 32(n - 2) - and an assertion
 * that cannot fail is not a test. An earlier version of this test made exactly
 * that mistake once the column cap it was written for stopped existing.
 */
void testSweepEmitsExactlyTheEarcutSeed() {

  printf( "\n== sweep emits exactly the earcut seed ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 4, 1.0, false );

  const size_t boundaryCount = mesh.vertices.size();

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  // A target of zero is unreachable, so the column search runs to its ceiling
  // - which is exactly the case the constant cap got wrong.
  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  check( built, "the small ribbon lofts" );

  if ( !built ) {
    return;
  }

  const size_t budget =
    static_cast< size_t >( conway::geometry::MAX_TRIANGLE_AMPLIFACTION ) *
    ( boundaryCount - 2 );

  printf( "  boundary %zu, emitted %zu triangles, expected %zu, budget %zu\n",
          boundaryCount, mesh.triangles.size(), boundaryCount - 2, budget );

  check( mesh.triangles.size() == boundaryCount - 2,
         "the sweep emits exactly boundary - 2 triangles, earcut's seed size" );

  check( mesh.triangles.size() <= budget,
         "which is inside MAX_TRIANGLE_AMPLIFACTION x that seed" );
}

/**
 * The seed introduces NO new boundary vertices, so it cannot create a
 * T-vertex against the neighbouring face.
 *
 * This is the topology finding, pinned directly. The construction this
 * replaced joined the two legs at one v, which forced an interpolated vertex
 * onto whichever leg had no sample there. The point is collinear, so the
 * solid stayed geometrically gap-free and every geometric check passed - but
 * Reify's weld does not subdivide the neighbour's edge, so the loft ended up
 * with two half-edges facing one. On Orbiter that took unpaired half-edges
 * from 1,448 to 70,698 (codex on conway-geom#190).
 *
 * The oracle is a two-face fixture: the ribbon plus a neighbour built from the
 * SAME leg-B polyline, which is what sharing a trim edge means. Every half
 * edge must pair.
 */
void testNoTVerticesAgainstANeighbour() {

  printf( "\n== no T-vertices against a neighbour ==\n" );

  // Legs deliberately sampled at DIFFERENT v - the motivating 2544/2545 case.
  WingedEdgeMesh< ParameterVertex > mesh;

  mesh.makeVertex( ribbonVertex( 0.0, 0.0 ) );

  for ( size_t at = 1; at < 9; ++at ) {
    mesh.makeVertex( ribbonVertex( 0.0, static_cast< double >( at ) / 9.0 ) );
  }

  mesh.makeVertex( ribbonVertex( 0.0, 1.0 ) );

  // Leg B: seven samples where leg A has nine, so their v values interleave
  // and almost never coincide.
  std::vector< uint32_t > legB;

  for ( size_t at = 1; at < 8; ++at ) {

    const double v = 1.0 - ( static_cast< double >( at ) / 8.0 );

    legB.push_back( mesh.makeVertex( ribbonVertex( 1.0, v ) ) );
  }

  const size_t boundaryCount = mesh.vertices.size();

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  check( built, "the unevenly sampled ribbon triangulates" );

  if ( !built ) {
    return;
  }

  printf( "  vertices %zu -> %zu\n", boundaryCount, mesh.vertices.size() );

  // The T-vertex property is not "no vertices were added" - interior
  // refinement adds plenty, and must be free to. It is that EVERY EDGE OF THE
  // ORIGINAL TRIM POLYLINE survives in the mesh unsplit, so the neighbouring
  // face, which builds that same polyline, meets it half-edge for half-edge.
  size_t missing = 0;

  for ( size_t at = 0; at < boundaryCount; ++at ) {

    const uint32_t from = static_cast< uint32_t >( at );
    const uint32_t to   = static_cast< uint32_t >( ( at + 1 ) % boundaryCount );

    bool found = false;

    for ( const conway::geometry::ConnectedTriangle& triangle : mesh.triangles ) {
      for ( uint32_t corner = 0; corner < 3; ++corner ) {

        const uint32_t a = triangle.vertices[ corner ];
        const uint32_t b = triangle.vertices[ ( corner + 1 ) % 3 ];

        if ( ( a == from && b == to ) || ( b == from && a == to ) ) {
          found = true;
        }
      }
    }

    if ( !found ) { ++missing; }
  }

  printf( "  trim-polyline edges missing from the mesh: %zu of %zu\n",
          missing, boundaryCount );

  check( missing == 0,
         "every shared boundary edge survives unsplit - no T-vertex" );

  check( mesh.vertices.size() == boundaryCount,
         "the seed itself adds no vertex at all" );
}

/**
 * A ribbon whose uv sits on a far-from-origin knot domain still passes the
 * validity gate.
 *
 * A b-spline's parameters come straight from the file's knot domain, which
 * carries whatever offset the exporter used. The gate compares the emitted
 * triangles' area against the boundary's shoelace area, and a shoelace term is
 * a difference of products: on a domain like [1e4, 1e4 + 0.03] the products
 * are ~1e8 and their difference ~1e-3, so the leading digits cancel and the
 * boundary area arrives as noise. The sweep is then rejected to the
 * ear-clipper and the long chords come back (codex on conway-geom#190).
 *
 * The failure mode is why this needs a test rather than an eyeball: the gate
 * exists in order to reject, so it rejecting a good face looks exactly like it
 * working. Nothing downstream would report anything wrong - the face would
 * just quietly go back to being badly triangulated.
 */
void testFarFromOriginKnotDomain() {

  printf( "\n== ribbon on a far-from-origin knot domain ==\n" );

  // The same ribbon as the coverage test, shifted onto an offset domain.
  //
  // The offset has to be chosen against the shape, not picked round: the
  // relative error of a raw shoelace goes as offset^2 * epsilon / area, so at
  // 1e4 on this 0.03-wide ribbon it is ~2e-7 and slips under the gate's 1e-6
  // unnoticed. 1e5 puts it at ~1e-4, comfortably over. The assertion below
  // measures the error rather than assuming it, so if this ever stops being
  // the bad case the test says so instead of passing quietly.
  constexpr double OFFSET = 1.0e5;
  constexpr double WIDTH  = 0.03;

  WingedEdgeMesh< ParameterVertex > mesh;

  const size_t steps = 24;

  for ( size_t at = 0; at <= steps; ++at ) {

    const double along =
      1.0 - ( static_cast< double >( at ) / static_cast< double >( steps ) );

    mesh.makeVertex(
      ParameterVertex{ ribbonPoint( 0.0, along ),
                       glm::dvec2( OFFSET, OFFSET + along ) } );
  }

  for ( size_t at = 1; at < steps; ++at ) {

    const double along =
      static_cast< double >( at ) / static_cast< double >( steps );

    mesh.makeVertex(
      ParameterVertex{ ribbonPoint( WIDTH, along ),
                       glm::dvec2( OFFSET + WIDTH, OFFSET + along ) } );
  }

  const size_t boundaryCount = mesh.vertices.size();

  // State the hazard as arithmetic, so this test cannot pass for the wrong
  // reason: on these coordinates the raw shoelace really does lose the area.
  {
    double raw = 0.0, shifted = 0.0;

    const glm::dvec2 reference = mesh.vertices[ 0 ].uv;

    for ( size_t at = 0; at < boundaryCount; ++at ) {

      const glm::dvec2& here = mesh.vertices[ at ].uv;
      const glm::dvec2& next = mesh.vertices[ ( at + 1 ) % boundaryCount ].uv;

      raw += ( here.x * next.y ) - ( next.x * here.y );

      const glm::dvec2 hereShifted = here - reference;
      const glm::dvec2 nextShifted = next - reference;

      shifted +=
        ( hereShifted.x * nextShifted.y ) - ( nextShifted.x * hereShifted.y );
    }

    printf( "  shoelace raw %.9e, relative to a boundary point %.9e\n",
            std::fabs( raw ) * 0.5, std::fabs( shifted ) * 0.5 );

    check( std::fabs( ( std::fabs( raw ) * 0.5 ) /
                      ( std::fabs( shifted ) * 0.5 ) - 1.0 ) > 1e-6,
           "the raw shoelace really is wrong here, so this test can fail" );
  }

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  printf( "  built=%d, triangles %zu, expected %zu\n",
          (int)built, mesh.triangles.size(), boundaryCount - 2 );

  check( built, "the offset ribbon is not rejected back to the ear-clipper" );

  check( built && mesh.triangles.size() == boundaryCount - 2,
         "and it triangulates completely" );
}

/**
 * An UNEVENLY SAMPLED ribbon, in both traversal directions, is tiled exactly
 * once and stays inside its own trim loop.
 *
 * This is the case the sweep's same-chain branch exists for and the one every
 * fixture above misses. When the two legs' v samples interleave one-for-one -
 * which is what `buildRibbon` arranges, and what the corpus's clean thread
 * flanks happen to do - consecutive sweep vertices always alternate chains, so
 * the same-chain branch NEVER RUNS and its predicate is never exercised. Give
 * one leg many samples where the other has none, as a near-horizontal cap on
 * a real trim does, and the branch fires on every one of them.
 *
 * With the chain-side sign taken from the LABEL rather than from the leg's
 * traversal direction, every triangle that branch emits lands outside the
 * polygon: 112 of 112 on the corpus face with 401 boundary points, 16 of 16 on
 * the one with 1,167 (bldrs-ai/conway#665). Both directions are built here
 * because the label happens to be right for one of them - a fix that merely
 * flipped the constant would pass a one-sided test and break the other side.
 *
 * The oracle is not the validity gate's own area sum: every emitted triangle's
 * uv centroid is tested against the boundary by ray casting, which is
 * independent of both the sweep and the gate. The area equality is checked too,
 * because overlap can hide from a centroid test and vice versa.
 */
void testUnevenlySampledRibbonStaysInsideItsTrim() {

  printf( "\n== unevenly sampled ribbon, both traversal directions ==\n" );

  // `descending` puts the TOP turning point at index 0, so leg A - the arc
  // between the two breaks in loop order - runs downward in v and has to be
  // reversed for the sweep. That is the case the label-based sign gets wrong.
  for ( int descending = 1; descending >= 0; --descending ) {

    printf( "  -- leg A runs %s in v --\n",
            descending != 0 ? "downward" : "upward" );

    constexpr size_t DENSE  = 40;
    constexpr size_t SPARSE = 5;
    constexpr double WIDTH  = 1.0;

    WingedEdgeMesh< ParameterVertex > mesh;

    // The densely sampled leg, from one turning point to the other.
    for ( size_t at = 0; at <= DENSE; ++at ) {

      const double along =
        static_cast< double >( at ) / static_cast< double >( DENSE );

      mesh.makeVertex(
        ribbonVertex( 0.0, descending != 0 ? 1.0 - along : along ) );
    }

    // The sparse leg back, STRICTLY BETWEEN the turning points, which the two
    // legs share as single vertices.
    for ( size_t at = 1; at < SPARSE; ++at ) {

      const double along =
        static_cast< double >( at ) / static_cast< double >( SPARSE );

      mesh.makeVertex(
        ribbonVertex( WIDTH, descending != 0 ? along : 1.0 - along ) );
    }

    const size_t boundaryCount = mesh.vertices.size();

    // PRECONDITIONS. Asserted, not assumed: a fixture that quietly stopped
    // being a two-break ribbon, or stopped running leg A the intended way
    // round, would make every check below pass without testing anything.
    const std::vector< size_t > breaks =
      uvMonotoneBreaks( mesh, boundaryCount );

    check( breaks.size() == 2, "the fixture is a two-break ribbon" );

    if ( breaks.size() != 2 ) {
      continue;
    }

    const double vFirst  = mesh.vertices[ breaks[ 0 ] ].uv.y;
    const double vSecond = mesh.vertices[ breaks[ 1 ] ].uv.y;

    printf( "    leg A spans v %.3f -> %.3f, %zu boundary points\n",
            vFirst, vSecond, boundaryCount );

    check( ( vFirst > vSecond ) == ( descending != 0 ),
           "leg A runs the direction this iteration is testing" );

    // And the legs really are unevenly sampled, which is what makes the
    // same-chain branch fire at all.
    check( DENSE > 4 * SPARSE,
           "one leg carries many more samples than the other" );

    double shoelace = 0.0;

    for ( size_t at = 0; at < boundaryCount; ++at ) {

      const glm::dvec2& p = mesh.vertices[ at ].uv;
      const glm::dvec2& q = mesh.vertices[ ( at + 1 ) % boundaryCount ].uv;

      shoelace += ( p.x * q.y ) - ( q.x * p.y );
    }

    shoelace = std::fabs( shoelace ) * 0.5;

    conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

    const bool built =
      conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

    check( built, "the unevenly sampled ribbon lofts" );

    if ( !built ) {
      continue;
    }

    check( mesh.triangles.size() == boundaryCount - 2,
           "it emits exactly boundary - 2 triangles" );

    // Ray-cast every triangle's uv centroid against the boundary. Independent
    // of the sweep and of the gate's area sum, which is the point.
    size_t outside = 0;
    double absolute = 0.0;

    for ( const conway::geometry::ConnectedTriangle& triangle :
            mesh.triangles ) {

      const glm::dvec2& a = mesh.vertices[ triangle.vertices[ 0 ] ].uv;
      const glm::dvec2& b = mesh.vertices[ triangle.vertices[ 1 ] ].uv;
      const glm::dvec2& c = mesh.vertices[ triangle.vertices[ 2 ] ].uv;

      absolute += std::fabs(
        0.5 * ( ( ( b.x - a.x ) * ( c.y - a.y ) ) -
                ( ( c.x - a.x ) * ( b.y - a.y ) ) ) );

      const glm::dvec2 centroid = ( a + b + c ) / 3.0;

      bool in = false;

      for ( size_t at = 0, previous = boundaryCount - 1;
            at < boundaryCount;
            previous = at++ ) {

        const glm::dvec2& here = mesh.vertices[ at ].uv;
        const glm::dvec2& last = mesh.vertices[ previous ].uv;

        if ( ( ( here.y > centroid.y ) != ( last.y > centroid.y ) ) &&
             ( centroid.x < ( ( last.x - here.x ) *
                              ( centroid.y - here.y ) /
                              ( last.y - here.y ) ) + here.x ) ) {
          in = !in;
        }
      }

      if ( !in ) {
        ++outside;
      }
    }

    printf( "    triangles %zu, outside the trim %zu, coverage %.9f\n",
            mesh.triangles.size(), outside, absolute / shoelace );

    check( outside == 0,
           "no emitted triangle has its centroid outside the trim loop" );

    check( std::fabs( ( absolute / shoelace ) - 1.0 ) < 1e-9,
           "and the emitted triangles cover the chart exactly once" );
  }
}

/**
 * A ribbon whose overlap is SMALL ENOUGH TO PASS the validity gate, so the
 * gate cannot be what pins this.
 *
 * The test above is refused by the gate when the sign is wrong, which is the
 * defect's usual outcome: the face falls back to the ear-clipper and ships the
 * chords. But the gate's tolerance is relative and loose (1e-6 of the chart's
 * area, for the accumulation reasons stated at the gate), so a face whose
 * same-chain branch fires over only a short stretch of its boundary overlaps by
 * LESS THAN THAT and is accepted with triangles lying outside itself. That is
 * not hypothetical: the corpus face with 1,112 boundary points sits at a
 * covered/shoelace ratio of 1.0000002 and ships 13 such triangles today
 * (bldrs-ai/conway#665).
 *
 * The fixture reproduces that regime deliberately. The legs interleave
 * one-for-one almost everywhere - the quiet case - with one short cluster of
 * extra samples on leg A, wandering in u by 1e-5. That amplitude was chosen by
 * measurement, not by taste: with the chain-side sign taken from the label, the
 * emitted mesh covers 1.000000546 of its chart, i.e. 5.5e-7 over - inside the
 * gate's 1e-6 - while four of its triangles have their centroid outside the
 * trim. The u wander matters; at exactly zero the mis-signed sweep fans over a
 * locally convex cap and the fan is accidentally valid, which is why the same
 * fixture with a straight leg proves nothing.
 *
 * So the only oracle that can fail here is the point-in-polygon census.
 */
void testOverlapBelowTheAreaGateIsStillCaught() {

  printf( "\n== overlap too small for the area gate ==\n" );

  constexpr size_t STEPS   = 24;
  constexpr size_t CLUSTER = 6;
  constexpr double WANDER  = 1e-5;
  constexpr double WIDTH   = 1.0;

  WingedEdgeMesh< ParameterVertex > mesh;

  // Leg A, downward from the top turning point - the reversed case - with a
  // short cluster of extra samples inside its first step, wandering slightly
  // in u the way a real trim curve does between its samples.
  mesh.makeVertex( ribbonVertex( 0.0, 1.0 ) );

  for ( size_t at = 1; at <= CLUSTER; ++at ) {

    const double within =
      static_cast< double >( at ) /
      ( static_cast< double >( CLUSTER + 1 ) * static_cast< double >( STEPS ) );

    mesh.makeVertex(
      ribbonVertex( WANDER * std::sin( 2.1 * static_cast< double >( at ) ),
                    1.0 - within ) );
  }

  for ( size_t at = 1; at <= STEPS; ++at ) {

    const double along =
      static_cast< double >( at ) / static_cast< double >( STEPS );

    mesh.makeVertex( ribbonVertex( 0.0, 1.0 - along ) );
  }

  // Leg B back up, interleaved half a step with leg A's regular samples.
  for ( size_t at = 1; at < STEPS; ++at ) {

    const double along =
      ( static_cast< double >( at ) - 0.5 ) / static_cast< double >( STEPS );

    mesh.makeVertex( ribbonVertex( WIDTH, along ) );
  }

  const size_t boundaryCount = mesh.vertices.size();

  const std::vector< size_t > breaks = uvMonotoneBreaks( mesh, boundaryCount );

  check( breaks.size() == 2, "the capped fixture is a two-break ribbon" );

  if ( breaks.size() != 2 ) {
    return;
  }

  check( mesh.vertices[ breaks[ 0 ] ].uv.y > mesh.vertices[ breaks[ 1 ] ].uv.y,
         "leg A runs downward, so it is the reversed case" );

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount );

  // Stated as a check rather than an early return: if this fixture ever stops
  // lofting, the census below silently stops running and the test goes quiet
  // instead of red.
  check( built,
         "the capped ribbon lofts - its overlap is under the gate's tolerance" );

  if ( !built ) {
    return;
  }

  size_t outside = 0;

  for ( const conway::geometry::ConnectedTriangle& triangle : mesh.triangles ) {

    const glm::dvec2& a = mesh.vertices[ triangle.vertices[ 0 ] ].uv;
    const glm::dvec2& b = mesh.vertices[ triangle.vertices[ 1 ] ].uv;
    const glm::dvec2& c = mesh.vertices[ triangle.vertices[ 2 ] ].uv;

    const glm::dvec2 centroid = ( a + b + c ) / 3.0;

    bool in = false;

    for ( size_t at = 0, previous = boundaryCount - 1;
          at < boundaryCount;
          previous = at++ ) {

      const glm::dvec2& here = mesh.vertices[ at ].uv;
      const glm::dvec2& last = mesh.vertices[ previous ].uv;

      if ( ( ( here.y > centroid.y ) != ( last.y > centroid.y ) ) &&
           ( centroid.x < ( ( last.x - here.x ) *
                            ( centroid.y - here.y ) /
                            ( last.y - here.y ) ) + here.x ) ) {
        in = !in;
      }
    }

    if ( !in ) {
      ++outside;
    }
  }

  printf( "  boundary %zu, triangles %zu, outside the trim %zu\n",
          boundaryCount, mesh.triangles.size(), outside );

  check( outside == 0,
         "no emitted triangle has its centroid outside the trim loop" );
}

}  // namespace

int main() {

  testPlateauAtWrapIsNotAReversal();
  testPlainRibbonDecomposes();
  testLoftCoversAndReusesTheBoundary();
  testNotchedBoundaryIsRejected();
  testConcaveRibbonTilesExactlyOnce();
  testSweepEmitsExactlyTheEarcutSeed();
  testNoTVerticesAgainstANeighbour();
  testFarFromOriginKnotDomain();
  testUnevenlySampledRibbonStaysInsideItsTrim();
  testOverlapBelowTheAreaGateIsStillCaught();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
