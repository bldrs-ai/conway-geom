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
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount, 1e-6 );

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
    conway::geometry::tryRibbonLoft( mesh, solve, mesh.vertices.size(), 1e-6 );

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
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount, 1e-6 );

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
 * A small ribbon at its column ceiling stays inside the amplification budget.
 *
 * The ceiling started as the constant 1 + MAX_TRIANGLE_AMPLIFACTION / 2, which
 * assumed rungs ~ boundary. That holds for a large ribbon and fails for a
 * small one, so the cap meant to DERIVE from MAX_TRIANGLE_AMPLIFACTION quietly
 * exceeded it (codex on conway-geom#190).
 */
void testSmallRibbonStaysInBudget() {

  printf( "\n== small ribbon stays in budget ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 4, 1.0, false );

  const size_t boundaryCount = mesh.vertices.size();

  conway::geometry::RationalNurbsInverseMethod solve( flatSurface() );

  // A target of zero is unreachable, so the column search runs to its ceiling
  // - which is exactly the case the constant cap got wrong.
  const bool built =
    conway::geometry::tryRibbonLoft( mesh, solve, boundaryCount, 0.0 );

  check( built, "the small ribbon lofts" );

  if ( !built ) {
    return;
  }

  const size_t budget =
    static_cast< size_t >( conway::geometry::MAX_TRIANGLE_AMPLIFACTION ) *
    ( boundaryCount - 2 );

  printf( "  boundary %zu, emitted %zu triangles, budget %zu\n",
          boundaryCount, mesh.triangles.size(), budget );

  check( mesh.triangles.size() <= budget,
         "emission stays within MAX_TRIANGLE_AMPLIFACTION x the earcut seed" );
}

}  // namespace

int main() {

  testPlateauAtWrapIsNotAReversal();
  testPlainRibbonDecomposes();
  testLoftCoversAndReusesTheBoundary();
  testNotchedBoundaryIsRejected();
  testConcaveRibbonTilesExactlyOnce();
  testSmallRibbonStaysInBudget();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
