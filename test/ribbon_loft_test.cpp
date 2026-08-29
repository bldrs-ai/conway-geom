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

  // Index 0 sits at the TOP turning point, so the wrap is where the plateau
  // goes - exactly the arrangement measured on the Orbiter faces.
  if ( plateauAtWrap ) {
    mesh.makeVertex( ribbonVertex( 0.0, 1.0 ) );
  }

  for ( size_t at = 0; at <= steps; ++at ) {

    const double v =
      1.0 - ( static_cast< double >( at ) / static_cast< double >( steps ) );

    mesh.makeVertex( ribbonVertex( 0.0, v ) );
  }

  for ( size_t at = 1; at <= steps; ++at ) {

    const double v = static_cast< double >( at ) / static_cast< double >( steps );

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

  WingedEdgeMesh< ParameterVertex > mesh = buildRibbon( 24, 1.0, true );

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

}  // namespace

int main() {

  testPlateauAtWrapIsNotAReversal();
  testPlainRibbonDecomposes();
  testLoftCoversAndReusesTheBoundary();
  testNotchedBoundaryIsRejected();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
