/*
 * Termination tests for `tesselate`'s ParameterVertex overload
 * (conway_geometry/operations/tesselation_utils.h).
 *
 * Splitting an edge replaces the two triangles on it with four, and queues
 * the four edges that meet the new vertex. Two of those are the bisection's
 * own halves - their uv span is half their parent's, so repeated subdivision
 * drives their deflection to zero. The other two are OFF-DIAGONALS: new edges
 * from a quad apex to the new point, which did not exist before the split and
 * which nothing makes finer than the edge that was split.
 *
 * On a quad whose four corners are strung out along one uv line, the
 * off-diagonal comes out spanning exactly as much uv as its parent did, so
 * refining it re-creates an off-diagonal congruent to the original and the
 * face subdivides for ever at a constant priority key. That is the runaway
 * behind bldrs-ai/conway-geom#208 and #210: on `Right_Hand.step`'s
 * `ADVANCED_FACE #18856` it spent the whole 9,664-triangle budget stacking
 * 4,510 copies of two points, and the 637 triangles that survived
 * `Geometry::Reify` were the earcut seed.
 *
 * The spine below is that configuration in the small: four vertices at
 * u = 0, 0.4, 0.8, 1.2 on a unit circle, triangulated as ( v0, v2, v1 ) and
 * ( v0, v2, v3 ), so the shared edge ( v0, v2 ) has v1 - its own uv midpoint -
 * as one apex and v3 as the other. Splitting it puts a vertex on v1, and the
 * off-diagonal ( v3, new ) spans u 0.4 to 1.2: the same 0.8 the parent
 * spanned, one step further along. Unbounded, that is a period-2 orbit.
 *
 * Verified by breaking rather than by reading, in a scratch copy of
 * tesselation_utils.h reached by an -I override - the working tree is never
 * reverted:
 *
 *   - the two off-diagonal `addCandidate` calls handed
 *     `std::numeric_limits< double >::infinity()` instead of `parentKey`:
 *     both spine assertions go red, at 1,000 triangles (its whole budget) and
 *     503 vertices over a largest coincident cluster of 251;
 *   - `if ( key >= keyCeiling )` made permanently false: the same two, the
 *     same numbers;
 *   - `if ( key >= keyCeiling )` made permanently true, so the ceiling
 *     refuses everything: the CONTROL goes red instead, falling from 48
 *     triangles to its 18-triangle seed.
 *
 * That last one is why the control is here. A bound that starved ordinary
 * refinement would fix the spine and facet every curved face in the corpus,
 * so a patch with real interior edges is refined here too and its triangle
 * count pinned exactly - and the first two variants leave it at 48, which is
 * the claim that this bound costs an ordinary face nothing.
 *
 * Standalone by design: includes the header directly and links nothing but
 * the Logger stubs below, matching outer_bound_order_test.cpp.
 */
#include "conway_geometry/operations/tesselation_utils.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

// The header's error paths call these; the rest is header-only.
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

/**
 * A unit cylinder, parameterised by angle in u. Curved enough that every
 * chord in these meshes reads a deflection well above the target, and cheap
 * enough to be exact: the surface point depends on nothing but the uv, which
 * is what makes the orbit reproduce bit-for-bit.
 */
glm::dvec3 cylinder( const glm::dvec3&, const glm::dvec2& uv ) {

  return glm::dvec3( cos( uv.x ), sin( uv.x ), uv.y );
}

ParameterVertex onCylinder( double u, double v ) {

  return ParameterVertex{ cylinder( glm::dvec3( 0.0 ), glm::dvec2( u, v ) ),
                          glm::dvec2( u, v ) };
}

/** How many vertices share their exact position with another vertex. */
size_t largestCoincidentCluster( const WingedEdgeMesh< ParameterVertex >& mesh ) {

  std::map< std::tuple< double, double, double >, size_t > counts;

  for ( const ParameterVertex& vertex : mesh.vertices ) {

    ++counts[ std::make_tuple(
      vertex.point.x, vertex.point.y, vertex.point.z ) ];
  }

  size_t largest = 0;

  for ( const auto& [ position, count ] : counts ) {
    largest = std::max( largest, count );
  }

  return largest;
}

// A target of 0.02 linear on a unit cylinder: a 0.8-radian chord sags 0.08
// and a 0.2-radian one sags 0.005, so every edge of the seeds below is well
// above it and refinement stops after a couple of halvings.
constexpr double TARGET_DEFLECTION         = 0.02;
constexpr double TARGET_DEFLECTION_SQUARED = TARGET_DEFLECTION * TARGET_DEFLECTION;

// Far above anything either seed needs, so that reaching it means the loop
// did not terminate rather than that it ran out of room.
constexpr int32_t GENEROUS_BUDGET = 1000;

// What the refinement loop reaches on the 4x4 grid below - with the
// off-diagonal bound and, measured, without it.
constexpr size_t ORDINARY_PATCH_TRIANGLES = 48;

void spineDoesNotRunAway() {

  WingedEdgeMesh< ParameterVertex > mesh;

  mesh.makeVertex( onCylinder( 0.0, 0.0 ) );
  mesh.makeVertex( onCylinder( 0.4, 0.0 ) );
  mesh.makeVertex( onCylinder( 0.8, 0.0 ) );
  mesh.makeVertex( onCylinder( 1.2, 0.0 ) );

  mesh.makeTriangle( 0, 2, 1 );
  mesh.makeTriangle( 0, 2, 3 );

  conway::geometry::tesselate(
    mesh, cylinder, GENEROUS_BUDGET, TARGET_DEFLECTION_SQUARED );

  check(
    static_cast< int32_t >( mesh.triangles.size() ) < GENEROUS_BUDGET,
    "the spine stops on its own rather than on its triangle budget" );

  check(
    largestCoincidentCluster( mesh ) <= 2,
    "and does so without stacking vertices on one point" );

  printf(
    "      spine: %zu triangles, %zu vertices, largest coincident cluster %zu\n",
    mesh.triangles.size(), mesh.vertices.size(),
    largestCoincidentCluster( mesh ) );
}

void ordinaryPatchIsRefinedAsBefore() {

  // A 4x4 grid rather than one quad: `tesselate` only splits INTERIOR edges,
  // so a single quad offers it exactly one and says nothing about whether the
  // bound costs an ordinary face any refinement.
  constexpr uint32_t SIDE = 4;

  WingedEdgeMesh< ParameterVertex > mesh;

  for ( uint32_t row = 0; row < SIDE; ++row ) {
    for ( uint32_t column = 0; column < SIDE; ++column ) {

      mesh.makeVertex( onCylinder( 0.8 * column, 0.8 * row ) );
    }
  }

  for ( uint32_t row = 0; row + 1 < SIDE; ++row ) {
    for ( uint32_t column = 0; column + 1 < SIDE; ++column ) {

      uint32_t corner = row * SIDE + column;

      mesh.makeTriangle( corner, corner + 1, corner + SIDE );
      mesh.makeTriangle( corner + 1, corner + SIDE + 1, corner + SIDE );
    }
  }

  conway::geometry::tesselate(
    mesh, cylinder, GENEROUS_BUDGET, TARGET_DEFLECTION_SQUARED );

  printf(
    "      patch: %zu triangles, %zu vertices\n",
    mesh.triangles.size(), mesh.vertices.size() );

  // Pinned exactly, because "the bound refuses nothing an ordinary face
  // needed" is the half of the claim a looser assertion would not hold. This
  // is the count the unbounded loop reaches, measured on both.
  check(
    mesh.triangles.size() == ORDINARY_PATCH_TRIANGLES,
    "an ordinary patch is refined to the same triangle count it always was" );

  check(
    largestCoincidentCluster( mesh ) == 1,
    "with every vertex on a position of its own" );
}

}  // namespace

int main() {

  printf( "\n=== a uv-collinear spine terminates ===\n" );
  spineDoesNotRunAway();

  printf( "\n=== an ordinary patch is refined exactly as before ===\n" );
  ordinaryPatchIsRefinedAsBefore();

  if ( failures != 0 ) {

    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
