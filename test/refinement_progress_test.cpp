/*
 * Subdivision-progress tests for `tesselate`'s ParameterVertex overload
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
 * `tesselate` stops it by bounding an off-diagonal at the key of the edge
 * whose split produced it, so a chain of them is strictly decreasing and
 * cannot close. THE TWO HALVES OF THIS FILE ARE THE TWO WAYS THAT BOUND CAN
 * BE WRONG:
 *
 *   - `spineDoesNotRunAway()` is the runaway in the small. Without the bound
 *     it does not terminate.
 *   - `codexAnisotropicCase()` is the bound applied where it means nothing.
 *     The key is deflection times chord and deflection is DIRECTIONAL, so
 *     across a face curved harder in one parameter than the other a genuine
 *     off-diagonal reads a higher key than the edge that spawned it and a
 *     blanket ceiling discards it for good. That case is why the ceiling is
 *     lifted for an off-diagonal nearly orthogonal to its parent, and this is
 *     the review finding on bldrs-ai/conway-geom#211.
 *
 * `ordinaryPatchIsRefinedAsBefore()` is the control for both: a patch with
 * real interior edges, refined to a pinned triangle count, so that a bound
 * which starved ordinary refinement could not pass this file either.
 *
 * Standalone by design: includes the header directly and links nothing but
 * the Logger stubs below, matching outer_bound_order_test.cpp.
 */
#include "conway_geometry/operations/tesselation_utils.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <tuple>
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

using conway::geometry::Edge;
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

/**
 * z = 0.1u^2 + v^2 - a regular surface, curved TEN TIMES harder across v than
 * along u. Nothing about it is degenerate; the anisotropy is the whole point.
 * See `codexAnisotropicCase()`.
 */
glm::dvec3 anisotropicParaboloid( const glm::dvec3&, const glm::dvec2& uv ) {

  return glm::dvec3( uv.x, uv.y, ( 0.1 * uv.x * uv.x ) + ( uv.y * uv.y ) );
}

template< typename SurfacePointFunction >
ParameterVertex on( SurfacePointFunction surface, double u, double v ) {

  return ParameterVertex{ surface( glm::dvec3( 0.0 ), glm::dvec2( u, v ) ),
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

/**
 * The worst squared deflection left on any live interior edge - what the
 * refinement still owes when it returns, in the units `tesselate` compares
 * against `minimumDeflection`.
 */
template< typename SurfacePointFunction >
double worstRemainingDeflection(
  const WingedEdgeMesh< ParameterVertex >& mesh,
  SurfacePointFunction                     surface ) {

  double worst = 0.0;

  for ( const Edge& edge : mesh.edges ) {

    if ( edge.border() ) {
      continue;
    }

    const ParameterVertex& v0 = mesh.vertices[ edge.vertices[ 0 ] ];
    const ParameterVertex& v1 = mesh.vertices[ edge.vertices[ 1 ] ];

    glm::dvec3 average = ( v0.point + v1.point ) * 0.5;
    glm::dvec3 delta   = surface( average, ( v0.uv + v1.uv ) * 0.5 ) - average;

    worst = std::max( worst, glm::dot( delta, delta ) );
  }

  return worst;
}

// A target of 0.02 linear on a unit cylinder: a 0.8-radian chord sags 0.08
// and a 0.2-radian one sags 0.005, so every edge of the seeds below is well
// above it and refinement stops after a couple of halvings.
constexpr double TARGET_DEFLECTION         = 0.02;
constexpr double TARGET_DEFLECTION_SQUARED = TARGET_DEFLECTION * TARGET_DEFLECTION;

// Far above anything either cylinder seed needs, so that reaching it means
// the loop did not terminate rather than that it ran out of room.
constexpr int32_t GENEROUS_BUDGET = 1000;

// What the refinement loop reaches on the 4x4 grid below - with the bound
// and, measured, without it.
constexpr size_t ORDINARY_PATCH_TRIANGLES = 48;

void spineDoesNotRunAway() {

  // Four vertices at u = 0, 0.4, 0.8, 1.2 on a unit circle, triangulated as
  // ( v0, v2, v1 ) and ( v0, v2, v3 ), so the shared edge ( v0, v2 ) has v1 -
  // its own uv midpoint - as one apex and v3 as the other. Splitting it puts
  // a vertex on v1, and the off-diagonal ( v3, new ) spans u 0.4 to 1.2: the
  // same 0.8 the parent spanned, one step further along, and PARALLEL to it,
  // which is what keeps it inside the ceiling's reach. Unbounded, that is a
  // period-2 orbit - measured, 1,000 triangles and 251 vertices on one point.
  WingedEdgeMesh< ParameterVertex > mesh;

  mesh.makeVertex( on( cylinder, 0.0, 0.0 ) );
  mesh.makeVertex( on( cylinder, 0.4, 0.0 ) );
  mesh.makeVertex( on( cylinder, 0.8, 0.0 ) );
  mesh.makeVertex( on( cylinder, 1.2, 0.0 ) );

  mesh.makeTriangle( 0, 2, 1 );
  mesh.makeTriangle( 0, 2, 3 );

  conway::geometry::tesselate(
    mesh, cylinder, GENEROUS_BUDGET, TARGET_DEFLECTION_SQUARED );

  printf(
    "      spine: %zu triangles, %zu vertices, largest coincident cluster %zu\n",
    mesh.triangles.size(), mesh.vertices.size(),
    largestCoincidentCluster( mesh ) );

  check(
    static_cast< int32_t >( mesh.triangles.size() ) < GENEROUS_BUDGET,
    "the spine stops on its own rather than on its triangle budget" );

  check(
    largestCoincidentCluster( mesh ) <= 2,
    "and does so without stacking vertices on one point" );
}

void codexAnisotropicCase() {

  // THE CASE THAT SAYS THE CEILING MAY NOT BE APPLIED BLANKET, raised in
  // review on bldrs-ai/conway-geom#211.
  //
  // Two triangles sharing ( -1, 0 ) - ( 1, 0 ) with apices ( 0, +-1 ), on a
  // surface ten times more curved across v than along u. The shared edge runs
  // the flat way: squared deflection 0.01 over a chord of 2, so key 0.02.
  // Each apex-to-midpoint off-diagonal runs the curved way: squared
  // deflection 0.0625 over a chord of sqrt( 2 ), so key 0.088 - FOUR TIMES
  // the key of the edge that spawned it, though it is a genuine halving of
  // the quad and the surface is regular everywhere. Bounded at the parent's
  // key both are discarded and never reconsidered, and since the two true
  // halves are already under target the queue empties four triangles in,
  // leaving interior edges at 0.0625 against a 0.005 target - TWELVE TIMES
  // over tolerance, and worse than the unbounded loop's own 0.0149.
  //
  // The off-diagonals are exactly orthogonal to the edge that made them, so
  // the ceiling is lifted for them and the refinement proceeds. The two
  // assertions are the two halves of the finding: that these edges are not
  // thrown away, and that what the loop settles at is far below what the
  // blanket ceiling settled at.
  constexpr double  TARGET = 0.005;
  constexpr int32_t BUDGET = 4000;

  // What an UNREFINED first-split off-diagonal reads on this seed, and so
  // exactly what the blanket ceiling leaves behind: 4 triangles and this.
  constexpr double DISCARDED_OFF_DIAGONAL = 0.0625;

  // What this loop settles at instead, measured. The unbounded loop reaches
  // 0.0149383 before its own orbit takes the rest of the budget, so this is
  // most of the way back to it from 0.0625.
  constexpr double SETTLES_AT = 0.0215724;

  WingedEdgeMesh< ParameterVertex > mesh;

  mesh.makeVertex( on( anisotropicParaboloid, -1.0,  0.0 ) );
  mesh.makeVertex( on( anisotropicParaboloid,  1.0,  0.0 ) );
  mesh.makeVertex( on( anisotropicParaboloid,  0.0,  1.0 ) );
  mesh.makeVertex( on( anisotropicParaboloid,  0.0, -1.0 ) );

  mesh.makeTriangle( 0, 1, 2 );
  mesh.makeTriangle( 1, 0, 3 );

  conway::geometry::tesselate( mesh, anisotropicParaboloid, BUDGET, TARGET );

  const double worst = worstRemainingDeflection( mesh, anisotropicParaboloid );

  printf(
    "      codex: %zu triangles, %zu vertices, largest coincident cluster %zu,"
    " worst remaining %.6g\n",
    mesh.triangles.size(), mesh.vertices.size(),
    largestCoincidentCluster( mesh ), worst );

  check(
    worst < DISCARDED_OFF_DIAGONAL,
    "the off-diagonals across the curved direction are refined, not discarded" );

  check(
    worst <= SETTLES_AT,
    "and the loop settles within a third of what discarding them settles at" );

  check(
    static_cast< int32_t >( mesh.triangles.size() ) < BUDGET,
    "while still stopping on its own rather than on its triangle budget" );

  check(
    largestCoincidentCluster( mesh ) == 1,
    "with every vertex on a position of its own" );
}

void ordinaryPatchIsRefinedAsBefore() {

  // A 4x4 grid rather than one quad: `tesselate` only splits INTERIOR edges,
  // so a single quad offers it exactly one and says nothing about whether the
  // bound costs an ordinary face any refinement.
  constexpr uint32_t SIDE = 4;

  WingedEdgeMesh< ParameterVertex > mesh;

  for ( uint32_t row = 0; row < SIDE; ++row ) {
    for ( uint32_t column = 0; column < SIDE; ++column ) {

      mesh.makeVertex( on( cylinder, 0.8 * column, 0.8 * row ) );
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

  printf( "\n=== an anisotropic face keeps its off-diagonals ===\n" );
  codexAnisotropicCase();

  printf( "\n=== an ordinary patch is refined exactly as before ===\n" );
  ordinaryPatchIsRefinedAsBefore();

  if ( failures != 0 ) {

    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
