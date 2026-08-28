/*
 * Termination test for tesselate()'s refinement loop on an exactly
 * representable uv lattice.
 *
 * tesselate() splits an edge at the midpoint of its two PARAMETERS and takes
 * the surface point there. Where a trim lands on a uv lattice whose spacing is
 * exact in doubles, that midpoint is a FIXED POINT of the split: the new
 * parameter is bitwise the parameter of a vertex the mesh already has, so the
 * surface hands back that vertex's point bitwise, and the "new" vertex is a
 * duplicate. The four triangles that replace the two include two exactly
 * degenerate ones, and one of the edges created carries the deflection the
 * split was supposed to remove - so the queue is handed back its own work and
 * the whole triangle budget goes into producing nothing. On face `#51059` of
 * Orbiter_v1.1_Gear_7.5.step that was 78,610 of 78,833 splits.
 *
 * Nothing about it is thread- or helix-specific, which is why this test does
 * not need the STEP file: a periodic surface plus an integer lattice is the
 * whole recipe. See bldrs-ai/conway#625.
 *
 * Both tests below fail if the guard in tesselate() is removed; that was
 * verified by removing it, not by reading. Without the guard the lattice case
 * runs its budget to zero (10,002 triangles here rather than 2).
 *
 * Standalone by design: it includes the header under test and links nothing,
 * so a plain compiler and the repo's include tree are all it needs, and it
 * does not depend on the conway_geom_native_tests target (which does not
 * currently build).
 */
#include "conway_geometry/operations/tesselation_utils.h"

#include <cmath>
#include <cstdio>
#include <string>

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

/**
 * A helical band about the z axis: v is the turn parameter and u the offset
 * across the band. The period is EIGHT units of v, a power of two, so a
 * lattice of integer v is exactly representable and closed under the
 * midpoints tesselate() takes - which is the whole hazard, in the form a
 * thread's trim hands it.
 *
 * @param uv The parameter to evaluate.
 * @return The point on the band.
 */
glm::dvec3 helicalBand( const glm::dvec2& uv ) {

  constexpr double RADIUS = 10.0;
  constexpr double PITCH  = 3.0;
  constexpr double PERIOD = 8.0;

  const double angle = ( 2.0 * PI * uv.y ) / PERIOD;

  return glm::dvec3(
    ( RADIUS + uv.x ) * std::cos( angle ),
    ( RADIUS + uv.x ) * std::sin( angle ),
    PITCH * uv.y );
}

/** A ParameterVertex on `helicalBand` at the given parameter. */
ParameterVertex bandVertex( double u, double v ) {

  const glm::dvec2 uv( u, v );

  return ParameterVertex{ helicalBand( uv ), uv };
}

/**
 * The two triangles that reproduce the measured trap, at its smallest.
 *
 * Four vertices one lattice step apart in v, and the two triangles that share
 * the edge spanning TWO steps. The apex of one of them sits exactly at that
 * edge's parametric midpoint, so the split lands bitwise on it. The
 * configuration is self-reproducing: the split's own output contains the same
 * arrangement one step along, which is why the measured face ping-ponged
 * 19,653 times between two endpoint pairs rather than stopping.
 *
 * @param mesh Mesh to build into.
 */
void buildLatticeTrap( WingedEdgeMesh< ParameterVertex >& mesh ) {

  const uint32_t a = mesh.makeVertex( bandVertex( 0.0, 0.0 ) );
  const uint32_t b = mesh.makeVertex( bandVertex( 0.0, 1.0 ) );
  const uint32_t c = mesh.makeVertex( bandVertex( 0.0, 2.0 ) );
  const uint32_t d = mesh.makeVertex( bandVertex( 0.0, 3.0 ) );

  mesh.makeTriangle( d, b, c );
  mesh.makeTriangle( b, d, a );
}

/**
 * Refinement on the exact lattice terminates rather than spending its budget.
 */
void testLatticeFixedPointTerminates() {

  printf( "\n== lattice fixed point ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh;

  buildLatticeTrap( mesh );

  // The trap, stated as arithmetic rather than assumed: the midpoint of the
  // spanning edge's two parameters IS the third vertex's parameter, bitwise,
  // so the surface returns that vertex's point bitwise.
  const glm::dvec2 midUV =
    ( mesh.vertices[ 3 ].uv + mesh.vertices[ 1 ].uv ) * 0.5;

  check( midUV == mesh.vertices[ 2 ].uv,
         "the split parameter is bitwise an existing vertex's parameter" );
  check( helicalBand( midUV ) == mesh.vertices[ 2 ].point,
         "the split point is bitwise an existing vertex's point" );

  // Deflection the refinement cannot meet, and a budget it would take 5,000
  // splits to exhaust. A target of zero is the honest statement of "this
  // never converges" - the same regime as the measured face, whose edges sat
  // at 122x their target for every one of the 19,653 iterations.
  const int32_t budget = 2 + ( 2 * 5000 );

  tesselate(
    mesh,
    []( const glm::dvec3&, const glm::dvec2& uv ) { return helicalBand( uv ); },
    budget,
    0.0 );

  printf( "  triangles after refinement: %zu (budget allowed %d)\n",
          mesh.triangles.size(), budget );

  check( mesh.triangles.size() == 2,
         "no triangle was spent on a split that adds no position" );
  check( mesh.vertices.size() == 4,
         "no duplicate vertex was added" );
}

/**
 * A surface off the lattice still refines, so the guard is not simply
 * stopping refinement.
 *
 * The same band, seeded a third of a lattice step off it, so no midpoint is
 * ever an existing parameter. Refinement must run to its budget here - if the
 * guard were rejecting on anything but exact duplication, this is where it
 * would show.
 */
void testOffLatticeStillRefines() {

  printf( "\n== off the lattice ==\n" );

  WingedEdgeMesh< ParameterVertex > mesh;

  const uint32_t a = mesh.makeVertex( bandVertex( 0.0, 0.0 ) );
  const uint32_t b = mesh.makeVertex( bandVertex( 0.0, 1.0 / 3.0 ) );
  const uint32_t c = mesh.makeVertex( bandVertex( 1.0, 2.0 / 3.0 ) );
  const uint32_t d = mesh.makeVertex( bandVertex( 1.0, 1.0 ) );

  mesh.makeTriangle( d, b, c );
  mesh.makeTriangle( b, d, a );

  const int32_t budget = 2 + ( 2 * 64 );

  tesselate(
    mesh,
    []( const glm::dvec3&, const glm::dvec2& uv ) { return helicalBand( uv ); },
    budget,
    0.0 );

  printf( "  triangles after refinement: %zu (budget allowed %d)\n",
          mesh.triangles.size(), budget );

  check( static_cast< int32_t >( mesh.triangles.size() ) == budget,
         "refinement off the lattice still spends its whole budget" );
}

}  // namespace

int main() {

  testLatticeFixedPointTerminates();
  testOffLatticeStillRefines();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
