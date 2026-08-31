/*
 * Pinning tests for Geometry::Normalize() (bldrs-ai/conway#680).
 *
 * The defect these pin was two bugs that CANCELLED, which is why nothing ever
 * went red on it: the function returned a never-written member instead of the
 * centre it had just computed (so callers were told the shift was zero), AND
 * it never dropped the reification (so the float32 buffer callers actually
 * upload still held the un-shifted positions). Composed world position was
 * therefore correct, and only the PRECISION of the f32 buffer was wrong — on a
 * Swiss LV95 model with the national-grid coordinates baked into the geometry,
 * vertices sit at ~2.6e6 m where an f32 ULP is 0.25 m, and the model visibly
 * jitters (bldrs-ai/test-models-private#97, bldrs-ai/Share#1634).
 *
 * So the assertions below are deliberately on the two halves SEPARATELY —
 * returned centre, and reified f32 content — because any test that only checks
 * `transform * vertex` passes against the bug. Each was verified to fail with
 * its own half of the fix reverted; the per-test notes say what the failure
 * looks like.
 *
 * Standalone by design: includes Geometry.cpp directly and defines the two
 * AABBTree symbols it references, matching corner_normals_test.cpp and
 * spherical_trim_test.cpp.
 */
#include "conway_geometry/representation/Geometry.cpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace conway::geometry {

// Referenced by ApplyRescale, which none of these tests call.
void AABBTree::applyRescale( const glm::dvec3&, const glm::dvec3& ) {}

/**
 * A root-box-only stand-in for the real BVH build, so this stays a
 * link-nothing translation unit (the real one lives in aabb_tree.cpp, which
 * pulls in the CSG headers).
 *
 * Only bounds() is exercised here, and only for what GetAABB() does with it:
 * serve the tree's root box in preference to scanning the vertices. That makes
 * a BVH which SURVIVED the normalize shift observable — it keeps reporting the
 * pre-shift box while the vertices have moved.
 *
 * The root unions TRIANGLE-REFERENCED vertices only, which is what the real
 * builder does (it boxes triangles, not the vertex array) and is not a detail
 * that can be skipped here: it is the only thing that makes the tree's answer
 * DIFFER from GetAABB()'s vertex-scan fallback, and a test where the two agree
 * cannot pin the ordering inside Normalize().
 */
AABBTree::AABBTree( const Geometry& mesh, double ) {

  box3 root;

  for ( const Triangle& triangle : mesh.triangles ) {

    for ( uint32_t corner = 0; corner < 3; ++corner ) {

      root.merge( mesh.vertices[ triangle.vertices[ corner ] ] );
    }
  }

  boxes_.push_back( root );
  triangles_.resize( mesh.triangles.size() );
}

}  // namespace conway::geometry

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

using conway::geometry::Geometry;

/** Floats per vertex in the reified stream: position then normal. */
constexpr size_t STRIDE = 6;

/**
 * Swiss-LV95-scale offset, standing in for Ecobau.ifc's baked national-grid
 * coordinates. 2.6e6 is in the binade where float32's ULP is 0.25 m, which is
 * the whole reason Normalize() exists; a test at a small offset would pass
 * against the unfixed code.
 *
 * The fractional parts are not decoration: a round 2646000 lands exactly on
 * the float32 grid at this magnitude (it is a multiple of 0.25), so the
 * quantization assertion below would read zero error and pass vacuously.
 */
const glm::dvec3 GEOREF( 2646123.37, 1249456.71, 412.53 );

/** Half-extent of the test box, in metres. */
constexpr double HALF_EXTENT = 10.0;

/**
 * An axis-aligned box of 8 vertices and 12 triangles, centred on `centre`.
 *
 * Its AABB centre is exactly `centre`, which is what every assertion about the
 * returned value is compared against.
 *
 * @param centre Where to put the box.
 * @return The geometry.
 */
Geometry makeBox( const glm::dvec3& centre ) {

  Geometry geometry;

  for ( int corner = 0; corner < 8; ++corner ) {

    geometry.MakeVertex(
      centre +
      glm::dvec3(
        ( corner & 1 ) ? HALF_EXTENT : -HALF_EXTENT,
        ( corner & 2 ) ? HALF_EXTENT : -HALF_EXTENT,
        ( corner & 4 ) ? HALF_EXTENT : -HALF_EXTENT ) );
  }

  // Two triangles per face, corners ordered so no face is degenerate.
  const uint32_t faces[ 6 ][ 4 ] = {
    { 0, 2, 3, 1 },  // -z
    { 4, 5, 7, 6 },  // +z
    { 0, 1, 5, 4 },  // -y
    { 2, 6, 7, 3 },  // +y
    { 0, 4, 6, 2 },  // -x
    { 1, 3, 7, 5 },  // +x
  };

  for ( const uint32_t( &face )[ 4 ] : faces ) {

    geometry.MakeTriangle( face[ 0 ], face[ 1 ], face[ 2 ] );
    geometry.MakeTriangle( face[ 0 ], face[ 2 ], face[ 3 ] );
  }

  return geometry;
}

/** The largest absolute position component in a reified vertex stream. */
double maxAbsPosition( const std::vector< float >& stream ) {

  double result = 0.0;

  for ( size_t where = 0; where + STRIDE <= stream.size(); where += STRIDE ) {

    for ( size_t axis = 0; axis < 3; ++axis ) {

      result = std::max( result, std::abs( static_cast< double >( stream[ where + axis ] ) ) );
    }
  }

  return result;
}

bool near( const glm::dvec3& left, const glm::dvec3& right, double tolerance ) {

  return
    std::abs( left.x - right.x ) <= tolerance &&
    std::abs( left.y - right.y ) <= tolerance &&
    std::abs( left.z - right.z ) <= tolerance;
}

/**
 * The centre that comes back is the one that was subtracted.
 *
 * Fails with `return centre;` reverted to the dead member: the returned vector
 * is (0,0,0), 2.6e6 m from the box.
 */
void returnsTheSubtractedCentre() {

  printf( "\n=== the returned centre is the shift, not zero ===\n" );

  Geometry geometry = makeBox( GEOREF );

  glm::dvec3 before = geometry.GetAABB().centre();

  check( near( before, GEOREF, 1e-6 ), "the box starts at the georeferenced centre" );

  glm::dvec3 centre = geometry.Normalize();

  check( near( centre, GEOREF, 1e-6 ), "Normalize() reports that centre" );

  check(
    std::abs( centre.x ) > 1.0 && std::abs( centre.y ) > 1.0,
    "and it is emphatically not the zero a dead member would give" );
}

/**
 * The f64 vertices move by the centre, and the AABB follows them.
 */
void shiftsTheVertices() {

  printf( "\n=== the vertices are shifted onto their own centre ===\n" );

  Geometry original = makeBox( GEOREF );
  Geometry geometry = makeBox( GEOREF );

  glm::dvec3 centre = geometry.Normalize();

  bool allShifted = true;

  for ( uint32_t vertex = 0; vertex < geometry.GetVertexCount(); ++vertex ) {

    // 1e-9 m, nine orders below anything that matters, and comfortably above
    // the ~5e-10 rounding of subtracting 2.6e6 from 2.6e6 in f64.
    allShifted = allShifted &&
      near( geometry.GetPoint( vertex ), original.GetPoint( vertex ) - centre, 1e-9 );
  }

  check( allShifted, "every vertex moved by exactly the reported centre" );

  check(
    near( geometry.GetAABB().centre(), glm::dvec3( 0.0 ), 1e-9 ),
    "so the AABB is now centred on the origin" );

  check(
    std::abs( geometry.GetAABB().interval().x - ( HALF_EXTENT * 2.0 ) ) < 1e-9,
    "and the box is the same size it was" );
}

/**
 * THE REGRESSION. A geometry reified BEFORE the shift must not keep serving
 * the pre-shift float32 buffer afterwards.
 *
 * This is the half that reaches the GPU: Share uploads GetVertexData(), and
 * `Reify()` early-returns on `isReified_`, so without the ClearReification()
 * in Normalize() the buffer below still holds 2.6e6-metre coordinates no
 * matter how far the f64 vertices moved.
 *
 * Fails with that ClearReification() removed: `after` comes back at 2646010,
 * not 10.
 */
void dropsThePreShiftReification() {

  printf( "\n=== a reification built before the shift does not survive it ===\n" );

  Geometry geometry = makeBox( GEOREF );

  // The state the walk actually hands to Normalize(): already reified, because
  // the CSG paths go through GetVertexStream().
  double beforeMax = maxAbsPosition( geometry.GetVertexStream() );

  check(
    beforeMax > 2.0e6,
    "the pre-normalize stream holds the raw georeferenced coordinates" );

  // What that magnitude costs in float32, which is the symptom being fixed:
  // 0.25 m of quantization at this scale, versus exact at 10 m.
  double quantized =
    static_cast< double >( static_cast< float >( GEOREF.x + HALF_EXTENT ) ) -
    ( GEOREF.x + HALF_EXTENT );

  check(
    std::abs( quantized ) > 0.05,
    "which float32 cannot even represent to a tenth of a metre" );

  geometry.Normalize();

  double afterMax = maxAbsPosition( geometry.GetVertexStream() );

  check(
    afterMax < HALF_EXTENT + 1e-6,
    "the stream re-reified from the shifted vertices" );

  check(
    afterMax > HALF_EXTENT - 1e-6,
    "and it is the whole box, not a collapsed one" );

  check(
    static_cast< double >( static_cast< float >( HALF_EXTENT ) ) == HALF_EXTENT,
    "at which magnitude float32 is exact" );
}

/**
 * A BVH built before the shift is dropped, so GetAABB() stops reporting the
 * pre-shift box — and the centre is measured from that tree BEFORE it goes.
 *
 * The orphan vertex is what gives the ordering claim teeth. GetAABB() has two
 * answers: `bvh->bounds()` (triangle-referenced vertices only) and the
 * vertex-scan fallback (every vertex, orphan included). With a box alone the
 * two agree and a reordered Normalize() would measure the same centre either
 * way; one vertex outside the box that no triangle references separates them
 * by 495 m on x.
 *
 * Fails with the bvh.reset() removed: the post-normalize AABB centre is still
 * the pre-shift one. Fails with the measurement moved after the reset: the
 * centre comes back as the orphan-inclusive one.
 */
void dropsTheStaleBvh() {

  printf( "\n=== a BVH does not outlive the positions it indexes ===\n" );

  Geometry geometry = makeBox( GEOREF );

  // Referenced by no triangle, so it lands in the vertex scan and not in the
  // tree.
  constexpr double ORPHAN_X_OFFSET = 1000.0;

  const glm::dvec3 ORPHAN = GEOREF + glm::dvec3( ORPHAN_X_OFFSET, 0.0, 0.0 );

  geometry.MakeVertex( ORPHAN );

  // Where the vertex-scan fallback would put the centre: midway between the
  // box's -x face and the orphan. Derived from the two constants that decide
  // it rather than written out, so moving the orphan or resizing the box
  // breaks this loudly instead of leaving a stale number that still passes.
  const glm::dvec3 SCAN_CENTRE =
    GEOREF + glm::dvec3( ( ORPHAN_X_OFFSET - HALF_EXTENT ) * 0.5, 0.0, 0.0 );

  geometry.MakeBVH();

  check( geometry.bvh.has_value(), "the geometry has a BVH going in" );

  check(
    near( geometry.GetAABB().centre(), GEOREF, 1e-6 ),
    "which is what GetAABB() is reading, and it excludes the orphan vertex" );

  glm::dvec3 centre = geometry.Normalize();

  check(
    near( centre, GEOREF, 1e-6 ),
    "the centre was measured from the tree, before it was dropped" );

  check(
    !near( centre, SCAN_CENTRE, 1.0 ),
    "and not from the vertex scan the reset would have fallen back to" );

  check( !geometry.bvh.has_value(), "the stale tree is gone" );

  check(
    near( geometry.GetAABB().centre(), SCAN_CENTRE - GEOREF, 1e-9 ),
    "so GetAABB() reports the shifted box, scanned rather than served stale" );
}

/**
 * The second call is the one IFC mapped items make: one representation-map
 * body, walked once per instancing product. It has to answer with the centre
 * the first call subtracted — the current near-zero AABB centre would place
 * every instance after the first on the coordination origin (conway#308).
 */
void secondCallRepeatsTheFirstCentre() {

  printf( "\n=== normalizing twice places the geometry twice the same ===\n" );

  Geometry geometry = makeBox( GEOREF );

  glm::dvec3 first = geometry.Normalize();

  glm::dvec3 firstVertex = geometry.GetPoint( 0 );

  glm::dvec3 second = geometry.Normalize();

  check( near( second, first, 0.0 ), "the second call returns the first centre exactly" );

  // Not redundant with the line above: two zeros also agree exactly, which is
  // what the unfixed function returned from both calls.
  check( near( second, GEOREF, 1e-6 ), "and that centre is the real one" );

  check(
    near( geometry.GetPoint( 0 ), firstVertex, 0.0 ),
    "and shifts nothing further" );

  check(
    near( geometry.GetAABB().centre(), glm::dvec3( 0.0 ), 1e-9 ),
    "so the box stays where the first call put it" );
}

/**
 * Empty geometry has an inverted AABB, whose centre() is -inf. Zero is the
 * only answer a caller can compose into a transform.
 *
 * Fails with the vertices.empty() guard removed: the returned components are
 * -inf, and an Ecobau load ends with 20 empty geometries (aborted booleans) to
 * feed it.
 */
void emptyGeometryReportsZero() {

  printf( "\n=== an empty geometry reports a finite zero ===\n" );

  Geometry geometry;

  glm::dvec3 centre = geometry.Normalize();

  check( std::isfinite( centre.x ) && std::isfinite( centre.y ) && std::isfinite( centre.z ),
    "the centre is finite" );

  check( near( centre, glm::dvec3( 0.0 ), 0.0 ), "and zero" );
}

/**
 * Clear() drops the reification and the normalize state with the vertices
 * those describe — the same "mutated the positions, kept the cache" shape as
 * the Normalize() defect, one function over.
 */
void clearDropsDerivedState() {

  printf( "\n=== Clear() takes the derived state with it ===\n" );

  Geometry geometry = makeBox( GEOREF );

  geometry.Normalize();

  check( !geometry.GetVertexStream().empty(), "the geometry reifies to something" );

  geometry.Clear();

  check( geometry.GetVertexStream().empty(), "a cleared geometry reifies to nothing" );

  Geometry refilled = makeBox( GEOREF );

  refilled.Normalize();
  refilled.Clear();

  // Refilled at a DIFFERENT centre on purpose. Refilling at GEOREF would make
  // this pass with Clear()'s resets reverted — the short-circuit would return
  // the remembered GEOREF, which is also the right answer for the new
  // contents, so the assertion could not tell the two apart.
  const glm::dvec3 MOVED = GEOREF + glm::dvec3( 500.0, -500.0, 0.0 );

  for ( uint32_t vertex = 0; vertex < 8; ++vertex ) {

    refilled.MakeVertex(
      MOVED +
      glm::dvec3(
        ( vertex & 1 ) ? HALF_EXTENT : -HALF_EXTENT,
        ( vertex & 2 ) ? HALF_EXTENT : -HALF_EXTENT,
        ( vertex & 4 ) ? HALF_EXTENT : -HALF_EXTENT ) );
  }

  check(
    near( refilled.Normalize(), MOVED, 1e-6 ),
    "and normalizes again from scratch, reporting the new contents' centre "
    "rather than short-circuiting on the old flag" );
}

/**
 * Baking a transform into the positions ends the normalized state, so the next
 * Normalize() measures the geometry where it now is.
 *
 * `applyTransform` and `normalize` are both embind entry points, so a caller
 * can interleave them freely. Fails with the reset at the end of
 * ApplyTransform removed: the second call short-circuits and returns the
 * pre-transform centre — which, unlike the old (0,0,0), looks entirely
 * plausible at the call site.
 */
void transformEndsTheNormalizedState() {

  printf( "\n=== baking a transform ends the normalized state ===\n" );

  Geometry geometry = makeBox( GEOREF );

  glm::dvec3 first = geometry.Normalize();

  check( near( first, GEOREF, 1e-6 ), "the first call centres the box" );

  const glm::dvec3 OFFSET( 100.0, 200.0, 300.0 );

  glm::dmat4 translation( 1.0 );

  translation[ 3 ] = glm::dvec4( OFFSET, 1.0 );

  geometry.ApplyTransform( translation );

  check(
    near( geometry.GetAABB().centre(), OFFSET, 1e-9 ),
    "the transform moved it off its centre" );

  glm::dvec3 second = geometry.Normalize();

  check( near( second, OFFSET, 1e-9 ), "so the next call reports where it now is" );

  check(
    !near( second, first, 1.0 ),
    "and not the centre the pre-transform frame had" );

  check(
    near( geometry.GetAABB().centre(), glm::dvec3( 0.0 ), 1e-9 ),
    "and it actually shifted, rather than only reporting" );
}

}  // namespace

int main() {

  returnsTheSubtractedCentre();
  shiftsTheVertices();
  dropsThePreShiftReification();
  dropsTheStaleBvh();
  secondCallRepeatsTheFirstCentre();
  emptyGeometryReportsZero();
  clearDropsDerivedState();
  transformEndsTheNormalizedState();

  if ( failures > 0 ) {

    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
