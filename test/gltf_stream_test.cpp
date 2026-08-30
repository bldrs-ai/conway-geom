/*
 * Pinning tests for appendComponentToGltfStream — the per-component half of
 * ConwayGeometryProcessor::GeometryToGltf, extracted into a header exactly so
 * these three rules can be tested without linking the writer
 * (bldrs-ai/conway#667):
 *
 *   1. The position bias is subtracted ONCE. The node translation adds it back
 *      once, so a second subtraction anywhere displaces the primitive by -bias
 *      against everything else in the file. Review found precisely that in the
 *      Draco encoder, which re-subtracted a bias the shared array already
 *      carried.
 *   2. Normals go through the normal matrix, and a MIRRORING placement gets its
 *      winding reversed rather than its normals negated — the invariant stated
 *      on TransformMirrors, and the same one Geometry::ApplyTransform follows.
 *      Baking the transform into the positions is what makes this the writer's
 *      job: a consumer only reverses winding for a mirroring NODE transform,
 *      and by then this one is gone.
 *   3. The stream is read through typed vectors. There is no address to
 *      truncate in this signature, which is the point: the writer used to reach
 *      the same data through GetVertexData(), an embind entry point that
 *      narrows the buffer's address to uint32_t.
 *
 * The assertions are on the emitted attributes rather than on internals, so
 * each was checked to fail with its own fix reverted.
 *
 * Standalone by design: includes the header under test and links nothing.
 */
#include "conway_geometry/operations/gltf_stream.h"

#include <glm/gtc/matrix_transform.hpp>

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

using conway::geometry::appendComponentToGltfStream;
using conway::geometry::GLTF_STREAM_STRIDE;
using conway::geometry::TransformMirrors;

/**
 * One triangle in the z = 0 plane, wound counter-clockwise, with +Z normals —
 * so its surface normal and its winding normal agree before any transform, and
 * any disagreement afterwards is the transform's doing.
 *
 * @return The interleaved [px, py, pz, nx, ny, nz] stream.
 */
std::vector< float > flatTriangleStream() {

  return {
    0.0f, 0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
    1.0f, 0.0f, 0.0f,   0.0f, 0.0f, 1.0f,
    0.0f, 1.0f, 0.0f,   0.0f, 0.0f, 1.0f };
}

/**
 * The winding normal of the emitted triangle at `triangle`.
 *
 * `pointOffset` is subtracted back off, because the emitted indices are rebased
 * onto the whole collection while `positions` here holds one component.
 */
glm::dvec3 windingNormal(
  const std::vector< float >&    positions,
  const std::vector< uint32_t >& indexData,
  size_t                         triangle,
  uint32_t                       pointOffset = 0 ) {

  auto corner = [ & ]( size_t which ) {

    uint32_t index = indexData[ ( triangle * 3 ) + which ] - pointOffset;

    return glm::dvec3(
      positions[ index * 3 ], positions[ ( index * 3 ) + 1 ], positions[ ( index * 3 ) + 2 ] );
  };

  glm::dvec3 p0 = corner( 0 );

  return glm::cross( corner( 1 ) - p0, corner( 2 ) - p0 );
}

glm::dvec3 normalAt( const std::vector< float >& normals, uint32_t vertex ) {

  return glm::dvec3(
    normals[ vertex * 3 ], normals[ ( vertex * 3 ) + 1 ], normals[ ( vertex * 3 ) + 2 ] );
}

/*
 * The bias is applied once, by this function, and by nothing else. A caller
 * that subtracts it again — as the Draco encoder did, on the array this
 * function had already biased — moves the whole primitive by -bias while the
 * node translation still adds only +bias back.
 */
void testPositionBiasIsAppliedExactlyOnce() {

  printf( "\n=== the position bias is subtracted exactly once ===\n" );

  std::vector< float >    vertexStream = flatTriangleStream();
  std::vector< uint32_t > indexStream  = { 0, 1, 2 };

  glm::dvec3 bias( 100.0, 200.0, 300.0 );
  glm::dmat4 placement = glm::translate( glm::dmat4( 1.0 ), glm::dvec3( 10.0, 20.0, 30.0 ) );

  std::vector< float >    positions;
  std::vector< float >    normals;
  std::vector< uint32_t > indexData;

  appendComponentToGltfStream(
    vertexStream, indexStream, placement, bias, 0, positions, normals, indexData );

  // Corner 0 is the local origin: placed at (10, 20, 30), biased to (-90, -180, -270).
  glm::dvec3 emitted( positions[ 0 ], positions[ 1 ], positions[ 2 ] );
  glm::dvec3 expected( -90.0, -180.0, -270.0 );

  printf( "      emitted (%.2f %.2f %.2f), expected (%.2f %.2f %.2f)\n",
          emitted.x, emitted.y, emitted.z, expected.x, expected.y, expected.z );

  check( glm::length( emitted - expected ) < 1e-4,
         "a placed vertex is rebased by the bias once, not twice" );

  // A second subtraction would land on the bias again, so the test can tell the
  // two apart by a whole bias vector.
  glm::dvec3 twiceBiased = expected - bias;

  check( glm::length( emitted - twiceBiased ) > 1.0,
         "double subtraction is far enough away to be distinguishable" );
}

/*
 * A mirroring placement: the winding must come out agreeing with the emitted
 * normal. Reverting the index reversal leaves them 180 degrees apart, which is
 * a primitive lit from the inside.
 */
void testMirroredPlacementKeepsNormalsWithWinding() {

  printf( "\n=== a mirrored placement reverses winding, not normals ===\n" );

  std::vector< float >    vertexStream = flatTriangleStream();
  std::vector< uint32_t > indexStream  = { 0, 1, 2 };

  glm::dmat4 mirror = glm::scale( glm::dmat4( 1.0 ), glm::dvec3( 1.0, 1.0, -1.0 ) );

  check( TransformMirrors( mirror ), "the fixture transform really does mirror" );

  std::vector< float >    positions;
  std::vector< float >    normals;
  std::vector< uint32_t > indexData;

  appendComponentToGltfStream(
    vertexStream, indexStream, mirror, glm::dvec3( 0.0 ), 0, positions, normals, indexData );

  glm::dvec3 winding = windingNormal( positions, indexData, 0 );
  glm::dvec3 emitted = normalAt( normals, indexData[ 0 ] );

  double agreement = glm::dot( glm::normalize( winding ), glm::normalize( emitted ) );

  printf( "      winding (%.2f %.2f %.2f), normal (%.2f %.2f %.2f), dot %.4f\n",
          winding.x, winding.y, winding.z, emitted.x, emitted.y, emitted.z, agreement );

  check( agreement > 0.99,
         "the emitted NORMAL faces the same way as the emitted winding" );

  // The normal still describes the mirrored SURFACE, pointing -Z here: this is
  // what separates "reverse the winding" from "negate the normals", which would
  // also satisfy the check above while turning the solid inside out.
  check( emitted.z < -0.99,
         "and it is the mirrored surface normal, not a negated copy" );

  // The reversal is a reordering, so the same three vertices are named.
  check( indexData.size() == 3 &&
         indexData[ 0 ] == 2 && indexData[ 1 ] == 1 && indexData[ 2 ] == 0,
         "the index triple is reversed, not rewritten" );
}

/*
 * The control: a non-mirroring placement must not be reversed. Without this,
 * "reverse when mirrored" could be implemented as "always reverse" and the test
 * above would still pass.
 */
void testOrdinaryPlacementIsNotReversed() {

  printf( "\n=== an ordinary placement keeps its winding ===\n" );

  std::vector< float >    vertexStream = flatTriangleStream();
  std::vector< uint32_t > indexStream  = { 0, 1, 2 };

  glm::dmat4 rotation =
    glm::rotate( glm::dmat4( 1.0 ), 0.7, glm::dvec3( 0.0, 1.0, 0.0 ) );

  check( !TransformMirrors( rotation ), "a rotation does not mirror" );

  std::vector< float >    positions;
  std::vector< float >    normals;
  std::vector< uint32_t > indexData;

  appendComponentToGltfStream(
    vertexStream, indexStream, rotation, glm::dvec3( 0.0 ), 7, positions, normals, indexData );

  check( indexData.size() == 3 &&
         indexData[ 0 ] == 7 && indexData[ 1 ] == 8 && indexData[ 2 ] == 9,
         "corner order is preserved, and rebased by the point offset" );

  glm::dvec3 winding = windingNormal( positions, indexData, 0, 7 );
  glm::dvec3 emitted = normalAt( normals, indexData[ 0 ] - 7 );

  check( glm::dot( glm::normalize( winding ), glm::normalize( emitted ) ) > 0.99,
         "normal and winding still agree" );
}

/*
 * A non-uniform scale is the case that separates the normal matrix from the
 * placement: transforming this normal BY THE PLACEMENT would leave it off the
 * surface by 45 degrees.
 */
void testNormalsUseTheNormalMatrix() {

  printf( "\n=== normals go through the normal matrix ===\n" );

  // A plane at 45 degrees in x/z, with its true normal.
  double root = std::sqrt( 0.5 );

  std::vector< float > vertexStream = {
    0.0f, 0.0f, 0.0f,   (float)root, 0.0f, (float)root,
    1.0f, 0.0f, -1.0f,  (float)root, 0.0f, (float)root,
    0.0f, 1.0f, 0.0f,   (float)root, 0.0f, (float)root };

  std::vector< uint32_t > indexStream = { 0, 1, 2 };

  glm::dmat4 squash = glm::scale( glm::dmat4( 1.0 ), glm::dvec3( 4.0, 1.0, 1.0 ) );

  std::vector< float >    positions;
  std::vector< float >    normals;
  std::vector< uint32_t > indexData;

  appendComponentToGltfStream(
    vertexStream, indexStream, squash, glm::dvec3( 0.0 ), 0, positions, normals, indexData );

  glm::dvec3 emitted = normalAt( normals, 0 );
  glm::dvec3 winding = glm::normalize( windingNormal( positions, indexData, 0 ) );

  printf( "      normal (%.4f %.4f %.4f), winding (%.4f %.4f %.4f)\n",
          emitted.x, emitted.y, emitted.z, winding.x, winding.y, winding.z );

  check( std::abs( glm::length( emitted ) - 1.0 ) < 1e-5,
         "the emitted normal is unit length, as glTF requires" );

  check( glm::dot( emitted, winding ) > 0.999,
         "under a non-uniform scale the normal still lies on the scaled surface" );

  // Scaling the normal by the placement instead would give (4/sqrt2, 0,
  // 1/sqrt2) normalized — about 62 degrees from the true one.
  glm::dvec3 wrong = glm::normalize( glm::dvec3( squash * glm::dvec4( root, 0.0, root, 0.0 ) ) );

  check( glm::dot( emitted, wrong ) < 0.9,
         "and it is measurably not the placement-transformed normal" );
}

/*
 * The degenerate guard: a singular placement has no normal matrix, so the
 * normals are emitted as +Z rather than as NaN, which would fail glTF
 * validation and, through the accessor bounds, truncate the whole file.
 */
void testSingularPlacementEmitsValidNormals() {

  printf( "\n=== a singular placement emits valid normals ===\n" );

  std::vector< float >    vertexStream = flatTriangleStream();
  std::vector< uint32_t > indexStream  = { 0, 1, 2 };

  glm::dmat4 collapse = glm::scale( glm::dmat4( 1.0 ), glm::dvec3( 1.0, 1.0, 0.0 ) );

  std::vector< float >    positions;
  std::vector< float >    normals;
  std::vector< uint32_t > indexData;

  appendComponentToGltfStream(
    vertexStream, indexStream, collapse, glm::dvec3( 0.0 ), 0, positions, normals, indexData );

  bool allFinite = true;

  for ( float value : normals ) {

    allFinite = allFinite && std::isfinite( value );
  }

  check( allFinite, "no NaN reaches the normal buffer" );

  check( std::abs( glm::length( normalAt( normals, 0 ) ) - 1.0 ) < 1e-5,
         "the substituted normal is still unit length" );
}

}  // namespace

int main() {

  testPositionBiasIsAppliedExactlyOnce();
  testMirroredPlacementKeepsNormalsWithWinding();
  testOrdinaryPlacementIsNotReversed();
  testNormalsUseTheNormalMatrix();
  testSingularPlacementEmitsValidNormals();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
