/*
 * Pinning tests for the analytic corner-normal path in Geometry
 * (bldrs-ai/conway#667), and for the three ways it was found to break its own
 * invariants in review.
 *
 * What is under test here is Reify()'s shading output and the mutators that
 * have to keep `corner_normals` in step with `triangles`. The tessellators
 * that PRODUCE the analytic normals are covered by the conway-side fixture
 * tests; these build the arrays by hand, which is the only way to construct
 * the mixed case (an analytic face welded to an uncovered one) deterministically
 * and at two different scales.
 *
 * Each check below was verified to fail with its own fix reverted, not by
 * reading — see the per-test notes for what the failure looks like.
 *
 * Standalone by design: includes Geometry.cpp directly and stubs the one
 * symbol that would otherwise pull in the BVH translation unit, matching
 * spherical_trim_test.cpp.
 */
#include "conway_geometry/representation/Geometry.cpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// Referenced by ApplyRescale, which none of these tests call. Defined here so
// the test links nothing but itself.
namespace conway::geometry {

void AABBTree::applyRescale( const glm::dvec3&, const glm::dvec3& ) {}

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

/** Analytic normal used by the tests: 10 degrees off +Z, about the x axis. */
const glm::dvec3 TILTED_NORMAL =
  glm::dvec3( 0.0, std::sin( 0.17453292519943295 ), std::cos( 0.17453292519943295 ) );

/**
 * Two triangles sharing the edge v1-v2, meeting at about 16 degrees — inside
 * Reify()'s 40-degree smoothing cutoff, so the corners at the shared vertices
 * form one local group and the accumulation semantics are what decide the
 * answer.
 *
 * @param scale Uniform scale applied to every position. The point of the
 *   parameter is that a correct answer does not depend on it: face normals are
 *   area-weighted, so scaling multiplies every term by the same s^2 and leaves
 *   the direction alone. Mixing a unit vector into that sum is what makes the
 *   result move with s.
 * @param analyticOnFirst Give triangle 0 analytic corner normals.
 * @param analyticOnSecond Give triangle 1 analytic corner normals.
 * @return The geometry, ready to Reify.
 */
Geometry twoFaceFan( double scale, bool analyticOnFirst, bool analyticOnSecond ) {

  Geometry geometry;

  geometry.MakeVertex( glm::dvec3( 0, 0, 0 ) * scale );
  geometry.MakeVertex( glm::dvec3( 1, 0, 0 ) * scale );
  geometry.MakeVertex( glm::dvec3( 0, 1, 0 ) * scale );
  geometry.MakeVertex( glm::dvec3( 1, 1, 0.2 ) * scale );

  uint32_t first  = geometry.MakeTriangle( 0, 1, 2 );
  uint32_t second = geometry.MakeTriangle( 1, 3, 2 );

  glm::vec3 analytic( TILTED_NORMAL );

  if ( analyticOnFirst ) {

    geometry.SetCornerNormals( first, analytic, analytic, analytic );
  }

  if ( analyticOnSecond ) {

    geometry.SetCornerNormals( second, analytic, analytic, analytic );
  }

  return geometry;
}

/**
 * The shading normal Reify() emitted for the output vertex at `position`.
 *
 * @param geometry A reified geometry.
 * @param position The position to look for, in the geometry's own space.
 * @return The normal, or a zero vector when no output vertex sits there.
 */
glm::dvec3 normalAt( Geometry& geometry, const glm::dvec3& position ) {

  geometry.Reify();

  const std::vector< float >& stream = geometry.GetVertexStream();

  double tolerance = 1e-6 * std::max( 1.0, glm::length( position ) );

  for ( size_t vertex = 0, end = stream.size() / STRIDE; vertex < end; ++vertex ) {

    glm::dvec3 emitted(
      stream[ ( vertex * STRIDE ) ],
      stream[ ( vertex * STRIDE ) + 1 ],
      stream[ ( vertex * STRIDE ) + 2 ] );

    if ( glm::length( emitted - position ) > tolerance ) {

      continue;
    }

    return glm::dvec3(
      stream[ ( vertex * STRIDE ) + 3 ],
      stream[ ( vertex * STRIDE ) + 4 ],
      stream[ ( vertex * STRIDE ) + 5 ] );
  }

  return glm::dvec3( 0.0 );
}

double angleBetween( const glm::dvec3& a, const glm::dvec3& b ) {

  double cosine = glm::dot( glm::normalize( a ), glm::normalize( b ) );

  return std::acos( std::min( 1.0, std::max( -1.0, cosine ) ) ) * ( 180.0 / 3.14159265358979323846 );
}

/*
 * A smoothing group that mixes an analytic corner with an uncovered one must
 * fall back for the WHOLE group.
 *
 * With the per-corner choice, the accumulator at the shared vertex summed a
 * unit analytic vector with an unnormalized, area-weighted face normal. Only
 * one of the two terms carries area units, so the answer moved with the
 * model's scale: at scale 1 the two terms are comparable, at scale 10 the face
 * normal is 100x the analytic one and swamps it. Reverting the group-wide
 * fallback puts about 5 degrees between the two scales here, and the same
 * geometry in metres and millimetres shades differently.
 */
void testMixedGroupIsScaleInvariant() {

  printf( "\n=== mixed smoothing group falls back consistently ===\n" );

  Geometry unitScale = twoFaceFan( 1.0, true, false );
  Geometry tenScale  = twoFaceFan( 10.0, true, false );

  glm::dvec3 atUnit = normalAt( unitScale, glm::dvec3( 1, 0, 0 ) );
  glm::dvec3 atTen  = normalAt( tenScale, glm::dvec3( 10, 0, 0 ) );

  check( glm::length( atUnit ) > 0.5 && glm::length( atTen ) > 0.5,
         "both scales emit a vertex at the shared corner" );

  double drift = angleBetween( atUnit, atTen );

  printf( "      scale 1: (%.6f %.6f %.6f)   scale 10: (%.6f %.6f %.6f)   drift %.4f deg\n",
          atUnit.x, atUnit.y, atUnit.z, atTen.x, atTen.y, atTen.z, drift );

  check( drift < 1e-6,
         "the shared vertex's normal does not move with the model's scale" );

  // And it is specifically the pre-#667 answer: the area-weighted sum of the
  // two face normals, which is what "fall back for the whole group" means.
  glm::dvec3 faceA( 0, 0, 1 );
  glm::dvec3 faceB( -0.2, -0.2, 1 );
  glm::dvec3 expected = glm::normalize( faceA + faceB );

  printf( "      expected (%.6f %.6f %.6f)   got (%.6f %.6f %.6f)\n",
          expected.x, expected.y, expected.z, atUnit.x, atUnit.y, atUnit.z );

  check( angleBetween( atUnit, expected ) < 1e-4,
         "the mixed group's normal is the face-normal average, not a blend" );

  // The analytic normal is 10 degrees off +Z and the face average is about 8
  // degrees off it the other way, so a test that could not tell them apart
  // would be vacuous. It can: they are this far apart.
  printf( "      analytic answer would be %.2f deg away\n",
          angleBetween( atUnit, TILTED_NORMAL ) );

  check( angleBetween( atUnit, TILTED_NORMAL ) > 5.0,
         "the two candidate answers are far enough apart to distinguish" );
}

/*
 * The control: a group where EVERY corner is analytic still takes the analytic
 * path. Without this, "fall back when mixed" could be implemented as "always
 * fall back" and the mixed test above would still pass.
 */
void testFullyAnalyticGroupKeepsTheAnalyticNormal() {

  printf( "\n=== fully analytic group keeps the analytic normal ===\n" );

  Geometry geometry = twoFaceFan( 1.0, true, true );

  glm::dvec3 shared = normalAt( geometry, glm::dvec3( 1, 0, 0 ) );

  printf( "      got (%.6f %.6f %.6f), analytic (%.6f %.6f %.6f)\n",
          shared.x, shared.y, shared.z,
          TILTED_NORMAL.x, TILTED_NORMAL.y, TILTED_NORMAL.z );

  check( angleBetween( shared, TILTED_NORMAL ) < 1e-4,
         "both corners analytic: the emitted normal IS the analytic normal" );
}

/*
 * ReverseFace() inverts a face's outward direction. The analytic normals have
 * to be negated as well as reordered, or Reify() — which prefers the analytic
 * vector over the face normal — shades the reversed face against its own
 * winding. Reverting the negation leaves this check reporting 180 degrees.
 */
void testReverseFaceNegatesAnalyticNormals() {

  printf( "\n=== reversing a face reverses its analytic normals ===\n" );

  Geometry geometry;

  geometry.MakeVertex( glm::dvec3( 0, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 1, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 0, 1, 0 ) );

  uint32_t triangle = geometry.MakeTriangle( 0, 1, 2 );

  // Deliberately three DIFFERENT normals, so a reversal that negated without
  // reordering (or reordered without negating) is visible in the array, not
  // just in the shading.
  glm::vec3 n0( 0.0f, 0.0f, 1.0f );
  glm::vec3 n1( 0.1f, 0.0f, 0.99499f );
  glm::vec3 n2( 0.0f, 0.1f, 0.99499f );

  geometry.SetCornerNormals( triangle, n0, n1, n2 );

  glm::dvec3 before = normalAt( geometry, glm::dvec3( 0, 0, 0 ) );

  geometry.ClearReification();
  geometry.ReverseFace( triangle );

  check( geometry.corner_normals[ 0 ] == -n2 &&
         geometry.corner_normals[ 1 ] == -n1 &&
         geometry.corner_normals[ 2 ] == -n0,
         "corners 0 and 2 swap AND every corner is negated" );

  glm::dvec3 after = normalAt( geometry, glm::dvec3( 0, 0, 0 ) );

  printf( "      before (%.4f %.4f %.4f)   after (%.4f %.4f %.4f)   %.2f deg apart\n",
          before.x, before.y, before.z, after.x, after.y, after.z,
          angleBetween( before, after ) );

  check( angleBetween( before, after ) > 179.0,
         "the shading normal of a reversed analytic face is inverted with it" );

  // ReverseFaces() is the bulk form of the same operation and must agree.
  Geometry bulk;

  bulk.MakeVertex( glm::dvec3( 0, 0, 0 ) );
  bulk.MakeVertex( glm::dvec3( 1, 0, 0 ) );
  bulk.MakeVertex( glm::dvec3( 0, 1, 0 ) );

  uint32_t bulkTriangle = bulk.MakeTriangle( 0, 1, 2 );

  bulk.SetCornerNormals( bulkTriangle, n0, n1, n2 );
  bulk.ReverseFaces();

  check( bulk.corner_normals[ 0 ] == -n2 &&
         bulk.corner_normals[ 1 ] == -n1 &&
         bulk.corner_normals[ 2 ] == -n0,
         "ReverseFaces() does exactly what ReverseFace() does, per triangle" );
}

/*
 * The all-or-nothing rollback a throwing face triangulator takes
 * (ConwayGeometryProcessor::rollbackFace) truncates the triangle list without
 * going through DeleteTriangle. With a bare resize, `corner_normals` stayed at
 * the failed face's length, and the next face's MakeTriangle appended zeros
 * PAST the live range — so the next face's own corners read the failed face's
 * normals instead of its zeros. Reverting TruncateTriangles' sync leaves the
 * recovered face shaded with the dead face's analytic normals, which this
 * check sees as a non-zero corner.
 */
void testRollbackTruncatesCornerNormals() {

  printf( "\n=== rolling a face back truncates its corner normals ===\n" );

  Geometry geometry;

  geometry.MakeVertex( glm::dvec3( 0, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 1, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 0, 1, 0 ) );

  size_t trianglesBefore = geometry.triangles.size();

  // The face that "throws": one analytic triangle, then the rollback.
  uint32_t doomed = geometry.MakeTriangle( 0, 1, 2 );

  geometry.SetCornerNormals(
    doomed, glm::vec3( TILTED_NORMAL ), glm::vec3( TILTED_NORMAL ), glm::vec3( TILTED_NORMAL ) );

  geometry.TruncateTriangles( trianglesBefore );

  check( geometry.corner_normals.size() == geometry.triangles.size() * 3,
         "corner_normals is back to three per surviving triangle" );

  // The next face, with no analytic normals of its own.
  uint32_t recovered = geometry.MakeTriangle( 0, 1, 2 );

  glm::vec3 inherited =
    geometry.corner_normals.empty() ?
      glm::vec3( 0.0f ) : geometry.corner_normals[ recovered * 3 ];

  printf( "      recovered face's corner 0 reads (%.4f %.4f %.4f)\n",
          inherited.x, inherited.y, inherited.z );

  check( inherited == glm::vec3( 0.0f ),
         "a face triangulated after a rollback inherits no analytic normal" );

  // And the shading follows: the recovered face is flat, at its face normal.
  glm::dvec3 shading = normalAt( geometry, glm::dvec3( 0, 0, 0 ) );

  check( angleBetween( shading, glm::dvec3( 0, 0, 1 ) ) < 1e-6,
         "the recovered face shades on its face normal, not the dead one's" );
}

/*
 * The other path that rewrites the triangle list outside MakeTriangle /
 * DeleteTriangle. ExtractTriangles parses a whole index stream in and — note —
 * decrements EVERY triangle's vertex indices afterwards, its own record that
 * it expects to own the list; the underlying parse_vector appends, so the
 * count below is 2 + 1 rather than 1. Either way the array length stops
 * matching the triangle count, so the normals have to go.
 */
void testExtractTrianglesDropsCornerNormals() {

  printf( "\n=== re-extracting triangles drops the analytic normals ===\n" );

  Geometry geometry;

  geometry.MakeVertex( glm::dvec3( 0, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 1, 0, 0 ) );
  geometry.MakeVertex( glm::dvec3( 0, 1, 0 ) );
  geometry.MakeVertex( glm::dvec3( 1, 1, 0 ) );

  uint32_t first  = geometry.MakeTriangle( 0, 1, 2 );

  geometry.MakeTriangle( 1, 3, 2 );

  geometry.SetCornerNormals(
    first, glm::vec3( TILTED_NORMAL ), glm::vec3( TILTED_NORMAL ), glm::vec3( TILTED_NORMAL ) );

  // One triangle, in the 1-based text form the extractor's buffers carry.
  const std::string replacement = "((1,2,3))";

  conway::ParseBuffer buffer;

  buffer.resize( replacement.size() );

  memcpy( buffer.data(), replacement.data(), replacement.size() );

  geometry.ExtractTriangles( buffer );

  printf( "      %zu triangle(s), %zu corner normal(s)\n",
          geometry.triangles.size(), geometry.corner_normals.size() );

  check( geometry.triangles.size() == 3,
         "the parsed triangles landed in the list" );

  check( geometry.corner_normals.empty() ||
         geometry.corner_normals.size() == geometry.triangles.size() * 3,
         "the 3-per-triangle invariant survives a wholesale replacement" );

  check( geometry.corner_normals.empty(),
         "normals that described the old triangles are dropped, not reused" );
}

}  // namespace

int main() {

  testMixedGroupIsScaleInvariant();
  testFullyAnalyticGroupKeepsTheAnalyticNormal();
  testReverseFaceNegatesAnalyticNormals();
  testRollbackTruncatesCornerNormals();
  testExtractTrianglesDropsCornerNormals();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
