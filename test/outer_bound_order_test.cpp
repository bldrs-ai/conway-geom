/*
 * Ring-0 ordering tests for the earcut-based triangulators in mesh_utils.h /
 * geometry_utils.h.
 *
 * mapbox::earcut takes ring 0 of the polygon it is given as the OUTER
 * boundary and every ring after it as a hole. Nothing about the order a STEP
 * or IFC file lists a face's bounds in satisfies that contract - Onshape's
 * AP242 writer lists hole loops first on some faces - so a triangulator that
 * forwards `bounds` in arrival order ear-clips a hole as the outline and the
 * face ships a patch the size of that hole instead of its surface. Nothing is
 * logged; the only evidence is the geometry.
 *
 * The planar path (TriangulateBounds) has always swapped the OUTERBOUND-typed
 * bound to index 0. TriangulateBspline did not, which is what
 * bldrs-ai/test-models#65 is. What is pinned here:
 *
 *   1. The permutation itself - outerBoundIndex / outerBoundFirstOrder - is
 *      identity whenever the outer bound is already first. That is the claim
 *      the "no digest churn on models that do not have the defect" argument
 *      rests on, and it is cheap to state exactly.
 *
 *   2. TriangulateBspline covers its OUTER loop when that loop is listed
 *      LAST. This is the behavioural test, and it is the one that goes red on
 *      the unfixed header: measured there, the face emits 4.0 of 117.0 mm^2
 *      and a bounding box 2.828 across instead of 15.556 - i.e. exactly the
 *      hole, which is the defect's signature.
 *
 *   3. The same face with its bounds listed the other way round - outer
 *      first - emits the SAME area. The reorder is a permutation, not a
 *      repair: if the two orders disagreed, something other than ring
 *      selection would be moving.
 *
 * Assertions are on emitted area and extent rather than on deviation from the
 * surface, for the reason spherical_trim_test.cpp records: every triangle in
 * this test lies exactly on a planar patch whether the right ring was clipped
 * or not, so a deflection metric cannot see the defect at all.
 *
 * Verified by reverting rather than by reading: with mesh_utils.h alone rolled
 * back to the unfixed revision, three assertions go red and the run exits 1.
 * ONE ASSERTION HERE DOES NOT REDDEN - "the inner loop is still a hole, not
 * paved over" - and that is deliberate. Clipping the hole as the outline emits
 * 4 mm^2, which is below the un-holed 121 too, so no ring choice trips it. It
 * is kept as a characterisation of the opposite failure, the one #595's
 * full-coverage grid would cause if it were ever widened onto multi-bound
 * faces, and it would catch that change.
 *
 * Standalone by design: includes mesh_utils.h directly and links nothing but
 * the Logger stubs below, matching spherical_trim_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

// mesh_utils.h's error paths call these; the rest of the header is
// header-only. Defining them here keeps the test linking nothing.
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

using conway::geometry::Geometry;
using conway::geometry::IfcBound3D;
using conway::geometry::IfcBoundType;
using conway::geometry::IfcSurface;
using conway::geometry::TriangulateBspline;
using conway::geometry::outerBoundFirstOrder;
using conway::geometry::outerBoundIndex;

IfcBound3D makeBound(
    std::vector< glm::dvec3 > points, IfcBoundType type ) {

  IfcBound3D bound;

  bound.type         = type;
  bound.orientation  = true;
  bound.curve.points = std::move( points );

  return bound;
}

/**
 * An axis-aligned rectangle in the z = 0 plane, sampled with `perSide`
 * points on each side so the inverse evaluation has a real polyline to walk
 * rather than four corners.
 */
std::vector< glm::dvec3 > rectangle(
    double loX, double loY, double hiX, double hiY,
    size_t perSide, bool clockwise ) {

  const glm::dvec2 corners[ 4 ] = {
    { loX, loY }, { hiX, loY }, { hiX, hiY }, { loX, hiY } };

  std::vector< glm::dvec3 > points;

  points.reserve( perSide * 4 );

  for ( size_t side = 0; side < 4; ++side ) {

    const size_t from = clockwise ? ( 3 - side ) : side;
    const size_t to   = clockwise ? ( ( 3 - side ) + 3 ) % 4 : ( side + 1 ) % 4;

    for ( size_t step = 0; step < perSide; ++step ) {

      const double along = static_cast< double >( step ) / perSide;

      points.emplace_back(
        corners[ from ].x + ( corners[ to ].x - corners[ from ].x ) * along,
        corners[ from ].y + ( corners[ to ].y - corners[ from ].y ) * along,
        0.0 );
    }
  }

  return points;
}

/**
 * A bilinear (degree 1 x degree 1) rational patch spanning
 * [0, side] x [0, side] at z = 0.
 *
 * Deliberately the simplest surface tinynurbs will accept: the defect under
 * test is which RING earcut is handed, and a curved patch would only add
 * inverse-evaluation error on top of the thing being measured. The uv chart
 * of this patch is exactly (x / side, y / side), so any coverage shortfall in
 * the output is the triangulation's, not the solve's.
 */
IfcSurface makePlanarBsplinePatch( double side ) {

  IfcSurface surface;

  surface.transformation = glm::dmat4( 1.0 );
  surface.sameSense      = true;
  surface.sameSenseKnown = true;

  conway::geometry::BSpline& spline = surface.BSplineSurface;

  spline.Active  = true;
  spline.UDegree = 1;
  spline.VDegree = 1;

  spline.ControlPoints = {
    { glm::dvec3( 0.0, 0.0, 0.0 ), glm::dvec3( 0.0, side, 0.0 ) },
    { glm::dvec3( side, 0.0, 0.0 ), glm::dvec3( side, side, 0.0 ) } };

  spline.WeightPoints = { { 1.0, 1.0 }, { 1.0, 1.0 } };

  spline.UKnots        = { 0.0, 1.0 };
  spline.VKnots        = { 0.0, 1.0 };
  spline.UMultiplicity = { 2.0, 2.0 };
  spline.VMultiplicity = { 2.0, 2.0 };

  return surface;
}

struct Emitted {
  size_t triangles = 0;
  double area      = 0.0;
  double extent    = 0.0;
};

Emitted measure( const Geometry& geometry ) {

  Emitted emitted;

  emitted.triangles = geometry.triangles.size();

  glm::dvec3 lo(  std::numeric_limits< double >::max() );
  glm::dvec3 hi( -std::numeric_limits< double >::max() );

  for ( const auto& triangle : geometry.triangles ) {

    const glm::dvec3 a = geometry.vertices[ triangle.vertices[ 0 ] ];
    const glm::dvec3 b = geometry.vertices[ triangle.vertices[ 1 ] ];
    const glm::dvec3 c = geometry.vertices[ triangle.vertices[ 2 ] ];

    emitted.area += 0.5 * glm::length( glm::cross( b - a, c - a ) );

    for ( const glm::dvec3& point : { a, b, c } ) {
      lo = glm::min( lo, point );
      hi = glm::max( hi, point );
    }
  }

  emitted.extent =
    emitted.triangles > 0 ? glm::length( hi - lo ) : 0.0;

  return emitted;
}

// The patch, its outer trim and its inner trim. The outer ring encloses
// 11 x 11 and the hole 2 x 2, so a correct face carries 121 - 4 = 117 mm^2
// and spans sqrt( 11^2 + 11^2 ) = 15.556; a face that clipped the HOLE as its
// outline carries 4 and spans 2.828. A factor of 29 on area and 5.5 on
// extent, which is why no tolerance here has to be delicate.
constexpr double SIDE       = 13.0;
constexpr double OUTER_AREA = 11.0 * 11.0;
constexpr double HOLE_AREA  = 2.0 * 2.0;
constexpr double TRUE_AREA  = OUTER_AREA - HOLE_AREA;
constexpr double TRUE_EXTENT = 15.5563491861;

std::vector< IfcBound3D > holeFirstBounds() {

  std::vector< IfcBound3D > bounds;

  // The hole first, the outer bound second - the order Onshape writes and the
  // order the unfixed triangulator forwards verbatim.
  bounds.push_back(
    makeBound( rectangle( 5.0, 5.0, 7.0, 7.0, 6, true ), IfcBoundType::BOUND ) );
  bounds.push_back(
    makeBound( rectangle( 1.0, 1.0, 12.0, 12.0, 10, false ),
               IfcBoundType::OUTERBOUND ) );

  return bounds;
}

std::vector< IfcBound3D > outerFirstBounds() {

  std::vector< IfcBound3D > bounds;

  bounds.push_back(
    makeBound( rectangle( 1.0, 1.0, 12.0, 12.0, 10, false ),
               IfcBoundType::OUTERBOUND ) );
  bounds.push_back(
    makeBound( rectangle( 5.0, 5.0, 7.0, 7.0, 6, true ), IfcBoundType::BOUND ) );

  return bounds;
}

void testPermutationIsIdentityWhenOuterIsFirst() {

  printf( "=== the permutation ===\n" );

  const std::vector< IfcBound3D > outerFirst = outerFirstBounds();
  const std::vector< IfcBound3D > holeFirst  = holeFirstBounds();

  check( outerBoundIndex( outerFirst ) == 0,
         "outerBoundIndex finds the outer bound at 0 when it is listed first" );
  check( outerBoundIndex( holeFirst ) == 1,
         "outerBoundIndex finds the outer bound at 1 when a hole is listed first" );

  const std::vector< size_t > unchanged = outerBoundFirstOrder( outerFirst );

  check( unchanged.size() == 2 && unchanged[ 0 ] == 0 && unchanged[ 1 ] == 1,
         "outerBoundFirstOrder is the identity when the outer bound is first" );

  const std::vector< size_t > swapped = outerBoundFirstOrder( holeFirst );

  check( swapped.size() == 2 && swapped[ 0 ] == 1 && swapped[ 1 ] == 0,
         "outerBoundFirstOrder puts a late outer bound first" );

  // No bound typed OUTERBOUND at all: every triangulator keeps exactly the
  // order it was given, which is what the front ends that emit no
  // FACE_OUTER_BOUND relied on before this helper existed.
  std::vector< IfcBound3D > untyped = holeFirstBounds();

  untyped[ 1 ].type = IfcBoundType::BOUND;

  check( outerBoundIndex( untyped ) == untyped.size(),
         "outerBoundIndex reports none when no bound declares itself outer" );

  const std::vector< size_t > untouched = outerBoundFirstOrder( untyped );

  check( untouched.size() == 2 && untouched[ 0 ] == 0 && untouched[ 1 ] == 1,
         "outerBoundFirstOrder is the identity when no bound is typed outer" );

  check( outerBoundFirstOrder( {} ).empty(),
         "outerBoundFirstOrder handles a face with no bounds" );
}

void testBsplineCoversItsOuterLoopWhenListedLast() {

  printf( "=== TriangulateBspline, outer bound listed last ===\n" );

  IfcSurface surface = makePlanarBsplinePatch( SIDE );

  std::vector< IfcBound3D > bounds = holeFirstBounds();

  Geometry geometry;

  TriangulateBspline( geometry, bounds, surface, 1.0, SIDE );

  const Emitted emitted = measure( geometry );

  printf( "      %zu triangles, area %.3f (want %.3f), extent %.3f (want %.3f)\n",
          emitted.triangles, emitted.area, TRUE_AREA,
          emitted.extent, TRUE_EXTENT );

  check( emitted.triangles > 0,
         "the face emits geometry at all" );

  // 2% of the true area. The failing mode is not a near miss: the unfixed
  // header emits 4.0 here, which is 3.4% of the right answer.
  check( std::abs( emitted.area - TRUE_AREA ) < TRUE_AREA * 0.02,
         "the emitted area is the outer loop minus the hole" );

  check( std::abs( emitted.extent - TRUE_EXTENT ) < TRUE_EXTENT * 0.02,
         "the emitted geometry spans the outer loop, not the hole" );

  // The hole must still be a hole. Paving it over is the opposite failure and
  // the area test alone would only catch it at 4 mm^2 out of 117, so state it
  // as its own assertion against the un-holed area.
  check( emitted.area < OUTER_AREA - HOLE_AREA * 0.5,
         "the inner loop is still a hole, not paved over" );
}

void testBothBoundOrdersAgree() {

  printf( "=== the two bound orders agree ===\n" );

  IfcSurface surfaceA = makePlanarBsplinePatch( SIDE );
  IfcSurface surfaceB = makePlanarBsplinePatch( SIDE );

  std::vector< IfcBound3D > holeFirst  = holeFirstBounds();
  std::vector< IfcBound3D > outerFirst = outerFirstBounds();

  Geometry fromHoleFirst;
  Geometry fromOuterFirst;

  TriangulateBspline( fromHoleFirst, holeFirst, surfaceA, 1.0, SIDE );
  TriangulateBspline( fromOuterFirst, outerFirst, surfaceB, 1.0, SIDE );

  const Emitted holeFirstEmitted  = measure( fromHoleFirst );
  const Emitted outerFirstEmitted = measure( fromOuterFirst );

  printf( "      hole-first area %.3f, outer-first area %.3f\n",
          holeFirstEmitted.area, outerFirstEmitted.area );

  check( outerFirstEmitted.triangles > 0,
         "the outer-first face emits geometry (the control)" );

  check( std::abs( holeFirstEmitted.area - outerFirstEmitted.area ) <
           outerFirstEmitted.area * 0.02,
         "declaration order does not change the area the face carries" );
}

}  // namespace

int main() {

  testPermutationIsIdentityWhenOuterIsFirst();
  testBsplineCoversItsOuterLoopWhenListedLast();
  testBothBoundOrdersAgree();

  if ( failures > 0 ) {
    printf( "%d failure(s)\n", failures );
    return 1;
  }

  printf( "all ok\n" );
  return 0;
}
