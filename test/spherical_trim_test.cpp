/*
 * Coverage tests for TriangulateSphericalSurface's trimmed path in
 * mesh_utils.h.
 *
 * The property pinned here is one a deflection metric CANNOT see. Every
 * triangle these tests are about lies exactly on the sphere, so deviation is
 * zero whether the face is tessellated correctly or crumpled into a fifth of
 * its footprint. What goes wrong is WHERE the triangles are, so the assertions
 * are on coverage and area, not on closeness to the surface:
 *
 *   1. Azimuth coverage. A band that circles the sphere must put triangles at
 *      every azimuth. On Orbiter's button-head dome (ADVANCED_FACE #50714,
 *      express solid #964) 19 of 24 azimuth sectors ship empty while 111 of
 *      244 triangles pile into a single 15-degree sector.
 *
 *   2. Area against the analytic zone. The fold is what makes coverage and
 *      area disagree: that dome carries 25.196 mm^2 of a 27.95 mm^2 analytic
 *      zone - about 90% of the right area inside about 21% of the right
 *      footprint. Either number alone is passable; together they are the
 *      signature. A test that checked only area would have passed on the
 *      broken engine.
 *
 *   3. An inner loop stays a hole. #595's full-coverage grid is correct only
 *      when the face IS the whole sphere. A face with a real inner trim must
 *      never be paved over - doing so replaces a socket with surface, adding
 *      volume silently, which is the failure mode #595's own comment block
 *      warns about for the area predicate it rejected.
 *
 * The band is deliberately placed PERPENDICULAR to the surface's own
 * placement axis. That is not an awkward case invented for the test: a
 * SPHERICAL_SURFACE's axis is arbitrary in the file, and on #50714 it lands
 * roughly perpendicular to the screw, which is exactly why the fixed-equator
 * hemisphere split bisects an otherwise unremarkable dome lengthwise (its
 * boundary's normalized z spans [-0.8985, +0.8901] for a band 1.4mm deep on a
 * 3.17mm sphere). A test using an axis-aligned band would be green on the
 * broken engine and prove nothing.
 *
 * Each test was verified to fail on the pre-fix engine by running it there,
 * not by reading. See bldrs-ai/conway#644.
 *
 * Standalone by design: includes mesh_utils.h directly and links nothing but
 * the Logger stub below, matching ribbon_loft_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// mesh_utils.h's CDT error paths call this; the rest of the header is
// header-only. Defining it here keeps the test linking nothing, at the cost of
// swallowing the message - which is fine, because every assertion below is on
// the emitted geometry rather than on log output.
void Logger::logError( const char*, ... ) {}

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

using conway::geometry::Geometry;
using conway::geometry::IfcBound3D;
using conway::geometry::IfcBoundType;
using conway::geometry::IfcSurface;
using conway::geometry::TriangulateSphericalSurface;

/**
 * A circle of constant latitude about `axis`, sampled `count` times.
 *
 * `axis` is the band's OWN axis, which the callers deliberately set to
 * something other than the surface's placement z so the boundary straddles the
 * placement equator the way #50714's does.
 */
std::vector< glm::dvec3 > latitudeRing(
    const glm::dvec3& centre,
    double            radius,
    const glm::dvec3& axis,
    double            polarAngle,
    size_t            count,
    bool              reverse = false ) {

  const glm::dvec3 unitAxis = glm::normalize( axis );

  // Any vector not parallel to the axis gives a usable first basis vector.
  const glm::dvec3 seed =
    std::abs( unitAxis.x ) < 0.9 ? glm::dvec3( 1, 0, 0 ) : glm::dvec3( 0, 1, 0 );

  const glm::dvec3 basisU = glm::normalize( glm::cross( unitAxis, seed ) );
  const glm::dvec3 basisV = glm::cross( unitAxis, basisU );

  const double ringRadius = radius * std::sin( polarAngle );
  const double height     = radius * std::cos( polarAngle );

  std::vector< glm::dvec3 > points;

  points.reserve( count );

  for ( size_t step = 0; step < count; ++step ) {

    const size_t index = reverse ? ( count - 1 - step ) : step;
    const double angle = ( index * 2.0 * PI ) / count;

    points.push_back(
      centre + unitAxis * height +
      basisU * ( ringRadius * std::cos( angle ) ) +
      basisV * ( ringRadius * std::sin( angle ) ) );
  }

  return points;
}

IfcBound3D makeBound(
    std::vector< glm::dvec3 > points, IfcBoundType type ) {

  IfcBound3D bound;

  bound.type        = type;
  bound.orientation = true;
  bound.curve.points = std::move( points );

  return bound;
}

/** A sphere whose placement axis is +z, centred at the origin. */
IfcSurface makeSphere( double radius ) {

  IfcSurface surface;

  surface.transformation      = glm::dmat4( 1.0 );
  surface.SphericalSurface.Radius = radius;
  surface.SphericalSurface.Active = true;
  surface.sameSense           = true;
  surface.sameSenseKnown      = true;

  return surface;
}

/**
 * Unpaired half-edges, matching scripts/render-side halfedges.cjs exactly:
 * UNDIRECTED edges keyed on world position, unpaired = odd incidence count.
 * This is the health metric the Orbiter epic reports throughout, so the test
 * has to compute it the same way or it would pin a different property than the
 * one the corpus gate measures.
 */
size_t unpairedHalfEdges( const Geometry& geometry ) {

  std::map< std::pair< std::string, std::string >, size_t > counts;

  auto key = []( const glm::dvec3& p ) {

    char buffer[ 96 ];

    snprintf( buffer, sizeof( buffer ), "%.17g,%.17g,%.17g", p.x, p.y, p.z );

    return std::string( buffer );
  };

  for ( const auto& triangle : geometry.triangles ) {

    const glm::dvec3 v[ 3 ] = {
      geometry.vertices[ triangle.vertices[ 0 ] ],
      geometry.vertices[ triangle.vertices[ 1 ] ],
      geometry.vertices[ triangle.vertices[ 2 ] ] };

    for ( int edge = 0; edge < 3; ++edge ) {

      std::string a = key( v[ edge ] );
      std::string b = key( v[ ( edge + 1 ) % 3 ] );

      counts[ a < b ? std::make_pair( a, b ) : std::make_pair( b, a ) ] += 1;
    }
  }

  size_t unpaired = 0;

  for ( const auto& [ _, count ] : counts ) {

    if ( ( count % 2 ) == 1 ) {
      unpaired += count;
    }
  }

  return unpaired;
}

struct Coverage {
  size_t triangles       = 0;
  size_t emptySectors    = 0;
  double fullestSectorArea = 0.0;
  double area            = 0.0;
};

/**
 * Azimuth coverage and total area of the emitted geometry, measured about
 * `axis` through `centre` - the BAND's axis, so "every sector occupied" means
 * the band goes all the way round.
 *
 * Sector count is 24 (15 degrees), matching the histogram measured on the
 * shipped mesh in conway#644 so the numbers here are comparable to the ones on
 * that issue.
 */
Coverage measure(
    const Geometry&   geometry,
    const glm::dvec3& centre,
    const glm::dvec3& axis ) {

  constexpr size_t SECTORS = 24;

  const glm::dvec3 unitAxis = glm::normalize( axis );
  const glm::dvec3 seed =
    std::abs( unitAxis.x ) < 0.9 ? glm::dvec3( 1, 0, 0 ) : glm::dvec3( 0, 1, 0 );
  const glm::dvec3 basisU = glm::normalize( glm::cross( unitAxis, seed ) );
  const glm::dvec3 basisV = glm::cross( unitAxis, basisU );

  // Binned by AREA, not by triangle count. Adaptive refinement is entitled to
  // put more triangles where the chart is more distorted, so a count histogram
  // flags healthy meshes; the surface a sector actually carries is the thing
  // that must be even.
  std::vector< double > histogram( SECTORS, 0.0 );

  Coverage coverage;

  const size_t triangleCount = geometry.triangles.size();

  for ( size_t where = 0; where < triangleCount; ++where ) {

    const auto& triangle = geometry.triangles[ where ];

    const glm::dvec3 a = geometry.vertices[ triangle.vertices[ 0 ] ];
    const glm::dvec3 b = geometry.vertices[ triangle.vertices[ 1 ] ];
    const glm::dvec3 c = geometry.vertices[ triangle.vertices[ 2 ] ];

    const double triangleArea = 0.5 * glm::length( glm::cross( b - a, c - a ) );

    coverage.area += triangleArea;

    const glm::dvec3 centroid = ( a + b + c ) / 3.0 - centre;

    double angle =
      std::atan2( glm::dot( centroid, basisV ), glm::dot( centroid, basisU ) );

    if ( angle < 0.0 ) {
      angle += 2.0 * PI;
    }

    size_t sector = static_cast< size_t >( ( angle / ( 2.0 * PI ) ) * SECTORS );

    if ( sector >= SECTORS ) {
      sector = SECTORS - 1;
    }

    histogram[ sector ] += triangleArea;
    ++coverage.triangles;
  }

  for ( double bucket : histogram ) {

    if ( bucket <= 0.0 ) {
      ++coverage.emptySectors;
    }

    coverage.fullestSectorArea = std::max( coverage.fullestSectorArea, bucket );
  }

  return coverage;
}

/**
 * The dome: a band between two latitudes about an axis PERPENDICULAR to the
 * sphere's placement axis, i.e. #50714's geometry in miniature.
 *
 * Returned outer-then-inner so the caller can drop the inner loop to get the
 * single-bound variant.
 */
std::vector< IfcBound3D > domeBounds(
    double            radius,
    const glm::dvec3& bandAxis,
    double            outerPolar,
    double            innerPolar,
    size_t            samples ) {

  // Outer rim and inner opening, as two latitude circles about bandAxis.
  // Reversed winding on the inner one, as a real inner trim loop arrives.
  std::vector< IfcBound3D > bounds;

  bounds.push_back( makeBound(
    latitudeRing( glm::dvec3( 0 ), radius, bandAxis, outerPolar, samples, false ),
    IfcBoundType::OUTERBOUND ) );

  bounds.push_back( makeBound(
    latitudeRing( glm::dvec3( 0 ), radius, bandAxis, innerPolar, samples, true ),
    IfcBoundType::BOUND ) );

  return bounds;
}

/** Analytic area of the spherical zone between two polar angles. */
double zoneArea( double radius, double polarOuter, double polarInner ) {

  return 2.0 * PI * radius * radius *
         ( std::cos( polarInner ) - std::cos( polarOuter ) );
}

/**
 * Test 1 - a band perpendicular to the placement axis spreads its area evenly.
 *
 * A GUARD, not a reproduction, and the distinction was earned the hard way.
 * The first version of this test binned triangle COUNT per sector and looked
 * like a spectacular reproduction: 1160 of 1536 triangles in one 15-degree
 * sector on the shipped engine. Re-binned by AREA it is 1.06x an even share -
 * the surface was spread correctly all along and the count was measuring
 * adaptive refinement doing its job, which it is entitled to concentrate
 * wherever the chart is most distorted.
 *
 * So this configuration does not reproduce anything, and the assertion below
 * passes on both engines. It is kept because it is the property the migration
 * must not break, and because writing down why the obvious metric is wrong is
 * worth more than deleting it. The regime that DOES reproduce is test 2's.
 */
void testPerpendicularBandSpreadsAreaEvenly() {

  printf( "=== perpendicular band does not crumple ===\n" );

  constexpr double RADIUS = 3.17192226735973;
  constexpr double OUTER  = 1.00;   // normalized z reaches 0.8415
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );   // perpendicular to placement z

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  const std::vector< IfcBound3D > bounds =
    domeBounds( RADIUS, bandAxis, OUTER, INNER, 24 );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  printf( "      triangles=%zu emptySectors=%zu fullestSectorArea=%.4f "
          "area=%.4f evenShare=%.4f\n",
          coverage.triangles, coverage.emptySectors,
          coverage.fullestSectorArea, coverage.area, coverage.area / 24.0 );

  check( coverage.triangles > 0, "the band emits geometry at all" );

  // A band of constant width carries the same surface in every 15-degree
  // sector, so the fullest may exceed an even 1/24 share only by the sampling
  // margin. Threefold is far outside that and well inside the fold, which puts
  // essentially the whole face into a handful of sectors.
  check( coverage.fullestSectorArea < ( coverage.area / 24.0 ) * 3.0,
         "no sector carries more than 3x an even share of the area" );
}

/**
 * Test 2 - a deeper band is not dropped entirely. THE REPRODUCTION.
 *
 * Past about normalized z 0.95 the hemisphere-split path stops emitting
 * anything at all: 0 triangles, all 24 sectors empty, zero area, and no
 * diagnostic anywhere. Band geometry was swept to find this boundary rather
 * than guessed at.
 *
 * This is the same root cause as the real dome's fold - a chart seam placed
 * without reference to the boundary - in its extreme. It is not the identical
 * symptom: #50714 crumples (19 of 24 sectors empty, ~90% of the area in ~21%
 * of the footprint) rather than vanishing, and no configuration of concentric
 * circular rings reproduced that crumple. The dome's fold is verified against
 * the real model in the corpus gate, not here.
 */
void testDeepBandIsNotDropped() {

  printf( "=== deep band is not dropped ===\n" );

  constexpr double RADIUS = 3.17192226735973;
  constexpr double OUTER  = 1.25;   // normalized z reaches 0.9490
  constexpr double INNER  = 0.908;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  const std::vector< IfcBound3D > bounds =
    domeBounds( RADIUS, bandAxis, OUTER, INNER, 24 );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  const double expected = zoneArea( RADIUS, OUTER, INNER );

  printf( "      triangles=%zu emptySectors=%zu area=%.4f expected=%.4f\n",
          coverage.triangles, coverage.emptySectors, coverage.area, expected );

  check( coverage.triangles > 0, "the deep band emits geometry at all" );

  check( coverage.emptySectors == 0,
         "the deep band reaches every azimuth" );

  check( coverage.area > expected * 0.90,
         "the deep band carries its analytic area" );
}

/**
 * Test 3 - area stays bounded by the analytic zone.
 *
 * A companion guard rather than a reproduction: the shipped engine already
 * passes this on the test-1 band (ratio 0.998), because folding conserves area
 * while destroying the footprint. It is here so that a "fix" cannot buy
 * coverage by emitting overlapping surface, which is the failure mode the
 * intermediate designs on conway-geom#190 kept producing.
 */
void testAreaStaysBoundedByAnalyticZone() {

  printf( "=== area stays bounded by the analytic zone ===\n" );

  constexpr double RADIUS = 3.17192226735973;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  const std::vector< IfcBound3D > bounds =
    domeBounds( RADIUS, bandAxis, OUTER, INNER, 24 );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  const double expected = zoneArea( RADIUS, OUTER, INNER );

  printf( "      area=%.4f expected=%.4f ratio=%.4f\n",
          coverage.area, expected, coverage.area / expected );

  check( coverage.area < expected * 1.02,
         "area does not exceed the analytic zone (no overlap)" );

  check( coverage.area > expected * 0.90,
         "area is not far short of the analytic zone" );
}

/**
 * Test 4 - the inner loop stays a hole.
 *
 * Guards #595's full-coverage grid from being widened to reach this case: that
 * grid would pave the socket, and the area would jump toward the whole sphere.
 */
void testInnerLoopIsNotPavedOver() {

  printf( "=== inner trim loop stays a hole ===\n" );

  constexpr double RADIUS = 3.17192226735973;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  const std::vector< IfcBound3D > bounds =
    domeBounds( RADIUS, bandAxis, OUTER, INNER, 24 );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  const double wholeSphere   = 4.0 * PI * RADIUS * RADIUS;
  const double band          = zoneArea( RADIUS, OUTER, INNER );
  const double capInsideHole = zoneArea( RADIUS, INNER, 0.0 );

  printf( "      area=%.4f band=%.4f capInsideHole=%.4f wholeSphere=%.4f\n",
          coverage.area, band, capInsideHole, wholeSphere );

  check( coverage.area < band + ( capInsideHole * 0.5 ),
         "the hole's cap is not filled in" );

  check( coverage.area < wholeSphere * 0.5,
         "the face is not the whole sphere" );
}

/**
 * Test 5 - a patch whose boundary REACHES the pole triangulates cleanly.
 *
 * The pole is the one place the (theta, phi) chart is not injective: it is a
 * single 3D point but an entire line (every theta at |phi| = pi/2). Right_Hand
 * carries 8+ spherical faces whose boundary reaches it, so this is corpus
 * behaviour rather than a hypothetical.
 *
 * The shape is a lune from the equator up to just short of the pole - two
 * meridian arcs and an equatorial arc, which is the Right_Hand shape. The apex
 * is placed 1e-5 radians SHORT of the pole deliberately: a sample exactly on
 * the axis trips the existing degenerate-sample guard and never reaches the
 * chart, so an on-axis apex would test nothing.
 *
 * Two assertions, and they do different jobs - measured, not assumed:
 *
 *   - AREA is the red proof. The legacy dual-hemisphere path overshoots this
 *     lune by 2.3% (0.818 against an analytic 0.800), which is overlap; the
 *     chart lands at 0.6% over, inside the sampling margin.
 *
 *   - UNPAIRED HALF-EDGES is a guard, and passes on both engines (60 against a
 *     60-edge boundary). An earlier draft claimed the legacy path tore this
 *     lune to 1794 unpaired edges; that number came from an ad-hoc DIRECTED
 *     edge counter, and the undirected odd-incidence metric the epic actually
 *     reports says 60. The guard is kept because a chart that failed to stitch
 *     at the pole would break it, but it is not what makes this test red.
 */
void testPoleReachingLuneTriangulatesCleanly() {

  printf( "=== pole-reaching lune triangulates cleanly ===\n" );

  constexpr double RADIUS = 1.0;
  constexpr double LO_THETA = 0.3;
  constexpr double HI_THETA = 1.1;
  constexpr int    STEPS    = 20;

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  std::vector< glm::dvec3 > points;

  auto onSphere = []( double theta, double phi ) {

    return glm::dvec3(
      std::cos( phi ) * std::cos( theta ),
      std::cos( phi ) * std::sin( theta ),
      std::sin( phi ) );
  };

  for ( int step = 0; step <= STEPS; ++step ) {
    points.push_back(
      onSphere( LO_THETA + ( HI_THETA - LO_THETA ) * step / STEPS, 0.0 ) );
  }

  for ( int step = 1; step < STEPS; ++step ) {
    points.push_back( onSphere( HI_THETA, ( PI / 2 ) * step / STEPS ) );
  }

  // Just short of the pole - see the note above on why not exactly on it.
  points.push_back(
    onSphere( ( LO_THETA + HI_THETA ) / 2, ( PI / 2 ) - 1e-5 ) );

  for ( int step = STEPS - 1; step >= 1; --step ) {
    points.push_back( onSphere( LO_THETA, ( PI / 2 ) * step / STEPS ) );
  }

  std::vector< IfcBound3D > bounds;

  bounds.push_back( makeBound( points, IfcBoundType::OUTERBOUND ) );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const size_t unpaired  = unpairedHalfEdges( geometry );
  const size_t boundary  = points.size();

  // A spherical lune from equator to pole has area R^2 * (longitude span).
  const double expected = ( HI_THETA - LO_THETA ) * RADIUS * RADIUS;

  const Coverage coverage =
    measure( geometry, glm::dvec3( 0 ), glm::dvec3( 0, 0, 1 ) );

  printf( "      triangles=%zu unpaired=%zu boundaryEdges=%zu "
          "area=%.6f expected=%.6f\n",
          coverage.triangles, unpaired, boundary,
          coverage.area, expected );

  check( coverage.triangles > 0, "the pole-reaching lune emits geometry" );

  check( unpaired == boundary,
         "unpaired half-edges equal the boundary edge count - no tear" );

  check( coverage.area < expected * 1.02 && coverage.area > expected * 0.95,
         "the lune carries its analytic area" );
}

/**
 * Test 6 - a collapsed inner loop fails the unwrap to the legacy path.
 *
 * If dedup takes a supplied loop below three points, erasing it and
 * triangulating the survivors pavesover the hole and adds volume silently -
 * the exact failure that ruled #595's full-coverage grid out of this job. The
 * unwrap must decline the whole face instead, so the "can never make output
 * worse" contract is literally true rather than true-except-here.
 *
 * WHAT THIS PINS, precisely, because the obvious reading is wrong: it pins
 * that the FALLBACK ENGAGED, not that the hole survives. Measured, the legacy
 * dual-hemisphere path paves this hole too (area 2.8708 against a 1.2578 band
 * and a 2.8884 band-plus-hole). So a collapsed inner trim produces bad
 * geometry either way; what this change buys is that the new chart does not
 * own that failure and the face behaves exactly as it did before this path
 * existed. Paving on a collapsed trim is a pre-existing defect and is tracked
 * separately.
 *
 * The triangle count is the pin because it is the only thing that separates
 * the two paths here (704 for the chart, 865 for legacy) - the areas agree to
 * four decimals. It is a golden: if legacy's refinement changes, update it.
 *
 * Red-proven: without the collapse check this reports 704.
 */
void testCollapsedInnerLoopFallsBackToLegacy() {

  printf( "=== collapsed inner loop falls back to legacy ===\n" );

  constexpr double RADIUS = 1.0;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;
  constexpr size_t LEGACY_TRIANGLE_COUNT = 865;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  std::vector< IfcBound3D > bounds;

  bounds.push_back( makeBound(
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, OUTER, 24, false ),
    IfcBoundType::OUTERBOUND ) );

  // A real inner trim that has degenerated: three points that all collapse to
  // one under the unwrap's 1e-12 parameter dedup.
  const glm::dvec3 seed =
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, INNER, 24, true )[ 0 ];

  bounds.push_back( makeBound(
    { seed, seed + glm::dvec3( 1e-15, 0, 0 ), seed + glm::dvec3( 0, 1e-15, 0 ) },
    IfcBoundType::BOUND ) );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  printf( "      triangles=%zu legacyExpected=%zu\n",
          geometry.triangles.size(), LEGACY_TRIANGLE_COUNT );

  check( geometry.triangles.size() == LEGACY_TRIANGLE_COUNT,
         "the face took the legacy path, not the unwrap" );
}

/**
 * Tests 7 and 8 - a bound lost anywhere before the CDT fails the unwrap.
 *
 * Test 6 covers a loop collapsing under the unwrap's own 1e-12 dedup. Review
 * on bldrs-ai/conway-geom#191 found two more places a supplied bound can go
 * missing while the survivors still triangulate, both of which paved the hole:
 *
 *   7. AT COLLECTION - a non-finite sample is skipped, leaving the loop under
 *      three points, so it is never collected at all. Counting collapses among
 *      the loops that survived cannot see this; the casualty is already gone.
 *
 *   8. AT THE CONSTRAINT WELD - triangulateUnwrappedLoops welds at 1e-9 while
 *      the dedup runs at 1e-12, so points that legitimately survive collection
 *      can merge inside the helper. A loop merged onto one vertex emits no
 *      edges and vanishes from the constraints while `cdtEdges.size() < 3`
 *      still passes on the outer loop alone.
 *
 * Unlike test 6, the legacy path does NOT pave either of these - it drops the
 * face - so here the assertion can be the property itself rather than a
 * triangle-count golden. Measured before the fix: 704 triangles at area 2.8703
 * for test 7, 768 at 2.8688 for test 8, both of which are band-plus-hole.
 */
void testNonFiniteSampleInInnerLoopDoesNotPave() {

  printf( "=== non-finite sample in inner loop does not pave ===\n" );

  constexpr double RADIUS = 1.0;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  std::vector< IfcBound3D > bounds;

  bounds.push_back( makeBound(
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, OUTER, 24, false ),
    IfcBoundType::OUTERBOUND ) );

  const std::vector< glm::dvec3 > innerRing =
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, INNER, 24, true );

  bounds.push_back( makeBound(
    { innerRing[ 0 ], innerRing[ 8 ],
      glm::dvec3( std::numeric_limits< double >::quiet_NaN(), 0.0, 0.0 ) },
    IfcBoundType::BOUND ) );

  // The legacy path throws its own non-finite guard on this input, which
  // AddFaceToGeometry catches per face in production. Either outcome is fine
  // here; what must not happen is a paved hole.
  try {
    TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );
  } catch ( const std::exception& ) {
  }

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  const double band          = zoneArea( RADIUS, OUTER, INNER );
  const double capInsideHole = zoneArea( RADIUS, INNER, 0.0 );

  printf( "      triangles=%zu area=%.4f band=%.4f band+hole=%.4f\n",
          coverage.triangles, coverage.area, band, band + capInsideHole );

  check( coverage.area < band + ( capInsideHole * 0.5 ),
         "the hole is not paved over when a sample is non-finite" );
}

void testWeldMergedInnerLoopDoesNotPave() {

  printf( "=== weld-merged inner loop does not pave ===\n" );

  constexpr double RADIUS = 1.0;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  IfcSurface surface = makeSphere( RADIUS );
  Geometry   geometry;

  std::vector< IfcBound3D > bounds;

  bounds.push_back( makeBound(
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, OUTER, 24, false ),
    IfcBoundType::OUTERBOUND ) );

  // Six points 1e-10 apart: above the 1e-12 dedup, below the helper's 1e-9
  // weld. They survive collection and then merge inside the CDT setup.
  const glm::dvec3 seed =
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, INNER, 24, true )[ 0 ];

  std::vector< glm::dvec3 > inner;

  for ( int step = 0; step < 6; ++step ) {
    inner.push_back( seed + glm::dvec3( 1e-10 * step, 1e-10 * step, 0.0 ) );
  }

  bounds.push_back( makeBound( inner, IfcBoundType::BOUND ) );

  TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

  const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

  const double band          = zoneArea( RADIUS, OUTER, INNER );
  const double capInsideHole = zoneArea( RADIUS, INNER, 0.0 );

  printf( "      triangles=%zu area=%.4f band=%.4f band+hole=%.4f\n",
          coverage.triangles, coverage.area, band, band + capInsideHole );

  check( coverage.area < band + ( capInsideHole * 0.5 ),
         "the hole is not paved over when a loop merges in the weld" );
}

/**
 * Test 9 - an inner loop that welds to only TWO distinct vertices is a slit,
 * not a boundary, and must not be triangulated across.
 *
 * The fourth and last member of the loop-guard family, and the only one that
 * was a genuine REGRESSION rather than a defect shared with the legacy path.
 * The earlier guards asked whether a loop reached the CDT at all, which a
 * two-vertex loop does: it emits one open edge. The CDT then triangulates
 * straight across that slit and fills the hole.
 *
 * Measured before the fix: 832 triangles at area 2.8692 against a 1.2578 band
 * and a 2.8884 band-plus-hole - the hole filled. The pre-#191 legacy path
 * emits 0 triangles for the same input, so unlike tests 6 to 8 this one was
 * strictly worse than doing nothing.
 *
 * Both spellings are pinned because they enter the guard by different routes:
 * `A,A,B,B` has four collected points and dedups to two, while `A,B,B` has the
 * bare minimum three and dedups to two. An implementation that only counted
 * collected points would pass the first and fail the second.
 */
void testTwoVertexInnerLoopIsNotTriangulatedAcross() {

  printf( "=== two-vertex inner loop is not triangulated across ===\n" );

  constexpr double RADIUS = 1.0;
  constexpr double OUTER  = 1.00;
  constexpr double INNER  = 0.737;

  const glm::dvec3 bandAxis( 1.0, 0.0, 0.0 );

  const std::vector< glm::dvec3 > innerRing =
    latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, INNER, 24, true );

  const glm::dvec3 pointA = innerRing[ 0 ];
  const glm::dvec3 pointB = innerRing[ 12 ];

  const double band          = zoneArea( RADIUS, OUTER, INNER );
  const double capInsideHole = zoneArea( RADIUS, INNER, 0.0 );

  auto runCase =
    [ & ]( const char* what, const std::vector< glm::dvec3 >& inner ) {

      IfcSurface surface = makeSphere( RADIUS );
      Geometry   geometry;

      std::vector< IfcBound3D > bounds;

      bounds.push_back( makeBound(
        latitudeRing( glm::dvec3( 0 ), RADIUS, bandAxis, OUTER, 24, false ),
        IfcBoundType::OUTERBOUND ) );

      bounds.push_back( makeBound( inner, IfcBoundType::BOUND ) );

      TriangulateSphericalSurface( geometry, bounds, surface, RADIUS * 2.0 );

      const Coverage coverage = measure( geometry, glm::dvec3( 0 ), bandAxis );

      printf( "      %-12s triangles=%zu area=%.4f band=%.4f band+hole=%.4f\n",
              what, coverage.triangles, coverage.area, band,
              band + capInsideHole );

      check( coverage.area < band + ( capInsideHole * 0.5 ),
             std::string( "the hole is not filled across a two-vertex loop (" ) +
               what + ")" );
    };

  runCase( "A,A,B,B", { pointA, pointA, pointB, pointB } );
  runCase( "A,B,B",   { pointA, pointB, pointB } );
}

}  // namespace

int main() {

  testPerpendicularBandSpreadsAreaEvenly();
  testDeepBandIsNotDropped();
  testAreaStaysBoundedByAnalyticZone();
  testInnerLoopIsNotPavedOver();
  testPoleReachingLuneTriangulatesCleanly();
  testCollapsedInnerLoopFallsBackToLegacy();
  testNonFiniteSampleInInnerLoopDoesNotPave();
  testWeldMergedInnerLoopDoesNotPave();
  testTwoVertexInnerLoopIsNotTriangulatedAcross();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
