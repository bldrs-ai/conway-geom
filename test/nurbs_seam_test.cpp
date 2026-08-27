/*
 * Soundness tests for the two closed-surface paths in mesh_utils.h: the seam
 * crossing in RationalNurbsInverseMethod's descent, and the seam-axis gate on
 * tryFullCoverageSeamGrid.
 *
 * Both behaviours are reachable only on a surface CLOSED in a parameter, which
 * the STEP corpus barely exercises - Orbiter has one such face in 721, and
 * Right_Hand one in 93 - so neither is pinned by the regression digests. Each
 * test below fails if its fix is reverted; that was verified by reverting, not
 * by reading. See bldrs-ai/conway#611.
 *
 * Standalone by design: it includes mesh_utils.h directly and links nothing,
 * so it builds with a plain compiler and does not depend on the
 * conway_geom_native_tests target (which does not currently build).
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

/**
 * A tube of radius `radius` about the z axis, closed in u (around) and open in
 * v (along), as a degree-2 rational NURBS - the standard nine-control-point
 * circle, so the weights are genuinely non-uniform and the rational path is
 * exercised rather than the polynomial one.
 *
 * u runs [-pi, pi] so the seam sits at both domain ends, as it does on
 * `B_SPLINE_SURFACE_WITH_KNOTS #1553`.
 */
tinynurbs::RationalSurface3d makeClosedTube( double radius, double height ) {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 1;

  // Nine control points, first == last, is the rational circle.
  const double corner = radius;
  const std::vector< glm::dvec2 > ring = {
    {  corner,  0.0 }, {  corner,  corner }, {  0.0,  corner },
    { -corner,  corner }, { -corner,  0.0 }, { -corner, -corner },
    {  0.0, -corner }, {  corner, -corner }, {  corner,  0.0 } };

  const double diagonal = std::sqrt( 2.0 ) / 2.0;
  const std::vector< double > ringWeights = {
    1.0, diagonal, 1.0, diagonal, 1.0, diagonal, 1.0, diagonal, 1.0 };

  const size_t rows = ring.size();
  const size_t cols = 2;

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( size_t row = 0; row < rows; ++row ) {
    for ( size_t col = 0; col < cols; ++col ) {
      points.push_back( glm::dvec3( ring[ row ].x, ring[ row ].y,
                                    col == 0 ? 0.0 : height ) );
      weights.push_back( ringWeights[ row ] );
    }
  }

  surface.control_points = tinynurbs::array2( rows, cols, points );
  surface.weights        = tinynurbs::array2( rows, cols, weights );

  // Degree-2, four spans, clamped: knots -pi .. pi with interior multiplicity 2.
  surface.knots_u = { -PI, -PI, -PI, -PI / 2, -PI / 2, 0.0, 0.0,
                      PI / 2, PI / 2, PI, PI, PI };
  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  return surface;
}

/**
 * The same closed tube, but with the u knots crowded into two narrow spans at
 * the ends so one span covers a quarter of the circumference over most of the
 * domain.
 *
 * That makes |dS/du| small there - about 0.26 against 1.0 on the uniform tube
 * - and a Gauss-Newton step is |e| / |dS/du|, so a residual of order the
 * diameter produces a step several times the half period. On a uniformly
 * parameterised circle that cannot happen: the step is bounded by
 * 2r / r = 2 < pi, which is why the uniform tube alone does not exercise the
 * wide-step path at all.
 */
tinynurbs::RationalSurface3d makeSkewedClosedTube( double radius, double height ) {

  tinynurbs::RationalSurface3d surface = makeClosedTube( radius, height );

  surface.knots_u = { -PI, -PI, -PI, -PI + 0.05, -PI + 0.05, -PI + 0.10,
                      -PI + 0.10, PI - 0.10, PI - 0.10, PI, PI, PI };

  return surface;
}

/** The same tube, but also closed in v by making the two v ends coincide. */
tinynurbs::RationalSurface3d makeDoublyClosedTube( double radius ) {

  tinynurbs::RationalSurface3d surface = makeClosedTube( radius, 0.0 );

  // height 0 already makes column 0 and column 1 identical, so the surface is
  // closed in v as well as u - degenerate as geometry, but it is the CLOSURE
  // FLAGS the gate reads, and this sets both.
  return surface;
}

// ---------------------------------------------------------------------------

/**
 * Armijo must never accept a trial that INCREASES the residual, and on a
 * closed axis that turns entirely on which displacement the predicted decrease
 * is measured over.
 *
 * This tests `trialDisplacement` directly rather than through a solve. The
 * defect only appears on a step wider than half a period, and `operator()`
 * makes those hard to reach: it seeds from the nearer of the grid sample and
 * the previous solution, which bounds the residual, and a Gauss-Newton step is
 * |e| / |dS/du| - so on a uniformly parameterised closed axis it cannot exceed
 * 2r / r = 2 < pi. Verified by measurement, not assumed: 1,024 seed/query
 * pairs over a closed tube, including one with its u knots crowded so
 * |dS/du| falls to ~0.26, produce not one negative predicted decrease. An
 * end-to-end test would therefore have passed with the defect in place, which
 * is worse than no test. The predicate is what is unsound, so the predicate is
 * what is pinned.
 */
void testWideStepAcrossSeamKeepsArmijoSound() {

  printf( "Armijo displacement across a seam\n" );

  tinynurbs::RationalSurface3d surface = makeClosedTube( 1.0, 4.0 );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  check( solve.closedU_ && !solve.closedV_, "tube is closed in u, open in v" );

  const double period = solve.max_extent.x - solve.min_extent.x;

  check( std::abs( period - ( 2.0 * PI ) ) < 1e-9, "u period is 2 pi" );

  // A step of 0.6 of a period: wider than half, so the period reduction that
  // this replaced would map it to 0.6 - 1.0 = -0.4 of a period, i.e. REVERSE
  // its sign.
  const double wide = 0.6 * period;

  const glm::dvec2 from( solve.min_extent.x + ( 0.2 * period ), 2.0 );
  const glm::dvec2 requested( wide, 0.0 );

  // Where the descent would actually land, brought back into the domain.
  const glm::dvec2 to(
    solve.min_extent.x +
      std::fmod( ( from.x - wide ) - solve.min_extent.x + ( 2.0 * period ), period ),
    from.y );

  const glm::dvec2 displacement = solve.trialDisplacement( from, to, requested );

  check( std::abs( displacement.x - requested.x ) < 1e-12,
         "closed axis credits the requested step, not its shortest representative" );

  // The shortest representative, which is what the first version used.
  const double reduced =
    ( from.x - to.x ) - ( period * std::round( ( from.x - to.x ) / period ) );

  check( ( reduced * requested.x ) < 0.0,
         "and the shortest representative really does reverse the sign here" );

  // `jte` is J^T e. `deltaUV` solves an SPD system, so dot( requested, jte ) is
  // positive by construction; the predicted decrease has to stay that way.
  const glm::dvec2 jte( 1.0, 0.0 );

  check( glm::dot( requested, jte ) > 0.0, "requested step is a descent direction" );
  check( glm::dot( displacement, jte ) > 0.0,
         "so the predicted decrease stays positive, and the Armijo threshold "
         "stays BELOW the current residual" );
  check( glm::dot( glm::dvec2( reduced, 0.0 ), jte ) < 0.0,
         "whereas the reduced form makes it negative - the accepted-increase case" );

  // An OPEN axis must still credit only what the clamp left, which is what
  // stopped the descent stalling at a domain corner.
  const glm::dvec2 vFrom( 0.0, solve.max_extent.y );
  const glm::dvec2 vRequested( 0.0, -1.0 );
  const glm::dvec2 vTo( 0.0, solve.max_extent.y );

  check( solve.trialDisplacement( vFrom, vTo, vRequested ).y == 0.0,
         "open axis credits only the movement the clamp left" );
}

/**
 * A general net rather than a pin for the fix above: the descent should never
 * hand back a uv further from the query than the seed it started from. It
 * passes with or without the displacement fix, for the reachability reason
 * given above, and is kept because it would catch a regression that makes the
 * descent diverge for any other reason.
 */
void testDescentNeverWorsensFromItsSeed() {

  printf( "descent never returns worse than its seed\n" );

  const double radius = 1.0;
  const double height = 4.0;

  tinynurbs::RationalSurface3d surface = makeSkewedClosedTube( radius, height );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  check( solve.closedU_, "skewed tube is detected closed in u" );

  double worstExcess = 0.0;
  size_t violations  = 0;

  for ( int seedStep = 0; seedStep < 32; ++seedStep ) {

    const double seedU = -PI + ( ( 2.0 * PI * seedStep ) / 32.0 );

    for ( int queryStep = 0; queryStep < 32; ++queryStep ) {

      const double queryU = -PI + ( ( 2.0 * PI * queryStep ) / 32.0 );

      solve.resetContinuity();

      const glm::dvec3 seedPoint = solve.evaluator.point( seedU, height * 0.5 );

      const glm::dvec2 seedUV = solve( seedPoint );

      const glm::dvec3 queryPoint = solve.evaluator.point( queryU, height * 0.5 );

      const double fromSeed =
        glm::distance( solve.evaluator.point( seedUV.x, seedUV.y ), queryPoint );

      const glm::dvec2 solved = solve( queryPoint );

      const double fromSolved =
        glm::distance( solve.evaluator.point( solved.x, solved.y ), queryPoint );

      if ( fromSolved > fromSeed ) {

        ++violations;
        worstExcess = std::max( worstExcess, fromSolved - fromSeed );
      }
    }
  }

  check( violations == 0,
         "1024 seed/query pairs, none ending further away than its seed "
         "(violations " + std::to_string( violations ) + ", worst excess " +
           std::to_string( worstExcess ) + ")" );
}

/**
 * Coincident Cartesian control rows do NOT make a rational surface closed.
 *
 * The endpoint curves share the surface's knots and degree, so with control
 * points P and weights w each is
 *
 *   C( t ) = sum_j N_j( t ) w_j P_j / sum_j N_j( t ) w_j
 *
 * and two rows of identical P with different w are different curves - they
 * agree only where the basis is interpolating, which is to say at their ends.
 * Believing them closed is not cosmetic: closedU_ gates both the seam crossing
 * in the descent and the full-coverage grid, so a false positive makes the
 * solver wrap across a real geometric discontinuity.
 */
void testClosureRequiresMatchingWeights() {

  printf( "closure detection accounts for rational weights\n" );

  tinynurbs::RationalSurface3d honest = makeClosedTube( 1.0, 4.0 );

  check( conway::geometry::RationalNurbsInverseMethod( honest ).closedU_,
         "matching weights across the seam: closed in u" );

  // Same control points, one seam weight disturbed. The rows still coincide
  // in P, so a Cartesian-only test still calls this closed.
  tinynurbs::RationalSurface3d disturbed = makeClosedTube( 1.0, 4.0 );

  const size_t lastRow = disturbed.weights.rows() - 1;

  for ( size_t col = 0; col < disturbed.weights.cols(); ++col ) {
    disturbed.weights( lastRow, col ) = disturbed.weights( lastRow, col ) * 0.5;
  }

  // The control points are still identical row 0 vs row last - so the ONLY
  // thing that can distinguish these surfaces is the weights.
  bool cartesianIdentical = true;

  for ( size_t col = 0; col < disturbed.control_points.cols(); ++col ) {
    cartesianIdentical = cartesianIdentical &&
      ( disturbed.control_points( 0, col ) ==
        disturbed.control_points( lastRow, col ) );
  }

  check( cartesianIdentical,
         "control points still coincide, so only the weights differ" );

  check( !conway::geometry::RationalNurbsInverseMethod( disturbed ).closedU_,
         "mismatched weights across the seam: NOT closed in u" );
}

// ---------------------------------------------------------------------------

/**
 * Build the boundary a seam-wrapping face has: a seam curve walked down one
 * side and back up the other, joined by two rings. This is the shape
 * `EDGE_LOOP #9507` has, reduced to the structure tryFullCoverageSeamGrid
 * actually reads.
 *
 * The two seam legs carry IDENTICAL 3D points and DIFFERENT uv - u = uMax on
 * one side, u = uMin on the other - which is what makes them one seam rather
 * than two curves.
 */
size_t buildSeamLoop(
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex >& mesh,
    double radius,
    double height,
    size_t seamCount,
    size_t ringCount ) {

  const auto onSurface =
    [ & ]( double angle, double along ) {
      return glm::dvec3(
        radius * std::cos( angle ), radius * std::sin( angle ), along );
    };

  // seamA: u = +pi, v ascending.
  for ( size_t at = 0; at < seamCount; ++at ) {

    const double along = ( height * at ) / double( seamCount - 1 );

    mesh.makeVertex( { onSurface( PI, along ), glm::dvec2( PI, along ) } );
  }

  // ringA: the far ring, interior points only, u descending from +pi to -pi.
  for ( size_t at = 1; at <= ringCount; ++at ) {

    const double angle = PI - ( ( 2.0 * PI * at ) / double( ringCount + 1 ) );

    mesh.makeVertex( { onSurface( angle, height ), glm::dvec2( angle, height ) } );
  }

  // seamB: u = -pi, the SAME points as seamA in reverse.
  for ( size_t at = 0; at < seamCount; ++at ) {

    const double along =
      ( height * ( seamCount - 1 - at ) ) / double( seamCount - 1 );

    mesh.makeVertex( { onSurface( PI, along ), glm::dvec2( -PI, along ) } );
  }

  // ringB: the near ring, interior points only.
  for ( size_t at = 1; at <= ringCount; ++at ) {

    const double angle = -PI + ( ( 2.0 * PI * at ) / double( ringCount + 1 ) );

    mesh.makeVertex( { onSurface( angle, 0.0 ), glm::dvec2( angle, 0.0 ) } );
  }

  return mesh.vertices.size();
}

/**
 * `seamPair` is true for a loop that wraps EITHER parameter, so the grid has to
 * establish for itself that the seam it is about to read is the u one.
 * Everything downstream assumes it: columns come from the ring, rows from the
 * leg. On a surface closed in both u and v that assumption is a coin flip, and
 * reading a v seam as a u seam paves the chart with collapsed geometry.
 *
 * Same boundary, same everything, two surfaces: closed in u only, and closed in
 * both. The first must be accepted and the second rejected.
 */
void testGridRejectsAmbiguousSeamAxis() {

  printf( "grid gate requires an unambiguous u seam\n" );

  const double radius = 1.0;
  const double height = 4.0;

  tinynurbs::RationalSurface3d uOnly  = makeClosedTube( radius, height );
  tinynurbs::RationalSurface3d both   = makeDoublyClosedTube( radius );

  conway::geometry::RationalNurbsInverseMethod uOnlySolve( uOnly );
  conway::geometry::RationalNurbsInverseMethod bothSolve( both );

  check( uOnlySolve.closedU_ && !uOnlySolve.closedV_,
         "single-closed tube reports closed in u only" );
  check( bothSolve.closedU_ && bothSolve.closedV_,
         "doubly-closed tube reports closed in u AND v" );

  const double deflection = 1e-4;

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop( mesh, radius, height, 6, 5 );

    check( conway::geometry::tryFullCoverageSeamGrid(
             mesh, uOnlySolve, boundaryCount, deflection ),
           "accepted on a surface closed in u only" );
  }

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop( mesh, radius, height, 6, 5 );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, bothSolve, boundaryCount, deflection ),
           "REJECTED on a surface closed in u and v - the ambiguous case" );
  }
}

}  // namespace

int main() {

  testWideStepAcrossSeamKeepsArmijoSound();
  testClosureRequiresMatchingWeights();
  testDescentNeverWorsensFromItsSeed();
  testGridRejectsAmbiguousSeamAxis();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
