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

#include <cstdarg>

// TriangulateBspline's CDT error paths call this; the rest of the header is
// header-only. Defining it here keeps the test linking nothing, matching
// spherical_trim_test.cpp.
void Logger::logError( const char*, ... ) {}

#include <algorithm>
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

/**
 * A torus: the same rational circle swept around a second one, so it is
 * genuinely closed in BOTH u and v.
 *
 * A degenerate stand-in - a tube of zero height, whose v ends coincide - sets
 * both closure flags too, but it is rejected by the chart contract for its
 * degeneracy rather than for its ambiguity, so it cannot show that the
 * seam-axis clause is doing anything. This surface is well formed in every
 * other respect, which is what makes it a test of that clause alone.
 */
tinynurbs::RationalSurface3d makeTorus( double tubeRadius, double ringRadius ) {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 2;

  const double diagonal = std::sqrt( 2.0 ) / 2.0;

  // The nine-point rational circle polygon, as offsets and weights.
  const std::vector< glm::dvec2 > circle = {
    {  1.0,  0.0 }, {  1.0,  1.0 }, {  0.0,  1.0 }, { -1.0,  1.0 }, { -1.0, 0.0 },
    { -1.0, -1.0 }, {  0.0, -1.0 }, {  1.0, -1.0 }, {  1.0,  0.0 } };

  const std::vector< double > circleWeights = {
    1.0, diagonal, 1.0, diagonal, 1.0, diagonal, 1.0, diagonal, 1.0 };

  const size_t side = circle.size();

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( size_t row = 0; row < side; ++row ) {
    for ( size_t col = 0; col < side; ++col ) {

      const double radial = ringRadius + ( tubeRadius * circle[ row ].x );

      points.push_back( glm::dvec3( radial * circle[ col ].x,
                                    radial * circle[ col ].y,
                                    tubeRadius * circle[ row ].y ) );

      weights.push_back( circleWeights[ row ] * circleWeights[ col ] );
    }
  }

  surface.control_points = tinynurbs::array2( side, side, points );
  surface.weights        = tinynurbs::array2( side, side, weights );

  const std::vector< double > knots = {
    -PI, -PI, -PI, -PI / 2, -PI / 2, 0.0, 0.0, PI / 2, PI / 2, PI, PI, PI };

  surface.knots_u = knots;
  surface.knots_v = knots;

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
 *
 * Points come from EVALUATING the surface, not from trigonometry. The u
 * parameter of a rational quadratic circle is not the geometric angle - it is
 * a rational reparameterisation of it - so placing points by cos/sin puts them
 * on the right circle at the wrong parameter, and the chart assertion would
 * correctly reject them.
 *
 * `ringVWander` displaces the INTERIOR ring points along v. At zero the rings
 * are isoparametric, which is the shape the grid requires. Non-zero makes them
 * wander off constant v while staying exactly ON the surface and leaving every
 * run count, seam-point equality and column ordering untouched - which is
 * precisely the face the structural checks alone cannot tell apart.
 */
struct SeamLoopShape {

  double height      = 4.0;
  size_t seamCount   = 6;
  size_t ringCount   = 5;

  /** Displace interior ring points along v, off their constant-v line. */
  double ringVWander = 0.0;

  /** Transpose two interior columns, so columnU stops being monotone. */
  bool scrambleColumns = false;

  /** Transpose two interior seam rows, so the v partition stops being monotone. */
  bool scrambleRows = false;

  /**
   * Give the FAR ring a different u partition from the near one, at the same
   * point count. Two differently shaped trim rings tessellated adaptively land
   * on the same count all the time; it does not make their partitions
   * correspond.
   */
  double farRingUSkew = 0.0;

  /**
   * Nudge the far seam leg off the u domain bound, keeping its POINTS bitwise
   * equal to the near leg's. The columns then span slightly less than the full
   * period while every point-equality check still passes.
   */
  double seamBUOffset = 0.0;

  /**
   * Displace INTERIOR seam points off the u domain bound, leaving the leg's
   * two endpoints - which become the first and last columns - where they were.
   */
  double seamUWander = 0.0;
};

size_t buildSeamLoop(
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex >& mesh,
    const conway::geometry::RationalNurbsInverseMethod& solve,
    const SeamLoopShape& shape ) {

  const auto place =
    [ & ]( double u, double v ) {
      mesh.makeVertex( { solve.evaluator.point( u, v ), glm::dvec2( u, v ) } );
    };

  // The v samples the seam is walked at, and the u samples the rings are.
  // Held as sequences so a transposition applies CONSISTENTLY to both legs and
  // both rings - otherwise the legs stop being exact reverses, or the far
  // ring stops aligning with the columns, and the loop would be rejected by an
  // earlier check instead of the one under test.
  std::vector< double > seamV;
  std::vector< double > ringU;

  for ( size_t at = 0; at < shape.seamCount; ++at ) {
    seamV.push_back( ( shape.height * at ) / double( shape.seamCount - 1 ) );
  }

  for ( size_t at = 1; at <= shape.ringCount; ++at ) {
    ringU.push_back( PI - ( ( 2.0 * PI * at ) / double( shape.ringCount + 1 ) ) );
  }

  if ( shape.scrambleRows && seamV.size() >= 4 ) {
    std::swap( seamV[ 1 ], seamV[ 2 ] );
  }

  if ( shape.scrambleColumns && ringU.size() >= 3 ) {
    std::swap( ringU[ 0 ], ringU[ 1 ] );
  }

  // seamA: u = +pi, except where the interior is deliberately pulled off it.
  for ( size_t at = 0; at < seamV.size(); ++at ) {

    const bool endpoint = at == 0 || at + 1 == seamV.size();

    place( endpoint ? PI : PI - shape.seamUWander, seamV[ at ] );
  }

  // ringA: the far ring, interior points only.
  for ( size_t at = 0; at < ringU.size(); ++at ) {

    place( ringU[ at ],
           shape.height -
             ( shape.ringVWander * ( double( at + 1 ) / double( ringU.size() ) ) ) );
  }

  // seamB: u = -pi, the SAME points as seamA in reverse. Placed from seamA's
  // own evaluation so the two legs are bitwise equal, which is what the run
  // analysis keys on.
  for ( size_t at = 0; at < shape.seamCount; ++at ) {

    const size_t mirror = shape.seamCount - 1 - at;

    mesh.makeVertex(
      { mesh.vertices[ mirror ].point,
        glm::dvec2( -PI + shape.seamBUOffset, seamV[ mirror ] ) } );
  }

  // ringB: the near ring. Walked in the opposite sense, so reversing it
  // reproduces ringA's column order - which is what lets the grid take its
  // columns from one ring and use them for the other.
  for ( size_t at = 0; at < ringU.size(); ++at ) {

    // Reversed, not negated: reverse( ringB ) has to reproduce ringA's column
    // order, because the grid takes its columns from the near ring and applies
    // them to the far one.
    const size_t mirror = ringU.size() - 1 - at;

    // Skew pulls the far ring's samples toward one side while keeping their
    // count, their order and their v identical.
    const double fraction = double( mirror + 1 ) / double( ringU.size() + 1 );
    const double skewed =
      ringU[ mirror ] + ( shape.farRingUSkew * std::sin( PI * fraction ) );

    place( skewed,
           shape.ringVWander * ( double( mirror + 1 ) / double( ringU.size() ) ) );
  }

  return mesh.vertices.size();
}

/**
 * The grid builds the tensor product columnU x rowV, which is a RECTANGULAR
 * ISOPARAMETRIC chart. The run analysis cannot establish the face is one: it
 * counts runs, checks they alternate and are pairwise equal in length, and
 * compares seam points. A loop whose other two runs wander in v satisfies all
 * of that and is not a rectangle, so paving one over it replaces the real trim
 * with constant-v boundaries that can leave the face entirely.
 *
 * Same surface, same run structure, same column ordering; the only difference
 * is whether the rings follow constant v.
 */
void testGridRequiresIsoparametricChart() {

  printf( "grid gate requires a rectangular isoparametric chart\n" );

  const double height = 4.0;

  tinynurbs::RationalSurface3d surface = makeClosedTube( 1.0, height );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  // sqrt( deflection ) = 0.01; |dS/dv| is height, so a v wander of 0.2 puts
  // the ring roughly 0.8 off its isocurve - decisively outside, not marginal.
  const double deflection = 1e-4;

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop( mesh, solve, SeamLoopShape{} );

    check( conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "accepted when the rings are isoparametric" );
  }

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .ringVWander = 0.2 } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when the rings wander off constant v - same run structure" );
  }

  // The grid indexes columns by a u partition and rows by a v partition, and
  // both have to be ORDERED for that to be a bijection onto the chart. A
  // transposed pair leaves every point exactly on the surface and exactly on
  // its own isoparametric line, so the checks above all pass; what it breaks
  // is the partition, and cells then fold back over one another.
  // Two end rings tessellated to the same NUMBER of points, at different u
  // positions. `ringA.size() == ringB.size()` passes; the partitions do not
  // correspond, so every interior column would come from the near ring while
  // the far ring's vertices were inserted unmatched, and the refinement would
  // test the far side at the wrong u.
  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .farRingUSkew = 0.35 } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when the far ring's u partition differs at equal count" );
  }

  // Columns that stop short of the u period. The offset is small enough that
  // the leg still lies on its constant-u line within the deflection tolerance,
  // so check 7 passes and only the domain-bound check can catch it - the chart
  // would be a strip, not the whole surface, while the function claims full
  // coverage.
  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .seamBUOffset = 0.005 } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when the columns stop short of the full u period" );
  }

  // Seam legs that leave their constant-u line while the columns, the rows and
  // both rings stay correct. Only the endpoints of the leg become columns, so
  // wandering its interior is invisible to every other clause.
  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .seamUWander = 0.2 } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when the seam legs wander off constant u" );
  }

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .scrambleColumns = true } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when two columns are transposed - u partition not monotone" );
  }

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount = buildSeamLoop(
        mesh, solve, SeamLoopShape{ .scrambleRows = true } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, solve, boundaryCount, deflection ),
           "REJECTED when two rows are transposed - v partition not monotone" );
  }
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
  tinynurbs::RationalSurface3d both   = makeTorus( radius, 4.0 );

  conway::geometry::RationalNurbsInverseMethod uOnlySolve( uOnly );
  conway::geometry::RationalNurbsInverseMethod bothSolve( both );

  check( uOnlySolve.closedU_ && !uOnlySolve.closedV_,
         "single-closed tube reports closed in u only" );
  check( bothSolve.closedU_ && bothSolve.closedV_,
         "doubly-closed tube reports closed in u AND v" );

  const double deflection = 1e-4;

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount =
      buildSeamLoop( mesh, uOnlySolve, SeamLoopShape{} );

    check( conway::geometry::tryFullCoverageSeamGrid(
             mesh, uOnlySolve, boundaryCount, deflection ),
           "accepted on a surface closed in u only" );
  }

  {
    conway::geometry::WingedEdgeMesh< conway::geometry::ParameterVertex > mesh;

    const size_t boundaryCount =
      buildSeamLoop( mesh, bothSolve, SeamLoopShape{ .height = 2.0 } );

    check( !conway::geometry::tryFullCoverageSeamGrid(
             mesh, bothSolve, boundaryCount, deflection ),
           "REJECTED on a surface closed in u and v - the ambiguous case" );
  }
}


/**
 * A multi-turn helical ribbon: v runs along the coil, u across its narrow
 * width. Consecutive turns pass close to each other in 3D, which is what gives
 * the inverse solve several near-equal basins for one query - the conway#594
 * configuration, and the one every cold-start failure in the corpus sits on.
 */
tinynurbs::RationalSurface3d makeHelicalRibbon(
    double coilRadius, double pitch, double turns, double width,
    int samples ) {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 3;

  const int alongCount  = samples;
  const int acrossCount = 2;

  surface.knots_u = { 0, 0, 1, 1 };

  surface.knots_v.clear();

  for ( int i = 0; i < surface.degree_v + 1; ++i ) {
    surface.knots_v.push_back( 0.0 );
  }

  for ( int i = 1; i < alongCount - surface.degree_v; ++i ) {
    surface.knots_v.push_back(
      static_cast< double >( i ) / ( alongCount - surface.degree_v ) );
  }

  for ( int i = 0; i < surface.degree_v + 1; ++i ) {
    surface.knots_v.push_back( 1.0 );
  }

  surface.control_points.resize( acrossCount, alongCount );
  surface.weights.resize( acrossCount, alongCount );

  for ( int j = 0; j < alongCount; ++j ) {

    const double t = turns * 2.0 * PI * j / ( alongCount - 1 );

    const glm::dvec3 centre(
      coilRadius * std::cos( t ),
      coilRadius * std::sin( t ),
      pitch * t / ( 2.0 * PI ) );

    const glm::dvec3 outward( std::cos( t ), std::sin( t ), 0.0 );

    for ( int i = 0; i < acrossCount; ++i ) {
      surface.control_points( i, j ) = centre + outward * ( width * ( i - 0.5 ) );
      surface.weights( i, j )        = 1.0;
    }
  }

  return surface;
}


/**
 * The first point of a closed trim must not cold-start into the wrong basin.
 *
 * resetContinuity() scopes the continuity seed to one bound, which is correct
 * and deliberate - see its own comment on why a cross-loop seed can converge
 * silently onto a different sheet. But it leaves the FIRST point of every
 * bound with no seed at all, cold-starting from the grid. On a surface that
 * passes near itself the grid can pick the wrong basin, and the bad uv then
 * seeds its successors until the descent recovers.
 *
 * Measured across Orbiter, Right_Hand and nist_ctc_02: 10 of 848 b-spline
 * faces cold-start into a residual more than 100x their own loop median and
 * above 0.1% of the face extent, the worst 26% of the face away from the
 * surface. Downstream it shows as a spurious v-monotone reversal where the
 * good run meets the bad head - the defect behind conway#647.
 *
 * The loop is closed, so its first point neighbours its last; re-solving the
 * head seeded from that last solution removes the cold start. This pins the
 * property directly on the solver rather than through triangulation, because
 * the triangulated consequence depends on the ribbon gate as well.
 *
 * Red-proven: with the second pass removed the first point's residual here is
 * 1.2e+00 against a loop median of 1.3e-03.
 */
void testClosedTrimHeadIsNotColdStarted() {

  printf( "=== closed trim head is not cold-started ===\n" );

  constexpr double COIL_RADIUS = 3.0;
  constexpr double PITCH       = 0.60;
  constexpr double TURNS       = 8.0;
  constexpr double WIDTH       = 0.30;
  constexpr int    SAMPLES     = 200;
  constexpr int    ALONG       = 400;

  // Start the loop 500 points in. A real STEP edge loop begins wherever the
  // file put it; starting at a domain corner is the case the seed grid always
  // gets right, so it would prove nothing.
  constexpr int    ROTATION    = 500;

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( COIL_RADIUS, PITCH, TURNS, WIDTH, SAMPLES );

  conway::geometry::RationalSurfaceEvaluator evaluator( surface );

  std::vector< glm::dvec3 > loop;

  for ( int i = 0; i < ALONG; ++i ) {
    loop.push_back( evaluator.point( 0.02, static_cast< double >( i ) / ALONG ) );
  }

  for ( int i = ALONG; i > 0; --i ) {
    loop.push_back( evaluator.point( 0.98, static_cast< double >( i ) / ALONG ) );
  }

  std::vector< glm::dvec3 > rotated;

  for ( size_t i = 0; i < loop.size(); ++i ) {
    rotated.push_back( loop[ ( i + ROTATION ) % loop.size() ] );
  }

  // Close it by DUPLICATION, as all 859 corpus trim bounds do. The gate accepts
  // explicit closure only - see reSolveClosedTrimHead on why proximity cannot
  // be made to work - so a fixture closing by adjacency would now be declined
  // and would test nothing.
  rotated.push_back( rotated.front() );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  solve.resetContinuity();

  std::vector< glm::dvec2 > solved;

  for ( const glm::dvec3& query : rotated ) {
    solved.push_back( solve( query ) );
  }

  // THE PRODUCTION HELPER, not a copy of it. An earlier version of this test
  // reimplemented the loop inline, which meant deleting the production call
  // site left the test green - caught in review on conway-geom#194.
  conway::geometry::reSolveClosedTrimHead( solve, rotated, solved, false );

  auto residual = [ & ]( size_t i ) {

    return glm::distance(
      solve.evaluator.point( solved[ i ].x, solved[ i ].y ), rotated[ i ] );
  };

  std::vector< double > rest;

  for ( size_t i = 1; i < rotated.size(); ++i ) {
    rest.push_back( residual( i ) );
  }

  std::sort( rest.begin(), rest.end() );

  const double head   = residual( 0 );
  const double median = rest[ rest.size() / 2 ];

  printf( "      headResidual=%.4e loopMedian=%.4e ratio=%.2e\n",
          head, median, head / median );

  // The head must be in the same league as the rest of its own loop. A
  // cold-start failure here is four to five orders out, so the bar does not
  // need to be tight to be decisive.
  check( head < median * 100.0,
         "the loop head solves to within 100x its own loop median" );
}


/**
 * The closure gate declines OPEN bounds at every sampling density, and declines
 * adjacency closure too.
 *
 * Everything reSolveClosedTrimHead does rests on the tail being the head. An
 * open bound's tail is somewhere unrelated, and seeding the head from it is the
 * cross-sheet failure resetContinuity exists to prevent, arriving by another
 * door.
 *
 * WHY THE GATE IS EXPLICIT CLOSURE AND NOT PROXIMITY. The first version
 * compared the closing gap to the loop's own median segment. That is unusable,
 * not merely mis-tuned: the gap is fixed by geometry while the median scales
 * with sampling. The same open one-turn arc measures
 *
 *     samples   120    60    40    30    20
 *     ratio    4.05  2.02  1.35  1.01  0.68
 *
 * against an adjacency-closed loop's ~1.0 - indistinguishable at 30 samples,
 * and at 20 it looks MORE closed than a closed loop. Every density here is
 * asserted, because a single-density fixture is exactly what let the earlier
 * 2.0 threshold look sound. Found on bldrs-ai/conway#655.
 *
 * Red-proven, and precisely: under the previous proximity gate at 2.0 the
 * 40-interval case is accepted (rewritten=1) and this test fails. The 30- and
 * 20-interval cases pass under BOTH gates - their re-solve happens to return
 * the head unchanged - so they are coverage of the gate's behaviour across the
 * density range, not proof of it. Only the 40 case discriminates. Saying so
 * because a row of five green checks reads like five proofs, and three of them
 * are not.
 */
void testOpenBoundIsNotReSeeded() {

  printf( "=== open bounds are declined at every sampling density ===\n" );

  constexpr double COIL_RADIUS = 3.0;
  constexpr double PITCH       = 0.60;
  constexpr double TURNS       = 8.0;
  constexpr double WIDTH       = 0.30;
  constexpr int    SAMPLES     = 200;

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( COIL_RADIUS, PITCH, TURNS, WIDTH, SAMPLES );

  conway::geometry::RationalSurfaceEvaluator evaluator( surface );

  // An open arc spanning exactly one coil turn: its tail sits one pitch from
  // its head - spatially adjacent, a whole turn away in v, a different sheet.
  for ( int intervals : { 120, 60, 40, 30, 20 } ) {

    std::vector< glm::dvec3 > open;

    for ( int i = 0; i <= intervals; ++i ) {
      open.push_back( evaluator.point(
        0.5, 0.30 + ( 1.0 / TURNS ) * i / static_cast< double >( intervals ) ) );
    }

    conway::geometry::RationalNurbsInverseMethod solve( surface );

    solve.resetContinuity();

    std::vector< glm::dvec2 > solved;

    for ( const glm::dvec3& query : open ) {
      solved.push_back( solve( query ) );
    }

    const std::vector< glm::dvec2 > before = solved;

    const size_t rewritten =
      conway::geometry::reSolveClosedTrimHead( solve, open, solved, false );

    std::vector< double > segments;

    for ( size_t i = 0; i + 1 < open.size(); ++i ) {
      segments.push_back( glm::distance( open[ i + 1 ], open[ i ] ) );
    }

    std::sort( segments.begin(), segments.end() );

    const double gap    = glm::distance( open.back(), open.front() );
    const double median = segments[ segments.size() / 2 ];

    printf( "      intervals=%-4d gap/median=%.2f rewritten=%zu\n",
            intervals, gap / median, rewritten );

    check( rewritten == 0 && solved == before,
           std::string( "open arc declined at " ) +
             std::to_string( intervals ) + " intervals" );
  }

  // ADJACENCY closure: a loop that genuinely returns to its start, with the
  // tail ONE sample increment from the head and no duplicate.
  //
  // Two earlier attempts at this fixture asserted nothing. The first pushed a
  // closing point and immediately popped it; the second used a long open arc
  // whose ends are 23.9 sample steps apart. Both were declined for reasons
  // unrelated to adjacency, so both passed with the gate removed. The guard
  // below - that the gap really is about one step - is what makes this one
  // real, and it is asserted rather than assumed.
  //
  // The geometry here is identical to testClosedPointListLoopIsAccepted's. The
  // ONLY difference is closedByConstruction, so the pair is a clean A/B on the
  // producer's closure signal: unflagged it is declined, flagged it is
  // accepted. That is the whole design - closure is a fact the producer states,
  // never something measured from the points.
  {
    constexpr int ALONG = 400;

    std::vector< glm::dvec3 > adjacent;

    for ( int i = 0; i < ALONG; ++i ) {
      adjacent.push_back(
        evaluator.point( 0.02, static_cast< double >( i ) / ALONG ) );
    }

    for ( int i = ALONG; i > 0; --i ) {
      adjacent.push_back(
        evaluator.point( 0.98, static_cast< double >( i ) / ALONG ) );
    }

    const double step = glm::distance( adjacent[ 1 ], adjacent[ 0 ] );
    const double gap  = glm::distance( adjacent.back(), adjacent.front() );

    conway::geometry::RationalNurbsInverseMethod solve( surface );

    solve.resetContinuity();

    std::vector< glm::dvec2 > solved;

    for ( const glm::dvec3& query : adjacent ) {
      solved.push_back( solve( query ) );
    }

    const std::vector< glm::dvec2 > before = solved;

    const size_t rewritten =
      conway::geometry::reSolveClosedTrimHead( solve, adjacent, solved, false );

    printf( "      adjacency: gap=%.4e step=%.4e gap/step=%.2f rewritten=%zu\n",
            gap, step, gap / step, rewritten );

    check( gap < step * 2.0,
           "the adjacency fixture really does close to within a sample step" );

    check( rewritten == 0 && solved == before,
           "adjacency closure alone is declined - the producer must say so" );
  }
}

/**
 * A closed POINT-LIST loop is accepted, on the producer's own say-so.
 *
 * GetLoop's point-list branch - IFCPOLYLOOP and the AP214 poly-loop path -
 * emits a closed polygon that does NOT repeat its head. An exact-duplication
 * gate declines that entire producer class and silently restores the
 * cold-start error for it. Found by review on conway-geom#195, and it is the
 * negative answer to "is there a front-end closure signal": the producer knows,
 * it just was not encoding it.
 *
 * WHAT THIS COVERS AND WHAT IT DOES NOT, stated because the gap is real:
 * `pointListLoopIsClosedPolygon` is the rule GetLoop applies, and the gate
 * honouring its result is asserted end to end below. What is NOT covered is
 * GetLoop's one-line CALL to it - ConwayGeometryProcessor.h reaches for draco
 * and the GLTF SDK, so these standalone tests cannot include the producer at
 * all (which is why run_native_tests.sh excludes the linked tests). Deleting
 * that call would leave this green. Closing that needs the linked test target
 * repaired, which is out of scope here.
 *
 * Red-proven: with the gate's closedByConstruction branch removed, the
 * un-duplicated loop is declined and rewritten drops to 0.
 */
void testClosedPointListLoopIsAccepted() {

  printf( "=== closed point-list loop is accepted ===\n" );

  constexpr double COIL_RADIUS = 3.0;
  constexpr double PITCH       = 0.60;
  constexpr double TURNS       = 8.0;
  constexpr double WIDTH       = 0.30;
  constexpr int    SAMPLES     = 200;
  constexpr int    ALONG       = 400;
  constexpr int    ROTATION    = 500;

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( COIL_RADIUS, PITCH, TURNS, WIDTH, SAMPLES );

  conway::geometry::RationalSurfaceEvaluator evaluator( surface );

  std::vector< glm::dvec3 > loop;

  for ( int i = 0; i < ALONG; ++i ) {
    loop.push_back( evaluator.point( 0.02, static_cast< double >( i ) / ALONG ) );
  }

  for ( int i = ALONG; i > 0; --i ) {
    loop.push_back( evaluator.point( 0.98, static_cast< double >( i ) / ALONG ) );
  }

  // The point-list form: rotated so the head is not at a domain corner, and
  // with NO head repeated - which is exactly what the producer emits.
  std::vector< glm::dvec3 > pointList;

  for ( size_t i = 0; i < loop.size(); ++i ) {
    pointList.push_back( loop[ ( i + ROTATION ) % loop.size() ] );
  }

  const bool closed =
    conway::geometry::pointListLoopIsClosedPolygon( pointList );

  printf( "      points=%zu headRepeated=%d closedByConstruction=%d\n",
          pointList.size(),
          ( pointList.front() == pointList.back() ) ? 1 : 0,
          closed ? 1 : 0 );

  check( pointList.front() != pointList.back(),
         "the producer form does NOT repeat the head - the case under test" );

  check( closed, "and the producer's rule calls it closed" );

  // A VERTEX_LOOP must not be flagged: one point encloses nothing.
  check( !conway::geometry::pointListLoopIsClosedPolygon(
           { glm::dvec3( 0.0 ) } ),
         "a one-point vertex loop is not a closed polygon" );

  conway::geometry::RationalNurbsInverseMethod solve( surface );

  solve.resetContinuity();

  std::vector< glm::dvec2 > solved;

  for ( const glm::dvec3& query : pointList ) {
    solved.push_back( solve( query ) );
  }

  const size_t rewritten =
    conway::geometry::reSolveClosedTrimHead( solve, pointList, solved, closed );

  auto residual = [ & ]( size_t i ) {

    return glm::distance(
      solve.evaluator.point( solved[ i ].x, solved[ i ].y ), pointList[ i ] );
  };

  std::vector< double > rest;

  for ( size_t i = 1; i < pointList.size(); ++i ) {
    rest.push_back( residual( i ) );
  }

  std::sort( rest.begin(), rest.end() );

  printf( "      rewritten=%zu headResidual=%.4e loopMedian=%.4e\n",
          rewritten, residual( 0 ), rest[ rest.size() / 2 ] );

  check( rewritten > 0,
         "an un-duplicated closed loop is accepted, not declined" );

  check( residual( 0 ) < rest[ rest.size() / 2 ] * 100.0,
         "and its head is solved to within 100x its own loop median" );
}

/**
 * The raw-array producer's closure composition.
 *
 * createSimpleBound3D (conway-api.cpp) is the optimised path extractFace uses
 * for point-list loops. Its closure is `pointListLoopClosure( callerSaysClosed,
 * points )`, and BOTH halves matter for different reasons:
 *
 *   - the caller's word is REQUIRED because a raw vertex array carries no
 *     entity type. Deriving closure inside would make the function assert
 *     something about every future caller it cannot know.
 *
 *   - the predicate is a DEGENERACY FILTER, not a check on the caller. A
 *     one-point VERTEX_LOOP is a perfectly legitimate thing to hand over and
 *     to call closed - it is the degenerate loop at a pole - but it encloses
 *     nothing, so it must not reach a consumer that treats it as a polygon.
 *
 * Missed in the first audit of this family, which fixed GetLoop and stopped -
 * review on conway-geom#195 found this second producer. A later pass of that
 * review then found this test VACUOUS: it re-implemented the AND in a local
 * lambda, so deleting the composition from production left it green. The
 * composition now lives in pointListLoopClosure (ConwayCurve.h) and the
 * producer calls it, so this exercises the shipped code rather than a copy.
 *
 * Red-proven by weakening the real helper, not by reading: dropping its
 * `callerDeclaresClosed &&` fails check 2, and dropping the predicate call
 * fails check 3.
 *
 * COVERAGE LIMIT, same as the GetLoop case: conway-api.cpp needs embind and the
 * standalone tests link nothing, so the entry point itself is not exercised -
 * only the helper it delegates to and the gate that consumes the result. A
 * caller passing the wrong argument, or none, is not catchable here at all;
 * that is what the required parameter on the TypeScript wrapper is for.
 */
void testRawArrayProducerClosureComposition() {

  printf( "=== raw-array producer closure composition ===\n" );

  using conway::geometry::pointListLoopClosure;

  const std::vector< glm::dvec3 > polygon = {
    glm::dvec3( 0.0, 0.0, 0.0 ),
    glm::dvec3( 1.0, 0.0, 0.0 ),
    glm::dvec3( 1.0, 1.0, 0.0 ),
    glm::dvec3( 0.0, 1.0, 0.0 ) };

  const std::vector< glm::dvec3 > vertexLoop = { glm::dvec3( 0.0, 0.0, 0.0 ) };

  check( pointListLoopClosure( true, polygon ),
         "a poly loop the caller calls closed is flagged closed" );

  check( !pointListLoopClosure( false, polygon ),
         "the caller's word is required - geometry alone does not flag it" );

  check( !pointListLoopClosure( true, vertexLoop ),
         "a one-point vertex loop is filtered out as degenerate" );
}

/**
 * End-to-end through TriangulateBspline: the head fix reaches the mesh.
 *
 * The other tests here exercise the solver. This one drives the production
 * triangulation entry point and asserts on the geometry it emits, so that
 * deleting the call to reSolveClosedTrimHead inside TriangulateBspline turns
 * this test red. Review on conway-geom#194 found that the solver-level tests
 * alone did not have that property.
 *
 * The assertion is emitted AREA against the analytic ribbon area. A wrong-basin
 * head puts the loop's first vertices in the wrong place in the chart, and the
 * triangulation smears across the strip - the same failure that made
 * conway#647's flank ship 328 mm2 where its symmetric sibling ships 49.
 *
 * Red-proven by disabling the production pass, not by reading.
 */
void testHeadFixReachesTheEmittedMesh() {

  printf( "=== head fix reaches the emitted mesh ===\n" );

  constexpr double COIL_RADIUS = 3.0;
  constexpr double PITCH       = 0.60;
  constexpr double TURNS       = 8.0;
  constexpr double WIDTH       = 0.30;
  constexpr int    SAMPLES     = 200;
  constexpr int    ALONG       = 400;
  constexpr int    ROTATION    = 500;

  const tinynurbs::RationalSurface3d nurbs =
    makeHelicalRibbon( COIL_RADIUS, PITCH, TURNS, WIDTH, SAMPLES );

  conway::geometry::RationalSurfaceEvaluator evaluator( nurbs );

  // The same surface, expressed the way the extractor hands it over.
  conway::geometry::IfcSurface surface;

  surface.transformation           = glm::dmat4( 1.0 );
  surface.BSplineSurface.Active    = true;
  surface.BSplineSurface.UDegree   = nurbs.degree_u;
  surface.BSplineSurface.VDegree   = nurbs.degree_v;

  for ( int i = 0; i < 2; ++i ) {

    std::vector< glm::dvec3 > row;
    std::vector< double >     weightRow;

    for ( int j = 0; j < SAMPLES; ++j ) {
      row.push_back( nurbs.control_points( i, j ) );
      weightRow.push_back( 1.0 );
    }

    surface.BSplineSurface.ControlPoints.push_back( row );
    surface.BSplineSurface.WeightPoints.push_back( weightRow );
  }

  // Knots arrive as distinct values plus multiplicities.
  auto pushKnots = []( const std::vector< double >& flat,
                       std::vector< glm::f64 >& values,
                       std::vector< glm::f64 >& multiplicities ) {

    for ( size_t i = 0; i < flat.size(); ) {

      size_t run = 1;

      while ( i + run < flat.size() && flat[ i + run ] == flat[ i ] ) {
        ++run;
      }

      values.push_back( flat[ i ] );
      multiplicities.push_back( static_cast< double >( run ) );

      i += run;
    }
  };

  pushKnots( nurbs.knots_u, surface.BSplineSurface.UKnots,
             surface.BSplineSurface.UMultiplicity );
  pushKnots( nurbs.knots_v, surface.BSplineSurface.VKnots,
             surface.BSplineSurface.VMultiplicity );

  // A closed ribbon trim, rotated so the head does not start at a domain
  // corner - the case the seed grid always gets right.
  std::vector< glm::dvec3 > loop;

  for ( int i = 0; i < ALONG; ++i ) {
    loop.push_back( evaluator.point( 0.02, static_cast< double >( i ) / ALONG ) );
  }

  for ( int i = ALONG; i > 0; --i ) {
    loop.push_back( evaluator.point( 0.98, static_cast< double >( i ) / ALONG ) );
  }

  conway::geometry::IfcBound3D bound;

  bound.type        = conway::geometry::IfcBoundType::OUTERBOUND;
  bound.orientation = true;

  for ( size_t i = 0; i < loop.size(); ++i ) {
    bound.curve.points.push_back( loop[ ( i + ROTATION ) % loop.size() ] );
  }

  // Close it explicitly, as every trim polyline in the corpus does: all 859
  // measured bounds repeat their first point as their last.
  bound.curve.points.push_back( bound.curve.points.front() );

  conway::geometry::Geometry geometry;

  conway::geometry::TriangulateBspline(
    geometry, { bound }, surface, 1.0, COIL_RADIUS * 2.0 );

  double area = 0.0;

  for ( const auto& triangle : geometry.triangles ) {

    const glm::dvec3 a = geometry.vertices[ triangle.vertices[ 0 ] ];
    const glm::dvec3 b = geometry.vertices[ triangle.vertices[ 1 ] ];
    const glm::dvec3 c = geometry.vertices[ triangle.vertices[ 2 ] ];

    area += 0.5 * glm::length( glm::cross( b - a, c - a ) );
  }

  // The ribbon is 0.96 of its width across and runs the full coil length.
  // Helix length for 8 turns at radius 3 and pitch 0.6.
  const double coilLength =
    TURNS * std::sqrt( std::pow( 2.0 * PI * COIL_RADIUS, 2.0 ) +
                       std::pow( PITCH, 2.0 ) );
  const double expected = coilLength * WIDTH * 0.96;

  printf( "      triangles=%zu area=%.4f expected=%.4f ratio=%.4f\n",
          geometry.triangles.size(), area, expected, area / expected );

  check( geometry.triangles.size() > 0, "the ribbon emits geometry" );

  // A smeared head inflates the area far past the strip; an inscribed
  // tessellation of the true strip lands just under it.
  check( area < expected * 1.10,
         "the emitted area does not overshoot the ribbon - head not smeared" );

  check( area > expected * 0.80,
         "the emitted area covers the ribbon" );
}

}  // namespace

int main() {

  testClosedTrimHeadIsNotColdStarted();
  testOpenBoundIsNotReSeeded();
  testClosedPointListLoopIsAccepted();
  testRawArrayProducerClosureComposition();
  testHeadFixReachesTheEmittedMesh();

  testWideStepAcrossSeamKeepsArmijoSound();
  testClosureRequiresMatchingWeights();
  testDescentNeverWorsensFromItsSeed();
  testGridRejectsAmbiguousSeamAxis();
  testGridRequiresIsoparametricChart();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
