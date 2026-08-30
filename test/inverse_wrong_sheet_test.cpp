/*
 * Soundness tests for the two wrong-sheet recovery tiers in
 * RationalNurbsInverseMethod: the symmetric retry from a REJECTED continuity
 * candidate, and the residual-triggered local grid refinement armed by
 * boundSeverityToChord.
 *
 * Both exist for one situation the 8x8 seed grid cannot handle: a support
 * surface that passes near ITSELF, where the nearest grid sample belongs to a
 * neighbouring sheet and the descent - being local - converges inside it. The
 * result is a converged-looking uv that is a whole turn wrong, which is what
 * turns a trim polyline into a self-overlapping uv chart (bldrs-ai/conway#642).
 *
 * The regression digests cannot pin either tier. They pin the shipped mesh of
 * whole models, so a solve that is wrong on 14 boundary points of 74,838 moves
 * a digest by an amount indistinguishable from noise, and the seed comparison
 * that makes the failure reachable at all is internal to one call. Each test
 * below fails if its tier is reverted; that was verified by reverting, not by
 * reading.
 *
 * Standalone by design: it includes mesh_utils.h directly and links nothing,
 * matching nurbs_seam_test.cpp and spherical_trim_test.cpp.
 */
#include "conway_geometry/operations/mesh_utils.h"

#include <cstdarg>

// TriangulateBspline's CDT error paths call this; the rest of the header is
// header-only. Defining it here keeps the test linking nothing.
void Logger::logError( const char*, ... ) {}

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

using namespace conway::geometry;

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
 * A helical ribbon: `turns` turns of a rational quadratic circle chained down
 * v, lifted by `pitch` per turn, with a linear radial extent of `width`
 * across u.
 *
 * This is the shape the defect needs, reduced to its essentials. Two things
 * matter and both are properties of a thread flank rather than of this
 * construction: consecutive turns lie `pitch` apart in 3D while being four
 * units apart in v, so the surface passes near itself; and the v domain is
 * long enough (24 units over an 8-sample grid) that one grid cell spans most
 * of a turn.
 *
 * Rational, not polynomial: the quarter-arc control net alternates on-curve
 * points (weight 1) with corner points at radius r/cos(45 deg) (weight
 * cos(45 deg)), so the weights are genuinely non-uniform and the rational
 * evaluation path is the one under test.
 */
tinynurbs::RationalSurface3d makeHelicalRibbon(
    double radius, double width, double pitch, int turns ) {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 2;

  const int spans   = 4 * turns;
  const int columns = 2 * spans + 1;

  const double cornerWeight = std::cos( PI / 4.0 );

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  points.reserve( 2 * size_t( columns ) );
  weights.reserve( 2 * size_t( columns ) );

  for ( int row = 0; row < 2; ++row ) {

    const double r = radius + ( row == 0 ? 0.0 : width );

    for ( int j = 0; j < columns; ++j ) {

      const double theta  = 0.25 * PI * double( j );
      const double z      = pitch * double( j ) / 8.0;
      const bool   corner = ( j % 2 ) == 1;
      const double cr     = corner ? r / cornerWeight : r;

      points.push_back(
        glm::dvec3( cr * std::cos( theta ), cr * std::sin( theta ), z ) );

      weights.push_back( corner ? cornerWeight : 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 2, size_t( columns ), points );
  surface.weights        = tinynurbs::array2( 2, size_t( columns ), weights );

  surface.knots_u = { 0.0, 0.0, 1.0, 1.0 };

  surface.knots_v.push_back( 0.0 );
  surface.knots_v.push_back( 0.0 );
  surface.knots_v.push_back( 0.0 );

  for ( int s = 1; s < spans; ++s ) {
    surface.knots_v.push_back( double( s ) );
    surface.knots_v.push_back( double( s ) );
  }

  surface.knots_v.push_back( double( spans ) );
  surface.knots_v.push_back( double( spans ) );
  surface.knots_v.push_back( double( spans ) );

  return surface;
}

double residual( const RationalNurbsInverseMethod& solve,
                 const glm::dvec2&                 uv,
                 const glm::dvec3&                 point ) {

  return glm::distance( solve.evaluator.point( uv.x, uv.y ), point );
}

/** The nearest of the 8x8 seed grid's samples, and how far it is. */
double bestGridDistance( const RationalNurbsInverseMethod& solve,
                         const glm::dvec3&                 point ) {

  double best = DBL_MAX;

  for ( size_t i = 0; i < INVERSE_GRID_SIDE; ++i ) {
    for ( size_t j = 0; j < INVERSE_GRID_SIDE; ++j ) {

      const glm::dvec2 uv =
        solve.min_extent + ( solve.max_extent - solve.min_extent ) *
        glm::dvec2( double( i ) * INVERSE_GRID_FACTOR,
                    double( j ) * INVERSE_GRID_FACTOR );

      best = std::min( best, glm::distance(
        solve.evaluator.point( uv.x, uv.y ), point ) );
    }
  }

  return best;
}

// The query, and the continuity candidate one chord back along the ribbon.
// Chosen so the seed comparison is a NEAR TIE that the grid wins - 1.154
// against 1.191 - which is the situation on the four #642 faces, where the
// local trim chord (0.97-2.75) slightly exceeds the distance to the
// neighbouring turn (0.86-1.39). A fixture where the grid wins comfortably
// would test the same code and a different geometry.
constexpr double QUERY_U = 0.5;
constexpr double QUERY_V = 2.75;
constexpr double CHORD_V = 0.25;

/**
 * The defect: the correct seed is computed, loses the nearest-in-3D
 * comparison to a wrong-sheet grid sample, and is discarded.
 *
 * Reverting the symmetric retry leaves the grid's answer, which is one full
 * turn of the ribbon away.
 */
void testRejectedContinuitySeedIsStillTried() {

  printf( "== rejected continuity candidate is still descended from\n" );

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( 3.0, 0.35, 1.0, 6 );

  check( tinynurbs::surfaceIsValid( surface ), "fixture surface is valid" );

  RationalNurbsInverseMethod solve( surface );

  const glm::dvec3 query = solve.evaluator.point( QUERY_U, QUERY_V );

  const glm::dvec2 continuityUv( QUERY_U, QUERY_V + CHORD_V );

  const double gridDistance = bestGridDistance( solve, query );

  const double continuityDistance =
    glm::distance( solve.evaluator.point( continuityUv.x, continuityUv.y ),
                   query );

  // The precondition the whole test rests on. Asserted rather than assumed:
  // if a change to the grid or the fixture ever made the continuity candidate
  // win, the seed comparison would take it, the symmetric retry would never
  // run, and the assertions below would pass while testing nothing.
  check( gridDistance < continuityDistance,
         "precondition: the grid sample is nearer than the continuity "
         "candidate, so the seed comparison discards the latter" );

  // What the solve does with no continuity candidate at all: the grid seed's
  // own basin, a whole turn away from the answer.
  solve.resetContinuity();

  const glm::dvec2 gridOnly       = solve( query );
  const double     gridOnlyResult = residual( solve, gridOnly, query );

  check( gridOnlyResult > 0.5,
         "precondition: the grid seed alone lands on the wrong sheet "
         "(residual " + std::to_string( gridOnlyResult ) + ")" );

  // The same query with the previous boundary point's uv available.
  solve.resetContinuity();
  solve.previousSolution_    = continuityUv;
  solve.hasPreviousSolution_ = true;

  const glm::dvec2 solved = solve( query );

  check( std::abs( solved.y - QUERY_V ) < 1.0e-3,
         "solved v returns to the correct turn" );

  check( residual( solve, solved, query ) < 1.0e-4,
         "solved uv reproduces the query point" );
}

/**
 * The accept rule, which is what makes a wrong continuity candidate cost time
 * rather than correctness.
 *
 * If the previous boundary point itself solved badly, its uv is a bad seed.
 * The retry cannot know that up front, so it accepts only on a strictly lower
 * FINAL residual. Here the candidate is deliberately on a far turn: the
 * retry fires, descends somewhere worse, and the answer must be unchanged -
 * bit for bit, not merely close.
 */
void testWorseContinuityCandidateIsRejected() {

  printf( "== a continuity candidate that descends worse is rejected\n" );

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( 3.0, 0.35, 1.0, 6 );

  RationalNurbsInverseMethod solve( surface );

  const glm::dvec3 query = solve.evaluator.point( QUERY_U, QUERY_V );

  solve.resetContinuity();

  const glm::dvec2 gridOnly = solve( query );

  // Far up the ribbon and across it, so its descent cannot reach the query's
  // neighbourhood and cannot beat the grid seed's own residual.
  const glm::dvec2 badCandidate( 0.05, 21.5 );

  const double badDistance =
    glm::distance( solve.evaluator.point( badCandidate.x, badCandidate.y ),
                   query );

  check( badDistance > bestGridDistance( solve, query ),
         "precondition: the bad candidate also loses the seed comparison, so "
         "it reaches the symmetric retry rather than the seed" );

  solve.resetContinuity();
  solve.previousSolution_    = badCandidate;
  solve.hasPreviousSolution_ = true;

  const glm::dvec2 withBadCandidate = solve( query );

  check( withBadCandidate == gridOnly,
         "the answer is unchanged when the retry ends at a higher residual" );
}

/**
 * The local grid refinement, and the same accept rule on it.
 *
 * boundSeverityToChord arms a residual-triggered re-scan of a small
 * neighbourhood of the best coarse sample. It is the tier that handles a
 * wrong-sheet seed with NO continuity candidate to fall back on - the head of
 * a trim loop, or a point whose predecessor is equally wrong.
 */
void testSeverityArmedRefinementRecoversTheSheet() {

  printf( "== severity-armed local refinement recovers the sheet\n" );

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( 3.0, 0.35, 1.0, 6 );

  RationalNurbsInverseMethod solve( surface );

  const glm::dvec3 query = solve.evaluator.point( QUERY_U, QUERY_V );

  solve.resetContinuity();

  const glm::dvec2 unarmed = solve( query );

  check( residual( solve, unarmed, query ) > 0.5,
         "precondition: unarmed, the grid seed alone lands on the wrong "
         "sheet" );

  // A chord of the order of the ribbon's own sampling. The gate is
  // 0.3 * chord, so a residual near one pitch (1.0) is far above it.
  solve.resetContinuity();
  solve.boundSeverityToChord( 0.25 );

  const glm::dvec2 armed = solve( query );

  check( std::abs( armed.y - QUERY_V ) < 1.0e-3,
         "armed, the refinement returns the correct turn" );

  check( residual( solve, armed, query ) < 1.0e-4,
         "armed, the solved uv reproduces the query point" );
}

/**
 * The refinement's false-trigger cost: a face whose points are already solved
 * cannot be moved by arming it.
 *
 * The severity gate is deliberately loose - it decides where to SPEND WORK,
 * not what to believe - so clean faces do trip it. This pins the claim that
 * doing so costs evaluations and nothing else.
 */
void testFalseTriggerLeavesTheAnswerUnchanged() {

  printf( "== a false severity trigger changes no answer\n" );

  const tinynurbs::RationalSurface3d surface =
    makeHelicalRibbon( 3.0, 0.35, 1.0, 6 );

  RationalNurbsInverseMethod solve( surface );

  // A point that solves correctly on the ordinary path: a continuity
  // candidate close enough to WIN the seed comparison, which is what a
  // finely-sampled trim loop gives every point after its first. Deliberately
  // not the near-tie above - this test is about the refinement's behaviour on
  // an answer that is already right, so it must not depend on either retry.
  const glm::dvec3 query = solve.evaluator.point( QUERY_U, QUERY_V );

  check( glm::distance( solve.evaluator.point( 0.5, 2.65 ), query ) <
           bestGridDistance( solve, query ),
         "precondition: the continuity candidate WINS the seed comparison, so "
         "this query solves on the ordinary path" );

  solve.resetContinuity();
  solve.previousSolution_    = glm::dvec2( 0.5, 2.65 );
  solve.hasPreviousSolution_ = true;

  const glm::dvec2 clean = solve( query );

  check( residual( solve, clean, query ) < 1.0e-4,
         "precondition: this query solves cleanly without the refinement" );

  // An absurdly small chord, so the gate trips on a residual that is already
  // at the convergence target.
  solve.resetContinuity();
  solve.previousSolution_    = glm::dvec2( 0.5, 2.65 );
  solve.hasPreviousSolution_ = true;
  solve.boundSeverityToChord( 1.0e-9 );

  const glm::dvec2 armed = solve( query );

  check( armed == clean,
         "the answer is unchanged when the refinement finds nothing better" );
}

}  // namespace

int main() {

  testRejectedContinuitySeedIsStillTried();
  testWorseContinuityCandidateIsRejected();
  testSeverityArmedRefinementRecoversTheSheet();
  testFalseTriggerLeavesTheAnswerUnchanged();

  if ( failures != 0 ) {
    printf( "\n%d check(s) failed\n", failures );
    return 1;
  }

  printf( "\nall checks passed\n" );
  return 0;
}
