/*
 * THE CASE THE SAMPLED DEFLECTION TEST CANNOT SEE.
 *
 * `tesselate`'s ParameterVertex overload decided whether to subdivide an edge
 * by reading the surface's departure from the chord at t = 1/4, 1/2 and 3/4
 * and keeping the worst. Inside one knot box of a bicubic tensor-product
 * surface the departure along a straight uv segment is a polynomial of degree
 * at most 6 which already vanishes at t = 0 and t = 1, so forcing it to
 * vanish at those three more points leaves
 *
 *     d( t ) = A t ( t - 1/4 )( t - 1/2 )( t - 3/4 )( t - 1 )( t - s )
 *
 * with A FREE. Every member of that family reads EXACTLY ZERO deflection and
 * is accepted, however large A is - and the family is not exotic, it is the
 * function class this code path produces.
 *
 * This file builds one such surface as a genuine Bezier patch, drives it
 * through the real `tesselate`, and shows the two answers:
 *
 *   - the sampled test accepts the chord ( no subdivision, 2 triangles out ),
 *   - the certificate bounds the departure ABOVE tolerance and refuses.
 *
 * Standalone, like refinement_progress_test.cpp: includes the headers and
 * links nothing but the Logger stubs.
 */
#include "conway_geometry/operations/deflection_certificate.h"
#include "conway_geometry/operations/tesselation_utils.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

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

using conway::geometry::CertificateOutcome;
using conway::geometry::NurbsDeflectionCertificate;
using conway::geometry::ParameterVertex;
using conway::geometry::RationalSurfaceEvaluator;
using conway::geometry::WingedEdgeMesh;

/** The adversarial departure, as a function of the diagonal parameter. */
double adversarialDeviation( double t, double amplitude, double extraRoot ) {

  return amplitude * t * ( t - 0.25 ) * ( t - 0.5 ) * ( t - 0.75 ) *
         ( t - 1.0 ) * ( t - extraRoot );
}

/**
 * A bicubic Bezier patch over the unit uv square, x = u, y = v, and z chosen
 * so that ALONG THE DIAGONAL u = v = t the height is exactly
 * `adversarialDeviation( t )`.
 *
 * The construction: write the target as monomials q_k t^k, place q_k on the
 * single (u^i, v^j) term with i + j = k listed below - each is within
 * bidegree ( 3, 3 ) - and convert that polynomial to the Bernstein grid with
 * u^i = sum_j ( C( j, i ) / C( 3, i ) ) B_j( u ), which is exact.
 *
 * The patch's four corners sit at z = 0 ( t = 0 and t = 1 are roots ), so the
 * chord across the diagonal is the straight line z = 0 and the departure from
 * it IS the height.
 */
tinynurbs::RationalSurface3d adversarialSurface(
  double amplitude,
  double extraRoot ) {

  // Monomial coefficients of the target, by expanding the roots.
  double monomial[ 7 ] = { 0.0 };

  {
    // Multiply out ( t - r ) for each root, starting from the constant 1.
    const double roots[ 6 ] = { 0.0, 0.25, 0.5, 0.75, 1.0, extraRoot };

    double accumulator[ 7 ] = { 1.0 };

    uint32_t degree = 0;

    for ( double root : roots ) {

      for ( uint32_t at = degree + 1; at > 0; --at ) {
        accumulator[ at ] =
          accumulator[ at - 1 ] - ( root * accumulator[ at ] );
      }

      accumulator[ 0 ] = -root * accumulator[ 0 ];

      ++degree;
    }

    for ( uint32_t at = 0; at < 7; ++at ) {
      monomial[ at ] = amplitude * accumulator[ at ];
    }
  }

  // ( i, j ) carrying t^k on the diagonal, one per k, all within ( 3, 3 ).
  const uint32_t powerU[ 7 ] = { 0, 1, 2, 3, 3, 3, 3 };
  const uint32_t powerV[ 7 ] = { 0, 0, 0, 0, 1, 2, 3 };

  // Degree-3 monomial -> Bernstein: row i holds u^i in the B_j basis.
  double toBernstein[ 4 ][ 4 ] = { { 0.0 } };

  const double binomial3[ 4 ] = { 1.0, 3.0, 3.0, 1.0 };

  for ( uint32_t i = 0; i < 4; ++i ) {

    for ( uint32_t j = i; j < 4; ++j ) {

      double numerator = 1.0;

      for ( uint32_t k = 0; k < i; ++k ) {
        numerator = numerator * ( j - k ) / ( k + 1 );
      }

      toBernstein[ i ][ j ] = numerator / binomial3[ i ];
    }
  }

  double height[ 4 ][ 4 ] = { { 0.0 } };

  for ( uint32_t k = 0; k < 7; ++k ) {

    for ( uint32_t a = 0; a < 4; ++a ) {

      for ( uint32_t b = 0; b < 4; ++b ) {

        height[ a ][ b ] +=
          monomial[ k ] *
          toBernstein[ powerU[ k ] ][ a ] *
          toBernstein[ powerV[ k ] ][ b ];
      }
    }
  }

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 3;
  surface.degree_v = 3;
  surface.knots_u  = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
  surface.knots_v  = surface.knots_u;

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 4; ++a ) {

    for ( uint32_t b = 0; b < 4; ++b ) {

      // x = u and y = v are themselves bicubic: their Bernstein grids are the
      // Greville abscissae, a / 3 and b / 3.
      points.push_back(
        glm::dvec3( a / 3.0, b / 3.0, height[ a ][ b ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 4, 4, points );
  surface.weights        = tinynurbs::array2( 4, 4, weights );

  return surface;
}

/** The mesh the refinement runs on: the unit uv square cut by its diagonal. */
WingedEdgeMesh< ParameterVertex > diagonalQuad(
  const RationalSurfaceEvaluator& evaluator ) {

  WingedEdgeMesh< ParameterVertex > mesh;

  const glm::dvec2 corners[ 4 ] = {
    { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

  for ( const glm::dvec2& uv : corners ) {
    mesh.makeVertex(
      ParameterVertex{ evaluator.point( uv.x, uv.y ), uv } );
  }

  // Both triangles share 0-2, the diagonal, which is therefore the only
  // INTERIOR edge - the other four are borders and `tesselate` skips those.
  mesh.makeTriangle( 0, 1, 2 );
  mesh.makeTriangle( 0, 2, 3 );

  return mesh;
}

/** Densely sampled truth, for comparing the bound against. */
double trueMaximumDeviation(
  const RationalSurfaceEvaluator& evaluator,
  uint32_t                        samples ) {

  const glm::dvec3 at0 = evaluator.point( 0.0, 0.0 );
  const glm::dvec3 at1 = evaluator.point( 1.0, 1.0 );

  double worst = 0.0;

  for ( uint32_t i = 0; i <= samples; ++i ) {

    const double t = static_cast< double >( i ) / samples;

    worst =
      std::max(
        worst,
        glm::length( evaluator.point( t, t ) -
                     ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  return worst;
}

/**
 * THE ADVERSARIAL CASE.
 *
 * `amplitude` is set so the surface departs from the diagonal chord by 0.1 -
 * a HUNDRED times the tolerance - while reading exactly zero at all three
 * sampled parameters.
 */
void certificateCatchesWhatSamplingMisses() {

  printf( "certificateCatchesWhatSamplingMisses\n" );

  constexpr double EXTRA_ROOT = 0.9;
  constexpr double TOLERANCE  = 1.0e-3;
  constexpr double PEAK       = 0.1;

  // Generous on purpose: the point of the last check is that the loop stops
  // because the BOUND came under tolerance, so the budget must not be what
  // stops it.
  constexpr int32_t CERTIFIED_BUDGET = 100000;

  // max | t ( t - 1/4 )( t - 1/2 )( t - 3/4 )( t - 1 )( t - 9/10 ) | over
  // [ 0, 1 ], from a 200,001-point sweep ( scratchpad ).
  constexpr double UNIT_PEAK = 0.0028885094447218523;

  const double amplitude = PEAK / UNIT_PEAK;

  const tinynurbs::RationalSurface3d surface =
    adversarialSurface( amplitude, EXTRA_ROOT );

  const RationalSurfaceEvaluator evaluator( surface );

  // The patch really is the polynomial it was built to be - if this fails,
  // nothing below means anything.
  double construction = 0.0;

  for ( uint32_t i = 0; i <= 1000; ++i ) {

    const double t = i / 1000.0;

    construction =
      std::max( construction,
                std::abs( evaluator.point( t, t ).z -
                          adversarialDeviation( t, amplitude, EXTRA_ROOT ) ) );
  }

  check( construction < 1.0e-12,
         "the patch's diagonal IS the adversarial polynomial (max error " +
           std::to_string( construction ) + ")" );

  const double truth = trueMaximumDeviation( evaluator, 100000 );

  check( truth > 0.09 && truth < 0.11,
         "the surface really does depart from the chord by ~0.1 (" +
           std::to_string( truth ) + ")" );

  // WHAT THE SAMPLED TEST READS. Exactly the three parameters `tesselate`
  // used, against exactly the chord it used.
  const glm::dvec3 at0 = evaluator.point( 0.0, 0.0 );
  const glm::dvec3 at1 = evaluator.point( 1.0, 1.0 );

  double sampled = 0.0;

  for ( double t : { 0.25, 0.5, 0.75 } ) {

    sampled =
      std::max( sampled,
                glm::length( evaluator.point( t, t ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( sampled < TOLERANCE,
         "the SAMPLED reading is under tolerance - it accepts this chord (" +
           std::to_string( sampled ) + " < " + std::to_string( TOLERANCE ) +
           ")" );

  check( sampled * 1000.0 < truth,
         "and it is wrong by more than three orders of magnitude" );

  // WHAT THE CERTIFICATE RETURNS.
  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, TOLERANCE * TOLERANCE );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound(
      glm::dvec2( 0.0, 0.0 ), glm::dvec2( 1.0, 1.0 ), at0, at1, bound );

  check( outcome == CertificateOutcome::Certified,
         "the certificate reaches a verdict on this chord" );

  check( bound >= truth,
         "the bound is an UPPER bound on the true departure (" +
           std::to_string( bound ) + " >= " + std::to_string( truth ) + ")" );

  check( bound > TOLERANCE,
         "and it is over tolerance - the certificate REFUSES this chord" );

  check( bound < 10.0 * truth,
         "the bound is not so loose as to be useless (" +
           std::to_string( bound / truth ) + "x the truth)" );

  // END TO END, through the real refinement loop.
  WingedEdgeMesh< ParameterVertex > sampledMesh = diagonalQuad( evaluator );

  conway::geometry::tesselate(
    sampledMesh,
    [ &evaluator ]( const glm::dvec3&, const glm::dvec2& uv ) {
      return evaluator.point( uv.x, uv.y );
    },
    CERTIFIED_BUDGET,
    TOLERANCE * TOLERANCE );

  check( sampledMesh.triangles.size() == 2,
         "WITHOUT the certificate `tesselate` accepts the diagonal as it is (" +
           std::to_string( sampledMesh.triangles.size() ) + " triangles)" );

  WingedEdgeMesh< ParameterVertex > certifiedMesh = diagonalQuad( evaluator );

  const NurbsDeflectionCertificate meshCertificate(
    evaluator, false, 0.0, 0.0, TOLERANCE * TOLERANCE );

  conway::geometry::tesselate(
    certifiedMesh,
    [ &evaluator ]( const glm::dvec3&, const glm::dvec2& uv ) {
      return evaluator.point( uv.x, uv.y );
    },
    CERTIFIED_BUDGET,
    TOLERANCE * TOLERANCE,
    meshCertificate );

  check( certifiedMesh.triangles.size() > 2,
         "WITH it the diagonal is subdivided (" +
           std::to_string( certifiedMesh.triangles.size() ) + " triangles)" );

  // And the refinement it produces actually converges: every surviving
  // interior edge is within tolerance by the same bound that forced the
  // splits, so the loop stopped because it was DONE, not because it ran out.
  check( certifiedMesh.triangles.size() <
           static_cast< size_t >( CERTIFIED_BUDGET ),
         "and it stops on the BOUND, not on the budget (" +
           std::to_string( certifiedMesh.triangles.size() ) + " of " +
           std::to_string( CERTIFIED_BUDGET ) + " triangles)" );
}

/**
 * The bound has to be an upper bound on chords that are NOT adversarial too,
 * or it would be an assertion that only ever fires on one contrived input.
 *
 * A sphere-like rational patch ( non-unit weights ) exercises the homogeneous
 * path, which is a different arm of `boundPiece` from the one above.
 */
void boundHoldsOnARationalPatch() {

  printf( "boundHoldsOnARationalPatch\n" );

  // A quarter cylinder as an exact rational quadratic in u, linear in v: the
  // classic weight-1, 1/sqrt(2), 1 construction, which is NOT polynomial.
  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  const double root = std::sqrt( 0.5 );

  std::vector< glm::dvec3 > points = {
    { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
    { 1.0, 1.0, 0.0 }, { 1.0, 1.0, 1.0 },
    { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 1.0 } };

  std::vector< double > weights = {
    1.0, 1.0, root, root, 1.0, 1.0 };

  surface.control_points = tinynurbs::array2( 3, 2, points );
  surface.weights        = tinynurbs::array2( 3, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  check( !evaluator.isPolynomial(),
         "the patch is recognised as RATIONAL, so the homogeneous arm runs" );

  // It is a unit quarter circle, so the surface point is at radius one.
  double radialError = 0.0;

  for ( uint32_t i = 0; i <= 200; ++i ) {

    const glm::dvec3 point = evaluator.point( i / 200.0, 0.5 );

    radialError =
      std::max( radialError,
                std::abs( std::hypot( point.x, point.y ) - 1.0 ) );
  }

  check( radialError < 1.0e-12,
         "and it is the exact quarter circle it was built as" );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-6 );

  // INSTRUMENT CHECK against an answer known in closed form: the mid-arc of a
  // chord subtending 90 degrees on a unit circle stands 1 - cos( 45 degrees )
  // = 0.2928932... off that chord. If the dense sweep below does not find
  // that number on the full span, the sweep is not measuring what it says.
  {
    const glm::dvec3 at0 = evaluator.point( 0.0, 0.5 );
    const glm::dvec3 at1 = evaluator.point( 1.0, 0.5 );

    const double middle =
      glm::length( evaluator.point( 0.5, 0.5 ) - ( ( at0 + at1 ) * 0.5 ) );

    check( std::abs( middle - ( 1.0 - std::cos( 0.25 * 3.14159265358979323846 ) ) ) <
             1.0e-12,
           "the mid-arc stands 1 - cos( 45 deg ) off the chord, as it must (" +
             std::to_string( middle ) + ")" );
  }

  // Chords of a range of lengths, each checked against dense truth.
  bool alwaysAbove = true;
  bool everTight   = false;

  for ( double span : { 1.0, 0.5, 0.25, 0.125, 0.0625 } ) {

    const glm::dvec2 uv0( 0.0, 0.5 );
    const glm::dvec2 uv1( span, 0.5 );

    const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
    const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

    double truth = 0.0;

    for ( uint32_t i = 0; i <= 20000; ++i ) {

      const double t = static_cast< double >( i ) / 20000.0;

      truth =
        std::max( truth,
                  glm::length(
                    evaluator.point( uv0.x + ( t * span ), 0.5 ) -
                    ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
    }

    double bound = 0.0;

    const CertificateOutcome outcome =
      certificate.bound( uv0, uv1, at0, at1, bound );

    printf( "    span %-8g truth %-14g bound %-14g ratio %g\n",
            span, truth, bound, bound / truth );

    alwaysAbove = alwaysAbove && ( outcome == CertificateOutcome::Certified ) &&
                  ( bound >= truth );

    everTight = everTight || ( bound < 4.0 * truth );
  }

  check( alwaysAbove,
         "the rational bound is above the dense-sampled truth at every span" );

  check( everTight, "and within 4x of it somewhere" );
}

/**
 * The bound must SHRINK as the chord does, or the refinement loop would never
 * terminate - which is the failure mode the whole surfaceAt cache exists to
 * remove, and it would be a poor trade to reintroduce it here.
 */
void boundConvergesUnderSubdivision() {

  printf( "boundConvergesUnderSubdivision\n" );

  const tinynurbs::RationalSurface3d surface =
    adversarialSurface( 34.61993180694947, 0.9 );

  const RationalSurfaceEvaluator evaluator( surface );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-12 );

  double previous   = std::numeric_limits< double >::infinity();
  double lastRatio  = 0.0;
  bool   monotone   = true;

  for ( uint32_t level = 0; level < 12; ++level ) {

    const double span = std::ldexp( 1.0, -static_cast< int >( level ) );

    const glm::dvec2 uv0( 0.0, 0.0 );
    const glm::dvec2 uv1( span, span );

    double bound = 0.0;

    certificate.bound(
      uv0, uv1,
      evaluator.point( uv0.x, uv0.y ),
      evaluator.point( uv1.x, uv1.y ),
      bound );

    printf( "    span %-12g bound %-16g ratio %g\n", span, bound,
            previous / bound );

    if ( level > 0 ) {

      monotone  = monotone && ( bound < previous );
      lastRatio = previous / bound;
    }

    previous = bound;
  }

  check( monotone, "the bound falls at every halving" );

  // ORDER TWO is what termination rests on, and it is what the bound
  // approaches rather than what it starts at: a long chord spanning several
  // sign changes of the departure has a loose coefficient hull, so the first
  // halvings gain less than the asymptotic factor of four. The measured
  // sequence over twelve levels runs 3.56, 1.89, 2.00, 2.80, 3.36, 3.68,
  // 3.84, 3.92, 3.96, ... -> 4.
  check( lastRatio > 3.9,
         "and approaches the factor of four an O( h^2 ) departure gives (" +
           std::to_string( lastRatio ) + ")" );

  check( previous < 3.0e-6,
         "reaching 3e-6 within twelve halvings (" +
           std::to_string( previous ) + ")" );
}


/**
 * KNOT-BOX CLIPPING, and what it is for.
 *
 * The polynomial argument holds INSIDE one knot box and nowhere else. A chord
 * crossing knot lines meets a different polynomial on each side, and
 * interpolating across the join fits a curve to a function that is not there.
 *
 * The cheapest surface that shows it: degree ONE in u, so the surface is a
 * fold-line in u, with knots at 1/4, 1/2 and 3/4 and heights alternating
 * 0, +h, 0, -h, 0. Total degree is two, so the uncut reading would interpolate
 * at the three Chebyshev-Lobatto nodes 0, 1/2 and 1 - and the zig-zag is zero
 * at all three. Cut at the knots it is bounded span by span and the peaks are
 * found.
 */
void knotClippingCatchesAZigZag() {

  printf( "knotClippingCatchesAZigZag\n" );

  constexpr double HEIGHT = 0.25;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  // The bump sits on an OUTER span and the span containing u = 1/2 is FLAT.
  // That matters for what the uncut reading does: every node of a piece is
  // evaluated as the polynomial of the span the piece's MIDPOINT falls in
  // (see boundPiece), so an uncut chord is read entirely as the flat middle
  // span and comes out at zero - an under-estimate, which is the direction a
  // bound may not err in. A symmetric zig-zag instead makes the uncut
  // reading EXTRAPOLATE a steep span and over-estimate, which hides the
  // defect; measured, it returned 0.5 against a truth of 0.25.
  const double heights[ 5 ] = { 0.0, HEIGHT, 0.0, 0.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 5; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      points.push_back(
        glm::dvec3( a / 4.0, static_cast< double >( b ), heights[ a ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 5, 2, points );
  surface.weights        = tinynurbs::array2( 5, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const glm::dvec2 uv0( 0.0, 0.5 );
  const glm::dvec2 uv1( 1.0, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 40000; ++i ) {

    const double t = static_cast< double >( i ) / 40000.0;

    truth =
      std::max( truth,
                glm::length( evaluator.point( t, 0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( std::abs( truth - HEIGHT ) < 1.0e-9,
         "the fold really does stand " + std::to_string( HEIGHT ) +
           " off its chord (" + std::to_string( truth ) + ")" );

  check( std::abs( evaluator.point( 0.625, 0.5 ).z ) < 1.0e-12,
         "and the span around u = 1/2 is flat, so an uncut reading of the "
         "whole chord lands at zero" );

  // The three parameters an UNCUT degree-two reading would interpolate at.
  double atNodes = 0.0;

  for ( double t : { 0.0, 0.5, 1.0 } ) {

    atNodes =
      std::max( atNodes,
                glm::length( evaluator.point( t, 0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( atNodes < 1.0e-12,
         "and it is invisible at the nodes a single uncut span would use (" +
           std::to_string( atNodes ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-8 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome == CertificateOutcome::Certified,
         "the certificate reaches a verdict across four knot spans" );

  check( bound >= truth,
         "and bounds the fold from above (" + std::to_string( bound ) +
           " >= " + std::to_string( truth ) + ")" );
}

/**
 * THE PERIODIC CHART'S CUT.
 *
 * The b-spline front end evaluates `point( wrapChartU( u ), v )`, so a chord
 * whose two ends sit on different sheets of the strip is a DISCONTINUOUS
 * function of t however smooth the surface is. Cutting at the sheet boundary
 * restores the affine map the bound needs; without the cut the certificate
 * can only decline, which on `Right_Hand.step` was 22,634 of 216,757
 * candidates.
 */
void wrappedChordIsStillCertified() {

  printf( "wrappedChordIsStillCertified\n" );

  constexpr double PERIOD = 1.0;

  // Bicubic Bezier over [ 0, 1 ]^2 with a plain bump, so the wrapped
  // evaluation has something to depart from its chord WITH.
  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 3;
  surface.degree_v = 3;
  surface.knots_u  = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
  surface.knots_v  = surface.knots_u;

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 4; ++a ) {

    for ( uint32_t b = 0; b < 4; ++b ) {

      // CLOSED in u: the first and last control rows coincide, so
      // S( 0, v ) == S( 1, v ) and the wrapped evaluation is continuous.
      // What is discontinuous is t -> u, which is what the cut is for.
      const double across[ 4 ] = { 0.0,  1.0, -1.0, 0.0 };
      const double height[ 4 ] = { 0.0, 0.35, 0.35, 0.0 };

      points.push_back(
        glm::dvec3( across[ a ], b / 3.0, height[ a ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 4, 4, points );
  surface.weights        = tinynurbs::array2( 4, 4, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = [ PERIOD ]( double u ) {

    const double offset = std::fmod( u, PERIOD );

    return offset < 0.0 ? offset + PERIOD : offset;
  };

  // Straddles the cut: 0.8 is on sheet 0, 1.3 on sheet 1.
  const glm::dvec2 uv0( 0.8, 0.5 );
  const glm::dvec2 uv1( 1.3, 0.5 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 40000; ++i ) {

    const double t = static_cast< double >( i ) / 40000.0;

    const double u = uv0.x + ( t * ( uv1.x - uv0.x ) );

    truth =
      std::max( truth,
                glm::length( evaluator.point( wrap( u ), uv0.y ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( truth > 0.05,
         "the wrapped chord really does depart from its chord (" +
           std::to_string( truth ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, true, 0.0, PERIOD, 1.0e-8 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  printf( "    wrap truth %g bound %g ratio %g\n", truth, bound, bound / truth );

  check( outcome == CertificateOutcome::Certified,
         "a chord across the chart cut is still CERTIFIED, not declined" );

  check( bound >= truth,
         "and bounded from above (" + std::to_string( bound ) + " >= " +
           std::to_string( truth ) + ")" );

  // TIGHTNESS is what pins the cut, not the inequality above. Read as one
  // affine piece the certificate evaluates the surface OUTSIDE its knot
  // domain for the part of the chord on the far sheet - it extrapolates the
  // end span - and the bound it returns is about a different function that
  // happens to sit above this one. Measured: 1.11x the truth with the cut,
  // 3.99x without it.
  check( bound < 1.5 * truth,
         "and the bound is tight, which reading it as one piece is not (" +
           std::to_string( bound / truth ) + "x)" );
}


/**
 * A FULL-MULTIPLICITY INTERIOR KNOT IS SAMPLED FROM THE PIECE'S OWN SIDE.
 *
 * A valid B-spline may repeat an interior knot `degree + 1` times, and there
 * the two spans are different polynomials with a STEP between them.
 * `findSpan` resolves a parameter sitting exactly on the knot to the span on
 * its right, and the Chebyshev-Lobatto nodes always include both endpoints of
 * a piece - so the piece to the LEFT of such a knot would be interpolated
 * through one sample taken from the polynomial to the RIGHT. The interpolant
 * is then fitted through a point that is not on the curve it is certifying,
 * and its coefficient hull bounds neither.
 *
 * This is the smallest case that shows it, and the numbers are the ones a
 * reviewer of bldrs-ai/conway-geom#214 worked out by hand: degrees ( 1, 1 ),
 * an interior knot of multiplicity 2, a left-hand deviation rising linearly
 * to 1 and a right-hand value at the knot of 2/3. The samples come out
 * [ 0, 1/2, 2/3 ], the degree-2 Bernstein coefficients [ 0, 2/3, 2/3 ], and
 * the bound 2/3 - so any tolerance between 2/3 and 1 was accepted when it
 * should have been refused.
 *
 * No surface in the smoke corpus has such a knot ( measured: 0 of 127
 * b-spline faces across `Right_Hand.step` and `nist_ctc_02_asme1_rc.stp` ),
 * which is exactly why this needs a test rather than a corpus run.
 */
void fullMultiplicityKnotIsSampledOneSided() {

  printf( "fullMultiplicityKnotIsSampledOneSided\n" );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  // 0.5 twice: multiplicity degree + 1, so the surface STEPS there.
  surface.knots_u = { 0.0, 0.0, 0.5, 0.5, 1.0, 1.0 };
  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  // Greville abscissae for degree 1 are the knots themselves, so x stays
  // continuous across the break and the step is purely in z.
  const double across[ 4 ] = { 0.0, 0.5, 0.5, 1.0 };
  const double height[ 4 ] = { 0.0, 1.0, 2.0 / 3.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 4; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      points.push_back(
        glm::dvec3( across[ a ], static_cast< double >( b ), height[ a ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 4, 2, points );
  surface.weights        = tinynurbs::array2( 4, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const double leftLimit =
    evaluator.point( std::nextafter( 0.5, 0.0 ), 0.5 ).z;

  const double rightValue = evaluator.point( 0.5, 0.5 ).z;

  check( std::abs( leftLimit - 1.0 ) < 1.0e-9 &&
           std::abs( rightValue - ( 2.0 / 3.0 ) ) < 1.0e-12,
         "the surface really does step at the knot, 1 to 2/3 (" +
           std::to_string( leftLimit ) + " -> " +
           std::to_string( rightValue ) + ")" );

  const glm::dvec2 uv0( 0.0, 0.5 );
  const glm::dvec2 uv1( 1.0, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 2000000; ++i ) {

    const double t = static_cast< double >( i ) / 2000000.0;

    truth =
      std::max( truth,
                glm::length( evaluator.point( t, 0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( truth > 0.999,
         "and the true departure from the chord approaches 1 (" +
           std::to_string( truth ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-12 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome == CertificateOutcome::Certified ||
           outcome == CertificateOutcome::Inconclusive,
         "the certificate reaches a verdict" );

  // THE ASSERTION. Reading the last node from the right-hand polynomial
  // returns 2/3 here, which is below the truth - a bound that is not one.
  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "a bound across a discontinuous knot is still an UPPER bound (" +
           std::to_string( bound ) + " vs " + std::to_string( truth ) + ")" );
}

/**
 * THE ERROR TERM IS GOVERNED BY THE CONTROL VALUES, NOT THE CHORD ENDPOINTS.
 *
 * The first version of this file scaled its floating-point inflation by
 * `max( |S0|, |S1| )`. That is not an error bound for what `boundPiece`
 * returns: an evaluation accumulates `basis * controlPointW` over
 * ( dU + 1 )( dV + 1 ) terms, so its absolute error scales with the largest
 * CONTROL value, and a patch can have control values far above the surface
 * values at the two ends of one chord. The reviewer of
 * bldrs-ai/conway-geom#214 raised it, and it is right as an analysis even
 * though no falsifying input could be built ( see the report ).
 *
 * What is checkable is the behaviour the corrected term produces: on a patch
 * whose control values sit 1e12 above its chord endpoints, and with a
 * tolerance fine enough that the evaluation error is a large fraction of it,
 * the certificate must DECLINE rather than certify. The endpoint-scaled
 * reading certifies it, because 1e-13 is under any gate.
 */
void errorTermTracksTheControlValues() {

  printf( "errorTermTracksTheControlValues\n" );

  constexpr double HUGE_CONTROL = 1.0e12;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 3;
  surface.degree_v = 3;
  surface.knots_u  = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
  surface.knots_v  = surface.knots_u;

  // Corners at zero - so both chord endpoints read O( 1 ) - with the two
  // interior rows at +/- 1e12, which cancel in the middle.
  const double height[ 4 ] = { 0.0, HUGE_CONTROL, -HUGE_CONTROL, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 4; ++a ) {

    for ( uint32_t b = 0; b < 4; ++b ) {

      points.push_back( glm::dvec3( a / 3.0, b / 3.0, height[ a ] ) );
      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 4, 4, points );
  surface.weights        = tinynurbs::array2( 4, 4, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  // A SHORT chord in the corner. Short, because that is what pulls the chord
  // endpoints far below the control values: on a long chord the surface
  // reaches the control magnitude somewhere and the endpoint-scaled reading
  // is accidentally the right size, which is why this defect needs a
  // deliberately local chord to show at all.
  const glm::dvec2 uv0( 0.0, 0.0 );
  const glm::dvec2 uv1( 1.0e-6, 1.0e-6 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  const double endpointScale =
    std::max( glm::length( at0 ), glm::length( at1 ) );

  check( endpointScale < evaluator.maxHomogeneousNorm() * 1.0e-5,
         "the chord endpoints sit far below the control values (" +
           std::to_string( endpointScale ) + " against " +
           std::to_string( evaluator.maxHomogeneousNorm() ) + ")" );

  check( evaluator.maxHomogeneousNorm() > 1.0e11,
         "while the control values are 1e12 (" +
           std::to_string( evaluator.maxHomogeneousNorm() ) + ")" );

  // A tolerance fine enough that the real evaluation error - about
  // 8 * 16 * eps * 1e12, amplified by ||M||inf = 46.2 - is a large fraction
  // of it. The endpoint-scaled reading computes 1e-13 here and sails through.
  const double tolerance = 1.0;

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome == CertificateOutcome::Inconclusive,
         "the certificate DECLINES rather than certifying at an accuracy it "
         "cannot deliver" );
}

}  // namespace

int main() {

  certificateCatchesWhatSamplingMisses();
  fullMultiplicityKnotIsSampledOneSided();
  errorTermTracksTheControlValues();
  knotClippingCatchesAZigZag();
  wrappedChordIsStillCertified();
  boundHoldsOnARationalPatch();
  boundConvergesUnderSubdivision();

  printf( "%s\n", failures == 0 ? "PASS" : "FAIL" );

  return failures == 0 ? 0 : 1;
}
