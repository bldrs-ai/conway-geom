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


/*
 * THE NEXT TWO TESTS ARE A PAIR, AND THE CONTRAST BETWEEN THEM IS THE POINT.
 *
 * Both build the same shape: a degree-1 surface with full-multiplicity breaks
 * at two nearby knots, flat at zero on the spans either side, and standing 5
 * off the chord only on the narrow span between them. Both are driven by a
 * chord from -1e9 to 1e9, long enough that the two knots round to the SAME
 * chord parameter. They differ in one number - how far apart the two knots
 * are - and that number decides which side of a real boundary the case falls
 * on:
 *
 *   span 1e-9 wide,       the collapse is a property of the CHORD being
 *   flat across itself    long, and nothing else. The box is still
 *                         enumerated, still carries its own knot-derived
 *                         corners, and is CERTIFIED. Nothing is lost and
 *                         nothing has to be declined.
 *
 *   span ONE ULP wide,    there is no representable parameter strictly
 *   surface CLIMBS        inside the span, so no chord - however short - can
 *   across it             put an interpolation node where the climb happens.
 *                         DECLINED, on the long chord and on every shorter
 *                         one.
 *
 * THE LINE BETWEEN THEM IS NOT WIDTH. A span an ulp wide that the surface is
 * FLAT across is bounded perfectly well: nodes that land in the wrong places
 * still read the right value, and the certificate's own gradient term says
 * so, because the term is built from the span's control points and they do
 * not differ. What cannot be bounded is a span too narrow to sample that the
 * surface also MOVES across - narrowness and variation together, neither
 * alone.
 *
 * That matters for how the refusal behaves. The first case never needs one.
 * The second can never be cleared by subdividing, because subdividing
 * shortens the chord and does not widen the span, so it is a stated
 * limitation rather than a guard that resolves itself. An earlier version of
 * this file asserted the second case WAS self-clearing; it passed only
 * because the narrow span was built flat, which is precisely the
 * distinction above - and giving that span a varying height is
 * `collapsedNodeParametersAreDeclined` below.
 */

/**
 * The shared shape - see the note above. `gap` is the only difference.
 *
 * Returned by value and BOUND TO A NAMED LOCAL at both call sites, never
 * handed straight to an evaluator: `RationalSurfaceEvaluator` keeps a
 * REFERENCE to the surface it is built from, so a temporary would dangle.
 */
tinynurbs::RationalSurface3d twoBreakSurface(
  double gap, double heightIn, double heightOut ) {

  const double low  = -1.0e9;
  const double high =  1.0e9;

  const double firstKnot  = 1.0;
  const double secondKnot =
    ( gap > 0.0 ) ?
      ( firstKnot + gap ) :
      std::nextafter( firstKnot, std::numeric_limits< double >::infinity() );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  // Multiplicity 2 = degree + 1 at BOTH knots, so the span between them is
  // independent of its neighbours.
  surface.knots_u = { low, low, firstKnot, firstKnot,
                      secondKnot, secondKnot, high, high };
  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  // Greville abscissae for degree 1 are the knots themselves, so x tracks u.
  const double across[ 6 ] = { low, firstKnot, firstKnot,
                               secondKnot, secondKnot, high };

  // The narrow span's two control heights. Equal makes it a flat shelf
  // standing off the chord; different makes the surface CLIMB across a span
  // that may be too narrow to sample.
  const double tall[ 6 ] = { 0.0, 0.0, heightIn, heightOut, 0.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 6; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      points.push_back(
        glm::dvec3( across[ a ], static_cast< double >( b ), tall[ a ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 6, 2, points );
  surface.weights        = tinynurbs::array2( 6, 2, weights );

  return surface;
}

/**
 * A NARROW SPAN WHOSE CHORD PARAMETERS COLLAPSE IS STILL CERTIFIED.
 *
 * The boxes used to be recovered by computing a cut parameter for every knot
 * the chord crossed, sorting them and stepping between consecutive ones - so
 * two distinct knots that rounded to one parameter produced one cut instead
 * of two, the zero-length step between them was skipped, and the span between
 * them was never certified while the result still claimed to cover the whole
 * chord. That was the third finding on bldrs-ai/conway-geom#214.
 *
 * The walk enumerates boxes by SPAN INDEX instead, and a box carries the
 * exact knot values that bound it, so a box whose chord parameters round
 * together is still a box with a real uv interval and is bounded like any
 * other. The defect is not guarded here - it is unrepresentable.
 *
 * Note what the sampled reading does, because it is the point of the whole
 * file: a MILLION-point sweep of this chord returns 2.4e-7. The span that
 * carries the departure is 1e-9 of `u` wide against a chord of 2e9, so no
 * sweep of any density lands in it. Only walking the knots finds it at all.
 */
void narrowSpanWithCollapsedParametersIsCertified() {

  printf( "narrowSpanWithCollapsedParametersIsCertified\n" );

  constexpr double LOW    = -1.0e9;
  constexpr double HIGH   =  1.0e9;
  constexpr double HEIGHT =  5.0;
  constexpr double GAP    =  1.0e-9;

  const double firstKnot  = 1.0;
  const double secondKnot = firstKnot + GAP;

  check( ( ( firstKnot - LOW ) / ( HIGH - LOW ) ) ==
           ( ( secondKnot - LOW ) / ( HIGH - LOW ) ),
         "the two knots map to the SAME chord parameter (" +
           std::to_string( ( firstKnot - LOW ) / ( HIGH - LOW ) ) + ")" );

  // Named, not a temporary - the evaluator keeps a reference to it.
  const tinynurbs::RationalSurface3d surface =
    twoBreakSurface( GAP, HEIGHT, HEIGHT );

  const RationalSurfaceEvaluator evaluator( surface );

  const glm::dvec2 uv0( LOW, 0.5 );
  const glm::dvec2 uv1( HIGH, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  const auto departureAt = [ & ]( double u ) {

    const double t = ( u - LOW ) / ( HIGH - LOW );

    return glm::length( evaluator.point( u, 0.5 ) -
                        ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) );
  };

  check( std::abs( departureAt( firstKnot ) - HEIGHT ) < 1.0e-6,
         "the narrow span stands " + std::to_string( HEIGHT ) +
           " off the chord (" +
           std::to_string( departureAt( firstKnot ) ) + ")" );

  double swept = 0.0;

  for ( uint32_t i = 0; i <= 1000000; ++i ) {

    swept =
      std::max( swept,
                departureAt( LOW + ( ( HIGH - LOW ) *
                                     ( static_cast< double >( i ) /
                                       1000000.0 ) ) ) );
  }

  check( swept < 1.0e-6,
         "and a million-point sweep of the chord cannot find it (" +
           std::to_string( swept ) + ")" );

  // Loose enough that the ROUNDING term does not decline this chord on its
  // own - the chord is 2e9 long, so that term is about 2.7e-5 here. Without
  // the gap the assertion below would fail for a reason that has nothing to
  // do with what it is pinning, which is the opposite mistake to the one
  // this file's history is full of.
  constexpr double TOLERANCE = 1.0e-2;

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, TOLERANCE * TOLERANCE );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome == CertificateOutcome::Certified,
         "the certificate reaches a verdict on the whole chord" );

  check( bound >= HEIGHT,
         "and its bound covers the span whose parameters collapsed (" +
           std::to_string( bound ) + " >= " +
           std::to_string( HEIGHT ) + ")" );
}

/**
 * A KNOT SPAN BELOW PARAMETER RESOLUTION IS DECLINED, PERMANENTLY.
 *
 * This is the other side of the pair. The span is ONE ULP of `u` wide AND
 * the surface climbs across it, from 0 at one end to 5 at the other. There
 * is no representable parameter strictly inside, so no chord can put a node
 * where the climb happens, and the gradient term - built from this span's
 * own control points, which now differ - is enormous. The certificate
 * declines.
 *
 * Unlike every other refusal in this file, THIS ONE DOES NOT CLEAR.
 * Subdividing shortens the chord; it does not widen the knot span. The face
 * will spend its triangle budget and stop, which is the correct behaviour
 * for a span nobody can bound, but it is a limitation to know about rather
 * than a guard that resolves itself.
 */
void subUlpKnotSpanIsDeclined() {

  printf( "subUlpKnotSpanIsDeclined\n" );

  constexpr double HEIGHT = 5.0;

  const double firstKnot  = 1.0;
  const double secondKnot =
    std::nextafter( firstKnot, std::numeric_limits< double >::infinity() );

  check( std::nextafter( firstKnot,
                         std::numeric_limits< double >::infinity() ) ==
           secondKnot,
         "the span is exactly one ulp wide - nothing lies strictly inside" );

  // Named, not a temporary - the evaluator keeps a reference to it.
  const tinynurbs::RationalSurface3d surface =
    twoBreakSurface( 0.0, 0.0, HEIGHT );

  const RationalSurfaceEvaluator evaluator( surface );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-4 );

  // The long chord it starts on, then the ones subdividing would produce.
  // Every one of them is declined: the refusal is about the SPAN, not the
  // chord, so making the chord shorter cannot reach it.
  const double reaches[ 5 ] =
    { 1.0e9, 1.0e3, 1.0, 1.0e-3, 1.0e-6 };

  bool declinedEverywhere = true;

  for ( double reach : reaches ) {

    const glm::dvec2 uv0( firstKnot - reach, 0.5 );
    const glm::dvec2 uv1( secondKnot + reach, 0.5 );

    const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
    const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

    double bound = 0.0;

    const CertificateOutcome outcome =
      certificate.bound( uv0, uv1, at0, at1, bound );

    printf( "    chord reach %-10g -> %s\n", reach,
            outcome == CertificateOutcome::Certified ? "Certified" :
              ( outcome == CertificateOutcome::Inconclusive ?
                  "Inconclusive" : "Unsupported" ) );

    declinedEverywhere =
      declinedEverywhere &&
      ( outcome == CertificateOutcome::Inconclusive );
  }

  check( declinedEverywhere,
         "a span an ulp wide is declined on every chord, however short - "
         "this refusal is NOT self-clearing" );
}

/**
 * NODES THAT COLLAPSE ONTO THE BOX'S CORNERS ARE DECLINED.
 *
 * A node's uv is computed as `from + s * ( to - from )`, and that sum rounds
 * to an ulp of the PARAMETER MAGNITUDE, not of the box width. On a box narrow
 * relative to where it sits, distinct Chebyshev-Lobatto parameters land on
 * the same representable uv - so the values handed to the interpolation are
 * of the right function at the wrong places, and the coefficient hull bounds
 * a curve that is not the one being certified.
 *
 * This is the fourth finding on bldrs-ai/conway-geom#214, built to its
 * author's recipe: a quadratic u span from 1e9 to the SECOND nextafter value,
 * linear in v, Bezier heights [ 0, 2, 0 ]. The span is two ulps wide, the one
 * representable interior parameter carries height 1, and the degree-3 nodes
 * round onto the zero-height ends. Measured before this was detected: bound
 * 6.2e-5 against a true departure of 1, certified at a tolerance of 0.01.
 *
 * IT IS NOT THE INDEX WALK THAT CLOSES THIS, and the distinction matters.
 * The walk guarantees the box EXISTS and is named by its span indices; it
 * says nothing about whether the nodes inside the box can be placed. What
 * closes it is the gradient term in boundPiece's error: the surface's
 * sensitivity per unit parameter, bounded from the control net, times the
 * parameter's own rounding. On this patch that product is enormous and the
 * error gate declines the chord.
 */
void collapsedNodeParametersAreDeclined() {

  printf( "collapsedNodeParametersAreDeclined\n" );

  const double first  = 1.0e9;
  const double second =
    std::nextafter( first, std::numeric_limits< double >::infinity() );
  const double third =
    std::nextafter( second, std::numeric_limits< double >::infinity() );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 1;
  surface.knots_u  = { first, first, first, third, third, third };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  const double across[ 3 ] = { first, second, third };
  const double height[ 3 ] = { 0.0, 2.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 3; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      points.push_back(
        glm::dvec3( across[ a ], static_cast< double >( b ), height[ a ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 3, 2, points );
  surface.weights        = tinynurbs::array2( 3, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  check( std::abs( evaluator.point( second, 0.5 ).z - 1.0 ) < 1.0e-9,
         "the one representable interior parameter carries height 1 (" +
           std::to_string( evaluator.point( second, 0.5 ).z ) + ")" );

  check( evaluator.point( first, 0.5 ).z == 0.0 &&
           evaluator.point( third, 0.5 ).z == 0.0,
         "and both ends of the span are zero" );

  const glm::dvec2 uv0( first, 0.5 );
  const glm::dvec2 uv1( third, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  const double middleT = ( second - first ) / ( third - first );

  const double truth =
    glm::length( evaluator.point( second, 0.5 ) -
                 ( ( at0 * ( 1.0 - middleT ) ) + ( at1 * middleT ) ) );

  check( truth > 0.99,
         "so the true departure from the chord is about 1 (" +
           std::to_string( truth ) + ")" );

  constexpr double TOLERANCE = 1.0e-2;

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, TOLERANCE * TOLERANCE );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "a box whose nodes cannot be placed is never certified below the "
         "truth (" + std::to_string( bound ) + " vs " +
           std::to_string( truth ) + ")" );

  check( outcome == CertificateOutcome::Inconclusive,
         "and it is declined outright" );
}


/**
 * A CHORD RUNNING DOWNWARD FROM A DISCONTINUOUS KNOT STARTS IN THE SPAN IT
 * IS ABOUT TO TRAVERSE.
 *
 * `directionalSpan` is the ONE span lookup left in the certificate - every
 * other span is reached by incrementing it - so it is the one place where
 * "which polynomial is this" can still be answered wrongly. `findSpan`
 * resolves a parameter sitting exactly on a knot to the span on its RIGHT,
 * which is what a chord climbing away from that knot wants and the opposite
 * of what a chord descending from it wants.
 *
 * Same surface as `fullMultiplicityKnotIsSampledOneSided`, walked downward
 * from the knot: the chord leaves u = 1/2 heading for u = 0, so every point
 * on it belongs to the LEFT polynomial, while the endpoint the caller cached
 * came from the right one. The departure is the difference between them.
 */
void descendingChordStartsInTheSpanItTraverses() {

  printf( "descendingChordStartsInTheSpanItTraverses\n" );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { 0.0, 0.0, 0.5, 0.5, 1.0, 1.0 };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  const double across[ 4 ] = { 0.0, 0.5, 0.5, 1.0 };
  const double height[ 4 ] = { 0.0, 1.0, 2.0 / 3.0, 2.0 / 3.0 };

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

  // DOWNWARD: starts exactly on the knot, heads away from it.
  const glm::dvec2 uv0( 0.5, 0.5 );
  const glm::dvec2 uv1( 0.0, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  check( std::abs( at0.z - ( 2.0 / 3.0 ) ) < 1.0e-12,
         "the cached start point came from the RIGHT polynomial (" +
           std::to_string( at0.z ) + ")" );

  double truth = 0.0;

  for ( uint32_t i = 1; i <= 200000; ++i ) {

    const double t = static_cast< double >( i ) / 200000.0;

    truth =
      std::max( truth,
                glm::length( evaluator.point( uv0.x + ( t * ( uv1.x - uv0.x ) ),
                                              0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( truth > 0.33,
         "and the chord departs from it by about a third (" +
           std::to_string( truth ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-12 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "the bound covers it (" + std::to_string( bound ) + " vs " +
           std::to_string( truth ) + ")" );

  check( outcome == CertificateOutcome::Certified,
         "and the chord is certified rather than declined" );
}


/**
 * WHEN A u BOUNDARY AND A v BOUNDARY TIE, BOTH ORDERS ARE BOUNDED.
 *
 * `tU` and `tV` are computed from exact knots but carried in doubles, so a
 * difference below their own rounding says nothing about which boundary the
 * chord reaches first. The walk used to always take u, which steps
 * `( oldU, oldV ) -> ( newU, oldV ) -> ( newU, newV )` and never bounds
 * `( oldU, newV )`.
 *
 * With full multiplicity on BOTH axes that skipped patch is an independent
 * polynomial - not a near-point sliver - so it can carry any departure at
 * all while every patch that IS bounded reads zero. This file previously
 * argued no such case could be built, and named the condition under which
 * the argument would fail; this is that condition. It is the fifth finding
 * on bldrs-ai/conway-geom#214, built to its author's recipe: a chord from
 * ( -1e9, -1e9 ) to ( 1e9, 1e9 ), a v knot at 1.0 and a u knot at
 * nextafter( 1.0 ), with the height living only on the omitted combination.
 *
 * Measured before the fix: bound 4.8e-5 against a true departure of 7, and
 * a two-million-point sweep of the chord finds only 5.1e-7.
 */
void tiedCrossAxisBoundariesCoverBothOrders() {

  printf( "tiedCrossAxisBoundariesCoverBothOrders\n" );

  constexpr double LOW    = -1.0e9;
  constexpr double HIGH   =  1.0e9;
  constexpr double HEIGHT =  7.0;

  const double vKnot = 1.0;
  const double uKnot =
    std::nextafter( vKnot, std::numeric_limits< double >::infinity() );

  check( ( ( vKnot - LOW ) / ( HIGH - LOW ) ) ==
           ( ( uKnot - LOW ) / ( HIGH - LOW ) ),
         "the u and v crossings land on the SAME chord parameter (" +
           std::to_string( ( vKnot - LOW ) / ( HIGH - LOW ) ) + ")" );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;
  surface.knots_u  = { LOW, LOW, uKnot, uKnot, HIGH, HIGH };
  surface.knots_v  = { LOW, LOW, vKnot, vKnot, HIGH, HIGH };

  const double alongU[ 4 ] = { LOW, uKnot, uKnot, HIGH };
  const double alongV[ 4 ] = { LOW, vKnot, vKnot, HIGH };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t i = 0; i < 4; ++i ) {

    for ( uint32_t j = 0; j < 4; ++j ) {

      // Height ONLY on ( left u span, right v span ) - the combination the
      // walk skips when it takes u first.
      const bool omitted = ( i <= 1 ) && ( j >= 2 );

      points.push_back(
        glm::dvec3( alongU[ i ], alongV[ j ], omitted ? HEIGHT : 0.0 ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 4, 4, points );
  surface.weights        = tinynurbs::array2( 4, 4, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  check( std::abs( evaluator.point( vKnot, vKnot ).z - HEIGHT ) < 1.0e-9,
         "the omitted patch really does stand " + std::to_string( HEIGHT ) +
           " up (" + std::to_string( evaluator.point( vKnot, vKnot ).z ) +
           ")" );

  const glm::dvec2 uv0( LOW, LOW );
  const glm::dvec2 uv1( HIGH, HIGH );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  const auto departure = [ & ]( double u, double v ) {

    const double t = ( u - LOW ) / ( HIGH - LOW );

    return glm::length( evaluator.point( u, v ) -
                        ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) );
  };

  double swept = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    const double t = static_cast< double >( i ) / 200000.0;
    const double at = LOW + ( t * ( HIGH - LOW ) );

    swept = std::max( swept, departure( at, at ) );
  }

  check( swept < 1.0e-5,
         "and no sweep of the chord finds it (" +
           std::to_string( swept ) + ")" );

  const double truth = departure( vKnot, vKnot );

  constexpr double TOLERANCE = 1.0e-2;

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, TOLERANCE * TOLERANCE );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "the bound covers the patch whose ordering was ambiguous (" +
           std::to_string( bound ) + " vs " + std::to_string( truth ) + ")" );

  // Covering both orders rather than declining is the point: an ambiguous
  // ordering is common enough on real geometry that refusing it would cost
  // subdivisions, and both boxes together are sound because each parameter
  // in the ambiguous stretch genuinely belongs to one of them.
  check( outcome == CertificateOutcome::Certified,
         "and it is CERTIFIED rather than declined" );
}

/**
 * THE WEIGHT'S PLACEMENT PERTURBATION REACHES THE WEIGHT FLOOR.
 *
 * `placementWeightError` bounds how far the sampled weights sit from the
 * intended ones when a box is narrow enough for the nodes to move in uv. It
 * was propagated into the numerator but not into the weight coefficients -
 * so the reconstructed weight hull omitted it, `weightFloor` could sit above
 * the true minimum denominator, and a denominator that is too large makes
 * the quotient too small. The sixth finding on bldrs-ai/conway-geom#214.
 *
 * The assertion is the invariant rather than a single number: on a rational
 * patch whose weights vary sharply across a narrow span, whatever the
 * certificate returns must still be an upper bound.
 */
void weightPlacementErrorReachesTheFloor() {

  printf( "weightPlacementErrorReachesTheFloor\n" );

  const double first  = 1.0e8;
  const double second =
    std::nextafter( std::nextafter( first,
                                    std::numeric_limits< double >::infinity() ),
                    std::numeric_limits< double >::infinity() );

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 1;
  surface.knots_u  = { first, first, first, second, second, second };
  surface.knots_v  = { 0.0, 0.0, 1.0, 1.0 };

  const double alongU[ 3 ] =
    { first, std::nextafter( first,
                             std::numeric_limits< double >::infinity() ),
      second };

  // Weights that swing hard within the span - which is what makes the weight
  // hull, and therefore the floor, sensitive to where the nodes land.
  const double weight[ 3 ] = { 1.0, 0.01, 1.0 };
  const double height[ 3 ] = { 0.0, 3.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 3; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      points.push_back(
        glm::dvec3( alongU[ a ], static_cast< double >( b ), height[ a ] ) );

      weights.push_back( weight[ a ] );
    }
  }

  surface.control_points = tinynurbs::array2( 3, 2, points );
  surface.weights        = tinynurbs::array2( 3, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  check( !evaluator.isPolynomial(),
         "the patch is rational, so the weight chain runs" );

  const glm::dvec2 uv0( first, 0.5 );
  const glm::dvec2 uv1( second, 0.5 );

  const glm::dvec3 at0 = evaluator.point( uv0.x, uv0.y );
  const glm::dvec3 at1 = evaluator.point( uv1.x, uv1.y );

  double truth = 0.0;

  for ( double u = first; u <= second;
        u = std::nextafter( u, std::numeric_limits< double >::infinity() ) ) {

    const double t = ( u - first ) / ( second - first );

    truth =
      std::max( truth,
                glm::length( evaluator.point( u, 0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  // Small, because the 0.01 middle weight pulls the curve off the control
  // point that carries the height - which is exactly the configuration that
  // makes the weight HULL, and so the floor, sensitive to where the nodes
  // land. The magnitude is not the point; the invariant below is.
  check( truth > 0.02,
         "and it really does depart from its chord (" +
           std::to_string( truth ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, false, 0.0, 0.0, 1.0e-4 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "whatever is certified is still an upper bound (" +
           std::to_string( bound ) + " vs " + std::to_string( truth ) + ")" );
}


/**
 * A PERIODIC CHORD THAT WRAPS IMMEDIATELY STILL WALKS ITS SPANS.
 *
 * `clampedU` records that the walk has run out of knot domain in the
 * direction it is travelling, and it is a LATCH. A strip restart puts the
 * chord back at the far end of the domain with the whole of it still to
 * cross, so the latch has to be cleared there - otherwise the walk stops
 * looking for u boundaries and everything after the wrap comes out as one
 * box carrying whatever span it happened to start in.
 *
 * This is the eighth finding on bldrs-ai/conway-geom#214, found by the
 * differential fuzzer rather than by review, and only once the fuzzer was
 * taught to drive periodic charts - a third of `walkPieces` had never been
 * executed by any test. Measured before the fix: a chord from u = 1 to
 * u = -1 on a strip of exactly that width starts clamped ( it begins on the
 * domain edge, descending ), wraps at once, and then emitted a SINGLE box
 * with spanU = 4 covering u from 1 down to -0.9999999979 - three spans read
 * as one. Bound 0.1618 against a true 0.1809.
 */
void periodicChordThatWrapsStillWalksItsSpans() {

  printf( "periodicChordThatWrapsStillWalksItsSpans\n" );

  constexpr double PERIOD = 2.0;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  // Three spans in u, so a walk that fails to advance is visible.
  surface.knots_u = { -1.0, -1.0, -0.4, -0.4, 0.3, 1.0, 1.0 };
  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  const double across[ 5 ] = { -1.0, -0.4, -0.4, 0.3, 1.0 };
  const double height[ 5 ] = { 0.0, 0.0, 0.9, 0.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 5; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      // CLOSED in u - a periodic chart is only built on a closed surface,
      // so the first and last control rows have to agree.
      const uint32_t row = ( a == 4 ) ? 0 : a;

      points.push_back(
        glm::dvec3( across[ a ], static_cast< double >( b ),
                    height[ row ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 5, 2, points );
  surface.weights        = tinynurbs::array2( 5, 2, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - ( -1.0 ), PERIOD );

    return -1.0 + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  // Starts exactly on the domain edge, descending - so the walk begins
  // clamped and wraps on its first step.
  const glm::dvec2 uv0( 1.0, 0.5 );
  const glm::dvec2 uv1( -1.0, 0.5 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 400000; ++i ) {

    const double t = static_cast< double >( i ) / 400000.0;
    const double u = uv0.x + ( t * ( uv1.x - uv0.x ) );

    truth =
      std::max( truth,
                glm::length( evaluator.point( wrap( u ), 0.5 ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( truth > 0.4,
         "the chord really does depart from its chord after the wrap (" +
           std::to_string( truth ) + ")" );

  const NurbsDeflectionCertificate certificate(
    evaluator, true, -1.0, PERIOD, 1.0e-8 );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and whatever is certified covers it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );

  check( outcome == CertificateOutcome::Certified,
         "the wrapped chord is certified rather than declined" );

  check( certificate.counters().strips > 0,
         "the strip restart really did fire (" +
           std::to_string( certificate.counters().strips ) + ")" );

  check( certificate.counters().clampsClear > 0,
         "and the clamp really was cleared by it (" +
           std::to_string( certificate.counters().clampsClear ) + ")" );
}

}  // namespace

int main() {

  certificateCatchesWhatSamplingMisses();
  fullMultiplicityKnotIsSampledOneSided();
  narrowSpanWithCollapsedParametersIsCertified();
  subUlpKnotSpanIsDeclined();
  collapsedNodeParametersAreDeclined();
  tiedCrossAxisBoundariesCoverBothOrders();
  weightPlacementErrorReachesTheFloor();
  periodicChordThatWrapsStillWalksItsSpans();
  descendingChordStartsInTheSpanItTraverses();
  errorTermTracksTheControlValues();
  knotClippingCatchesAZigZag();
  wrappedChordIsStillCertified();
  boundHoldsOnARationalPatch();
  boundConvergesUnderSubdivision();

  printf( "%s\n", failures == 0 ? "PASS" : "FAIL" );

  return failures == 0 ? 0 : 1;
}
