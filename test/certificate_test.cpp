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
// THE INSTRUMENT TESTS THE CODE, NOT THE CONFIGURATION.
//
// `CERTIFICATE_CERTIFIES_PERIODIC_CHARTS` ships at 0, so a periodic chord in
// a release build takes the sampled test. That is a decision about what to
// expose, and following it here would point the only thing that has ever
// found a periodic defect away from the periodic path - which is the one
// place we already know we cannot size the risk. Everything below therefore
// runs with the periodic path ON.
#define CERTIFICATE_CERTIFIES_PERIODIC_CHARTS 1

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


/**
 * FUZZ SEED 51, minimised: a periodic chord whose u-domain end and strip top
 * are the same point, with a v boundary tied to them.
 *
 * This is the ninth finding on bldrs-ai/conway-geom#214 and it lived in code
 * added for the FIFTH. `endsOnStrip` read `!ambiguous && ...`, so a tie
 * between the two knot axes - an ordering question inside one sheet - also
 * suppressed the strip crossing, which is not an ordering question but a
 * change of frame. Every box after the tie stayed on the sheet the chord had
 * already left.
 *
 * The arrangement is the ordinary one for a periodic chart rather than a
 * contrivance: the knot domain ends where the strip ends, so `tU` and
 * `tStrip` are the same expression and agree bit for bit. All the fuzzer had
 * to add was a v boundary within `separation` of them.
 *
 * The three sub-ulp u spans at the bottom of the domain are what makes the
 * miss visible: the chord's last 1.4e-15 of parameter belongs to them and
 * carries a different control row, so certifying it against span 4 is off by
 * 1.2x rather than by a rounding.
 */
void periodicKnotTieStillCrossesTheStrip() {

  printf( "periodicKnotTieStillCrossesTheStrip\n" );

  constexpr double STRIP_MIN = -1.0;
  constexpr double PERIOD    =  2.0;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 3;

  // Three sub-ulp spans at the bottom of the u domain. They are what the
  // walk has to reach after the wrap, and what it certified span 4 against
  // instead.
  surface.knots_u = { -1.0, -1.0,
                      -0.99999999999999989, -0.99999999999999978,
                      -0.99999999799999972, 1.0, 1.0 };

  surface.knots_v = { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 };

  const double net[ 60 ] = {
    555.30102133302501, 65.140426731815637, -671.69463827467939,
    -387.62508136699967, 693.93097078287371, 148.272119698653,
    969.85539435141084, 414.65515633865573, -253.16965167964986,
    516.7112896770094, 612.21330529454815, -572.52643293363735,
    379.06830333418486, 314.85759129876635, -971.05294396675197,
    124.58739535255211, -110.33469481842383, -330.62818235244981,
    264.39031433017692, -588.93013781199875, 185.92106082125122,
    -692.8471379204841, 637.04875813570641, 922.38010021307844,
    -96.699201316808171, -882.27702333225568, 170.32169744249097,
    666.33689723755754, -674.21753955852512, 337.57164129784007,
    863.99718498663844, 808.71405221529892, -14.887922113306331,
    926.08941688566335, -496.35912436497136, -380.62148294540771,
    115.39132720436066, 56.893751438588083, -120.24050461411821,
    -822.46410300602804, -328.35688412680247, -450.7302651795809,
    507.58412215198848, -256.85685290765218, 780.40904163375274,
    -917.0618030691254, -140.33937138685212, 949.16949352532072,
    // CLOSED in u - the last row repeats the first, as a periodic chart
    // requires.
    555.30102133302501, 65.140426731815637, -671.69463827467939,
    -387.62508136699967, 693.93097078287371, 148.272119698653,
    969.85539435141084, 414.65515633865573, -253.16965167964986,
    516.7112896770094, 612.21330529454815, -572.52643293363735 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 20; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 5, 4, points );
  surface.weights        = tinynurbs::array2( 5, 4, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec2 uv0( -1.153887887186434, 0.54600406465490914 );
  const glm::dvec2 uv1( -0.99999999999999978, -1.0 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  const auto departure = [ & ]( long double t ) {

    const double u =
      static_cast< double >( (long double)uv0.x + ( t * ( (long double)uv1.x - uv0.x ) ) );

    const double v =
      static_cast< double >( (long double)uv0.y + ( t * ( (long double)uv1.y - uv0.y ) ) );

    const glm::dvec3 q = evaluator.point( wrap( u ), v );

    const long double dx = (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
    const long double dy = (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
    const long double dz = (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

    return static_cast< double >( sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
  };

  double sweep = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    sweep = std::max( sweep, departure( (long double)i / 200000.0L ) );
  }

  // THE KNOT PROBE IS THE POINT. No uniform sweep lands in a sub-ulp span,
  // so the departure that the miss hides is only visible by evaluating at
  // the knot itself - which is where the b-spline call site's own subdivision
  // would land too.
  const long double tKnot =
    ( (long double)surface.knots_u[ 2 ] - uv0.x ) /
    ( (long double)uv1.x - uv0.x );

  const double atKnot = departure( tKnot );

  const double truth = std::max( sweep, atKnot );

  check( tKnot > 0.0L && tKnot < 1.0L,
         "the tied knot is interior to the chord (" +
           std::to_string( (double)tKnot ) + ")" );

  check( atKnot > 1700.0,
         "the sub-ulp span really does carry a large departure (" +
           std::to_string( atKnot ) + ")" );

  check( atKnot > sweep * 1.2,
         "and one a uniform sweep does not see (" + std::to_string( atKnot ) +
           " vs " + std::to_string( sweep ) + ")" );

  const double tolerance = 5.4642684472996372e-07;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD, tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  // The tie has to actually be a tie, or this case stops testing the
  // interaction it was minimised for.
  check( certificate.counters().ambiguous > 0,
         "the two knot axes really do tie here (" +
           std::to_string( certificate.counters().ambiguous ) + ")" );

  check( certificate.counters().strips > 0,
         "the strip crossing fires in spite of the tie (" +
           std::to_string( certificate.counters().strips ) + ")" );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below the truth (" +
           std::to_string( bound ) + " vs " + std::to_string( truth ) + ")" );
}


/**
 * FUZZ SEED 1 OF THE FAR-SHEET GENERATOR, minimised: a periodic chord a
 * million and a half sheets outside its strip.
 *
 * This is the case the `| shift |` term in `roundingU` was written for, and
 * until this test it had no red proof - not because it was inert, but because
 * the fuzzer could not reach the regime. The stock generator draws u from
 * [ lowU * 1.5, highU * 1.5 ], so the walk's shift never exceeded ONE period
 * and an ulp of the raw parameter was within a factor of three of an ulp of
 * the wrapped one. Here the raw parameter is 1.58e6 and the wrapped one is
 * O(1), so the two differ by a factor of about 7e6.
 *
 * The callback evaluates `point( wrapChartU( u ), v )` - one exact `fmod` -
 * and the walk evaluates `point( u - shift, v )` - three rounded operations.
 * The gap between them is an ulp of the RAW parameter, and a bound that
 * prices node placement at an ulp of the WRAPPED parameter is under by that
 * same factor: measured 6.437016711e-13 against a true departure of
 * 4.148887769e-11.
 */
void chordManySheetsOutsideTheStripIsPricedByItsRawParameter() {

  printf( "chordManySheetsOutsideTheStripIsPricedByItsRawParameter\n" );

  constexpr double STRIP_MIN = -1.0;
  constexpr double PERIOD    =  2.0;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  surface.knots_u = { -1.0, -1.0, -0.40266291262089327, 1.0, 1.0 };

  surface.knots_v = { -1.0, -1.0,
                      -0.61210848541668317, -0.61210848341668311,
                      -0.612108483416683, 1.0, 1.0 };

  const double net[ 45 ] = {
    0.052911171003463126, 0.16590149205574273, -0.045694485954048736,
    0.036517219798542222, 0.13548374136138874, 0.13569607054083888,
    -0.049908481206195704, -0.092117586814640315, -0.077188447506570293,
    0.04456451635024361, -0.12116806075050654, -0.094116024429135442,
    -0.014455023125218208, -0.14657714361014568, -0.045219622990830655,
    -0.069632949116618117, 0.067747662297891864, -0.11804404130381119,
    0.12357586605036991, 0.12773533670499271, -0.15194700979759118,
    0.037444559395139743, -0.097194160558549347, -0.10370575225061768,
    0.099712115874956486, -0.13928271018335184, -0.14529671712880021,
    0.089543831029257503, -0.087115120885457129, 0.11874110783743064,
    // CLOSED in u.
    0.052911171003463126, 0.16590149205574273, -0.045694485954048736,
    0.036517219798542222, 0.13548374136138874, 0.13569607054083888,
    -0.049908481206195704, -0.092117586814640315, -0.077188447506570293,
    0.04456451635024361, -0.12116806075050654, -0.094116024429135442,
    -0.014455023125218208, -0.14657714361014568, -0.045219622990830655 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 15; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 3, 5, points );
  surface.weights        = tinynurbs::array2( 3, 5, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const double tolerance = 1.7830893502681724e-06;

  // The same chord twice: once where the fuzzer found it, and once slid back
  // into the strip. The GEOMETRY is identical - only the raw magnitude of the
  // parameter differs - so anything that separates the two answers is the
  // distance and nothing else.
  const double SHEETS = 1583448.0;

  const auto run = [ & ]( double offset, double& bound, double& truth ) {

    const glm::dvec2 uv0( 1583447.5973370874 - offset, 1.0 );
    const glm::dvec2 uv1( 1583447.0 - offset, 1.0 );

    const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
    const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

    const auto departure = [ & ]( long double t ) {

      const double u =
        static_cast< double >( (long double)uv0.x +
                               ( t * ( (long double)uv1.x - uv0.x ) ) );

      const glm::dvec3 q = evaluator.point( wrap( u ), uv0.y );

      const long double dx =
        (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
      const long double dy =
        (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
      const long double dz =
        (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

      return static_cast< double >(
        sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
    };

    truth = 0.0;

    for ( uint32_t i = 0; i <= 200000; ++i ) {

      truth = std::max( truth, departure( (long double)i / 200000.0L ) );
    }

    const NurbsDeflectionCertificate certificate(
      evaluator, true, STRIP_MIN, PERIOD, tolerance * tolerance );

    bound = 0.0;

    const CertificateOutcome outcome =
      certificate.bound( uv0, uv1, at0, at1, bound );

    return outcome;
  };

  double farBound  = 0.0;
  double farTruth  = 0.0;
  double nearBound = 0.0;
  double nearTruth = 0.0;

  const CertificateOutcome farOutcome  = run( 0.0, farBound, farTruth );
  const CertificateOutcome nearOutcome = run( SHEETS, nearBound, nearTruth );

  // `std::to_string` is %f, and everything here is at 1e-11. Print in a
  // form that still says something when the assertion goes red.
  const auto show = []( double value ) {

    char text[ 32 ];

    snprintf( text, sizeof( text ), "%.10g", value );

    return std::string( text );
  };

  check( farTruth > 4.0e-11,
         "the far chord really does depart from its chord (" +
           show( farTruth ) + ")" );

  check( nearOutcome == CertificateOutcome::Certified,
         "the same geometry inside the strip is certified, so nothing here "
         "is about the surface" );

  check( farOutcome != CertificateOutcome::Certified || farBound >= farTruth,
         "and nothing is certified below the truth a million sheets out (" +
           show( farBound ) + " vs " + show( farTruth ) + ")" );
}


/**
 * FUZZ SEED 6, minimised: a descending periodic chord whose u-domain edge and
 * strip edge are 8e-9 apart on a chart whose knots are at 2.4e7 - so their
 * two cut parameters DIVIDE TO THE SAME DOUBLE.
 *
 * The tenth finding on bldrs-ai/conway-geom#214, and the last comparison in
 * the walk still made in parameter space. `tU` and `tStrip` are both
 * `( edge - uv0.x ) / du`; at this magnitude an ulp of the numerator is
 * 7.5e-9, so two edges 8e-9 apart round together and `tStrip <= tTo` decides
 * the order by which way the `<=` happens to fall. It fell to the strip, the
 * walk jumped a sheet, and the three u spans still between it and the strip
 * edge - which is where the chord's last 8e-9 of parameter actually lies -
 * were never entered. Bound 0.04044181424 against a true departure of
 * 0.05351126679.
 *
 * The fix compares the two edges in EVALUATED u, where both are knot-or-strip
 * values and the comparison is exact. This test pins the outcome, not the
 * route: the walk may certify or decline, but it may not certify short.
 */
void tiedStripAndKnotEdgesAreOrderedInParameterSpaceNotChordSpace() {

  printf( "tiedStripAndKnotEdgesAreOrderedInParameterSpaceNotChordSpace\n" );

  constexpr double STRIP_MIN = -24461668.268462215;
  constexpr double PERIOD    =  48923336.536924429;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 3;
  surface.degree_v = 1;

  // Three micro-spans at the bottom of the u domain, 4e-9 apart, on knots of
  // magnitude 2.4e7 - a ratio of 1.6e-16, just under one ulp.
  surface.knots_u = { -24461668.268462215, -24461668.268462215,
                      -24461668.268462215, -24461668.268462215,
                      -24461668.268462211, -24461668.268462211,
                      -24461668.268462207,
                      24461668.268462215, 24461668.268462215,
                      24461668.268462215, 24461668.268462215 };

  surface.knots_v = { -24461668.268462215, -24461668.268462215,
                      -11751770.912085051,
                      24461668.268462215, 24461668.268462215 };

  const double net[ 63 ] = {
    0.0066021523315360764, 0.0087403881186380247, -0.017838196908007746,
    -0.00036635123264366909, -0.037925092929460218, 0.034415441423718503,
    -0.046240897473654216, -0.015784834020217187, -0.039309082375280364,
    0.017570293017203982, -0.011130315439521088, -0.027377173131882249,
    -0.018902762490813903, -0.018974040382710486, -0.010325786262893509,
    0.038114051640173625, 0.014375747837490439, 0.021334372231788613,
    -0.0052054615923575395, -0.022852232836698082, -0.036244851294394462,
    -0.0017487195966010075, 0.025877725551360135, 0.010921565304606598,
    0.026140166429110261, -0.036557824616160549, -0.0080551014071843299,
    0.039606786982324983, 0.025214589202848449, 0.015162067819354232,
    -0.034982981155498999, -0.047897659130310792, 0.03628890198449803,
    -0.016782497112719605, 0.011908324705000958, -0.03839697114252083,
    0.016485717532553075, 0.0078283831644718263, 0.030964196182101521,
    0.020850631199250146, 0.019178761521681055, 0.013136798616221912,
    -0.027548854220578334, -0.048154409062731718, -0.014318669272199297,
    -0.013689020277261689, -0.038135412920725822, 0.022894971297073068,
    0.041231321101795747, 0.043183506054162484, -0.037909026968017318,
    0.010092325060880133, 0.024474695099902764, 0.013528708125882444,
    // CLOSED in u.
    0.0066021523315360764, 0.0087403881186380247, -0.017838196908007746,
    -0.00036635123264366909, -0.037925092929460218, 0.034415441423718503,
    -0.046240897473654216, -0.015784834020217187, -0.039309082375280364 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 21; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 7, 3, points );
  surface.weights        = tinynurbs::array2( 7, 3, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec2 uv0( 10211859.265748158, -8091389.4596511573 );
  const glm::dvec2 uv1( -24461668.268462215, 24461668.268462215 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  const auto departure = [ & ]( long double t ) {

    const double u =
      static_cast< double >( (long double)uv0.x + ( t * ( (long double)uv1.x - uv0.x ) ) );

    const double v =
      static_cast< double >( (long double)uv0.y + ( t * ( (long double)uv1.y - uv0.y ) ) );

    const glm::dvec3 q = evaluator.point( wrap( u ), v );

    const long double dx =
      (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
    const long double dy =
      (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
    const long double dz =
      (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

    return static_cast< double >( sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
  };

  double sweep = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    sweep = std::max( sweep, departure( (long double)i / 200000.0L ) );
  }

  // The departure that the skipped spans carry sits two ulps of t below the
  // chord's end, where no uniform sweep lands.
  long double tKnot =
    ( (long double)surface.knots_u[ 4 ] - uv0.x ) /
    ( (long double)uv1.x - uv0.x );

  tKnot = std::nextafterl( std::nextafterl( tKnot, 0.0L ), 0.0L );

  const double atKnot = departure( tKnot );

  const double truth = std::max( sweep, atKnot );

  check( atKnot > 0.05,
         "the skipped micro-spans really do carry a departure (" +
           std::to_string( atKnot ) + ")" );

  check( atKnot > sweep * 1.2,
         "and one a uniform sweep does not see (" + std::to_string( atKnot ) +
           " vs " + std::to_string( sweep ) + ")" );

  const double tolerance = 4.6917256171007448e-07;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD, tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );
}


/**
 * FUZZ SEED 3 AT TRIAL 32,525, minimised: the walk's periodic reduction and
 * the callback's landing on OPPOSITE ENDS OF THE STRIP.
 *
 * The reviewer on bldrs-ai/conway-geom#214 named this mechanism - "use the
 * callback's reduction for the initial periodic point" - while attributing
 * it to a case where it was not the cause (on fuzz seed 51 the two reductions
 * agree bit for bit; what went wrong there was the strip crossing being
 * suppressed by a knot tie). It is a real defect all the same, and this is
 * the case that shows it.
 *
 * `wrapChartU` is `uMin + fmod( u - uMin, P )`: IEEE `fmod` is exact, so the
 * callback's reduction has no error at all. The walk used
 * `u - floor( ( u - uMin ) / P ) * P`: four rounded operations, and near a
 * sheet boundary the rounded quotient crosses an integer where the exact one
 * does not. Here u = 79627595.655702367 with P = 114.74700735250697: `floor`
 * puts the start 6.2e-9 above the BOTTOM of the strip, `fmod` puts it 1.3e-9
 * below the TOP - a full period apart.
 *
 * `du` is zero, so the walk never moves in u and never gets a second chance:
 * every box is bounded at a u the callback does not evaluate.
 */
void periodicStartIsReducedTheWayTheCallbackReducesIt() {

  printf( "periodicStartIsReducedTheWayTheCallbackReducesIt\n" );

  constexpr double STRIP_MIN = -57.373503676253485;
  constexpr double PERIOD    = 114.74700735250697;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 3;

  surface.knots_u = { -57.373503676253485, -57.373503676253485,
                      -57.373503676253485,
                      -57.373503676253478, -57.373503676253478,
                      57.373503676253485, 57.373503676253485,
                      57.373503676253485 };

  surface.knots_v = { -57.373503676253485, -57.373503676253485,
                      -57.373503676253485, -57.373503676253485,
                      -57.373503561506475,
                      57.373503676253485, 57.373503676253485,
                      57.373503676253485, 57.373503676253485 };

  // THE FUZZER'S OWN NET, not a substitute. A first attempt at this test
  // used a tidy closed net of its own and it was VACUOUS - the case declined
  // either way, so the assertion passed with the defect restored. The
  // numbers below are what trial 32,525 generated, and the last row repeats
  // the first, which is what makes the chart periodic.
  const double net[ 75 ] = {
    -0.0017517532670473183, 0.0015746121395702139, 0.0015821301951952165,
    0.0011702115466867172, -0.00053782024544115161, -0.0014410850813068557,
    -0.00096316625058488455, -0.0022013435825237121, 0.00059372150553900667,
    -0.0015870419184512066, 0.001234968195446543, -0.001454928435764211,
    0.0010687300938512258, 0.00010077528341219952, 0.0014430209278028286,
    0.0019772823348125474, 0.00070006936510851395, 0.0020098337655775548,
    0.00071094844842253552, -0.0010006343639879726, 0.0021567154530632046,
    -0.0014721930330668129, -0.002108944335109299, 0.0018325025641636188,
    -0.00090043907719711291, 0.0019139541911185752, 0.0016533185517177778,
    -0.0023864771573792005, 0.0020152529889563707, 0.002076481882942544,
    -0.00085896604098827925, 0.0018886305100227567, -0.0012434574668067522,
    -0.00018737484334745711, -0.00064339061424294511, -0.0012541330485418374,
    -0.0022992604075862096, -0.00033161013839504106, 0.0020483258091436294,
    0.00079329795713613, 0.0020689672792076597, 0.0010065252781206817,
    0.0022849961782679408, -0.0019619367574640304, 0.0015740266774232257,
    -0.00085014514229452866, -0.0014687989490391508, -0.00054461444967614193,
    -0.00023425483018200792, 0.001847974228920695, -0.00068300935272389246,
    0.0014788345878406182, 0.002266438613054565, 0.0023472733559065655,
    0.0013740424836580166, 0.00054786877420530183, -0.0006150195648733012,
    -0.0021484157315325232, 0.00089823235151480912, 0.0018398522840405679,
    -0.0017517532670473183, 0.0015746121395702139, 0.0015821301951952165,
    0.0011702115466867172, -0.00053782024544115161, -0.0014410850813068557,
    -0.00096316625058488455, -0.0022013435825237121, 0.00059372150553900667,
    -0.0015870419184512066, 0.001234968195446543, -0.001454928435764211,
    0.0010687300938512258, 0.00010077528341219952, 0.0014430209278028286 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 25; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 5, 5, points );
  surface.weights        = tinynurbs::array2( 5, 5, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  // The two reductions have to actually disagree, or this case stops being
  // the case it was minimised for.
  const double raw = 79627595.655702367;

  const double byFloor =
    raw - ( std::floor( ( raw - STRIP_MIN ) / PERIOD ) * PERIOD );

  check( byFloor != wrap( raw ),
         "the floor reduction and the callback's fmod really do disagree" );

  check( std::abs( byFloor - wrap( raw ) ) > 0.5 * PERIOD,
         "and by most of a period, not by a rounding (" +
           std::to_string( std::abs( byFloor - wrap( raw ) ) ) + " of " +
           std::to_string( PERIOD ) + ")" );

  // A chord that moves only in v, so the walk never steps in u and the start
  // point is the only u it ever uses.
  const glm::dvec2 uv0( raw, -57.373503676253485 );
  const glm::dvec2 uv1( raw,  57.373503676253485 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    const double t = static_cast< double >( i ) / 200000.0;
    const double v = uv0.y + ( t * ( uv1.y - uv0.y ) );

    truth =
      std::max( truth,
                glm::length( evaluator.point( wrap( uv0.x ), v ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  check( truth > 0.0,
         "the chord really does depart from its chord (" +
           std::to_string( truth ) + ")" );

  // The v knot the chord crosses is where the departure peaks, and it is two
  // ulps of t from the chord's start - no uniform sweep lands there.
  {
    const long double tKnot =
      ( (long double)surface.knots_v[ 4 ] - uv0.y ) /
      ( (long double)uv1.y - uv0.y );

    const double v =
      static_cast< double >(
        (long double)uv0.y +
        std::nextafterl( std::nextafterl( tKnot, 0.0L ), 0.0L ) *
          ( (long double)uv1.y - uv0.y ) );

    const double t = ( v - uv0.y ) / ( uv1.y - uv0.y );

    truth =
      std::max( truth,
                glm::length( evaluator.point( wrap( uv0.x ), v ) -
                             ( ( at0 * ( 1.0 - t ) ) + ( at1 * t ) ) ) );
  }

  const double tolerance = 4.5874072345266852e-06;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD, tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome == CertificateOutcome::Certified,
         "the case is certified rather than declined, so the assertion below "
         "is about the bound and not about a refusal" );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );
}


/**
 * FUZZ SEED 7 AT TRIAL 173,559, minimised: a chord whose v moves 2e-9 in
 * total and whose far end sits exactly on a v knot of multiplicity
 * degree + 1.
 *
 * THIS IS THE CASE THAT RED-PROVES A GUARD DELETED EARLIER IN THIS PR. The
 * second arm of boundPiece's unplaceable test - "a box whose own extent on an
 * axis is below the parameter's rounding AND that sits on a span boundary the
 * surface steps across" - was removed on the reasoning that it declined
 * 46,000 of 250,000 fuzz cases while catching nothing. Nothing could break it
 * and see a test go red, so the only discriminator left was cost. That
 * reasoning was wrong in the way it is usually wrong: the generator had not
 * reached the shape it defends.
 *
 * The shape: `findSpan` resolves a parameter sitting exactly ON a knot to the
 * span on its RIGHT. The walk assigns spans by exact knot arithmetic and
 * evaluates one-sided, which is what closed the second finding - but THE CALL
 * SITE does not. It calls `point( u, v )`, and `point` uses `findSpan` on the
 * rounded double. Normally that matters at a single parameter, which is one
 * point of a continuous piece. Here it matters over a stretch: v moves 2e-9
 * over the whole chord, so every t above 1 - 2.8e-8 has a rounded v that IS
 * the knot. Twenty-eight million ulps of t are evaluated on the span to the
 * knot's right while every box says the span to its left, and the two
 * polynomials differ by 248.57 at the parameter that matters.
 *
 * Certified at 335.968276 against a true departure of 408.4816988 - an
 * unsound answer, not merely a loose one.
 */
void steppingAxisBoxBelowResolutionIsDeclined() {

  printf( "steppingAxisBoxBelowResolutionIsDeclined\n" );

  constexpr double STRIP_MIN = -1.0;
  constexpr double PERIOD    =  2.0;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  surface.knots_u = { -1.0, -1.0,
                      -0.99999999799999995,
                      -0.99999999599999989, -0.99999999599999989,
                      1.0, 1.0 };

  // MULTIPLICITY 2 AT DEGREE 1 - the surface steps here, and the chord's far
  // end sits exactly on it.
  surface.knots_v = { -1.0, -1.0,
                      -0.81127201727053222,
                      -0.81127201527053217, -0.81127201527053217,
                      1.0, 1.0 };

  const double net[ 75 ] = {
    -169.9957553632207, -41.132363392136909, -145.87756369007079,
    139.80372154338238, 12.506191339502578, -50.707897174967272,
    -175.05345266721827, 27.168742267987881, 86.624815630241358,
    -22.117766762763551, -176.07785699109371, -132.86391356717436,
    -1.5696927129000642, -165.53884904146835, -8.9130757260750784,
    18.699445888633903, 144.64542753602871, 41.530717206185869,
    71.574650064996547, -125.54396950023363, 90.867903285862667,
    -15.908345442602638, -66.867651452344816, -13.863102913457084,
    38.89792383958121, 103.88908784372774, 158.25960720542184,
    -168.75481017041835, 115.72812312156969, -146.25775659258218,
    -117.13578217143505, 126.60730397272948, -177.30463605391586,
    167.04055469283242, 25.65881483694946, -143.25310443513311,
    -69.339173731589341, 17.079050772636293, -97.620718302579107,
    -101.33055154393435, -156.89376430033667, 110.26488451361746,
    -117.02486187873714, 138.61061564601948, -166.86941065159783,
    95.187055101484077, 174.33909506806532, -171.42430851085646,
    37.162837216709534, 56.459273897988879, 115.57593631492381,
    157.63172340659526, 42.339565775580837, -2.0623714671639277,
    -2.0810360459615538, -4.352442881652081, 30.797692461413845,
    141.57775191022643, 91.376348188891029, -142.44427486761816,
    -169.9957553632207, -41.132363392136909, -145.87756369007079,
    139.80372154338238, 12.506191339502578, -50.707897174967272,
    -175.05345266721827, 27.168742267987881, 86.624815630241358,
    -22.117766762763551, -176.07785699109371, -132.86391356717436,
    -1.5696927129000642, -165.53884904146835, -8.9130757260750784 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 25; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 5, 5, points );
  surface.weights        = tinynurbs::array2( 5, 5, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec2 uv0( 1.0, -0.81127201727053222 );
  const glm::dvec2 uv1( -1.0, -0.81127201527053217 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  // The rounded v really does sit ON the knot over a macroscopic stretch of
  // the chord, which is why a single-point argument does not cover this.
  double firstOnKnot = 1.0;

  {
    double low  = 0.0;
    double high = 1.0;

    for ( uint32_t i = 0; i < 100; ++i ) {

      const double mid = 0.5 * ( low + high );
      const double v   = uv0.y + ( mid * ( uv1.y - uv0.y ) );

      if ( v >= surface.knots_v[ 3 ] ) { high = mid; } else { low = mid; }
    }

    firstOnKnot = high;
  }

  char stretch[ 32 ];

  snprintf( stretch, sizeof( stretch ), "%.4g", 1.0 - firstOnKnot );

  check( 1.0 - firstOnKnot > 1.0e-9,
         "the rounded v sits on the stepping knot over a macroscopic stretch "
         "of the chord (" + std::string( stretch ) + " of it)" );

  check( RationalSurfaceEvaluator::findSpan(
           surface.degree_v, surface.knots_v, surface.knots_v[ 3 ] ) != 2,
         "and findSpan puts that v on the span to the knot's RIGHT, not the "
         "one the walk walks" );

  const auto departure = [ & ]( long double t ) {

    const double u =
      static_cast< double >( (long double)uv0.x + ( t * ( (long double)uv1.x - uv0.x ) ) );

    const double v =
      static_cast< double >( (long double)uv0.y + ( t * ( (long double)uv1.y - uv0.y ) ) );

    const glm::dvec3 q = evaluator.point( wrap( u ), v );

    const long double dx =
      (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
    const long double dy =
      (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
    const long double dz =
      (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

    return static_cast< double >( sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
  };

  double truth = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    truth = std::max( truth, departure( (long double)i / 200000.0L ) );
  }

  {
    long double tKnot =
      ( (long double)surface.knots_u[ 2 ] - uv0.x ) /
      ( (long double)uv1.x - uv0.x );

    tKnot = std::nextafterl( std::nextafterl( tKnot, 0.0L ), 0.0L );

    truth = std::max( truth, departure( tKnot ) );
  }

  check( truth > 400.0,
         "the mismatched polynomial really does carry the departure (" +
           std::to_string( truth ) + ")" );

  const double tolerance = 1.114693126170528;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD, tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );
}


/**
 * GRAZING GENERATOR, TRIAL 272: a box that sits a few ulps short of the knot
 * the surface steps at, and so passes an EQUALITY test for "touches the
 * boundary" while the caller evaluates the whole of it on the far side.
 *
 * The thirteenth finding on bldrs-ai/conway-geom#214, and the first one this
 * file's generator was BUILT to produce rather than waited for - the shape
 * is the eleventh finding's, and the free generator reached it once in
 * 2,800,000 trials.
 *
 * The box's u corner is 41.412005379020393; the knot is 41.4120053790204.
 * They differ by 7.1e-15, which is the periodic shift the walk removed, and
 * `roundingU` - the quantity the rest of boundPiece uses to price exactly
 * that displacement - is 7.4e-14, ten times larger. So the corner is inside
 * the rounding of the boundary and the equality test said it was not.
 * Meanwhile the caller's own parameter over that whole box is
 * 41.41200537902045653, above the knot, on the other polynomial.
 *
 * Bound 1.639305956 against a true departure of 2.013100244.
 */
void grazingBoxShortOfAStepKnotIsDeclined() {

  printf( "grazingBoxShortOfAStepKnotIsDeclined\n" );

  constexpr double STRIP_MIN = -649.28405088089676;
  constexpr double PERIOD    = 1298.5681017617935;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 2;

  surface.knots_u = {
    -649.28405088089676, -649.28405088089676, -649.28405088089676,
    41.4120053790204, 41.4120053790204, 41.4120053790204,
    649.28405088089676, 649.28405088089676, 649.28405088089676 };

  surface.knots_v = {
    -649.28405088089676, -649.28405088089676, -649.28405088089676,
    -216.42801696029892, -216.42801696029892, -216.42801696029892,
    216.42801696029892, 216.42801696029892, 649.27781367627631,
    649.27781367627631, 649.28405088089676, 649.28405088089676,
    649.28405088089676 };

  // THE FUZZER'S OWN NET. A synthetic stand-in was tried first for an
  // earlier case in this family and came out VACUOUS - it declined with the
  // defect restored, so the assertion proved nothing. These are the numbers
  // the generator produced.
  const double net[ 180 ] = {
    0.6955265944896869, 0.4195818939617868, 0.23649686284885174,
    0.43506442593133676, 0.35349079661584981, 0.24507864971235727,
    1.482435777576337, 0.51140376283918543, -0.093813248560431717,
    3.0715151123892643, 0.070186737338705063, -0.09494420126431663,
    4.126487336092814, -0.71348450135283503, -0.14499906947274765,
    4.1091813922114184, -0.067044286904324624, 0.10208992476339468,
    6.2420222988543728, -0.97879065202761639, -0.011724978349541948,
    6.4859114917319562, -0.32132108199680953, -0.047449718841585942,
    8.850372349347154, -0.5387577305095832, -0.14864147209030787,
    9.9638367204209484, -0.55681999747149624, 0.23324048210634496,
    -0.23989928207377065, 0.48456392004468363, -0.10850431639749408,
    1.3589609814097299, 0.15442991925399663, -0.031262400324812978,
    1.2107131495574148, 0.24677633536057808, 0.13043359521834985,
    2.988537491183286, -0.90512159570873529, -0.027038313803487757,
    4.1761591743625344, -0.81640372389125315, -0.11160497950537446,
    4.387754541351585, 0.6768346595081316, 0.20188321897383504,
    6.602796901892531, -0.13290421978248923, -0.18711567188316958,
    6.2600465976595574, 0.12187437927312117, 0.0048292243015378156,
    8.524106602124526, 0.92269790688049835, 0.1571596557312907,
    9.3630074053930219, -0.83017773915200155, -0.1651702911476719,
    -0.24685021626395542, -0.67010256586881467, 0.13050054303307057,
    1.5653427678358236, -0.94367890023758272, -0.13579327411574066,
    1.5109457067530849, -0.83598041172761128, -0.19622596015997851,
    3.0710489211629115, 0.38177229059754181, 0.18087397648697151,
    3.9028292038003229, 0.092179503479963953, -0.23447521456791881,
    5.5390267104448947, 0.98129259939590718, 0.065282022777113147,
    6.8518726760440956, 0.54572820006592049, -0.038350037180193863,
    6.7661701531672565, -0.33196561876674369, -0.20806323365521595,
    7.4926404710672019, 0.38286235981459704, -0.072417462803534616,
    8.6844168050590635, 0.97060894862993297, -0.0124975491737635,
    -0.56869291957948498, 0.68657017261064146, 0.017606502390226986,
    1.6513828580308323, 0.77061243809262736, 0.093264912038198922,
    1.5806252723709937, -0.99787560993785274, -0.010074038390889384,
    2.2375421029582769, -0.30436525819028826, 0.022082181214219121,
    4.3285786233286592, 0.83950710652961003, -0.084678663139051591,
    4.7987229762055907, 0.31333012286022055, 0.033536574205692116,
    5.4786653152798568, 0.14824188345882283, 0.23505708283341914,
    6.5009154189360849, -0.9008563743945468, 0.22937890894804491,
    8.4014953110408932, 0.91437835082919094, -0.14643755537291914,
    8.3941140927363396, -0.26604684867930062, -0.026645643756667392,
    -0.73820238363975021, -0.88750723274595655, 14.639762184519677,
    1.154800489251494, 0.8825031822839664, 14.797947804105135,
    2.5488791991337116, -0.86969832943806336, 14.888347127871308,
    3.9663849170579715, 0.7264588325171708, 14.72453115190781,
    3.8631494803974702, -0.78482172588137944, 14.862338645927965,
    4.5014423534938608, -0.26878504074880905, 15.080944681142245,
    6.2722896760119626, -0.95422213355722385, 14.63522670700293,
    7.8265476506490073, -0.5102222353164545, 14.651967255366719,
    7.2376431288224286, -0.47120671976867357, 14.896362803206376,
    9.9907446073782591, -0.15269085086784862, 15.007982684967594,
    0.6955265944896869, 0.4195818939617868, 0.23649686284885174,
    0.43506442593133676, 0.35349079661584981, 0.24507864971235727,
    1.482435777576337, 0.51140376283918543, -0.093813248560431717,
    3.0715151123892643, 0.070186737338705063, -0.09494420126431663,
    4.126487336092814, -0.71348450135283503, -0.14499906947274765,
    4.1091813922114184, -0.067044286904324624, 0.10208992476339468,
    6.2420222988543728, -0.97879065202761639, -0.011724978349541948,
    6.4859114917319562, -0.32132108199680953, -0.047449718841585942,
    8.850372349347154, -0.5387577305095832, -0.14864147209030787,
    9.9638367204209484, -0.55681999747149624, 0.23324048210634496 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 60; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 6, 10, points );
  surface.weights        = tinynurbs::array2( 6, 10, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec2 uv0( 41.412005379012506,-649.28405088089676 );
  const glm::dvec2 uv1( 41.4120053790204,649.28405088089676 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  const auto departure = [ & ]( long double t ) {

    const double u =
      static_cast< double >( (long double)uv0.x + ( t * ( (long double)uv1.x - uv0.x ) ) );

    const double v =
      static_cast< double >( (long double)uv0.y + ( t * ( (long double)uv1.y - uv0.y ) ) );

    const glm::dvec3 q = evaluator.point( wrap( u ), v );

    const long double dx =
      (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
    const long double dy =
      (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
    const long double dz =
      (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

    return static_cast< double >( sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
  };

  double sweep = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    sweep = std::max( sweep, departure( (long double)i / 200000.0L ) );
  }

  // The departure lives at the v knots the chord crosses, two ulps of t
  // below each - which is inside the grazing window and nowhere a uniform
  // sweep lands.
  double truth = sweep;

  for ( double knot : surface.knots_v ) {

    long double t =
      ( (long double)knot - uv0.y ) / ( (long double)uv1.y - uv0.y );

    if ( !( t > 0.0L ) || !( t < 1.0L ) ) continue;

    for ( int nudge = -2; nudge <= 2; ++nudge ) {

      long double tt = t;

      for ( int k = 0; k < std::abs( nudge ); ++k ) {
        tt = std::nextafterl( tt, nudge < 0 ? 0.0L : 1.0L );
      }

      truth = std::max( truth, departure( tt ) );
    }
  }

  check( truth > sweep,
         "the grazing window carries a departure a uniform sweep does not "
         "see (" + std::to_string( truth ) + " vs " +
           std::to_string( sweep ) + ")" );

  const double tolerance = 0.20308719174249568;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD,
    tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );
}



/**
 * GRAZING GENERATOR, TRIAL 2,967: the same shape again, but with the
 * displacement LARGER than `roundingU` rather than smaller.
 *
 * The fourteenth finding on bldrs-ai/conway-geom#214. `wrapChartU( x )` is
 * `uMin + fmod( x - uMin, P )`, and while `fmod` is exact the subtraction
 * and the addition each round AT THE MAGNITUDE OF `uMin`, not at the
 * magnitude of x. The walk's own interior points are `x - shift` with
 * `shift` fixed from the chord's start, so the two reductions drift apart by
 * a few ulps of the STRIP as the chord runs.
 *
 * Here the strip origin is -27637.36 while the parameter is -399.37, so ulps
 * of the strip are seventy times ulps of the parameter. The box's corner is
 * 1.6e-12 from the step knot and `roundingU` was 7.1e-13 - the box was
 * judged not to touch the boundary, and every parameter the caller evaluates
 * inside it lands on the far side.
 *
 * Bound 2.181686448 against a true departure of 2.563162239. This case
 * red-proves BOTH the rounding-aware boundary test and the strip term in
 * `roundingU`; trial 272 red-proves only the first.
 */
void reductionDriftAcrossAStepKnotIsPriced() {

  printf( "reductionDriftAcrossAStepKnotIsPriced\n" );

  constexpr double STRIP_MIN = -27637.360607875835;
  constexpr double PERIOD    = 55274.72121575167;

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 2;
  surface.degree_v = 1;

  surface.knots_u = {
    -27637.360607875835, -27637.360607875835, -27637.360607875835,
    -399.37070255423532, -399.37070255423532, -399.37070255423532,
    27637.360607875835, 27637.360607875835, 27637.360607875835 };

  surface.knots_v = {
    -27637.360607875835, -27637.360607875835, -9212.4535359586116,
    9212.4535359586116, 27570.420060896493, 27637.360607875835,
    27637.360607875835 };

  // THE FUZZER'S OWN NET. A synthetic stand-in was tried first for an
  // earlier case in this family and came out VACUOUS - it declined with the
  // defect restored, so the assertion proved nothing. These are the numbers
  // the generator produced.
  const double net[ 90 ] = {
    -0.10637260889229627, -0.80506420120276667, -0.24765280605537204,
    0.98463971838264741, 0.99172047752049197, 0.11647206058754805,
    2.1835887004361521, -0.63327470946302045, -0.19835542808978235,
    2.6691221737777013, -0.4410476128906391, 0.071430659019328879,
    3.0513059480333982, 0.88726749799577753, 0.064484775276447692,
    0.78511784534872886, 0.92382349200403269, 0.057331797302886678,
    1.2352289347475816, 0.9906114940427373, 0.24054541991464978,
    2.8724926199273368, 0.30218567376989047, 0.20546407757403229,
    3.5698529452138383, 0.87099645825953487, 0.20642441299446562,
    4.3041424943265003, 0.043582547934907012, -0.14130170125062019,
    0.83982062358012177, 0.75066005931565494, -0.02107276667048813,
    0.66376246072443457, -0.88486588603542171, 0.13651035432610525,
    2.2338456749256581, 0.61660382940451597, -0.21321876778327886,
    3.9526730214869357, -0.15225280378600226, 0.080363187549146842,
    3.4717995597296558, -0.42930476694252828, -0.24353180854833673,
    0.51427662144011999, 0.42957847705122632, 0.12455023851094976,
    1.5581623964857281, 0.17412151270780818, -0.21874665061059273,
    2.0502443051206969, 0.56193305292029105, 0.11254665766060451,
    2.0873407819724754, -0.24344124118713495, 0.21221126001482599,
    4.3312049067885603, 0.99544719650764901, 0.071692380666408373,
    -0.33931832901448034, 0.85893528658806173, 271.71326891478748,
    0.93050701563107241, 0.94466649292117144, 271.66915794704204,
    1.590620405079854, 0.99989348676669865, 271.80179625231739,
    3.360409062753658, -0.92194044308811329, 271.6817916104477,
    3.6351709573108231, 0.99084865809803091, 271.86662880808427,
    -0.10637260889229627, -0.80506420120276667, -0.24765280605537204,
    0.98463971838264741, 0.99172047752049197, 0.11647206058754805,
    2.1835887004361521, -0.63327470946302045, -0.19835542808978235,
    2.6691221737777013, -0.4410476128906391, 0.071430659019328879,
    3.0513059480333982, 0.88726749799577753, 0.064484775276447692 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t at = 0; at < 30; ++at ) {

    points.push_back(
      glm::dvec3( net[ ( 3 * at ) ], net[ ( 3 * at ) + 1 ],
                  net[ ( 3 * at ) + 2 ] ) );

    weights.push_back( 1.0 );
  }

  surface.control_points = tinynurbs::array2( 6, 5, points );
  surface.weights        = tinynurbs::array2( 6, 5, weights );

  const RationalSurfaceEvaluator evaluator( surface );

  const auto wrap = []( double u ) {

    const double offset = std::fmod( u - STRIP_MIN, PERIOD );

    return STRIP_MIN + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec2 uv0( -399.37070255424464,-27637.360607875835 );
  const glm::dvec2 uv1( -399.37070255423532,27637.360607875835 );

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  const auto departure = [ & ]( long double t ) {

    const double u =
      static_cast< double >( (long double)uv0.x + ( t * ( (long double)uv1.x - uv0.x ) ) );

    const double v =
      static_cast< double >( (long double)uv0.y + ( t * ( (long double)uv1.y - uv0.y ) ) );

    const glm::dvec3 q = evaluator.point( wrap( u ), v );

    const long double dx =
      (long double)q.x - ( (long double)at0.x + ( t * ( (long double)at1.x - at0.x ) ) );
    const long double dy =
      (long double)q.y - ( (long double)at0.y + ( t * ( (long double)at1.y - at0.y ) ) );
    const long double dz =
      (long double)q.z - ( (long double)at0.z + ( t * ( (long double)at1.z - at0.z ) ) );

    return static_cast< double >( sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) ) );
  };

  double sweep = 0.0;

  for ( uint32_t i = 0; i <= 200000; ++i ) {

    sweep = std::max( sweep, departure( (long double)i / 200000.0L ) );
  }

  // The departure lives at the v knots the chord crosses, two ulps of t
  // below each - which is inside the grazing window and nowhere a uniform
  // sweep lands.
  double truth = sweep;

  for ( double knot : surface.knots_v ) {

    long double t =
      ( (long double)knot - uv0.y ) / ( (long double)uv1.y - uv0.y );

    if ( !( t > 0.0L ) || !( t < 1.0L ) ) continue;

    for ( int nudge = -2; nudge <= 2; ++nudge ) {

      long double tt = t;

      for ( int k = 0; k < std::abs( nudge ); ++k ) {
        tt = std::nextafterl( tt, nudge < 0 ? 0.0L : 1.0L );
      }

      truth = std::max( truth, departure( tt ) );
    }
  }

  check( truth > sweep,
         "the grazing window carries a departure a uniform sweep does not "
         "see (" + std::to_string( truth ) + " vs " +
           std::to_string( sweep ) + ")" );

  const double tolerance = 0.0001000306590146737;

  const NurbsDeflectionCertificate certificate(
    evaluator, true, STRIP_MIN, PERIOD,
    tolerance * tolerance );

  double bound = 0.0;

  const CertificateOutcome outcome =
    certificate.bound( uv0, uv1, at0, at1, bound );

  check( outcome != CertificateOutcome::Certified || bound >= truth,
         "and nothing is certified below it (" + std::to_string( bound ) +
           " vs " + std::to_string( truth ) + ")" );
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
  periodicKnotTieStillCrossesTheStrip();
  chordManySheetsOutsideTheStripIsPricedByItsRawParameter();
  tiedStripAndKnotEdgesAreOrderedInParameterSpaceNotChordSpace();
  periodicStartIsReducedTheWayTheCallbackReducesIt();
  steppingAxisBoxBelowResolutionIsDeclined();
  grazingBoxShortOfAStepKnotIsDeclined();
  reductionDriftAcrossAStepKnotIsPriced();
  descendingChordStartsInTheSpanItTraverses();
  errorTermTracksTheControlValues();
  knotClippingCatchesAZigZag();
  wrappedChordIsStillCertified();
  boundHoldsOnARationalPatch();
  boundConvergesUnderSubdivision();

  printf( "%s\n", failures == 0 ? "PASS" : "FAIL" );

  return failures == 0 ? 0 : 1;
}
