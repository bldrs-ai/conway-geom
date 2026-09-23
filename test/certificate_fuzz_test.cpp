/*
 * DIFFERENTIAL FUZZER for the deflection certificate.
 *
 * Six P1s in three rounds were all found by a reviewer reading the code, and
 * twice we convinced ourselves the class was closed by reasoning about it.
 * This checks the property directly instead: generate a random NURBS, a
 * random chord, ask the certificate for a bound, and compare it against a
 * dense sampling of the actual departure. A `Certified` verdict whose bound
 * sits below what sampling can already see is a definite defect.
 *
 * The generator is biased HARD toward the shapes that produced findings:
 * full-multiplicity knots, knot spacings down to one ulp, coordinate
 * magnitudes up to 1e9, weights spanning orders of magnitude, and chords
 * that start or end exactly on a knot.
 *
 * Sampling is a lower bound on the true supremum, so a miss proves nothing;
 * a hit proves a defect. Knot parameters are sampled from BOTH sides, since
 * a full-multiplicity knot is where the surface steps.
 *
 * THIS FUZZER HAS BEEN SEEN TO FAIL, which is the only thing that makes a
 * clean run worth reading. Validated against known-bad code before it was
 * trusted:
 *
 *   - against commit 4b070b9, which carried the fifth finding, it reported
 *     834 violations in 20,000 trials, rediscovering that family without
 *     being told what to look for;
 *   - against a build with the node-placement gradient term removed it
 *     reports violations immediately;
 *   - it found the seventh finding - a chord whose evaluated parameter is
 *     PINNED to a knot because the step per unit chord is below the
 *     parameter's resolution - with no reviewer involved;
 *   - it found the eighth, in the periodic path, within minutes of being
 *     taught to drive periodic charts at all.
 *
 * It also reports which BRANCHES of the certificate it reached, because the
 * eighth finding lived in one that nothing had ever executed. A zero in that
 * line is a hole in this file, not a clean bill of health for the code.
 *
 * HOW LONG IT TAKES TO REDISCOVER EACH FINDING, measured by reintroducing
 * each one and running until this file catches it again - with the
 * unmodified build's own violations subtracted, so that running into the
 * residual does not count as a rediscovery. The harness that produced these
 * reintroduces each finding into a scratch copy of the header reached by an
 * `-I` override, exactly as the red proofs do, and is not committed for the
 * same reason they are not:
 *
 *     finding                                        first caught at trial
 *     #2  node's span resolved by findSpan                              5
 *     #5  ambiguous cross-axis order always taking u                    5
 *     #8  clamp not cleared by a strip restart                      2,752
 *     #4  node-placement gradient term dropped                     17,575
 *     #9a wrap entry recomputed instead of pinned to the edge      45,004
 *     #8b extrapolated box bounded by an in-span gradient          51,982
 *     #6  weight placement error not carried to the floor      NOT CAUGHT
 *
 * #6 IS NOT REDISCOVERED AT ALL, in 200,000 trials across three seeds. It is
 * a real missing link - the weight hull is rebuilt from samples whose
 * perturbation was priced only into the numerator - but its effect is below
 * anything this generator produces, exactly as it was below anything a
 * hand-built case could produce. Nothing here defends it; it is carried on
 * the argument alone, and that is worth knowing.
 *
 * #3 and #7 cannot be reintroduced to measure, because the code that carried
 * them no longer exists: #3's cut-parameter sorting was replaced by the index
 * walk, and #7's stepping-axis guard was DELETED once the span-width test
 * made breaking it undetectable.
 *
 * The default budget below is set above the largest number in that table
 * that is not "NOT CAUGHT".
 */
#include "conway_geometry/operations/deflection_certificate.h"

#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
#include "logging/Logger.h"

void Logger::logError( const char*, ... ) {}
void Logger::logWarning( const char*, ... ) {}

using namespace conway::geometry;

namespace {

std::mt19937_64 rng;

double uniform( double a, double b ) {
  return a + ( ( b - a ) *
               ( static_cast< double >( rng() >> 11 ) / 9007199254740992.0 ) );
}

uint32_t pick( uint32_t n ) { return static_cast< uint32_t >( rng() % n ); }

/** A knot vector of `count` control points at `degree`, with random
 *  interior multiplicities and random - sometimes sub-ulp - spacing. */
std::vector< double > makeKnots(
  uint32_t degree, uint32_t count, double low, double high ) {

  std::vector< double > knots;

  for ( uint32_t i = 0; i <= degree; ++i ) knots.push_back( low );

  const uint32_t interior = count - degree - 1;

  double at = low;

  for ( uint32_t i = 0; i < interior; ) {

    // Spacing: usually ordinary, sometimes down at the representable floor.
    const uint32_t how = pick( 6 );

    double step;

    if ( how == 0 ) {
      step = std::nextafter( at, high ) - at;              // one ulp
    } else if ( how == 1 ) {
      step = ( high - low ) * 1.0e-9;
    } else {
      step = ( high - low ) * uniform( 0.02, 0.3 );
    }

    at = std::min( at + step, std::nextafter( high, low ) );

    // Multiplicity up to degree + 1, which makes the surface STEP.
    uint32_t multiplicity = 1 + pick( degree + 1 );

    multiplicity = std::min( multiplicity, interior - i );

    for ( uint32_t m = 0; m < multiplicity; ++m ) knots.push_back( at );

    i += multiplicity;
  }

  for ( uint32_t i = 0; i <= degree; ++i ) knots.push_back( high );

  return knots;
}

}  // namespace

int main( int argc, char** argv ) {

  // 100,000 by default, which is ~14 seconds. Not a round number picked for
  // comfort: the budget has to be larger than what it takes to rediscover
  // the findings this file exists to catch, and the largest of those that it
  // CAN rediscover needs 79,231 trials. See the table above.
  const uint64_t trials = ( argc > 1 ) ? strtoull( argv[ 1 ], nullptr, 10 ) : 100000;
  const uint64_t seed   = ( argc > 2 ) ? strtoull( argv[ 2 ], nullptr, 10 ) : 1;

  rng.seed( seed );

  uint64_t certified = 0, declined = 0, unsupported = 0, violations = 0;
  uint64_t skipped = 0;

  // Which BRANCHES of the certificate the generator has actually reached.
  // #8 lived in a branch nothing executed, so this is reported alongside the
  // violation count: a zero here is a hole in the instrument, not a clean
  // bill of health for the code.
  NurbsDeflectionCertificate::Counters total;

  const auto accumulate =
    []( NurbsDeflectionCertificate::Counters& into,
        const NurbsDeflectionCertificate::Counters& from ) {
      into.strips += from.strips; into.clampsClear += from.clampsClear;
      into.descending += from.descending; into.rationalPath += from.rationalPath;
      into.ambiguous += from.ambiguous;
      into.extrapolated += from.extrapolated; into.unplaceable += from.unplaceable;
      into.illConditioned += from.illConditioned; into.tooManySpans += from.tooManySpans;
      into.badPiece += from.badPiece;
    };

  for ( uint64_t trial = 0; trial < trials; ++trial ) {

    // TWO GENERATORS. The free one explores widely; the ALIGNED one builds
    // the shape that a purely random search will never hit - a u knot and a
    // v knot whose chord parameters round together, with full multiplicity
    // on both axes so the four span combinations are independent patches.
    // That is the family of the fifth finding, and a fuzzer that cannot
    // produce it is not an instrument worth trusting.
    const bool aligned = ( pick( 2 ) == 0 );

    if ( aligned ) {

      const double magnitude =
        std::pow( 10.0, uniform( 0.0, 9.0 ) );

      const double low = -magnitude, high = magnitude;

      // A knot somewhere inside, and its neighbour a few ulps away on the
      // other axis - so the two crossings land on the same parameter.
      const double base = uniform( low * 0.5, high * 0.5 );

      double other = base;

      for ( uint32_t step = 0, n = pick( 4 ); step <= n; ++step ) {
        other = std::nextafter( other, high );
      }

      const uint32_t degree = 1 + pick( 2 );

      tinynurbs::RationalSurface3d surface;

      surface.degree_u = degree;
      surface.degree_v = degree;

      surface.knots_u.clear();
      surface.knots_v.clear();

      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_u.push_back( low );
      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_u.push_back( other );
      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_u.push_back( high );

      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_v.push_back( low );
      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_v.push_back( base );
      for ( uint32_t i = 0; i <= degree; ++i ) surface.knots_v.push_back( high );

      const uint32_t count = 2 * ( degree + 1 );

      const bool rationalAligned = pick( 3 ) == 0;

      // Height concentrated on ONE span combination, so a patch that is
      // skipped shows up and one that is covered does not.
      const uint32_t hotU = pick( 2 ), hotV = pick( 2 );

      const double tall = std::pow( 10.0, uniform( -1.0, 2.0 ) );

      std::vector< glm::dvec3 > points;
      std::vector< double >     weights;

      for ( uint32_t i = 0; i < count; ++i ) {
        for ( uint32_t j = 0; j < count; ++j ) {
          const bool hot =
            ( ( i < degree + 1 ) == ( hotU == 0 ) ) &&
            ( ( j < degree + 1 ) == ( hotV == 0 ) );
          points.push_back(
            glm::dvec3( surface.knots_u[ i + 1 ], surface.knots_v[ j + 1 ],
                        hot ? tall : 0.0 ) );
          weights.push_back(
            rationalAligned ? std::pow( 10.0, uniform( -2.0, 2.0 ) ) : 1.0 );
        }
      }

      surface.control_points = tinynurbs::array2( count, count, points );
      surface.weights        = tinynurbs::array2( count, count, weights );

      const RationalSurfaceEvaluator alignedEvaluator( surface );

      if ( !alignedEvaluator.supportsFastPath() ) { ++skipped; continue; }

      const glm::dvec2 uv0( low, low );
      const glm::dvec2 uv1( high, high );

      const glm::dvec3 at0 = alignedEvaluator.point( uv0.x, uv0.y );
      const glm::dvec3 at1 = alignedEvaluator.point( uv1.x, uv1.y );

      if ( !std::isfinite( at0.x ) || !std::isfinite( at1.x ) ) { ++skipped; continue; }

      const double tolerance = std::pow( 10.0, uniform( -6.0, 0.0 ) );

      NurbsDeflectionCertificate certificate(
        alignedEvaluator, false, 0.0, 0.0, tolerance * tolerance );

      double bound = 0.0;

      const CertificateOutcome outcome =
        certificate.bound( uv0, uv1, at0, at1, bound );

      accumulate( total, certificate.counters() );

      if ( outcome != CertificateOutcome::Certified ) { ++declined; continue; }

      ++certified;

      double truth = 0.0;

      for ( uint32_t i = 0; i <= 4000; ++i ) {
        const long double t = (long double)i / 4000.0L;
        const double u = (double)( (long double)low + t * ( (long double)high - low ) );
        const glm::dvec3 q = alignedEvaluator.point( u, u );
        const long double cx = (long double)at0.x + t * ( (long double)at1.x - at0.x );
        const long double cy = (long double)at0.y + t * ( (long double)at1.y - at0.y );
        const long double cz = (long double)at0.z + t * ( (long double)at1.z - at0.z );
        truth = std::max( truth, (double)sqrtl(
          ( (long double)q.x - cx ) * ( (long double)q.x - cx ) +
          ( (long double)q.y - cy ) * ( (long double)q.y - cy ) +
          ( (long double)q.z - cz ) * ( (long double)q.z - cz ) ) );
      }

      // And the four corner combinations around the aligned knots, which is
      // where a skipped patch lives.
      for ( double u : { std::nextafter( base, low ), base, other,
                         std::nextafter( other, high ) } ) {
        for ( double v : { std::nextafter( base, low ), base, other,
                           std::nextafter( other, high ) } ) {
          const double t = ( u - low ) / ( high - low );
          const glm::dvec3 q = alignedEvaluator.point( u, v );
          truth = std::max( truth, glm::length(
            q - ( at0 * ( 1.0 - t ) + at1 * t ) ) );
        }
      }

      if ( std::isfinite( truth ) && bound < truth ) {
        ++violations;
        if ( violations <= 8 ) {
          printf( "VIOLATION(aligned) trial=%llu bound=%.10g < truth=%.10g "
                  "ratio=%.4g  degree=%u magnitude=%.3g rational=%d tol=%.3g\n",
                  (unsigned long long)trial, bound, truth,
                  truth / std::max( bound, 1e-300 ), degree, magnitude,
                  (int)rationalAligned, tolerance );
        }
      }

      continue;
    }

    const uint32_t degreeU = 1 + pick( 3 );
    const uint32_t degreeV = 1 + pick( 3 );

    const uint32_t countU = degreeU + 1 + pick( 4 );
    const uint32_t countV = degreeV + 1 + pick( 4 );

    // Parameter domains, sometimes far from the origin.
    const double magnitude =
      ( pick( 3 ) == 0 ) ? std::pow( 10.0, uniform( 0.0, 9.0 ) ) : 1.0;

    const double lowU = -magnitude, highU = magnitude;
    const double lowV = -magnitude, highV = magnitude;

    tinynurbs::RationalSurface3d surface;

    surface.degree_u = degreeU;
    surface.degree_v = degreeV;
    surface.knots_u  = makeKnots( degreeU, countU, lowU, highU );
    surface.knots_v  = makeKnots( degreeV, countV, lowV, highV );

    if ( surface.knots_u.size() != countU + degreeU + 1 ||
         surface.knots_v.size() != countV + degreeV + 1 ) { ++skipped; continue; }

    const double spread =
      std::pow( 10.0, uniform( -3.0, 3.0 ) );

    const bool rational = pick( 2 ) == 0;

    std::vector< glm::dvec3 > points;
    std::vector< double >     weights;

    for ( uint32_t i = 0; i < countU; ++i ) {
      for ( uint32_t j = 0; j < countV; ++j ) {
        points.push_back(
          glm::dvec3( uniform( -spread, spread ),
                      uniform( -spread, spread ),
                      uniform( -spread, spread ) ) );
        weights.push_back(
          rational ? std::pow( 10.0, uniform( -2.0, 2.0 ) ) : 1.0 );
      }
    }

    // A periodic chart is only ever built on a surface that is genuinely
    // CLOSED in u, so make it so - otherwise the wrap introduces a step that
    // no production surface has, and a violation would say nothing about
    // what ships.
    const bool wantPeriodic = pick( 3 ) == 0;

    if ( wantPeriodic ) {

      for ( uint32_t j = 0; j < countV; ++j ) {

        const uint32_t firstAt = ( 0 * countV ) + j;
        const uint32_t lastAt  = ( ( countU - 1 ) * countV ) + j;

        points[ lastAt ]  = points[ firstAt ];
        weights[ lastAt ] = weights[ firstAt ];
      }
    }

    surface.control_points = tinynurbs::array2( countU, countV, points );
    surface.weights        = tinynurbs::array2( countU, countV, weights );

    const RationalSurfaceEvaluator evaluator( surface );

    if ( !evaluator.supportsFastPath() ) { ++skipped; continue; }

    // A chord, sometimes with an end pinned exactly to a knot.
    const auto endpoint = [ & ]( bool onKnot ) {
      if ( onKnot ) {
        return glm::dvec2( surface.knots_u[ pick( (uint32_t)surface.knots_u.size() ) ],
                           surface.knots_v[ pick( (uint32_t)surface.knots_v.size() ) ] );
      }
      return glm::dvec2( uniform( lowU * 1.5, highU * 1.5 ),
                         uniform( lowV, highV ) );
    };

    const glm::dvec2 uv0 = endpoint( pick( 3 ) == 0 );
    const glm::dvec2 uv1 = endpoint( pick( 3 ) == 0 );

    if ( uv0 == uv1 ) { ++skipped; continue; }

    const auto preWrap = [ & ]( double u ) {
      if ( !wantPeriodic ) return u;
      const double offset = std::fmod( u - lowU, highU - lowU );
      return lowU + ( offset < 0.0 ? offset + ( highU - lowU ) : offset );
    };

    const glm::dvec3 at0 = evaluator.point( preWrap( uv0.x ), uv0.y );
    const glm::dvec3 at1 = evaluator.point( preWrap( uv1.x ), uv1.y );

    if ( !std::isfinite( at0.x ) || !std::isfinite( at1.x ) ) { ++skipped; continue; }

    const double tolerance = std::pow( 10.0, uniform( -9.0, 1.0 ) );

    // The b-spline call site drives a PERIODIC chart on a closed surface,
    // where the callback evaluates `point( wrap( u ), v )`. Exercise that
    // too - without it a third of the certificate's branches never run.
    const bool periodic = wantPeriodic;

    const double stripMin    = lowU;
    const double stripPeriod = highU - lowU;

    const auto wrap = [ & ]( double u ) {
      if ( !periodic ) return u;
      const double offset = std::fmod( u - stripMin, stripPeriod );
      return stripMin + ( offset < 0.0 ? offset + stripPeriod : offset );
    };

    NurbsDeflectionCertificate certificate(
      evaluator, periodic, stripMin, stripPeriod, tolerance * tolerance );

    double bound = 0.0;

    const CertificateOutcome outcome =
      certificate.bound( uv0, uv1, at0, at1, bound );

    accumulate( total, certificate.counters() );

    if ( outcome == CertificateOutcome::Unsupported ) { ++unsupported; continue; }
    if ( outcome == CertificateOutcome::Inconclusive ) { ++declined; continue; }

    ++certified;

    // Dense truth. Every sample is a LOWER bound on the true supremum.
    const auto departure = [ & ]( long double t ) {
      const double u = (double)( (long double)uv0.x + t * ( (long double)uv1.x - uv0.x ) );
      const double v = (double)( (long double)uv0.y + t * ( (long double)uv1.y - uv0.y ) );
      const glm::dvec3 q = evaluator.point( wrap( u ), v );
      const long double cx = (long double)at0.x + t * ( (long double)at1.x - at0.x );
      const long double cy = (long double)at0.y + t * ( (long double)at1.y - at0.y );
      const long double cz = (long double)at0.z + t * ( (long double)at1.z - at0.z );
      const long double dx = (long double)q.x - cx;
      const long double dy = (long double)q.y - cy;
      const long double dz = (long double)q.z - cz;
      return (double)sqrtl( dx * dx + dy * dy + dz * dz );
    };

    double truth = 0.0;

    for ( uint32_t i = 0; i <= 3000; ++i ) {
      truth = std::max( truth, departure( (long double)i / 3000.0L ) );
    }

    // And at every knot the chord crosses, from BOTH sides - that is where a
    // full-multiplicity surface steps, and no uniform sweep lands there.
    const auto probeKnots = [ & ]( const std::vector< double >& knots,
                                   double a, double b ) {
      if ( a == b ) return;
      for ( double knot : knots ) {
        const long double t = (long double)( knot - a ) / ( (long double)b - a );
        if ( t <= 0.0L || t >= 1.0L ) continue;
        for ( int step = -2; step <= 2; ++step ) {
          long double tt = t;
          for ( int k = 0; k < std::abs( step ); ++k ) {
            tt = std::nextafterl( tt, step < 0 ? 0.0L : 1.0L );
          }
          truth = std::max( truth, departure( tt ) );
        }
      }
    };

    probeKnots( surface.knots_u, uv0.x, uv1.x );
    probeKnots( surface.knots_v, uv0.y, uv1.y );

    if ( !std::isfinite( truth ) ) continue;

    if ( bound < truth ) {
      ++violations;
      if ( violations <= 2 ) {
        printf( "VIOLATION trial=%llu  bound=%.10g < truth=%.10g  ratio=%.4g\n",
                (unsigned long long)trial, bound, truth,
                truth / std::max( bound, 1e-300 ) );
        // Dump the whole case as C++, so it can be replayed standalone.
        printf( "  s.degree_u=%u; s.degree_v=%u;\n", degreeU, degreeV );
        printf( "  s.knots_u={" );
        for ( double k : surface.knots_u ) printf( "%.17g,", k );
        printf( "};\n  s.knots_v={" );
        for ( double k : surface.knots_v ) printf( "%.17g,", k );
        printf( "};\n  const double P[]={" );
        for ( const glm::dvec3& q : points ) printf( "%.17g,%.17g,%.17g,", q.x, q.y, q.z );
        printf( "};\n  const double W[]={" );
        for ( double q : weights ) printf( "%.17g,", q );
        printf( "};\n  const uint32_t NU=%u, NV=%u;\n", countU, countV );
        printf( "  glm::dvec2 uv0(%.17g,%.17g), uv1(%.17g,%.17g); double TOL=%.17g;\n",
                uv0.x, uv0.y, uv1.x, uv1.y, tolerance );
        printf( "  const bool PERIODIC=%d; const double SMIN=%.17g, SPER=%.17g;\n",
                (int)periodic, stripMin, stripPeriod );
      }
    }
  }

  printf( "BRANCHES reached: strips=%llu clampsClear=%llu descending=%llu "
          "rational=%llu ambiguous=%llu extrapolated=%llu "
          "unplaceable=%llu illCond=%llu manyPieces=%llu badPiece=%llu\n",
          (unsigned long long)total.strips, (unsigned long long)total.clampsClear,
          (unsigned long long)total.descending, (unsigned long long)total.rationalPath,
          (unsigned long long)total.ambiguous,
          (unsigned long long)total.extrapolated, (unsigned long long)total.unplaceable,
          (unsigned long long)total.illConditioned, (unsigned long long)total.tooManySpans,
          (unsigned long long)total.badPiece );
  printf( "trials=%llu certified=%llu declined=%llu unsupported=%llu skipped=%llu "
          "VIOLATIONS=%llu\n",
          (unsigned long long)trials, (unsigned long long)certified,
          (unsigned long long)declined, (unsigned long long)unsupported,
          (unsigned long long)skipped, (unsigned long long)violations );

  return violations == 0 ? 0 : 1;
}
