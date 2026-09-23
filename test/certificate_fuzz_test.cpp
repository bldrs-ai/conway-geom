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
 *     taught to drive periodic charts at all;
 *   - and it found FOUR MORE - #9, #10, #11 and #12 - in the round that was
 *     supposed to be the clean one. See the table below.
 *
 * It also reports which BRANCHES of the certificate it reached, because the
 * eighth finding lived in one that nothing had ever executed. A zero in that
 * line is a hole in this file, not a clean bill of health for the code.
 *
 * FAR-SHEET CHORDS. The generator used to draw u from
 * [ lowU * 1.5, highU * 1.5 ], so ( u - lowU ) / P landed in [ -0.25, 1.25 ]
 * and the walk's shift was never more than ONE period. A whole regime was
 * therefore unreachable: a chord many sheets outside its strip, where an ulp
 * of the RAW parameter is orders of magnitude larger than an ulp of the
 * wrapped one. That regime is what the `| shift |` term in `roundingU`
 * exists for, and it is why that term had no red proof for a round. It has
 * one now - see the table. Both endpoints are pushed out by up to 2^20
 * periods on half the periodic trials.
 *
 * THE SHAPE IS BUILT, NOT WAITED FOR. There are three generators. The free
 * one explores; the ALIGNED one constructs two knots whose chord parameters
 * round together, which is the fifth finding's family; the GRAZING one
 * constructs the eleventh's - a chord whose extent on one axis is below the
 * parameter's resolution, ending exactly on a knot the surface steps at,
 * with knots on the other axis inside that window so the walk emits a box
 * the caller resolves entirely elsewhere.
 *
 * That third one is the reason this round exists. The eleventh finding was
 * reached ONCE in 2,800,000 free trials, which is luck rather than
 * instrumentation - and the guard it red-proved had already been deleted for
 * want of a proof. With the shape constructed it is caught at trial 94, and
 * the generator immediately produced three more findings of the same family.
 *
 * THE PERIODIC PATH SHIPS DISABLED AND IS FUZZED ANYWAY.
 * `CERTIFICATE_CERTIFIES_PERIODIC_CHARTS` is 0 in a release build, so a
 * periodic chord takes the sampled test; this file defines it to 1. The
 * reason to ship it off is that the size of that family is unknown, which is
 * precisely the reason to keep looking - a fuzzer that followed the shipping
 * default would stop looking where the risk is concentrated and report
 * nothing, and the nothing would mean nothing. About a third of the trials
 * here drive periodic charts.
 *
 * AND THE SEAM IS WATCHED DIRECTLY. CERTIFICATE_AUDIT_CALLER_SPANS, defined
 * above, compares the span the walk assigned each box against the span the
 * CALL SITE's own `findSpan` names for the same doubles. `inside` - a box
 * placed on a span that does not contain it - is a hard failure and this
 * program exits non-zero on it. The fifteenth finding was found that way and
 * produced no violation at all in 800,000 trials, so nothing else would have
 * found it.
 *
 * HOW LONG IT TAKES TO REDISCOVER EACH FINDING, measured by reintroducing
 * each one and running until this file catches it again - with the
 * unmodified build's own violations subtracted, so that running into the
 * residual does not count as a rediscovery ( `test/certificate_rediscover.sh` ):
 *
 *     finding                                        first caught at trial
 *     #2  node's span resolved by findSpan                              5
 *     #5  ambiguous cross-axis order always taking u                    5
 *     #11 stepping-axis arm of the unplaceable test                    94
 *     #13 boundary test taken as exact equality                     1,553
 *     #4b node placement priced by the WRAPPED parameter             1,475
 *     #8  clamp not cleared by a strip restart                      2,752
 *     #14 strip magnitude left out of roundingU                     5,387
 *     #4  node-placement gradient term dropped                     17,575
 *     #9  strip boundary suppressed by a knot tie                  24,405
 *     #12 periodic start reduced with floor, not fmod               32,525
 *     #9a wrap entry recomputed instead of pinned to the edge      45,004
 *     #8b extrapolated box bounded by an in-span gradient          51,982
 *     #10 strip/knot order decided in chord space                 124,473
 *     #15 rounding priced at the box's magnitude, not the chord's
 *                                       AUDIT ONLY - 1 in ~600,000 trials
 *     #6  weight placement error not carried to the floor      NOT CAUGHT
 *
 * TWO ENTRIES ARE THE POINT OF THIS FILE. #15 is caught by the audit and by
 * nothing else - it never became a wrong number. #6 is caught by nothing at
 * all, in 200,000 trials across three seeds, and is carried on its argument.
 *
 * AND THE HISTORY OF "NO RED PROOF" IN THIS PR IS THREE FOR THREE. The
 * `| shift |` term ( #4b ), the stepping-axis arm ( #11 ) and the span-width
 * guard ( 3b8 in `certificate_redprove.sh` ) were each at some point
 * defensible only by argument. Two of the three turned out to be load-
 * bearing as soon as the generator reached their shape, and one of those had
 * already been DELETED on the reasoning that it cost declines and caught
 * nothing. 3b8 is the one still in that state; it now costs 4 declines in
 * 120,000 trials and none on the corpus, and it stays.
 *
 * #10 fired at 124,473 when it was written and NO LONGER DOES: #12 moves the
 * periodic start by 3.7e-9, which is enough to stop that chord's two cut
 * parameters rounding together. The fix is kept on its argument; the number
 * above is the measurement that existed, not one that reproduces today.
 *
 * #3 and #7 cannot be reintroduced to measure, because the code that carried
 * them no longer exists: #3's cut-parameter sorting was replaced by the index
 * walk. (#7's stepping-axis guard was deleted and is now BACK, as #11.)
 *
 * The default budget below is set above the largest number in that table
 * that is not "NOT CAUGHT".
 */
// Watch the walk/caller seam directly rather than by its consequences. This
// has to come before the include - see CERTIFICATE_AUDIT_CALLER_SPANS.
#define CERTIFICATE_AUDIT_CALLER_SPANS 1

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

#include <cstdio>
#include <cstdlib>
#include <algorithm>
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

  // 100,000 by default. Not a round number picked for comfort: the budget
  // has to be larger than what it takes to rediscover the findings this file
  // exists to catch, and the largest of those that it can rediscover by a
  // VIOLATION needs 124,473 - so the default no longer covers the table on
  // one seed, and the sweep that gates a change runs fourteen. See above.
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
      into.callerSpanAgree += from.callerSpanAgree;
      into.callerSpanEdge += from.callerSpanEdge;
      into.callerSpanInside += from.callerSpanInside;
      into.callerSpanWhole += from.callerSpanWhole;
      into.callerSpanPair += from.callerSpanPair;
      into.callerSpanWholeStep += from.callerSpanWholeStep;
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

    // GENERATOR THREE: THE GRAZING CHORD, built rather than waited for.
    //
    // This is the eleventh finding's shape. The free generator reached it
    // once, at trial 173,559 on one seed out of fourteen, which means it is
    // effectively unreachable at any budget this file would run by default -
    // and the finding was then found by its CONSEQUENCE, a bound under a
    // sampled truth, rather than by anything watching the seam. Both halves
    // of that are fixed here: the shape is constructed, and
    // CERTIFICATE_AUDIT_CALLER_SPANS watches the seam directly.
    //
    // The shape has three parts, and all three are needed:
    //
    //   1. a knot K of multiplicity degree + 1 on one axis, so the surface
    //      STEPS there rather than merely kinking;
    //   2. a chord whose extent on that axis is so small that a MACROSCOPIC
    //      range of t has a rounded parameter equal to K. The fraction is
    //      ulp( K ) / ( 2 |d| ), so `d` is chosen from a target fraction
    //      rather than the other way round;
    //   3. the chord ENDING on K, so there is no next box to own that
    //      double - which is what makes it unsound rather than merely
    //      covered by the neighbour.
    //
    // Multiplicity `degree` is generated as well as `degree + 1`. At
    // multiplicity `degree` the surface is C0 - the two polynomials AGREE at
    // K - so it should be harmless, and `stepsAnywhere`'s threshold says so.
    // That is an argument, and arguments are what this file exists to check.
    const bool grazing = !aligned && ( pick( 2 ) == 0 );

    if ( grazing ) {

      const uint32_t degreeSlow = 1 + pick( 3 );
      const uint32_t degreeFast = 1 + pick( 3 );

      // The axis the chord grazes along. Generated both ways round: the
      // walk, the strip and the periodic reduction are all u-specific, so a
      // v-only version of this would leave half the seam untested.
      const bool grazeU = ( pick( 2 ) == 0 );

      const double magnitude = std::pow( 10.0, uniform( -2.0, 6.0 ) );

      const double low = -magnitude, high = magnitude;

      // K, strictly inside, and the fraction of the chord that will round
      // onto it.
      const double K = uniform( low * 0.6, high * 0.6 );

      const double fraction = std::pow( 10.0, uniform( -9.0, -2.0 ) );

      const double ulpK =
        std::nextafter( std::abs( K ) + 1.0,
                        std::numeric_limits< double >::infinity() ) -
        ( std::abs( K ) + 1.0 );

      double step = ulpK / ( 2.0 * fraction );

      if ( !( step > 0.0 ) || !std::isfinite( step ) ) { ++skipped; continue; }

      // MULTIPLICITY degree + 1 STEPS, degree ONLY KINKS. Both are built.
      const uint32_t slowMultiplicity =
        ( pick( 4 ) == 0 ) ? degreeSlow : ( degreeSlow + 1 );

      std::vector< double > slowKnots;

      for ( uint32_t i = 0; i <= degreeSlow; ++i ) slowKnots.push_back( low );
      for ( uint32_t i = 0; i < slowMultiplicity; ++i ) slowKnots.push_back( K );
      for ( uint32_t i = 0; i <= degreeSlow; ++i ) slowKnots.push_back( high );

      // The axis the chord actually travels along gets several interior
      // knots, so the walk emits a run of boxes rather than one.
      // MORE SPANS THAN THE WALK WILL CARRY, one time in four. The walk
      // gives up at CERTIFICATE_MAX_PIECES, and until this that branch had
      // never executed in any run - the free generator's surfaces are not
      // deep enough and its chords do not cross every span of one.
      const uint32_t fastSpans =
        ( pick( 4 ) == 0 ) ? ( 8 + pick( 24 ) ) : ( 1 + pick( 6 ) );

      std::vector< double > fastKnots;

      for ( uint32_t i = 0; i <= degreeFast; ++i ) fastKnots.push_back( low );

      for ( uint32_t i = 1; i < fastSpans; ++i ) {

        const double at =
          low + ( ( high - low ) * static_cast< double >( i ) /
                  static_cast< double >( fastSpans ) );

        const uint32_t multiplicity = 1 + pick( degreeFast + 1 );

        for ( uint32_t m = 0; m < multiplicity; ++m ) fastKnots.push_back( at );
      }

      // AND SOME KNOTS INSIDE THE GRAZING WINDOW. Without these the walk
      // never emits a box whose extent on the slow axis is zero: every box
      // straddles the point where the rounded parameter reaches K, so its
      // two corners differ and only one node sits on the knot. Seed 7 had
      // them - its u knots were 1e-9 from the chord's end while the grazing
      // window was 2.8e-8 wide - and that is the difference between a box
      // the caller resolves elsewhere at ONE node and one it resolves
      // elsewhere at EVERY node.
      {
        const uint32_t inWindow = 1 + pick( 3 );

        std::vector< double > late;

        for ( uint32_t i = 0; i < inWindow; ++i ) {

          const double atT =
            1.0 - ( fraction * std::pow( 10.0, uniform( -2.0, 0.0 ) ) );

          const double value = low + ( ( high - low ) * atT );

          if ( value > fastKnots.back() && value < high ) {
            late.push_back( value );
          }
        }

        std::sort( late.begin(), late.end() );

        for ( double value : late ) {

          const uint32_t multiplicity = 1 + pick( degreeFast + 1 );

          for ( uint32_t m = 0; m < multiplicity; ++m ) {
            fastKnots.push_back( value );
          }
        }
      }

      for ( uint32_t i = 0; i <= degreeFast; ++i ) fastKnots.push_back( high );

      tinynurbs::RationalSurface3d surface;

      surface.degree_u = grazeU ? degreeSlow : degreeFast;
      surface.degree_v = grazeU ? degreeFast : degreeSlow;
      surface.knots_u  = grazeU ? slowKnots : fastKnots;
      surface.knots_v  = grazeU ? fastKnots : slowKnots;

      const uint32_t countU =
        static_cast< uint32_t >( surface.knots_u.size() ) - surface.degree_u - 1;

      const uint32_t countV =
        static_cast< uint32_t >( surface.knots_v.size() ) - surface.degree_v - 1;

      if ( countU < surface.degree_u + 1 || countV < surface.degree_v + 1 ) {
        ++skipped; continue;
      }

      // The two sides of K have to be FAR APART, or a step that is taken on
      // the wrong side of it costs nothing and the case proves nothing. The
      // index of the last control row before K is what divides them.
      const uint32_t before = degreeSlow + slowMultiplicity - 1;

      const double relief = std::pow( 10.0, uniform( 0.0, 3.0 ) );

      const bool rationalGrazing = ( pick( 4 ) == 0 );

      const bool wantPeriodicGrazing = ( pick( 3 ) == 0 );

      std::vector< glm::dvec3 > points;
      std::vector< double >     weights;

      for ( uint32_t i = 0; i < countU; ++i ) {

        for ( uint32_t j = 0; j < countV; ++j ) {

          const uint32_t slowAt = grazeU ? i : j;
          const uint32_t fastAt = grazeU ? j : i;

          const double side = ( slowAt < before ) ? 0.0 : relief;

          points.push_back(
            glm::dvec3( uniform( -1.0, 1.0 ) + static_cast< double >( fastAt ),
                        uniform( -1.0, 1.0 ),
                        side + uniform( -0.25, 0.25 ) ) );

          weights.push_back(
            rationalGrazing ? std::pow( 10.0, uniform( -2.0, 2.0 ) ) : 1.0 );
        }
      }

      if ( wantPeriodicGrazing ) {

        // CLOSED in u, which is what a periodic chart requires.
        for ( uint32_t j = 0; j < countV; ++j ) {
          points[ ( ( countU - 1 ) * countV ) + j ]  = points[ j ];
          weights[ ( ( countU - 1 ) * countV ) + j ] = weights[ j ];
        }
      }

      surface.control_points = tinynurbs::array2( countU, countV, points );
      surface.weights        = tinynurbs::array2( countU, countV, weights );

      const RationalSurfaceEvaluator grazeEvaluator( surface );

      if ( !grazeEvaluator.supportsFastPath() ) { ++skipped; continue; }

      // APPROACH K, AND END ON IT. Both directions of approach are built:
      // from below, `findSpan` at K names the span above and the walk names
      // the one below, which is the disagreement; from above they agree, and
      // that arm is generated so the difference between them is measured
      // rather than assumed.
      const bool fromBelow = ( pick( 4 ) != 0 );

      if ( !fromBelow ) step = -step;

      const double slowStart = K - step;
      const double slowEnd   = K;

      if ( slowStart == slowEnd ) { ++skipped; continue; }

      // ... and travel the whole of the other axis, so the boxes cover real
      // stretches of chord.
      const double fastStart = low;
      const double fastEnd   = high;

      const glm::dvec2 uv0(
        grazeU ? slowStart : fastStart, grazeU ? fastStart : slowStart );

      const glm::dvec2 uv1(
        grazeU ? slowEnd : fastEnd, grazeU ? fastEnd : slowEnd );

      const double stripMinG    = low;
      const double stripPeriodG = high - low;

      const auto wrapG = [ & ]( double u ) {
        if ( !wantPeriodicGrazing ) return u;
        const double offset = std::fmod( u - stripMinG, stripPeriodG );
        return stripMinG + ( offset < 0.0 ? offset + stripPeriodG : offset );
      };

      const glm::dvec3 at0 = grazeEvaluator.point( wrapG( uv0.x ), uv0.y );
      const glm::dvec3 at1 = grazeEvaluator.point( wrapG( uv1.x ), uv1.y );

      if ( !std::isfinite( at0.x ) || !std::isfinite( at1.x ) ) {
        ++skipped; continue;
      }

      const double tolerance = std::pow( 10.0, uniform( -9.0, 1.0 ) );

      NurbsDeflectionCertificate certificate(
        grazeEvaluator, wantPeriodicGrazing, stripMinG, stripPeriodG,
        tolerance * tolerance );

      double bound = 0.0;

      const CertificateOutcome outcome =
        certificate.bound( uv0, uv1, at0, at1, bound );

      accumulate( total, certificate.counters() );

      if ( outcome == CertificateOutcome::Unsupported ) { ++unsupported; continue; }
      if ( outcome == CertificateOutcome::Inconclusive ) { ++declined; continue; }

      ++certified;

      const auto departureG = [ & ]( long double t ) {
        const double u = (double)( (long double)uv0.x + t * ( (long double)uv1.x - uv0.x ) );
        const double v = (double)( (long double)uv0.y + t * ( (long double)uv1.y - uv0.y ) );
        const glm::dvec3 q = grazeEvaluator.point( wrapG( u ), v );
        const long double cx = (long double)at0.x + t * ( (long double)at1.x - at0.x );
        const long double cy = (long double)at0.y + t * ( (long double)at1.y - at0.y );
        const long double cz = (long double)at0.z + t * ( (long double)at1.z - at0.z );
        const long double dx = (long double)q.x - cx;
        const long double dy = (long double)q.y - cy;
        const long double dz = (long double)q.z - cz;
        return (double)sqrtl( ( dx * dx ) + ( dy * dy ) + ( dz * dz ) );
      };

      double truth = 0.0;

      for ( uint32_t i = 0; i <= 3000; ++i ) {
        truth = std::max( truth, departureG( (long double)i / 3000.0L ) );
      }

      // THE SWEEP ABOVE CANNOT SEE THIS CASE. The stretch whose rounded
      // parameter equals K is `fraction` of the chord, and `fraction` goes
      // down to 1e-9. Sample it directly, geometrically, from one part in
      // ten of it up to the chord's end.
      for ( uint32_t i = 0; i <= 400; ++i ) {

        const long double reach =
          (long double)fraction * 10.0L *
          powl( 1.0e-4L, (long double)i / 400.0L );

        truth = std::max( truth, departureG( 1.0L - reach ) );
      }

      for ( int nudge = 0; nudge <= 4; ++nudge ) {

        long double tt = 1.0L;

        for ( int k = 0; k < nudge; ++k ) tt = std::nextafterl( tt, 0.0L );

        truth = std::max( truth, departureG( tt ) );
      }

      if ( !std::isfinite( truth ) ) continue;

      if ( bound < truth ) {
        ++violations;
        if ( violations <= 4 ) {
          printf( "VIOLATION(grazing) trial=%llu bound=%.10g < truth=%.10g "
                  "ratio=%.4g  grazeU=%d mult=%u/%u fraction=%.3g "
                  "periodic=%d fromBelow=%d rational=%d\n",
                  (unsigned long long)trial, bound, truth,
                  truth / std::max( bound, 1e-300 ), (int)grazeU,
                  slowMultiplicity, degreeSlow, fraction,
                  (int)wantPeriodicGrazing, (int)fromBelow,
                  (int)rationalGrazing );
          printf( "  s.degree_u=%u; s.degree_v=%u;\n",
                  surface.degree_u, surface.degree_v );
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
                  (int)wantPeriodicGrazing, stripMinG, stripPeriodG );
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

    glm::dvec2 uv0 = endpoint( pick( 3 ) == 0 );
    glm::dvec2 uv1 = endpoint( pick( 3 ) == 0 );

    // FAR-SHEET VARIANT. The stock generator draws u from
    // [ lowU * 1.5, highU * 1.5 ], so ( u - lowU ) / P lands in [ -0.25, 1.25 ]
    // and the walk's shift is never more than ONE period. That is the whole
    // reason the `| shift |` term in roundingU has no red proof: the regime it
    // prices - a chord many sheets outside the strip, where an ulp of the raw
    // parameter is orders of magnitude larger than an ulp of the wrapped one -
    // is unreachable. This pushes both endpoints out by up to 2^20 periods.
    if ( wantPeriodic && pick( 2 ) == 0 ) {

      const double sheets =
        (double)( (int64_t)pick( 2097153 ) - 1048576 );

      const double offset = sheets * ( highU - lowU );

      uv0.x += offset;
      uv1.x += offset;
    }

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

  // THE WALK/CALLER SEAM. `inside` must be zero - anything else means a box
  // was placed on a span that does not contain it. `whole` is the eleventh
  // finding's signature, and a zero there means this generator is not
  // reaching that shape, NOT that the shape is impossible.
  printf( "CALLER SPANS: agree=%llu edge=%llu inside=%llu whole=%llu "
          "wholeOnStep=%llu pair=%llu\n",
          (unsigned long long)total.callerSpanAgree,
          (unsigned long long)total.callerSpanEdge,
          (unsigned long long)total.callerSpanInside,
          (unsigned long long)total.callerSpanWhole,
          (unsigned long long)total.callerSpanWholeStep,
          (unsigned long long)total.callerSpanPair );
  // THE AUDIT IS A GATE, NOT A READOUT. `inside` means a box was placed on a
  // span that does not contain it, which no amount of sampling is guaranteed
  // to turn into a violation - the fifteenth finding was found this way and
  // produced none in 800,000 trials. A run that reports it has failed.
  if ( total.callerSpanInside > 0 ) {
    printf( "AUDIT FAILURE: %llu box(es) placed on a span that does not "
            "contain them\n",
            (unsigned long long)total.callerSpanInside );
  }

  printf( "trials=%llu certified=%llu declined=%llu unsupported=%llu skipped=%llu "
          "VIOLATIONS=%llu\n",
          (unsigned long long)trials, (unsigned long long)certified,
          (unsigned long long)declined, (unsigned long long)unsupported,
          (unsigned long long)skipped, (unsigned long long)violations );

  return ( violations == 0 && total.callerSpanInside == 0 ) ? 0 : 1;
}
