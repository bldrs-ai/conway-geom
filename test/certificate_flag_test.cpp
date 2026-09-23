/*
 * THE SWITCH ITSELF, BUILT BOTH WAYS.
 *
 * `CERTIFICATE_CERTIFIES_PERIODIC_CHARTS` decides whether a periodic chord
 * is certified or handed back to the sampled test. A flag that does not
 * actually gate anything is the same failure as a test that does not run,
 * and this PR has already had one of those - so this file is compiled twice,
 * with the constant at 0 and at 1, and asserts the OPPOSITE outcome each
 * time. `run_native_tests.sh` builds it both ways.
 *
 * It also asserts the gate is NARROW: the same surface, walked by a chord
 * that does not use the periodic chart, is certified either way. Without
 * that half, a change that disabled the certificate outright would pass.
 */
#include "conway_geometry/operations/deflection_certificate.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "logging/Logger.h"

void Logger::logError( const char*, ... ) {}
void Logger::logWarning( const char*, ... ) {}

using namespace conway::geometry;

namespace {

int failures = 0;

void check( bool condition, const std::string& what ) {

  printf( "  %s  %s\n", condition ? "ok  " : "FAIL", what.c_str() );

  if ( !condition ) ++failures;
}

constexpr double PERIOD = 2.0;

tinynurbs::RationalSurface3d buildClosedSurface() {

  tinynurbs::RationalSurface3d surface;

  surface.degree_u = 1;
  surface.degree_v = 1;

  surface.knots_u = { -1.0, -1.0, -0.4, -0.4, 0.3, 1.0, 1.0 };
  surface.knots_v = { 0.0, 0.0, 1.0, 1.0 };

  const double across[ 5 ] = { -1.0, -0.4, -0.4, 0.3, 1.0 };
  const double height[ 5 ] = { 0.0, 0.0, 0.9, 0.0, 0.0 };

  std::vector< glm::dvec3 > points;
  std::vector< double >     weights;

  for ( uint32_t a = 0; a < 5; ++a ) {

    for ( uint32_t b = 0; b < 2; ++b ) {

      // CLOSED in u, which is what a periodic chart needs.
      const uint32_t row = ( a == 4 ) ? 0 : a;

      points.push_back(
        glm::dvec3( across[ a ], static_cast< double >( b ),
                    height[ row ] ) );

      weights.push_back( 1.0 );
    }
  }

  surface.control_points = tinynurbs::array2( 5, 2, points );
  surface.weights        = tinynurbs::array2( 5, 2, weights );

  return surface;
}

CertificateOutcome ask( const RationalSurfaceEvaluator& evaluator,
                        bool                            periodic,
                        const glm::dvec2&               uv0,
                        const glm::dvec2&               uv1,
                        double&                         bound ) {

  const auto wrap = [ & ]( double u ) {

    if ( !periodic ) return u;

    const double offset = std::fmod( u - ( -1.0 ), PERIOD );

    return -1.0 + ( offset < 0.0 ? offset + PERIOD : offset );
  };

  const glm::dvec3 at0 = evaluator.point( wrap( uv0.x ), uv0.y );
  const glm::dvec3 at1 = evaluator.point( wrap( uv1.x ), uv1.y );

  const NurbsDeflectionCertificate certificate(
    evaluator, periodic, -1.0, PERIOD, 1.0e-8 );

  bound = 0.0;

  return certificate.bound( uv0, uv1, at0, at1, bound );
}

}  // namespace

int main() {

  printf( "certificateFlagGatesThePeriodicPath "
          "(CERTIFICATE_CERTIFIES_PERIODIC_CHARTS=%d)\n",
          CERTIFICATE_CERTIFIES_PERIODIC_CHARTS );

  const tinynurbs::RationalSurface3d surface = buildClosedSurface();

  const RationalSurfaceEvaluator evaluator( surface );

  double periodicBound    = 0.0;
  double nonPeriodicBound = 0.0;

  // A chord that crosses the chart cut, so it can only be answered by the
  // periodic path.
  const CertificateOutcome periodic =
    ask( evaluator, true, glm::dvec2( 0.6, 0.5 ), glm::dvec2( 1.4, 0.5 ),
         periodicBound );

  // The same surface, the same span structure, no periodic chart.
  const CertificateOutcome plain =
    ask( evaluator, false, glm::dvec2( -0.8, 0.5 ), glm::dvec2( 0.2, 0.5 ),
         nonPeriodicBound );

#if CERTIFICATE_CERTIFIES_PERIODIC_CHARTS

  check( periodic != CertificateOutcome::Unsupported,
         "with the periodic path ON, a periodic chord reaches the "
         "certificate" );

  check( periodic == CertificateOutcome::Certified && periodicBound > 0.0,
         "and is bounded rather than declined (" +
           std::to_string( periodicBound ) + ")" );

#else

  check( periodic == CertificateOutcome::Unsupported,
         "with the periodic path OFF, a periodic chord is handed back to the "
         "sampled test" );

  check( periodicBound == 0.0,
         "and no bound is produced for it (" +
           std::to_string( periodicBound ) + ")" );

#endif

  // THE GATE IS NARROW. Without this, disabling the certificate outright
  // would satisfy the half above.
  check( plain == CertificateOutcome::Certified && nonPeriodicBound > 0.0,
         "a NON-periodic chord on the same surface is certified either way (" +
           std::to_string( nonPeriodicBound ) + ")" );

  printf( "%s\n", failures == 0 ? "PASS" : "FAIL" );

  return failures == 0 ? 0 : 1;
}
