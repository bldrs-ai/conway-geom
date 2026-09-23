#pragma once

/*
 * A CERTIFICATE for the chord-deflection test in `tesselate`'s
 * ParameterVertex overload - an UPPER BOUND on how far the surface departs
 * from a mesh edge's chord, over the WHOLE edge, rather than a reading at a
 * few sampled parameters.
 *
 * WHY A BOUND AND NOT SAMPLES.
 *
 * The refinement loop asks one question per edge: "is the surface within
 * `sqrt( minimumDeflection )` of this chord?". Sampling answers a different
 * question - "is it within tolerance AT THESE POINTS" - and the two come
 * apart inside the exact function class this code path produces. Restricted
 * to one knot box and to a straight line in uv, a tensor-product B-spline of
 * bidegree ( dU, dV ) is a polynomial of degree at most dU + dV in the line
 * parameter; the deviation from the chord is such a polynomial vanishing at
 * t = 0 and t = 1. Sampling at k interior points forces k more roots and
 * leaves
 *
 *     d( t ) = A * t * ( t - 1 ) * prod_i ( t - s_i )
 *
 * with A UNBOUNDED. For the bicubic case ( degree 6 ) three interior samples
 * leave a whole one-parameter family the test cannot see at all. See
 * `certificateCatchesWhatSamplingMisses` in the refinement tests, which
 * builds exactly that surface and shows the sampled reading accepting a chord
 * the bound refuses.
 *
 * WHAT IT RESTS ON.
 *
 * Write the deviation on a sub-interval in the Bernstein basis of degree n,
 *
 *     d( t ) = sum_j b_j * B_{ j, n }( t ),   B_{ j, n } >= 0,
 *                                             sum_j B_{ j, n }( t ) = 1.
 *
 * Non-negativity and the partition of unity give, by the triangle inequality,
 *
 *     || d( t ) || = || sum_j b_j B_j( t ) ||
 *                 <= sum_j B_j( t ) || b_j ||
 *                 <= ( max_j || b_j || ) * sum_j B_j( t )
 *                  = max_j || b_j ||.
 *
 * THAT INEQUALITY IS THE WHOLE CERTIFICATE. It holds for any norm, needs no
 * convexity argument beyond the two stated facts, and it is what lets a
 * finite set of coefficients bound a continuum of parameters. The price is
 * that it is an OVER-estimate: the hull of the coefficients contains the
 * curve strictly, so a chord can be refused whose true deviation is inside
 * tolerance. That direction is the safe one.
 *
 * RATIONAL SURFACES. A NURBS with non-unit weights is not polynomial, so
 * interpolating `d` and reading its coefficients would be bounding the wrong
 * function. This handles it in the homogeneous form instead. With
 * S = A / w and the chord C( t ) = ( 1 - t ) S0 + t S1,
 *
 *     d( t ) = S( P( t ) ) - C( t ) = N( t ) / W( t ),
 *     N( t ) = A( P( t ) ) - W( t ) C( t ),      deg N <= dU + dV + 1,
 *     W( t ) = w( P( t ) ),                      deg W <= dU + dV,
 *
 * both POLYNOMIAL. Bounding each by its own coefficients,
 *
 *     max_t || d( t ) || <= ( max_j || n_j || ) / ( min_j W_j ),
 *
 * valid whenever `min_j W_j > 0`, because W( t ) >= min_j W_j by the same
 * partition-of-unity argument applied to a scalar. A weight hull that reaches
 * zero makes the quotient unbounded, and that case is REFUSED rather than
 * guessed - see `Outcome::Inconclusive`.
 *
 * WHAT IS NOT COVERED, stated here so it is not mistaken for coverage:
 *
 *   - only EDGES are certified, not triangle interiors. A triangle whose
 *     three edges are each within tolerance can still bow in its middle.
 *   - BORDER edges are never subdivided by `tesselate` at all, so a trim
 *     boundary's chords carry no certificate.
 *   - the bound is about the NURBS this object holds. That it is the same
 *     surface the `tesselate` callback evaluates is a property of the call
 *     site, not something checked here - though a mismatch shows up in the
 *     bound, since the Bernstein endpoint coefficients ARE the endpoint
 *     deviations and a mismatched surface makes them non-zero.
 *   - floating point. The inflation term below is scaled by the interpolation
 *     matrix's infinity norm and an ulp estimate of the evaluation error; it
 *     is an ESTIMATE with a safety factor, not an interval-arithmetic proof.
 */

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "conway_geometry/operations/nurbs_utils.h"

namespace conway::geometry {

/**
 * Highest total degree ( dU + dV, plus one for a rational surface ) this will
 * certify. Above it the node->Bernstein matrix's infinity norm grows fast
 * enough that the inflation term below stops being small compared with a
 * tessellation tolerance, and the honest answer is "cannot certify" rather
 * than a number nobody should trust. Every B-spline surface in the STEP
 * corpus measured here is bicubic or lower, i.e. n <= 7.
 */
constexpr uint32_t CERTIFICATE_MAX_DEGREE = 13;

/**
 * Most knot boxes a single chord may pass through and still be certified.
 * The deviation is only polynomial WITHIN a box, so a chord crossing knot
 * lines is bounded box by box; this caps that work.
 *
 * Refusing above the cap is conservative - the edge is subdivided - and it
 * terminates, because halving an edge halves the boxes it crosses, so a chord
 * over the cap is under it within log2( boxes / cap ) levels.
 *
 * 16 rather than the 8 the cut-based arrangement used, because the two count
 * DIFFERENT THINGS and 8 is not the same budget here. That one capped the
 * cut parameters collected per axis, in two nested passes, so a chord could
 * reach roughly 8 x 8 boxes before being refused; this caps the boxes
 * themselves, once. Measured on `Right_Hand.step`: at 16 the walk refuses
 * NOTHING on the whole corpus and the model comes out at 201,223 triangles,
 * against 201,311 for the arrangement this replaces - the 88 chords that one
 * refused at its cap are now certified. At 8 it refuses enough to cost 200
 * triangles back, at 201,423. Degenerate and welded-open counts are
 * identical at all three.
 */
constexpr uint32_t CERTIFICATE_MAX_PIECES = 16;

/**
 * Multiple of the interpolation error the bound is inflated by before it is
 * compared with the tolerance. 8 rather than 1 because the per-node error is
 * itself an estimate ( a few ulps of the surface magnitude, dominated by the
 * cancellation in `surfacePoint - chord` ), not a proved envelope.
 */
constexpr double CERTIFICATE_ERROR_SAFETY = 8.0;

/**
 * If the inflation reaches this fraction of the tolerance the certificate
 * declines rather than returning a bound whose own error is comparable with
 * what it is bounding. Self-policing: it is what makes CERTIFICATE_MAX_DEGREE
 * a performance cut-off rather than the thing soundness hangs on.
 */
constexpr double CERTIFICATE_MAX_ERROR_FRACTION = 0.25;

/**
 * COVERAGE AUDIT OF THE WALK/CALLER SEAM. Off in every shipping build; the
 * differential fuzzer defines it.
 *
 * The walk owns each box's span as an INTEGER and evaluates one-sided through
 * `pointAtSpan`. The call site does not: it calls `point( u, v )`, which
 * resolves the span with `findSpan` on the ROUNDED double. Those two answers
 * are allowed to differ at a box's own bounding knot - the next box owns that
 * double and bounds it against the right polynomial - and are NOT allowed to
 * differ anywhere else.
 *
 * The eleventh finding on bldrs-ai/conway-geom#214 was a box where they
 * differed at EVERY node, because the chord's extent on that axis was below
 * the parameter's resolution and its far end sat on a knot. Nothing watched
 * that seam; it was found by its consequence - a bound below a sampled truth
 * - at fuzz trial 173,559, which is luck rather than instrumentation. These
 * counters watch it directly.
 */
/**
 * WHETHER THE CERTIFICATE IS CONSULTED ON A PERIODIC CHART. Off.
 *
 * A periodic chord takes the sampled deflection test on its own - which is
 * what shipped before this file existed - and every non-periodic chord keeps
 * the certificate.
 *
 * WHY IT IS OFF. Not because the periodic path is known broken; it is not,
 * and the whole of it is covered by tests and by the differential fuzzer,
 * which drives it at full strength whatever this constant says. It is off
 * because SEVEN OF THE LAST ELEVEN FINDINGS on bldrs-ai/conway-geom#214 were
 * reachable only through it - #4b, #8, #9, #9a, #10, #12 and #14 - and
 * because each of those was found when an instrument got better rather than
 * when someone reasoned harder. Four rounds of that in a row is not evidence
 * that the family is closed; it is evidence that its size is unknown. This
 * constant is how a change whose risk cannot be sized still ships: the 94% of
 * the work that never touched the family goes out, and the 6% that carried it
 * waits.
 *
 * SIX PER CENT is the measured share, not an estimate. Instrumented builds
 * over the smoke corpus: `Right_Hand.step` makes 110,000 certificate calls,
 * 6,967 of them periodic; `nist_ctc_02_asme1_rc.stp` makes 80,000, none of
 * them periodic; every other model in the corpus makes fewer than 200 calls
 * in total.
 *
 * WHAT WOULD JUSTIFY TURNING IT ON. A round of `test/certificate_fuzz.cpp`
 * that adds no new finding on the periodic path - meaning a round where the
 * GENERATOR was extended and found nothing, not one where it was run again
 * unchanged and found nothing. Both rounds that were declared clean on this
 * PR were clean in the second sense, and both were followed by a round that
 * found more as soon as the generator reached further. The bar is a generator
 * that reaches a shape it did not before and comes back empty.
 *
 * Turning it on is this one line. Nothing else in this file branches on it,
 * and nothing about the periodic code is compiled out - so the tests and the
 * fuzzer keep exercising the path, and the day it is turned on it is not
 * being run for the first time.
 */
#ifndef CERTIFICATE_CERTIFIES_PERIODIC_CHARTS
#define CERTIFICATE_CERTIFIES_PERIODIC_CHARTS 0
#endif

#ifndef CERTIFICATE_AUDIT_CALLER_SPANS
#define CERTIFICATE_AUDIT_CALLER_SPANS 0
#endif



/**
 * Values at n + 1 nodes -> the n + 1 Bernstein coefficients of the degree-n
 * polynomial through them, as a dense matrix, plus that matrix's infinity
 * norm for the error term.
 *
 * CHEBYSHEV-LOBATTO nodes, ( 1 - cos( pi i / n ) ) / 2, not equally spaced
 * ones. Measured infinity norms of the inverse ( scratchpad `cond.py` ):
 *
 *     n     equally spaced      Lobatto
 *     6            89.24          46.20
 *     7           210.23          85.80
 *    10          3648.0          733.16
 *
 * so the error term the bound carries is halved at the bicubic degree this
 * corpus actually uses and better than five times smaller by n = 10. Lobatto
 * also INCLUDES both endpoints, which matters here for a second reason: the
 * deviation is supposed to vanish at t = 0 and t = 1, the endpoint Bernstein
 * coefficients are exactly the endpoint values, so evaluating there costs
 * nothing extra and any disagreement between this object's surface and the
 * caller's lands in the bound instead of hiding under it.
 *
 * Built per `tesselate` call rather than in a shared static: a 14x14 inverse
 * is a few thousand flops against a face's tens of thousands of surface
 * evaluations, and a per-object table has no initialisation race to reason
 * about in the threaded wasm build.
 */
class BernsteinNodes {
 public:

  BernsteinNodes() = default;

  explicit BernsteinNodes( uint32_t degree ) { build( degree ); }

  void build( uint32_t degree ) {

    degree_ = degree;

    const uint32_t count = degree + 1;

    nodes_.resize( count );

    for ( uint32_t i = 0; i < count; ++i ) {
      nodes_[ i ] =
        0.5 * ( 1.0 - std::cos( 3.14159265358979323846 * i / degree ) );
    }

    // Exact at the ends: cos( 0 ) and cos( pi ) are exact, but the arithmetic
    // around them need not land on 0 and 1, and the endpoints being EXACTLY
    // the edge's own endpoints is what makes the endpoint coefficients the
    // endpoint deviations.
    nodes_.front() = 0.0;
    nodes_.back()  = 1.0;

    // Vandermonde in the Bernstein basis, inverted by Gauss-Jordan with
    // partial pivoting.
    std::vector< double > work( count * 2 * count, 0.0 );

    for ( uint32_t i = 0; i < count; ++i ) {

      const double t   = nodes_[ i ];
      const double oneT = 1.0 - t;

      double binomial = 1.0;

      for ( uint32_t j = 0; j < count; ++j ) {

        work[ ( i * 2 * count ) + j ] =
          binomial * std::pow( t, static_cast< double >( j ) ) *
          std::pow( oneT, static_cast< double >( degree - j ) );

        binomial = binomial * ( degree - j ) / ( j + 1 );
      }

      work[ ( i * 2 * count ) + count + i ] = 1.0;
    }

    for ( uint32_t column = 0; column < count; ++column ) {

      uint32_t pivot = column;

      for ( uint32_t row = column + 1; row < count; ++row ) {

        if ( std::abs( work[ ( row * 2 * count ) + column ] ) >
             std::abs( work[ ( pivot * 2 * count ) + column ] ) ) {
          pivot = row;
        }
      }

      for ( uint32_t k = 0; k < 2 * count; ++k ) {
        std::swap( work[ ( column * 2 * count ) + k ],
                   work[ ( pivot * 2 * count ) + k ] );
      }

      const double diagonal = work[ ( column * 2 * count ) + column ];

      if ( diagonal == 0.0 ) {
        degree_ = 0;
        return;
      }

      for ( uint32_t k = 0; k < 2 * count; ++k ) {
        work[ ( column * 2 * count ) + k ] /= diagonal;
      }

      for ( uint32_t row = 0; row < count; ++row ) {

        if ( row == column ) {
          continue;
        }

        const double factor = work[ ( row * 2 * count ) + column ];

        if ( factor == 0.0 ) {
          continue;
        }

        for ( uint32_t k = 0; k < 2 * count; ++k ) {
          work[ ( row * 2 * count ) + k ] -=
            factor * work[ ( column * 2 * count ) + k ];
        }
      }
    }

    inverse_.resize( count * count );

    normInfinity_ = 0.0;

    for ( uint32_t i = 0; i < count; ++i ) {

      double rowSum = 0.0;

      for ( uint32_t j = 0; j < count; ++j ) {

        const double value = work[ ( i * 2 * count ) + count + j ];

        inverse_[ ( i * count ) + j ] = value;

        rowSum += std::abs( value );
      }

      normInfinity_ = std::max( normInfinity_, rowSum );
    }

    measureInverseError( degree );
  }

  uint32_t degree() const { return degree_; }

  double node( uint32_t i ) const { return nodes_[ i ]; }

  double normInfinity() const { return normInfinity_; }

  /**
   * Bound on how far the STORED inverse is from the exact one, in the
   * infinity norm.
   *
   * `normInfinity()` amplifies perturbations of exact node values; it says
   * nothing about the inverse itself being inexact, and the Gauss-Jordan
   * above is ordinary floating point. Measured a posteriori rather than
   * estimated: with `R = I - M V`, `M = ( I - R ) V^-1` exactly, so
   * `V^-1 - M = R V^-1` and `|| V^-1 - M || <= ||R|| ||M|| / ( 1 - ||R|| )`.
   * The residual is itself computed in floating point, so it is widened by
   * the rounding of forming `V` and of the product before being used.
   *
   * Measured values, degree 2 to 15 ( scratchpad /tmp/err.cpp ): 1.7e-16 at
   * n = 2, 1.2e-13 at n = 6, 3.5e-13 at n = 7, 1.9e-9 at n = 13. Small
   * everywhere the degree cap allows, which is the point - it is carried so
   * the bound does not REST on it being small.
   */
  double inverseError() const { return inverseError_; }

  /** True when the table is usable at all. */
  bool valid() const { return degree_ > 0 && inverseError_ < 1.0; }

  /** Bernstein coefficient `i` from the values at the nodes. */
  template< typename Value >
  Value coefficient( uint32_t i, const Value* values ) const {

    const uint32_t count = degree_ + 1;

    Value result = values[ 0 ] * inverse_[ i * count ];

    for ( uint32_t j = 1; j < count; ++j ) {
      result += values[ j ] * inverse_[ ( i * count ) + j ];
    }

    return result;
  }

 private:

  /** See inverseError(). */
  void measureInverseError( uint32_t degree ) {

    const uint32_t count = degree + 1;

    double residual = 0.0;

    for ( uint32_t i = 0; i < count; ++i ) {

      double rowSum = 0.0;

      for ( uint32_t j = 0; j < count; ++j ) {

        double product = 0.0;

        for ( uint32_t k = 0; k < count; ++k ) {

          // V[ k ][ j ] = B_{ j, degree }( node_k ), rebuilt here rather than
          // kept, so the residual is not read off the same array the
          // elimination consumed.
          double binomial = 1.0;

          for ( uint32_t b = 0; b < j; ++b ) {
            binomial = binomial * ( degree - b ) / ( b + 1 );
          }

          const double t = nodes_[ k ];

          product +=
            inverse_[ ( i * count ) + k ] * binomial *
            std::pow( t, static_cast< double >( j ) ) *
            std::pow( 1.0 - t, static_cast< double >( degree - j ) );
        }

        rowSum += std::abs( product - ( i == j ? 1.0 : 0.0 ) );
      }

      residual = std::max( residual, rowSum );
    }

    // Widen by the rounding that forming the residual itself carries: the
    // Bernstein values are a few ulps off, and the row product is a sum of
    // `count` terms. Both are amplified by the inverse's own norm.
    const double slack =
      8.0 * count * std::numeric_limits< double >::epsilon() * normInfinity_;

    residual += slack;

    inverseError_ =
      residual < 0.5 ?
        ( normInfinity_ * residual / ( 1.0 - residual ) ) :
        std::numeric_limits< double >::infinity();
  }

  uint32_t              degree_        = 0;
  std::vector< double > nodes_;
  std::vector< double > inverse_;
  double                normInfinity_  = 0.0;
  double                inverseError_  = std::numeric_limits< double >::infinity();
};

/** What a certificate attempt concluded. */
enum class CertificateOutcome : uint8_t {

  /** `bound` holds an upper bound on the deviation over the whole edge. */
  Certified,

  /**
   * No bound could be established for THIS chord - it wrapped a periodic
   * chart, crossed too many knot spans, met a weight hull reaching zero, or
   * hit non-finite arithmetic. The caller must subdivide: the answer is
   * unknown, not small.
   */
  Inconclusive,

  /**
   * This surface is not one that can be certified at all ( not a NURBS, or
   * above the degree cap ). The caller keeps its sampled reading; the
   * guarantee does not extend to this surface.
   */
  Unsupported
};

/**
 * The null certificate - what a `tesselate` call that has no NURBS behind it
 * ( the cone and cylinder projections, which are functions of the incoming
 * POSITION and not of uv alone, so there is no `S( uv )` to bound ) gets by
 * default. It reports `Unsupported`, never `Inconclusive`, and that
 * distinction is load bearing: `Inconclusive` means subdivide, and a surface
 * that can never be certified would then subdivide until the triangle budget
 * ran out on every face.
 */
struct NoDeflectionCertificate {

  static constexpr bool CERTIFIES = false;

  CertificateOutcome bound(
    const glm::dvec2&,
    const glm::dvec2&,
    const glm::dvec3&,
    const glm::dvec3&,
    double& ) const {

    return CertificateOutcome::Unsupported;
  }
};

/**
 * Certificate for a rational B-spline surface evaluated as
 * `evaluator.point( wrap( u ), v )`.
 *
 * Holds the wrap's parameters rather than the wrapping lambda so it can ask
 * the question the lambda cannot answer: whether the chord CROSSES the wrap,
 * where t -> u is discontinuous and no polynomial argument survives.
 */
class NurbsDeflectionCertificate {
 public:

  NurbsDeflectionCertificate(
    const RationalSurfaceEvaluator& evaluator,
    bool                            periodic,
    double                          stripUMin,
    double                          stripPeriod,
    double                          toleranceSquared )
    : evaluator_( &evaluator ),
      periodic_( periodic ),
      stripUMin_( stripUMin ),
      stripPeriod_( stripPeriod ),
      tolerance_( std::sqrt( toleranceSquared ) ) {

    rational_ = !evaluator.isPolynomial();

    const uint32_t totalDegree =
      evaluator.degreeU() + evaluator.degreeV() + ( rational_ ? 1u : 0u );

    // A chart this build does not certify reports UNSUPPORTED, which the
    // call site already handles by taking the sampled reading as the whole
    // test - see CERTIFICATE_CERTIFIES_PERIODIC_CHARTS, and the
    // `case CertificateOutcome::Unsupported` arm in tesselation_utils.h. It
    // is deliberately the same answer as a cone or a cylinder, which have no
    // parameterisation to bound: "outside the guarantee", not "in tolerance".
    supported_ =
      evaluator.supportsFastPath() &&
      totalDegree >= 1 &&
      totalDegree <= CERTIFICATE_MAX_DEGREE &&
      ( !periodic || stripPeriod > 0.0 ) &&
      ( !periodic || CERTIFICATE_CERTIFIES_PERIODIC_CHARTS );

    if ( supported_ ) {

      interpolation_.build( totalDegree );

      supported_ =
        interpolation_.degree() == totalDegree && interpolation_.valid();

      // Does either axis carry an interior knot of multiplicity
      // degree + 1 - the multiplicity at which the surface STEPS? Only on
      // such an axis can a node nudged across a span boundary by parameter
      // rounding read a value that is not near the one intended, so only
      // there does boundPiece have to refuse a box it cannot separate from
      // its neighbour. Computed once, because it is a property of the knot
      // vector and not of any chord.
      const auto stepsAnywhere =
        []( uint32_t degree, const std::vector< double >& knots ) {

          for ( size_t at = degree + 1; at + degree + 1 < knots.size(); ) {

            size_t run = at;

            while ( run + 1 < knots.size() && knots[ run + 1 ] == knots[ at ] ) {
              ++run;
            }

            if ( ( run - at + 1 ) > degree ) {
              return true;
            }

            at = run + 1;
          }

          return false;
        };

      stepsU_ = stepsAnywhere( evaluator.degreeU(), evaluator.knotsU() );
      stepsV_ = stepsAnywhere( evaluator.degreeV(), evaluator.knotsV() );
    }
  }

  static constexpr bool CERTIFIES = true;

  /**
   * How many spans a chord crossed, how many chords were bounded and how many
   * were declined - read by the cost probe, and nothing else depends on it.
   */
  struct Counters {
    uint64_t certified    = 0;
    uint64_t inconclusive = 0;
    uint64_t evaluations  = 0;
    uint64_t spans        = 0;
    uint64_t wrapped      = 0;
    uint64_t tooManySpans = 0;
    uint64_t unplaceable  = 0;
    uint64_t ambiguous    = 0;
    uint64_t extrapolated = 0;

    /**
     * Branch counters, for the COVERAGE audit rather than for cost. A
     * certificate defect has twice now lived in a branch no test executed,
     * so the fuzzer reports which of these it has reached and the ones
     * still at zero are the list of what is not being exercised.
     */
    uint64_t strips       = 0;
    uint64_t clampsClear  = 0;
    uint64_t descending   = 0;
    uint64_t rationalPath = 0;
    uint64_t illConditioned = 0;
    uint64_t badPiece     = 0;

    /**
     * THE WALK/CALLER SEAM, under CERTIFICATE_AUDIT_CALLER_SPANS.
     *
     *   callerSpanAgree   every node of the box resolves, through the
     *                     caller's own `findSpan`, to the span the box
     *                     carries;
     *   callerSpanEdge    the only nodes that do not are sitting exactly on
     *                     a knot that BOUNDS the box - benign, because the
     *                     next box owns that double and bounds it against
     *                     the polynomial the caller will use;
     *   callerSpanInside  a node strictly inside the box's knot interval
     *                     resolves elsewhere. THIS MUST STAY ZERO; a
     *                     non-zero here means the walk put a box on a span
     *                     that does not contain it;
     *   callerSpanWhole   EVERY node resolves elsewhere, so the caller
     *                     evaluates the entire box on another polynomial.
     *                     This is the eleventh finding's signature.
     *   callerSpanPair    one of the two boxes emitted for an undecidable
     *                     cross-axis ordering. Excluded from the others
     *                     because their corners are the chord's rather than
     *                     a knot's; between them they cover both spans the
     *                     caller can resolve to.
     */
    uint64_t callerSpanAgree  = 0;
    uint64_t callerSpanEdge   = 0;
    uint64_t callerSpanInside = 0;
    uint64_t callerSpanWhole  = 0;
    uint64_t callerSpanPair   = 0;
    uint64_t callerSpanWholeStep = 0;
  };

  const Counters& counters() const { return counters_; }

  /** False when no chord on this surface can be certified at all. */
  bool supported() const { return supported_; }

  /**
   * Upper bound on `max_t || S( P( t ) ) - ( ( 1 - t ) S0 + t S1 ) ||` for
   * the straight uv segment P from `uv0` to `uv1`.
   *
   * `S0` and `S1` are the caller's own surface positions at the endpoints -
   * the `surfaceAt` cache - not re-evaluated here, so the bound is against
   * the chord the mesh will actually carry.
   */
  CertificateOutcome bound(
    const glm::dvec2& uv0,
    const glm::dvec2& uv1,
    const glm::dvec3& surface0,
    const glm::dvec3& surface1,
    double&           result ) const {

    if ( !supported_ ) {
      return CertificateOutcome::Unsupported;
    }

    // The chord's own magnitude, which is one of the things the error term
    // is built from. It is NOT on its own a scale for the evaluation error -
    // see maxHomogeneousNorm() in nurbs_utils.h, and the ERROR PROPAGATION
    // note on boundPiece.
    const double chordScale =
      std::max( glm::length( surface0 ), glm::length( surface1 ) );

    // THE PARAMETER SCALE THE CHORD IS READ AT, per axis.
    //
    // A box's corners are `uv0 + t * d` - or a knot - and the rounding in
    // that expression is at the magnitude of `uv0` and `d`, NOT at the
    // magnitude of the box. A box sitting at u = -350 on a chord that runs
    // to u = 67,727 carries the larger number's rounding, and pricing it at
    // the smaller one puts the box's own corner outside its own span by more
    // than the error term admits. Found by the CALLER-SPAN AUDIT rather than
    // by a violation: `callerSpanInside` went to 1 in 800,000 trials, which
    // is a box placed on a span that does not contain it.
    //
    // This is the first finding's shape a fourth time - an error term scaled
    // by the wrong quantity - and the pattern is worth stating plainly: every
    // rounding in this file has to be priced at the magnitude of the
    // EXPRESSION THAT PRODUCED the number, not of the number.
    const double parameterScaleU =
      std::max( std::abs( uv0.x ), std::abs( uv1.x ) ) +
      std::abs( uv1.x - uv0.x );

    const double parameterScaleV =
      std::max( std::abs( uv0.y ), std::abs( uv1.y ) ) +
      std::abs( uv1.y - uv0.y );

    Piece pieces[ CERTIFICATE_MAX_PIECES ];

    uint32_t pieceCount = 0;

    if ( !walkPieces( uv0, uv1, pieces, pieceCount ) ) {
      ++counters_.tooManySpans;
      ++counters_.inconclusive;
      return CertificateOutcome::Inconclusive;
    }

    double worst      = 0.0;
    double worstError = 0.0;

    for ( uint32_t at = 0; at < pieceCount; ++at ) {

      ++counters_.spans;

      double pieceBound = 0.0;
      double pieceError = 0.0;

      if ( !boundPiece(
             pieces[ at ], surface0, surface1, chordScale,
             parameterScaleU, parameterScaleV,
             pieceBound, pieceError ) ) {

        ++counters_.badPiece;
        ++counters_.inconclusive;
        return CertificateOutcome::Inconclusive;
      }

      worst      = std::max( worst, pieceBound );
      worstError = std::max( worstError, pieceError );
    }

    // THE BOUND'S OWN ERROR HAS TO BE SMALL COMPARED WITH WHAT IT IS
    // BOUNDING, or it is not a bound anybody should act on. `pieceBound`
    // already carries its error outward; this declines the cases where that
    // outward widening is itself a large fraction of the target, which is
    // what keeps CERTIFICATE_MAX_DEGREE a cost cut-off rather than the thing
    // soundness hangs on.
    if ( worstError >= tolerance_ * CERTIFICATE_MAX_ERROR_FRACTION ) {
      ++counters_.illConditioned;
      ++counters_.inconclusive;
      return CertificateOutcome::Inconclusive;
    }

    result = worst;

    ++counters_.certified;

    return CertificateOutcome::Certified;
  }

 private:

  /**
   * ONE KNOT BOX THE CHORD PASSES THROUGH, named by its span INDICES.
   *
   * The indices are the whole point. An earlier arrangement recovered the
   * pieces by computing a cut parameter for every knot the chord crossed,
   * sorting those parameters and stepping between consecutive ones - so a
   * piece existed only if its two bounding parameters came out distinct, and
   * the polynomial to certify it against was recovered afterwards by looking
   * a parameter back up with `findSpan`. Both steps put exact integer
   * structure - WHICH KNOT SPAN IS THIS - through a lossy round trip into
   * doubles, and both produced a P1 on bldrs-ai/conway-geom#214:
   *
   *   - a node parameter sitting exactly on a knot resolved to the span on
   *     its RIGHT, so a piece to the left of a discontinuous knot was
   *     certified against the polynomial it was not covering;
   *   - two distinct knots that rounded to one parameter produced one cut
   *     instead of two, the zero-length step between them was skipped, and
   *     the span between them was never certified at all while the result
   *     still claimed to cover the whole chord.
   *
   * Here a piece exists because the walk emitted it, and carries the span it
   * belongs to as an integer it was emitted with. `tFrom`/`tTo` are derived
   * FROM the boxes and are used only to place nodes and to read the chord;
   * nothing is recovered from them.
   *
   * `from` and `to` are in EVALUATED uv - the periodic shift is already
   * removed - and the axis that ended the piece is pinned to the exact knot
   * value that ended it, rather than recomputed from a parameter.
   */
  struct Piece {

    int        spanU = 0;
    int        spanV = 0;

    /**
     * The periodic shift already removed from `from` and `to`, kept only so
     * the error term can price it - see ROUNDING IS GOVERNED BY THE RAW
     * PARAMETER in boundPiece.
     */
    double     shift = 0.0;
    glm::dvec2 from  = glm::dvec2( 0.0 );
    glm::dvec2 to    = glm::dvec2( 0.0 );
    double     tFrom = 0.0;
    double     tTo   = 0.0;

    /**
     * True for the two boxes emitted when a cross-axis ordering is
     * undecidable. Their corners are the CHORD's, not a knot's, because
     * which knot bounds the stretch is exactly what is not known - so unlike
     * every other box they are not expected to sit inside their own span's
     * knot interval. Only the caller-span audit reads this.
     */
    bool       covering = false;
  };

  /**
   * The span the walk should START in at `at`, moving in the direction
   * `ascending`, skipping spans of zero width.
   *
   * THIS IS THE ONLY SPAN LOOKUP IN THE CERTIFICATE. Every other span comes
   * from incrementing this one, which is what makes "certified against the
   * wrong polynomial" a question that can only be asked once per chord
   * instead of once per node. `findSpan` resolves a parameter sitting
   * exactly on a knot to the span on its right - the span the walk wants
   * when it is moving up, and the wrong one when it is moving down - so the
   * direction is applied here explicitly.
   *
   * THE DIRECTION IS AN OPTIMISATION, NOT A SOUNDNESS REQUIREMENT, which is
   * worth writing down because it is the one lookup that could still answer
   * wrongly. If it returns the span on the wrong side, that span's exit
   * boundary lies AT or BEHIND the chord's start, so the walk's first box
   * comes out with zero extent - contributing a deviation of zero, since its
   * single point is the cached endpoint itself - and the very next step
   * lands in the right span. Measured: the `walk-ignores-direction` red
   * proof removes this block entirely and every assertion still passes,
   * including `descendingChordStartsInTheSpanItTraverses`, which is built
   * to catch exactly this. What the block buys is not correctness but one
   * fewer wasted box per descending chord.
   *
   * Returns -1 when the direction leaves the knot domain with no non-empty
   * span in it.
   */
  static int directionalSpan(
    uint32_t                     degree,
    const std::vector< double >& knots,
    double                       at,
    bool                         ascending ) {

    const int first = static_cast< int >( degree );
    const int last  =
      static_cast< int >( knots.size() ) - static_cast< int >( degree ) - 2;

    if ( last < first ) {
      return -1;
    }

    int span = RationalSurfaceEvaluator::findSpan( degree, knots, at );

    if ( !ascending ) {

      // Descending from exactly a knot: the span the chord is about to
      // traverse is the one to the LEFT, and a repeated knot means stepping
      // back over every copy of it.
      while ( span > first && knots[ span ] >= at ) {
        --span;
      }
    }

    while ( span >= first && span <= last &&
            !RationalSurfaceEvaluator::spanIsEvaluable( knots, span ) ) {
      span += ascending ? 1 : -1;
    }

    return ( span < first || span > last ) ? -1 : span;
  }

  /**
   * Step one span index in the direction of travel, skipping spans of zero
   * width, and report when the walk has run off the end of the knot domain.
   */
  static int advanceSpan(
    int                          span,
    bool                         ascending,
    uint32_t                     degree,
    const std::vector< double >& knots,
    bool&                        clamped ) {

    const int first = static_cast< int >( degree );
    const int last  =
      static_cast< int >( knots.size() ) - static_cast< int >( degree ) - 2;

    do {
      span += ascending ? 1 : -1;
    } while ( span >= first && span <= last &&
              !RationalSurfaceEvaluator::spanIsEvaluable( knots, span ) );

    if ( span < first || span > last ) {
      span    = ( span < first ) ? first : last;
      clamped = true;
    }

    return span;
  }

  /**
   * Enumerate the knot boxes the chord passes through, in order along it.
   *
   * The walk advances in INDEX space: at each step it asks which axis's box
   * ends first along the chord, emits the box it is in, steps that axis's
   * index, and repeats. Progress is therefore guaranteed by the indices
   * advancing, not by the chord parameter increasing - which is what lets a
   * box whose parameter interval collapses still be emitted and bounded
   * rather than vanishing.
   *
   * The periodic chart's cut is folded in as a RESTART of the u walk: when
   * the chord reaches the edge of the strip, the shift changes and the u
   * index resumes at the far end of the knot domain. It is not a separate
   * kind of cut, so there is no second source of boundaries that could
   * collide with the first.
   */
  bool walkPieces(
    const glm::dvec2& uv0,
    const glm::dvec2& uv1,
    Piece*            pieces,
    uint32_t&         count ) const {

    const std::vector< double >& knotsU = evaluator_->knotsU();
    const std::vector< double >& knotsV = evaluator_->knotsV();

    const uint32_t degreeU = evaluator_->degreeU();
    const uint32_t degreeV = evaluator_->degreeV();

    const int firstU = static_cast< int >( degreeU );
    const int firstV = static_cast< int >( degreeV );
    const int lastU  =
      static_cast< int >( knotsU.size() ) - static_cast< int >( degreeU ) - 2;
    const int lastV  =
      static_cast< int >( knotsV.size() ) - static_cast< int >( degreeV ) - 2;

    const double du = uv1.x - uv0.x;
    const double dv = uv1.y - uv0.y;

    const bool upU = du >= 0.0;
    const bool upV = dv >= 0.0;

    if ( !upU || !upV ) {
      ++counters_.descending;
    }

    if ( rational_ ) {
      ++counters_.rationalPath;
    }

    // THE WALK STARTS WHERE THE CALLBACK STARTS, BIT FOR BIT.
    //
    // The call site evaluates `point( wrapChartU( u ), v )`, and
    // `wrapChartU` is `uMin + fmod( u - uMin, P )` - one exact IEEE `fmod`.
    // This used to reduce with `u - floor( ( u - uMin ) / P ) * P`, which is
    // a divide, a floor, a multiply and a subtract, each rounded. The two
    // agree on most inputs and DISAGREE BY A WHOLE PERIOD near a sheet
    // boundary, because the rounded quotient crosses an integer where the
    // exact one does not.
    //
    // Measured on fuzz seed 3 at trial 32,525
    // (`periodicStartIsReducedTheWayTheCallbackReducesIt`): a chord at
    // u = 79627595.655702367 on a strip of period 114.74700735250697. The
    // `floor` form gives -57.3735036700963974, six nanometres above the
    // BOTTOM of the strip; `fmod` gives 57.373503674944139163, one nanometre
    // below the TOP. On a chord with du = 0 - the walk never moves in u -
    // every box is then bounded at a u the callback never evaluates.
    // Bound 0.003925974707 against a true departure of 0.00470532459.
    //
    // Reducing with `fmod` and DERIVING the shift from it, rather than the
    // other way round, makes the disagreement unrepresentable: `from.x` is
    // by construction the number the callback will pass to the evaluator.
    // `shift` keeps its meaning - raw = evaluated + shift - and is used only
    // to place boundaries and to price rounding, where one ulp is already
    // carried.
    glm::dvec2 from( uv0.x, uv0.y );

    if ( periodic_ ) {

      const double offset = std::fmod( uv0.x - stripUMin_, stripPeriod_ );

      from.x =
        stripUMin_ + ( offset < 0.0 ? offset + stripPeriod_ : offset );
    }

    double shift = uv0.x - from.x;

    int spanU = directionalSpan( degreeU, knotsU, from.x, upU );
    int spanV = directionalSpan( degreeV, knotsV, from.y, upV );

    if ( spanU < 0 || spanV < 0 ) {
      return false;
    }

    const double infinity = std::numeric_limits< double >::infinity();

    // Once the walk runs off the end of a knot domain there are no further
    // edges on that axis - the surface clamps there and the last span's
    // polynomial is extended, which is what the evaluator does too. Without
    // these the walk would keep finding the same domain edge and emit
    // zero-width boxes until it hit the cap.
    bool clampedU =
      ( du == 0.0 ) ||
      ( upU ? ( from.x >= knotsU[ lastU + 1 ] ) : ( from.x <= knotsU[ firstU ] ) );

    bool clampedV =
      ( dv == 0.0 ) ||
      ( upV ? ( from.y >= knotsV[ lastV + 1 ] ) : ( from.y <= knotsV[ firstV ] ) );

    double tFrom = 0.0;

    count = 0;

    while ( true ) {

      if ( count >= CERTIFICATE_MAX_PIECES ) {
        return false;
      }

      // Where this box ends on each axis, in that axis's own parameter.
      const double edgeU = upU ? knotsU[ spanU + 1 ] : knotsU[ spanU ];
      const double edgeV = upV ? knotsV[ spanV + 1 ] : knotsV[ spanV ];

      // And where the periodic strip ends - in EVALUATED u, so it can be
      // compared with `edgeU` without either of them going through a
      // division first. See WHICH BOUNDARY COMES FIRST below.
      const double edgeStripLocal =
        upU ? ( stripUMin_ + stripPeriod_ ) : stripUMin_;

      const double edgeStrip = edgeStripLocal + shift;

      const double tU =
        clampedU ? infinity : ( ( ( edgeU + shift ) - uv0.x ) / du );

      const double tV =
        clampedV ? infinity : ( ( edgeV - uv0.y ) / dv );

      const double tStrip =
        ( periodic_ && du != 0.0 ) ?
          ( ( edgeStrip - uv0.x ) / du ) : infinity;

      double tTo = 1.0;

      tTo = std::min( tTo, std::min( tU, std::min( tV, tStrip ) ) );

      // Never run backwards: a knot whose parameter rounds below where the
      // walk already is gives a box of zero parameter width, which is
      // emitted and bounded like any other - its uv interval comes from the
      // knots and is not degenerate.
      tTo = std::max( tTo, tFrom );

      // WHEN TWO AXES TIE, THEIR ORDER IS NOT DECIDABLE, SO COVER BOTH.
      //
      // `tU` and `tV` are computed from exact knots but carried in doubles,
      // so a difference below their own rounding says nothing about which
      // boundary the chord really reaches first. Picking one - this used to
      // always pick u - walks `( oldU, oldV ) -> ( newU, oldV ) ->
      // ( newU, newV )` and never bounds `( oldU, newV )`. With full
      // multiplicity on both axes that skipped patch is an INDEPENDENT
      // polynomial, not a near-point sliver, and it can carry any departure
      // at all while every emitted patch reads zero. Measured on the
      // reviewer's case: a bound of 4.8e-5 against a true departure of 7.
      //
      // That was the fifth finding on bldrs-ai/conway-geom#214, and it is
      // the one this file shipped on an ARGUMENT - that such a box had zero
      // extent on both axes so no departure could develop in it. The
      // argument named its own failure condition, full multiplicity in both
      // axes, and that is the condition the reviewer used.
      //
      // Both intermediate boxes are emitted instead of declining. Each t in
      // the ambiguous interval genuinely belongs to one of them, and the
      // bound is the maximum over boxes, so covering the interval twice is
      // sound and costs two boxes where declining would cost a subdivision.
      const double separation =
        8.0 * std::numeric_limits< double >::epsilon() *
        std::max( std::max( std::abs( tU ), std::abs( tV ) ), 1.0 );

      // THE STRIP BOUNDARY OUTRANKS A KNOT TIE, and the order of these two
      // declarations is the whole of that rule - it used to read
      // `endsOnStrip = !ambiguous && ...`, which is how the ninth finding on
      // bldrs-ai/conway-geom#214 got in.
      //
      // A tie between `tU` and `tV` is an ordering question INSIDE one sheet:
      // both orders are covered, both boxes are bounded, and the walk carries
      // on with the same `shift`. A strip crossing is not an ordering
      // question at all - it changes the frame the walk is in. Suppressing it
      // because two knot boundaries happened to tie leaves every later box on
      // the sheet the chord has already left, certified against a span the
      // callback will never evaluate.
      //
      // Measured on fuzz seed 51 (`periodicKnotTieStillCrossesTheStrip`): the
      // chord's u-domain end and its strip top are THE SAME POINT - the usual
      // arrangement for a periodic chart, so `tU == tStrip` bit for bit - and
      // a v boundary landed 1.4e-16 away, inside `separation`. The walk
      // emitted three boxes, all on sheet `shift = -2`, all on span 4, and
      // `strips` stayed at 0. The chord's last 1.4e-15 of parameter really
      // lies in spans 1-3, which are sub-ulp wide and carry a different
      // control row: bound 1425.375608 against a true departure of 1721.11071.
      //
      // Taking the strip first costs nothing, because the tie is re-examined
      // on the next iteration from the new sheet - on seed 51 it fires there
      // and emits its pair as usual.
      // WHICH BOUNDARY COMES FIRST IS A QUESTION IN u, NOT IN t.
      //
      // `tU` and `tStrip` are both `( edge - uv0.x ) / du`, and on a chart
      // whose knots are at 1e7 and whose micro-spans are 8e-9 apart, TWO
      // DISTINCT EDGES DIVIDE TO THE SAME DOUBLE. Ordering them by t then
      // picks whichever the `<=` happens to favour, and if that is the strip
      // the walk jumps sheets over every u span still between it and the
      // strip edge - which is the index-space failure the walk exists to
      // prevent, reintroduced through the one comparison still made in
      // parameter space.
      //
      // Both edges are knot-or-strip values in EVALUATED u, so compare them
      // there: exact, and no division. Equality goes to the strip, because
      // then there is no u interval between the two to skip, and the span
      // on the far side is re-derived by `directionalSpan` after the wrap.
      const bool stripBeforeU =
        clampedU ||
        ( upU ? ( edgeU >= edgeStripLocal ) : ( edgeU <= edgeStripLocal ) );

      const bool endsOnStrip =
        periodic_ && ( du != 0.0 ) && ( tStrip <= tTo ) && stripBeforeU;

      const bool ambiguous =
        !endsOnStrip &&
        !clampedU && !clampedV && std::isfinite( tU ) && std::isfinite( tV ) &&
        ( std::abs( tU - tV ) <= separation ) &&
        ( std::min( tU, tV ) <= 1.0 );
      const bool endsOnU     = !endsOnStrip && !clampedU && ( tU <= tTo );
      const bool endsOnV     =
        !endsOnStrip && !endsOnU && !clampedV && ( tV <= tTo );
      const bool endsChord   = !endsOnU && !endsOnV && !endsOnStrip;

      // The far corner of this box. The axis that ended it is PINNED to the
      // knot value that ended it; the other is read off the chord.
      glm::dvec2 to(
        endsOnStrip ? ( edgeStrip - shift )
                    : ( endsOnU ? edgeU
                                : ( ( uv0.x + ( tTo * du ) ) - shift ) ),
        endsOnV ? edgeV : ( uv0.y + ( tTo * dv ) ) );

      if ( endsChord ) {
        to = glm::dvec2( uv1.x - shift, uv1.y );
      }

      pieces[ count++ ] = Piece{ spanU, spanV, shift, from, to, tFrom, tTo };

      if ( endsChord && !ambiguous ) {
        return true;
      }

      from  = to;
      tFrom = tTo;

      if ( ambiguous ) {

        const double late =
          std::min( 1.0, std::max( std::max( tU, tV ), tFrom ) );

        const glm::dvec2 beyond(
          ( uv0.x + ( late * du ) ) - shift, uv0.y + ( late * dv ) );

        const int steppedU =
          advanceSpan( spanU, upU, degreeU, knotsU, clampedU );

        const int steppedV =
          advanceSpan( spanV, upV, degreeV, knotsV, clampedV );

        if ( count + 2 > CERTIFICATE_MAX_PIECES ) {
          return false;
        }

        // One box for each order the two boundaries could have come in. The
        // corners are the chord's own, not pinned to either knot, because
        // which knot bounds this stretch is exactly what is not known.
        pieces[ count++ ] =
          Piece{ steppedU, spanV, shift, from, beyond, tFrom, late, true };

        pieces[ count++ ] =
          Piece{ spanU, steppedV, shift, from, beyond, tFrom, late, true };

        ++counters_.ambiguous;

        spanU = steppedU;
        spanV = steppedV;
        from  = beyond;
        tFrom = late;

        if ( late >= 1.0 ) {
          return true;
        }

        continue;
      }

      if ( endsOnStrip ) {

        // The strip's edge is not a cut of its own - it is where the u walk
        // RESTARTS, at the far end of the knot domain, one period along.
        shift += upU ? stripPeriod_ : -stripPeriod_;

        // THE ENTRY POINT AFTER A WRAP IS THE STRIP EDGE, EXACTLY - not a
        // parameter recomputed from the chord. `( uv0.x + tTo * du ) - shift`
        // is three rounded operations, and it lands PAST the edge often
        // enough to matter: measured, it came out 1.4e-14 above the domain
        // start, which stepped the walk straight over two knot spans an ulp
        // wide sitting there and left them unbounded while the callback's
        // own wrap put parameters squarely inside them. Bound 0.01577
        // against a true 0.02063.
        //
        // The chord crosses the strip edge; where it resumes is that edge
        // and nothing else, so say so rather than deriving it.
        from.x = upU ? stripUMin_ : ( stripUMin_ + stripPeriod_ );

        spanU = directionalSpan( degreeU, knotsU, from.x, upU );

        if ( spanU < 0 ) {
          return false;
        }

        // AND CLEAR THE CLAMP. `clampedU` records that the walk has run out
        // of knot domain in the direction of travel, and it is a latch - but
        // a strip restart puts the chord back at the FAR end of the domain
        // with the whole of it still to cross. Leaving the latch set makes
        // the walk stop looking for u boundaries entirely, so everything
        // after the wrap comes out as ONE box carrying the span it happened
        // to start in.
        //
        // Measured on the case the fuzzer found: a chord from u = 1 to
        // u = -1 on a strip of exactly that width starts clamped ( it begins
        // on the domain edge, descending ), wraps immediately, and then
        // emitted a single box with spanU = 4 covering u from 1 down to
        // -0.9999999979 - three spans read as one. Bound 0.1618 against a
        // true 0.1809.
        //
        // The eighth finding on bldrs-ai/conway-geom#214, and the one the
        // fuzzer found only once it was taught to drive periodic charts at
        // all: a third of this function's branches had never been executed
        // by any test.
        const bool wasClamped = clampedU;

        clampedU =
          upU ? ( from.x >= knotsU[ lastU + 1 ] )
              : ( from.x <= knotsU[ firstU ] );

        ++counters_.strips;

        if ( wasClamped && !clampedU ) {
          ++counters_.clampsClear;
        }

      } else if ( endsOnU ) {

        // Off the end of the knot domain with no strip to wrap into: the
        // surface clamps there, so advanceSpan stays in the boundary span
        // and stops the walk looking for more u edges.
        spanU = advanceSpan( spanU, upU, degreeU, knotsU, clampedU );

      } else if ( endsOnV ) {

        spanV = advanceSpan( spanV, upV, degreeV, knotsV, clampedV );
      }
    }
  }

  /**
   * Bound the deviation on one knot box, and the error that bound carries.
   *
   * The deviation restricted to a sub-interval is the same polynomial
   * composed with an affine map, so it has the same degree and the same
   * treatment; what it is measured against stays the WHOLE chord, which is
   * what the mesh will carry.
   *
   * ONE-SIDED BY CONSTRUCTION. Every node is evaluated as the polynomial of
   * `piece.spanU` / `piece.spanV`, the indices the walk emitted this box
   * with. No parameter is ever resolved back to a span here, so a node
   * sitting exactly on a knot cannot be read from the far side of it - the
   * question is not asked. See the note on `Piece`.
   *
   * NODES ARE PLACED IN uv, NOT IN t. `piece.from` and `piece.to` are the
   * box's own corners, with the axis that ended the box pinned to the exact
   * knot that ended it, so the nodes span the box even when the box's chord
   * parameters `tFrom` and `tTo` round to the same double. `t` is read only
   * to place the chord, which is affine and therefore insensitive to it at
   * that scale.
   *
   * ERROR PROPAGATION. `error` is an outward bound on everything between the
   * exact deviation polynomial and the number returned, from four sources,
   * none of which the chord's endpoints alone govern:
   *
   *   1. THE NODE VALUES. `pointHomogeneousAtSpan` accumulates
   *      `basis * controlPointW` over ( dU + 1 )( dV + 1 ) terms whose basis
   *      factors are a partition of unity, so the absolute error scales with
   *      the largest CONTROL value, not with the surface value - a patch
   *      whose values are small can have larger control values that cancel.
   *      `maxHomogeneousNorm()` is that scale. The chord subtracted from it
   *      contributes its own magnitude.
   *   2. THE INTERPOLATION. `normInfinity()` covers only a perturbation of
   *      exact node values. The inverse is itself computed in floating point
   *      ( `inverseError()`, measured a posteriori ) and applying it is a
   *      dot product of `degree + 1` terms, and both are carried here.
   *   3. WHERE THE NODES ACTUALLY LANDED. A node's uv is computed as
   *      `from + s * ( to - from )`, and that sum rounds to an ulp of the
   *      PARAMETER MAGNITUDE - not of the box width. So on a box that is
   *      narrow relative to where it sits, the nodes are not at the
   *      parameters the interpolation assumes and the values handed to it
   *      are of the right function at the wrong places. That is the fourth
   *      finding on bldrs-ai/conway-geom#214: a quadratic span two ulps wide
   *      at u = 1e9, with a bump at the one representable interior
   *      parameter, returned 6.2e-5 against a true departure of 1.
   *
   *      It is carried as an error in the node VALUE, via a bound on the
   *      surface's own gradient taken from the CONTROL NET - a B-spline's
   *      derivative control points are `degree * dP / spanWidth`, so
   *      `| dS/du | <= degreeU * 2 * maxHomogeneousNorm / spanWidthU` - times
   *      the parameter's rounding. Not as a misplacement of the abscissa:
   *      the rounding moves a node OFF the chord in uv as well as along it,
   *      and only the along-chord part is expressible as a shifted abscissa.
   *      Bounding the value directly covers both.
   *   4. THE RATIONAL QUOTIENT. `numerator / weight` amplifies the
   *      numerator's error by `1 / weight` AND turns any over-estimate of
   *      the weight's lower bound into an under-estimate of the quotient. So
   *      the weight hull is taken at its own LOWER bound, and a lower bound
   *      that reaches zero is refused rather than divided by.
   */
  bool boundPiece(
    const Piece&      piece,
    const glm::dvec3& surface0,
    const glm::dvec3& surface1,
    double            chordScale,
    double            parameterScaleU,
    double            parameterScaleV,
    double&           result,
    double&           error ) const {

    const uint32_t count = interpolation_.degree() + 1;

    if ( !RationalSurfaceEvaluator::spanIsEvaluable(
           evaluator_->knotsU(), piece.spanU ) ||
         !RationalSurfaceEvaluator::spanIsEvaluable(
           evaluator_->knotsV(), piece.spanV ) ) {
      return false;
    }

    const glm::dvec2 across = piece.to - piece.from;

    glm::dvec3 numerators[ CERTIFICATE_MAX_DEGREE + 1 ];
    double     weights[ CERTIFICATE_MAX_DEGREE + 1 ];

    double largestNode   = 0.0;
    double largestWeight = 0.0;

    for ( uint32_t i = 0; i < count; ++i ) {

      const double s = interpolation_.node( i );

      const glm::dvec2 uv = piece.from + ( across * s );

      const double at = piece.tFrom + ( ( piece.tTo - piece.tFrom ) * s );

      const glm::dvec3 chord =
        ( surface0 * ( 1.0 - at ) ) + ( surface1 * at );

      ++counters_.evaluations;

      if ( rational_ ) {

        const glm::dvec4 homogeneous =
          evaluator_->pointHomogeneousAtSpan(
            piece.spanU, piece.spanV, uv.x, uv.y );

        // N( t ) = A( t ) - W( t ) * chord( t ), both polynomial.
        numerators[ i ] =
          glm::dvec3( homogeneous ) - ( homogeneous.w * chord );

        weights[ i ] = homogeneous.w;

      } else {

        numerators[ i ] =
          evaluator_->pointAtSpan(
            piece.spanU, piece.spanV, uv.x, uv.y ) - chord;

        weights[ i ] = 1.0;
      }

      if ( !std::isfinite( numerators[ i ].x ) ||
           !std::isfinite( numerators[ i ].y ) ||
           !std::isfinite( numerators[ i ].z ) ||
           !std::isfinite( weights[ i ] ) ) {
        return false;
      }

      largestNode   = std::max( largestNode, glm::length( numerators[ i ] ) );
      largestWeight = std::max( largestWeight, std::abs( weights[ i ] ) );
    }

    double largestNumerator     = 0.0;
    double largestWeightHull    = 0.0;
    double smallestWeight       = std::numeric_limits< double >::infinity();

    for ( uint32_t i = 0; i < count; ++i ) {

      largestNumerator =
        std::max(
          largestNumerator,
          glm::length( interpolation_.coefficient( i, numerators ) ) );

      const double weight = interpolation_.coefficient( i, weights );

      largestWeightHull = std::max( largestWeightHull, std::abs( weight ) );
      smallestWeight    = std::min( smallestWeight, weight );
    }

    if ( !std::isfinite( largestNumerator ) ||
         !std::isfinite( smallestWeight ) ) {
      return false;
    }

    const double epsilon = std::numeric_limits< double >::epsilon();

    const double terms =
      static_cast< double >( evaluator_->degreeU() + 1 ) *
      static_cast< double >( evaluator_->degreeV() + 1 );

    // 1. THE NODE VALUES.
    const double evaluationError =
      CERTIFICATE_ERROR_SAFETY * terms * epsilon *
      evaluator_->maxHomogeneousNorm();

    double weightError =
      rational_ ?
        ( CERTIFICATE_ERROR_SAFETY * terms * epsilon *
          evaluator_->maxWeight() ) :
        0.0;

    const double chordError =
      CERTIFICATE_ERROR_SAFETY * epsilon * chordScale;

    // 3. WHERE THE NODES ACTUALLY LANDED, as an error in the VALUE.
    // ROUNDING IS GOVERNED BY THE RAW PARAMETER, NOT THE WRAPPED ONE.
    //
    // On a periodic chart the callback evaluates `point( wrapChartU( u ), v )`
    // and the walk evaluates `point( u - shift, v )`, and THOSE ARE NOT THE
    // SAME NUMBER. IEEE `fmod` is exact; `u - floor( ( u - uMin ) / P ) * P`
    // is three rounded operations. Measured on the case the fuzzer found,
    // the two disagree on 98,591 of 200,001 sampled parameters along one
    // chord - by an ulp of the RAW parameter, which on a chart whose chords
    // run far outside the strip is much larger than an ulp of the wrapped
    // one this used to read.
    //
    // So the magnitude that prices a node's placement is the raw parameter's,
    // and `| shift |` is what carries it: raw = wrapped + shift, so
    // | raw | <= | wrapped | + | shift |. Zero on every non-periodic chord,
    // which is why nothing outside the periodic path moves.
    //
    // This is the eighth finding on bldrs-ai/conway-geom#214 and the same
    // shape as the first: an error term scaled by the wrong quantity.
    // AND BY THE STRIP'S OWN MAGNITUDE, because the reduction is done
    // THROUGH `stripUMin_`.
    //
    // `wrapChartU( x )` is `uMin + fmod( x - uMin, P )`: `fmod` is exact, but
    // the subtraction and the addition each round at the magnitude of `uMin`,
    // not at the magnitude of x. The walk's own points are `x - shift` with
    // `shift` fixed from the chord's start, so for every t but the first the
    // two disagree by a few ulps OF THE STRIP, which can be larger than an
    // ulp of the parameter.
    //
    // Measured on the grazing generator's trial 2,967: a box whose u corner
    // is 1.6e-12 from the knot the surface steps at, with `roundingU` at
    // 7.1e-13 before this term - so the box was judged not to touch the
    // boundary, and the caller evaluated the whole of it on the far side.
    // Bound 2.181686448 against a true departure of 2.563162239.
    const double roundingU =
      CERTIFICATE_ERROR_SAFETY * epsilon *
      ( parameterScaleU + std::abs( piece.shift ) +
        ( periodic_ ? ( std::abs( stripUMin_ ) + stripPeriod_ ) : 0.0 ) );

    const double roundingV =
      CERTIFICATE_ERROR_SAFETY * epsilon * parameterScaleV;

#if CERTIFICATE_AUDIT_CALLER_SPANS

    // Walk the same nodes again, asking the question the CALL SITE asks:
    // given this double, which span does `findSpan` name? See
    // CERTIFICATE_AUDIT_CALLER_SPANS.
    {
      const std::vector< double >& auditKnotsU = evaluator_->knotsU();
      const std::vector< double >& auditKnotsV = evaluator_->knotsV();

      const uint32_t auditDegreeU = evaluator_->degreeU();
      const uint32_t auditDegreeV = evaluator_->degreeV();

      // A node's disagreement is BENIGN when the node sits exactly on a knot
      // that bounds this box: the next box the walk emits carries that span
      // and bounds that double against the polynomial the caller will use.
      // A node is excused when it is within the placement rounding of a knot
      // that BOUNDS this box. Exact equality is not the right test: the
      // corner of the axis that did NOT end the piece is read off the chord,
      // so it lands up to an ulp outside its own span, and `findSpan` then
      // names the neighbour. That ulp is what `roundingU` / `roundingV`
      // already price, and the neighbouring box covers the double.
      const auto onOwnBoundary =
        []( const std::vector< double >& knots, int span, double at,
            double rounding ) {

          return std::abs( at - knots[ span ] ) <= rounding ||
                 std::abs( at - knots[ span + 1 ] ) <= rounding;
        };

      uint32_t agreeing = 0;
      uint32_t edge     = 0;
      uint32_t inside   = 0;

      for ( uint32_t i = 0; i < count; ++i ) {

        const double s = interpolation_.node( i );

        const glm::dvec2 uv = piece.from + ( across * s );

        const bool sameU =
          RationalSurfaceEvaluator::findSpan( auditDegreeU, auditKnotsU, uv.x )
            == piece.spanU;

        const bool sameV =
          RationalSurfaceEvaluator::findSpan( auditDegreeV, auditKnotsV, uv.y )
            == piece.spanV;

        if ( sameU && sameV ) {
          ++agreeing;
          continue;
        }

        const bool excusedU =
          sameU ||
          onOwnBoundary( auditKnotsU, piece.spanU, uv.x, roundingU );

        const bool excusedV =
          sameV ||
          onOwnBoundary( auditKnotsV, piece.spanV, uv.y, roundingV );

        if ( excusedU && excusedV ) { ++edge; } else { ++inside; }
      }

      if ( piece.covering ) {

        // Not an error and not a coverage signal: these two boxes are
        // deliberately not inside their own span, and between them they
        // cover both spans the caller can resolve to.
        ++counters_.callerSpanPair;

      } else if ( inside > 0 ) {
        ++counters_.callerSpanInside;
      } else if ( agreeing == count ) {
        ++counters_.callerSpanAgree;
      } else if ( agreeing == 0 && piece.tTo > piece.tFrom ) {

        // A box of zero parameter width is excluded: the walk emits those at
        // the chord's ends and wherever two boundaries coincide, and their
        // single point is the cached endpoint, whose value the caller
        // supplied. What is left is a box covering a REAL stretch of chord,
        // every node of which the caller resolves to another polynomial.
        ++counters_.callerSpanWhole;

        if ( stepsU_ || stepsV_ ) {
          ++counters_.callerSpanWholeStep;
        }

      } else {
        ++counters_.callerSpanEdge;
      }
    }

#endif

    double spreadPointU  = 0.0;
    double spreadWeightU = 0.0;
    double spreadPointV  = 0.0;
    double spreadWeightV = 0.0;

    evaluator_->controlSpread(
      piece.spanU, piece.spanV, true, spreadPointU, spreadWeightU );

    evaluator_->controlSpread(
      piece.spanU, piece.spanV, false, spreadPointV, spreadWeightV );

    const double gradientU =
      gradientBound( evaluator_->degreeU(), evaluator_->knotsU(),
                     piece.spanU, spreadPointU );

    const double gradientV =
      gradientBound( evaluator_->degreeV(), evaluator_->knotsV(),
                     piece.spanV, spreadPointV );

    const double placementError =
      ( gradientU * roundingU ) + ( gradientV * roundingV );

    const double placementWeightError =
      rational_ ?
        ( ( gradientBound( evaluator_->degreeU(), evaluator_->knotsU(),
                           piece.spanU, spreadWeightU ) * roundingU ) +
          ( gradientBound( evaluator_->degreeV(), evaluator_->knotsV(),
                           piece.spanV, spreadWeightV ) * roundingV ) ) :
        0.0;

    if ( !std::isfinite( placementError ) ||
         !std::isfinite( placementWeightError ) ) {
      ++counters_.unplaceable;
      return false;
    }

    // A BOX OUTSIDE THE KNOT DOMAIN IS EXTRAPOLATED, AND THE GRADIENT BOUND
    // DOES NOT HOLD THERE.
    //
    // `findSpan` clamps a parameter past either end of the domain into the
    // boundary span, and the basis is then evaluated outside that span - so
    // the surface is a polynomial EXTRAPOLATION, whose derivative grows
    // without any relation to `controlSpread / spanWidth`. The placement
    // term would be priced from a gradient that only describes the span's
    // interior.
    //
    // Found by the fuzzer: a chord reaching 8.4e6 beyond a domain of 2.6e8
    // returned a bound 0.1% below the truth, with both around 1.5e12 because
    // extrapolation had taken a surface of size 1e-3 there.
    //
    // Declining rather than pricing it, because an extrapolated box is
    // outside what the certificate is a statement about at all.
    {
      const std::vector< double >& knotsU = evaluator_->knotsU();
      const std::vector< double >& knotsV = evaluator_->knotsV();

      const double lowU  = knotsU[ evaluator_->degreeU() ];
      const double highU = knotsU[ knotsU.size() - evaluator_->degreeU() - 1 ];
      const double lowV  = knotsV[ evaluator_->degreeV() ];
      const double highV = knotsV[ knotsV.size() - evaluator_->degreeV() - 1 ];

      // Outside by a ROUNDING is not extrapolation, and refusing it would be
      // expensive: a trim solve lands parameters on the domain edge, where a
      // couple of ulps either way is the normal case rather than a
      // pathology. Measured on `Right_Hand.step`, refusing those too
      // declined 126 chords and cost 122 degenerate triangles.
      //
      // The tolerance is scaled by the DOMAIN, not by the parameter. A face
      // whose u domain is [ -1.1e-17, 0.0796 ] has boxes sitting at the
      // near-zero end whose own magnitude is 1e-17, so a parameter-scaled
      // tolerance comes out at 2e-32 and refuses an overshoot of 4.5e-32
      // that is two ulps of nothing. What decides whether a parameter is
      // outside the domain is the domain's own scale. Against the case the
      // fuzzer found - 8.4e6 beyond a domain of 2.6e8, four parts in a
      // thousand - this still declines by six orders of magnitude.
      const double reachU =
        roundingU + ( CERTIFICATE_ERROR_SAFETY * epsilon *
                      std::max( std::abs( lowU ), std::abs( highU ) ) );

      const double reachV =
        roundingV + ( CERTIFICATE_ERROR_SAFETY * epsilon *
                      std::max( std::abs( lowV ), std::abs( highV ) ) );

      if ( piece.from.x < lowU - reachU || piece.from.x > highU + reachU ||
           piece.to.x   < lowU - reachU || piece.to.x   > highU + reachU ||
           piece.from.y < lowV - reachV || piece.from.y > highV + reachV ||
           piece.to.y   < lowV - reachV || piece.to.y   > highV + reachV ) {
        ++counters_.extrapolated;
        return false;
      }
    }

    // A BOX TOO THIN TO SEPARATE FROM ITS NEIGHBOUR, ON AN AXIS THAT STEPS.
    //
    // `placementError` prices a node landing slightly off its intended
    // parameter using the surface's GRADIENT, which is the right price only
    // while the node stays inside this box. Where the box's extent on an
    // axis is no larger than that axis's own rounding, and the box touches a
    // span boundary, a node can land on the far side of the boundary - and
    // across a knot of multiplicity degree + 1 the two sides are different
    // polynomials with a STEP between them, so the cost is the step, not the
    // gradient.
    //
    // Found by the differential fuzzer rather than by review: a chord
    // descending in u by 2e-9 while v swept a whole domain, on a surface
    // with a multiplicity-3 knot in u, had its evaluated u PINNED to the
    // knot across a stretch of the chord - so the callback read one span for
    // that stretch while the walk bounded the other. Bound 0.0134 against a
    // true 0.0143.
    //
    // Declining rather than covering both sides, because unlike the
    // ambiguous-ordering case this is not a stretch of chord with two
    // candidate spans - it is a box whose own identity is not resolvable.
    //
    // THE TEST IS ON THE SPAN'S OWN WIDTH, NOT ON HOW MUCH OF IT THIS CHORD
    // CROSSES, and the difference is what makes it affordable. A chord
    // running almost parallel to v has boxes with no u extent at all, and
    // there is nothing wrong with those - u is constant, the span is
    // unambiguous, and the nodes are exactly where they should be. What
    // cannot be sampled is a KNOT SPAN narrower than the parameter's own
    // resolution, because then no node placed by
    // `from + s * ( to - from )` lands inside it whatever the chord does.
    //
    // Measured: testing the box's extent instead of the span's declines
    // ordinary v-parallel chords and costs 686 degenerate triangles on
    // `Right_Hand.step` ( 1,075 -> 1,761 ). Testing the span's width costs
    // nothing there and still refuses every case the fuzzer found.
    //
    // Only a box that covers a STRETCH of the chord can hide a mismatch:
    // the walk emits zero-length boxes at the chord's ends and wherever two
    // boundaries coincide, and those are single points whose value is the
    // cached endpoint.
    const auto spanIsBelowResolution =
      []( const std::vector< double >& knots, int span, double rounding ) {

        return ( knots[ span + 1 ] - knots[ span ] ) <= rounding;
      };

    // THIS ARM STILL HAS NO RED PROOF. The one that used to sit beside it was
    // deleted for want of one and HAS BEEN PUT BACK - see below - so the
    // reasoning that removed it is worth reading before it is reused here.
    //
    // At the time, 200,000 trials across three seeds found nothing either arm
    // caught, and the only discriminator left was cost: the stepping arm
    // declined 46,000 of 250,000 fuzz cases, this one 30 of 250,000 and none
    // at all on the corpus. So the expensive one went. Then the fuzzer was
    // taught to drive chords far outside the strip and to probe knots from
    // both sides, and it produced fuzz seed 7 at trial 173,559 - a case the
    // deleted arm catches and nothing else does. "Nothing can red-prove it"
    // had meant "the generator has not reached it", not "it is not carrying
    // anything".
    //
    // This arm is kept on its argument, with that correction attached: a knot
    // span narrower than the parameter's own resolution cannot be sampled,
    // whatever the chord does, because no node placed by
    // `from + s * ( to - from )` lands inside it. At 30 in 250,000 the price
    // of being wrong about it is low, and the price of being wrong the other
    // way has now been demonstrated once.
    // THE SECOND ARM IS BACK, AND THIS TIME IT HAS A RED PROOF.
    //
    // It catches a box whose own extent on an axis is below the parameter's
    // rounding AND that sits on a span boundary the surface STEPS across. It
    // was deleted earlier in this PR on the reasoning that it declined 46,000
    // of 250,000 fuzz cases while catching nothing - the discriminator being
    // COST, since nothing could red-prove it. The fuzzer then produced fuzz
    // seed 7 at trial 173,559, which it catches and nothing else does, so the
    // reasoning was wrong: it was catching something the generator had not
    // yet reached.
    //
    // What it catches there: the chord's v moves 2e-9 in total and its far
    // end sits exactly on a v knot of multiplicity degree + 1, so the ROUNDED
    // v equals that knot for every t above 1 - 2.8e-8. The call site resolves
    // that double with `findSpan` and gets the span on the knot's right; the
    // walk had the whole chord on the span to its left; the two polynomials
    // differ by 248.57 at the point that matters. Bound 335.968276 against a
    // true departure of 408.4816988.
    //
    // `stepsU_` / `stepsV_` keep it qualified to axes that actually step -
    // without that it also refuses every box of a chord running parallel to
    // the other axis, measured at 686 extra degenerate triangles on
    // `Right_Hand.step`.
    // WITHIN THE ROUNDING OF THE BOUNDARY, not exactly on it. This used to
    // test `from == knots[ span ]` and so on, and that is too narrow by
    // exactly the quantity the rest of this function is built to price.
    //
    // The box's corners are the WALK's numbers; the caller's parameter at
    // the same t is `wrapChartU( uv0.x + t * du )`, and the two differ by up
    // to `roundingU`. So a box can sit a few ulps short of a knot, pass an
    // equality test, and still have every parameter the caller evaluates
    // land on the far side of it. Measured on the grazing generator's trial
    // 272: the box's u corner is 41.412005379020393 against a knot at
    // 41.4120053790204 - short by 7.1e-15, which is the periodic shift -
    // while `roundingU` is 7.4e-14, and the caller evaluates the whole box
    // at 41.41200537902045653, on the other side. Bound 1.639305956 against
    // a true departure of 2.013100244.
    //
    // The error term prices that displacement as a GRADIENT times the
    // rounding. Across a knot the surface steps at, the cost is a JUMP, and
    // no gradient covers it - which is the whole reason this arm exists.
    const auto touchesBoundary =
      []( const std::vector< double >& knots, int span,
          double from, double to, double rounding ) {

        return std::abs( from - knots[ span ] ) <= rounding ||
               std::abs( from - knots[ span + 1 ] ) <= rounding ||
               std::abs( to - knots[ span ] ) <= rounding ||
               std::abs( to - knots[ span + 1 ] ) <= rounding;
      };

    if ( piece.tTo > piece.tFrom &&
         ( spanIsBelowResolution( evaluator_->knotsU(), piece.spanU,
                                  roundingU ) ||
           spanIsBelowResolution( evaluator_->knotsV(), piece.spanV,
                                  roundingV ) ||
           ( stepsU_ && std::abs( across.x ) <= roundingU &&
             touchesBoundary( evaluator_->knotsU(), piece.spanU,
                              piece.from.x, piece.to.x, roundingU ) ) ||
           ( stepsV_ && std::abs( across.y ) <= roundingV &&
             touchesBoundary( evaluator_->knotsV(), piece.spanV,
                              piece.from.y, piece.to.y, roundingV ) ) ) ) {
      ++counters_.unplaceable;
      return false;
    }

    const double nodeError =
      placementError +
      ( rational_ ?
        ( evaluationError + ( largestWeight * chordError ) +
          ( chordScale * ( weightError + placementWeightError ) ) +
          ( epsilon * largestWeight * chordScale ) ) :
        ( evaluationError + chordError ) );

    // 2. THE INTERPOLATION.
    const double applyError =
      interpolation_.inverseError() +
      ( count * epsilon * interpolation_.normInfinity() );

    const double numeratorError =
      ( interpolation_.normInfinity() * nodeError ) +
      ( applyError * largestNode );

    // `placementWeightError` belongs HERE as well as in `nodeError`. It
    // bounds how far the sampled weights are from the intended ones, and the
    // weight hull is reconstructed from those same samples - so leaving it
    // out lets `weightFloor` sit above the true minimum denominator, and a
    // denominator that is too large makes the quotient too small. That is
    // the one direction this bound may not err in. The sixth finding on
    // bldrs-ai/conway-geom#214, and the same shape as the first: a stage of
    // the chain not carrying an error the stage before it computed.
    const double weightCoefficientError =
      rational_ ?
        ( ( interpolation_.normInfinity() *
            ( weightError + placementWeightError ) ) +
          ( applyError * largestWeight ) ) :
        0.0;

    // 4. THE RATIONAL QUOTIENT. The weight hull is taken at its own LOWER
    // bound: over-estimating it would under-estimate the quotient, which is
    // the one direction a bound may not err in. A lower bound that reaches
    // zero makes N / W unbounded on the box and is refused - "the bound does
    // not exist" is not the same as "the bound is small".
    const double weightFloor = smallestWeight - weightCoefficientError;

    if ( !( weightFloor > 0.0 ) ) {
      return false;
    }

    // The un-widened reading, kept only so `error` can say how much of what
    // is returned is outward widening rather than measurement.
    const double raw = largestNumerator / smallestWeight;

    // One more rounding for the division and the sums above.
    result =
      ( ( largestNumerator + numeratorError ) / weightFloor ) *
      ( 1.0 + ( 4.0 * epsilon ) );

    error = result - raw;

    if ( !std::isfinite( result ) || !( error >= 0.0 ) ) {
      return false;
    }

    return true;
  }

  /**
   * Bound on how fast the surface can move per unit of one parameter, over
   * one span, read off the CONTROL NET rather than from a derivative
   * evaluation.
   *
   * A B-spline's derivative is itself a B-spline whose control points are
   * `degree * ( P[ i + 1 ] - P[ i ] ) / ( knot difference )`, and the basis
   * is a partition of unity, so the derivative is bounded by the largest of
   * them. `spread` is that largest control step over the span's OWN window
   * ( see controlSpread ), and the denominator here is the span's own width
   * - which is never larger than the knot differences the real formula uses,
   * so the result is never smaller than the real bound.
   *
   * On ordinary geometry the product with a few ulps of parameter rounding
   * is many orders below the tolerance; on a span too narrow to place nodes
   * in, it grows until the error gate declines the chord.
   */
  static double gradientBound(
    uint32_t                     degree,
    const std::vector< double >& knots,
    int                          span,
    double                       spread ) {

    const double width = knots[ span + 1 ] - knots[ span ];

    if ( !( width > 0.0 ) ) {
      return std::numeric_limits< double >::infinity();
    }

    return static_cast< double >( degree ) * spread / width;
  }

  const RationalSurfaceEvaluator* evaluator_   = nullptr;
  BernsteinNodes                  interpolation_;
  bool                            periodic_    = false;
  bool                            rational_    = false;
  bool                            supported_   = false;
  bool                            stepsU_      = false;
  bool                            stepsV_      = false;
  double                          stripUMin_   = 0.0;
  double                          stripPeriod_ = 0.0;
  double                          tolerance_   = 0.0;
  mutable Counters                counters_;
};

}  // namespace conway::geometry
