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
 * Most knot spans a single chord may cross and still be certified. The
 * deviation is only polynomial WITHIN a knot box, so a chord crossing knot
 * lines is cut at the crossings and bounded piecewise; this caps that work.
 *
 * Refusing above the cap is conservative - the edge is subdivided - and it
 * terminates, because halving an edge halves the spans it crosses, so a chord
 * over the cap is under it within log2( spans / cap ) levels.
 */
constexpr uint32_t CERTIFICATE_MAX_SPANS = 8;

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

    supported_ =
      evaluator.supportsFastPath() &&
      totalDegree >= 1 &&
      totalDegree <= CERTIFICATE_MAX_DEGREE &&
      ( !periodic || stripPeriod > 0.0 );

    if ( supported_ ) {

      interpolation_.build( totalDegree );

      supported_ =
        interpolation_.degree() == totalDegree && interpolation_.valid();
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
    uint64_t illConditioned = 0;
    uint64_t badPiece     = 0;
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

    // The chord's own magnitude, which is one of the three things the error
    // term below is built from. It is NOT on its own a scale for the
    // evaluation error - see maxHomogeneousNorm() in nurbs_utils.h, and the
    // ERROR PROPAGATION note on boundPiece.
    const double chordScale =
      std::max( glm::length( surface0 ), glm::length( surface1 ) );

    // STAGE 1: the cuts that are not u knots - the periodic chart's SHEET
    // boundaries and the v knot lines.
    //
    // The chart wrap makes t -> u( t ) discontinuous where the segment steps
    // across the cut, and a discontinuous function has no polynomial through
    // it. Cutting there rather than declining is worth the code: measured on
    // `Right_Hand.step`, declining a wrapped chord left 22,634 of 216,757
    // candidates ( 10.4% ) uncertifiable, all of which then had to be
    // subdivided on no evidence, and solid #19715 came out with 1,737
    // degenerate triangles against 164 before the certificate.
    double cuts[ CERTIFICATE_MAX_SPANS + 1 ];

    uint32_t cutCount = 0;

    cuts[ cutCount++ ] = 0.0;

    if ( periodic_ ) {

      const double sheet0 = std::floor( ( uv0.x - stripUMin_ ) / stripPeriod_ );
      const double sheet1 = std::floor( ( uv1.x - stripUMin_ ) / stripPeriod_ );

      const double lowSheet  = std::min( sheet0, sheet1 );
      const double highSheet = std::max( sheet0, sheet1 );

      for ( double sheet = lowSheet + 1.0; sheet <= highSheet; sheet += 1.0 ) {

        if ( !addCut( stripUMin_ + ( sheet * stripPeriod_ ),
                      uv0.x, uv1.x, cuts, cutCount ) ) {

          ++counters_.tooManySpans;
          ++counters_.inconclusive;
          return CertificateOutcome::Inconclusive;
        }
      }
    }

    if ( !collectCrossings(
           evaluator_->degreeV(), evaluator_->knotsV(), uv0.y, uv1.y,
           cuts, cutCount ) ) {

      ++counters_.tooManySpans;
      ++counters_.inconclusive;
      return CertificateOutcome::Inconclusive;
    }

    cuts[ cutCount++ ] = 1.0;

    std::sort( cuts, cuts + cutCount );

    double worst     = 0.0;
    double worstError = 0.0;

    for ( uint32_t piece = 0; piece + 1 < cutCount; ++piece ) {

      const double from = cuts[ piece ];
      const double to   = cuts[ piece + 1 ];

      if ( !( to > from ) ) {
        continue;
      }

      // STAGE 2: this piece lies on ONE sheet, so `shift` is constant across
      // it and t -> u is affine again. Read at the MIDPOINT rather than at an
      // end, because a cut sits exactly on a sheet boundary and rounding
      // there could put it on either side.
      double shift = 0.0;

      if ( periodic_ ) {

        const double middle = 0.5 * ( from + to );

        shift =
          std::floor(
            ( ( uv0.x + ( middle * ( uv1.x - uv0.x ) ) ) - stripUMin_ ) /
            stripPeriod_ ) * stripPeriod_;
      }

      double inner[ CERTIFICATE_MAX_SPANS + 1 ];

      uint32_t innerCount = 0;

      inner[ innerCount++ ] = from;

      const double pieceU0 = uv0.x + ( from * ( uv1.x - uv0.x ) ) - shift;
      const double pieceU1 = uv0.x + ( to * ( uv1.x - uv0.x ) ) - shift;

      if ( !collectCrossings(
             evaluator_->degreeU(), evaluator_->knotsU(),
             pieceU0, pieceU1, inner, innerCount, from, to ) ) {

        ++counters_.tooManySpans;
        ++counters_.inconclusive;
        return CertificateOutcome::Inconclusive;
      }

      inner[ innerCount++ ] = to;

      std::sort( inner, inner + innerCount );

      for ( uint32_t part = 0; part + 1 < innerCount; ++part ) {

        if ( !( inner[ part + 1 ] > inner[ part ] ) ) {
          continue;
        }

        ++counters_.spans;

        double pieceBound = 0.0;
        double pieceError = 0.0;

        if ( !boundPiece(
               uv0, uv1, surface0, surface1, chordScale,
               inner[ part ], inner[ part + 1 ], shift,
               pieceBound, pieceError ) ) {

          ++counters_.badPiece;
          ++counters_.inconclusive;
          return CertificateOutcome::Inconclusive;
        }

        worst      = std::max( worst, pieceBound );
        worstError = std::max( worstError, pieceError );
      }
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

  /** Record the parameter at which `a -> b` passes `value`, if it does. */
  static bool addCut(
    double    value,
    double    a,
    double    b,
    double*   cuts,
    uint32_t& cutCount ) {

    const double low  = std::min( a, b );
    const double high = std::max( a, b );

    if ( value <= low || value >= high ) {
      return true;
    }

    if ( cutCount + 1 >= CERTIFICATE_MAX_SPANS + 1 ) {
      return false;
    }

    cuts[ cutCount++ ] = ( value - a ) / ( b - a );

    return true;
  }

  /**
   * Append the parameters at which the segment `a -> b` crosses an INTERIOR
   * knot of `knots`. False if that would take the piece count over the cap.
   *
   * Interior only: the first and last `degree + 1` knots clamp the domain and
   * are not places the polynomial changes.
   */
  bool collectCrossings(
    uint32_t                     degree,
    const std::vector< double >& knots,
    double                       a,
    double                       b,
    double*                      cuts,
    uint32_t&                    cutCount,
    double                       tFrom = 0.0,
    double                       tTo   = 1.0 ) const {

    if ( knots.size() < ( 2 * degree ) + 2 ) {
      return true;
    }

    const double low  = std::min( a, b );
    const double high = std::max( a, b );

    if ( !( high > low ) ) {
      return true;
    }

    const double span = b - a;

    for ( size_t at = degree + 1, end = knots.size() - degree - 1;
          at < end;
          ++at ) {

      const double knot = knots[ at ];

      if ( knot <= low || knot >= high ) {
        continue;
      }

      if ( cutCount + 1 >= CERTIFICATE_MAX_SPANS + 1 ) {
        return false;
      }

      cuts[ cutCount++ ] =
        tFrom + ( ( ( knot - a ) / span ) * ( tTo - tFrom ) );
    }

    return true;
  }

  /**
   * Bound the deviation on the sub-interval [ `from`, `to` ] of the chord,
   * and the error that bound carries.
   *
   * The deviation restricted to a sub-interval is the same polynomial
   * composed with an affine map, so it has the same degree and the same
   * treatment; what it is measured against stays the WHOLE chord, which is
   * what the mesh will carry.
   *
   * ONE-SIDED EVALUATION. Every node is evaluated as the polynomial of the
   * spans this piece lies in, named from the piece's MIDPOINT, not of the
   * spans its parameters fall in. At an interior knot of multiplicity
   * `degree + 1` the two sides are different polynomials with a step between
   * them, and `findSpan` resolves a parameter sitting exactly on the knot to
   * the right-hand one - so a piece whose last node is that knot would be
   * given one sample from the polynomial it is NOT certifying, and the
   * interpolant would be fitted through a point that is not on the curve.
   * The Chebyshev-Lobatto nodes make this certain rather than unlikely:
   * they always include both endpoints. Measured on the smallest case that
   * shows it - degrees ( 1, 1 ), a knot of multiplicity 2, left-hand
   * deviation rising to 1 and right-hand knot value 2/3 - the samples come
   * out [ 0, 1/2, 2/3 ], the coefficients [ 0, 2/3, 2/3 ], and the bound
   * 0.667 against a true 1.0. See `fullMultiplicityKnotIsSampledOneSided`.
   *
   * ERROR PROPAGATION. `error` is an outward bound on everything between the
   * exact deviation polynomial and the number returned, and it is built from
   * three sources, none of which the chord's endpoints alone govern:
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
   *   3. THE RATIONAL QUOTIENT. `numerator / weight` amplifies the
   *      numerator's error by `1 / weight` AND turns any over-estimate of
   *      the weight's lower bound into an under-estimate of the quotient. So
   *      the weight hull is taken at its own LOWER bound, and a lower bound
   *      that reaches zero is refused rather than divided by.
   */
  bool boundPiece(
    const glm::dvec2& uv0,
    const glm::dvec2& uv1,
    const glm::dvec3& surface0,
    const glm::dvec3& surface1,
    double            chordScale,
    double            from,
    double            to,
    double            uShift,
    double&           result,
    double&           error ) const {

    const uint32_t count = interpolation_.degree() + 1;

    // The spans this piece lies in, read at the midpoint so that a cut
    // sitting exactly on a knot resolves to the piece's own side of it.
    const double middle = 0.5 * ( from + to );

    const glm::dvec2 midUV =
      ( uv0 * ( 1.0 - middle ) ) + ( uv1 * middle ) -
      glm::dvec2( uShift, 0.0 );

    const int spanU =
      RationalSurfaceEvaluator::findSpan(
        evaluator_->degreeU(), evaluator_->knotsU(), midUV.x );

    const int spanV =
      RationalSurfaceEvaluator::findSpan(
        evaluator_->degreeV(), evaluator_->knotsV(), midUV.y );

    // A span with an empty knot interval has no polynomial to name, and the
    // Cox-de-Boor denominators are only bounded away from zero on a
    // non-empty one.
    if ( !RationalSurfaceEvaluator::spanIsEvaluable(
           evaluator_->knotsU(), spanU ) ||
         !RationalSurfaceEvaluator::spanIsEvaluable(
           evaluator_->knotsV(), spanV ) ) {
      return false;
    }

    glm::dvec3 numerators[ CERTIFICATE_MAX_DEGREE + 1 ];
    double     weights[ CERTIFICATE_MAX_DEGREE + 1 ];

    double largestNode   = 0.0;
    double largestWeight = 0.0;

    for ( uint32_t i = 0; i < count; ++i ) {

      const double at = from + ( ( to - from ) * interpolation_.node( i ) );

      const glm::dvec2 uv =
        ( uv0 * ( 1.0 - at ) ) + ( uv1 * at ) - glm::dvec2( uShift, 0.0 );

      const glm::dvec3 chord =
        ( surface0 * ( 1.0 - at ) ) + ( surface1 * at );

      ++counters_.evaluations;

      if ( rational_ ) {

        const glm::dvec4 homogeneous =
          evaluator_->pointHomogeneousAtSpan( spanU, spanV, uv.x, uv.y );

        // N( t ) = A( t ) - W( t ) * chord( t ), both polynomial.
        numerators[ i ] =
          glm::dvec3( homogeneous ) - ( homogeneous.w * chord );

        weights[ i ] = homogeneous.w;

      } else {

        numerators[ i ] =
          evaluator_->pointAtSpan( spanU, spanV, uv.x, uv.y ) - chord;

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

    const double epsilon = std::numeric_limits< double >::epsilon();

    const double terms =
      static_cast< double >( evaluator_->degreeU() + 1 ) *
      static_cast< double >( evaluator_->degreeV() + 1 );

    // 1. THE NODE VALUES.
    const double evaluationError =
      CERTIFICATE_ERROR_SAFETY * terms * epsilon *
      evaluator_->maxHomogeneousNorm();

    const double weightError =
      rational_ ?
        ( CERTIFICATE_ERROR_SAFETY * terms * epsilon *
          evaluator_->maxWeight() ) :
        0.0;

    const double chordError =
      CERTIFICATE_ERROR_SAFETY * epsilon * chordScale;

    // N = A - W * chord, so the chord's error enters scaled by the weight
    // and the weight's error scaled by the chord.
    const double nodeError =
      rational_ ?
        ( evaluationError + ( largestWeight * chordError ) +
          ( chordScale * weightError ) +
          ( epsilon * largestWeight * chordScale ) ) :
        ( evaluationError + chordError );

    // 2. THE INTERPOLATION: exact inverse on perturbed values, plus the
    // inverse's own error, plus the rounding of applying it.
    const double applyError =
      interpolation_.inverseError() +
      ( count * epsilon * interpolation_.normInfinity() );

    const double numeratorError =
      ( interpolation_.normInfinity() * nodeError ) +
      ( applyError * largestNode );

    const double weightCoefficientError =
      rational_ ?
        ( ( interpolation_.normInfinity() * weightError ) +
          ( applyError * largestWeight ) ) :
        0.0;

    double largestNumerator = 0.0;
    double smallestWeight   = std::numeric_limits< double >::infinity();

    for ( uint32_t i = 0; i < count; ++i ) {

      largestNumerator =
        std::max(
          largestNumerator,
          glm::length( interpolation_.coefficient( i, numerators ) ) );

      smallestWeight =
        std::min( smallestWeight,
                  interpolation_.coefficient( i, weights ) );
    }

    if ( !std::isfinite( largestNumerator ) ||
         !std::isfinite( smallestWeight ) ) {
      return false;
    }

    // 3. THE RATIONAL QUOTIENT. The weight hull is taken at its own LOWER
    // bound: over-estimating it would under-estimate the quotient, which is
    // the one direction a bound may not err in. A lower bound that reaches
    // zero makes N / W unbounded on the piece and is refused - "the bound
    // does not exist" is not the same as "the bound is small".
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

  const RationalSurfaceEvaluator* evaluator_   = nullptr;
  BernsteinNodes                  interpolation_;
  bool                            periodic_    = false;
  bool                            rational_    = false;
  bool                            supported_   = false;
  double                          stripUMin_   = 0.0;
  double                          stripPeriod_ = 0.0;
  double                          tolerance_   = 0.0;
  mutable Counters                counters_;
};

}  // namespace conway::geometry
