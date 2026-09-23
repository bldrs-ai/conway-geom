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
    uint64_t unplaceable  = 0;
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

    // The chord's own magnitude, which is one of the things the error term
    // is built from. It is NOT on its own a scale for the evaluation error -
    // see maxHomogeneousNorm() in nurbs_utils.h, and the ERROR PROPAGATION
    // note on boundPiece.
    const double chordScale =
      std::max( glm::length( surface0 ), glm::length( surface1 ) );

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
    glm::dvec2 from  = glm::dvec2( 0.0 );
    glm::dvec2 to    = glm::dvec2( 0.0 );
    double     tFrom = 0.0;
    double     tTo   = 0.0;
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

    double shift =
      periodic_ ?
        ( std::floor( ( uv0.x - stripUMin_ ) / stripPeriod_ ) *
          stripPeriod_ ) :
        0.0;

    glm::dvec2 from( uv0.x - shift, uv0.y );

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

      // And where the periodic strip ends, in RAW u.
      const double edgeStrip =
        upU ? ( stripUMin_ + shift + stripPeriod_ )
            : ( stripUMin_ + shift );

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

      // ONE AXIS PER STEP, even when two boundaries land on the same
      // parameter. Advancing both at once would skip the box between them -
      // and two boundaries tie here whenever their parameters round
      // together, which is the same rounding that produced the third
      // finding. Stepping one at a time emits that box with its own span
      // indices and its own knot-derived corners; it is thin, and it is
      // bounded rather than assumed away.
      const bool endsOnStrip = periodic_ && ( du != 0.0 ) && ( tStrip <= tTo );
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

      pieces[ count++ ] = Piece{ spanU, spanV, from, to, tFrom, tTo };

      if ( endsChord ) {
        return true;
      }

      from  = to;
      tFrom = tTo;

      if ( endsOnStrip ) {

        // The strip's edge is not a cut of its own - it is where the u walk
        // RESTARTS, at the far end of the knot domain, one period along.
        shift += upU ? stripPeriod_ : -stripPeriod_;

        from.x = ( uv0.x + ( tTo * du ) ) - shift;

        spanU = directionalSpan( degreeU, knotsU, from.x, upU );

        if ( spanU < 0 ) {
          return false;
        }

      } else if ( endsOnU ) {

        do {
          spanU += upU ? 1 : -1;
        } while ( spanU >= firstU && spanU <= lastU &&
                  !RationalSurfaceEvaluator::spanIsEvaluable(
                    knotsU, spanU ) );

        // Off the end of the knot domain with no strip to wrap into: the
        // surface clamps there, so stay in the boundary span and stop
        // looking for more u edges.
        if ( spanU < firstU || spanU > lastU ) {
          spanU   = ( spanU < firstU ) ? firstU : lastU;
          clampedU = true;
        }
      }

      else if ( endsOnV ) {

        do {
          spanV += upV ? 1 : -1;
        } while ( spanV >= firstV && spanV <= lastV &&
                  !RationalSurfaceEvaluator::spanIsEvaluable(
                    knotsV, spanV ) );

        if ( spanV < firstV || spanV > lastV ) {
          spanV    = ( spanV < firstV ) ? firstV : lastV;
          clampedV = true;
        }
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
    const double roundingU =
      CERTIFICATE_ERROR_SAFETY * epsilon *
      std::max( std::abs( piece.from.x ), std::abs( piece.to.x ) );

    const double roundingV =
      CERTIFICATE_ERROR_SAFETY * epsilon *
      std::max( std::abs( piece.from.y ), std::abs( piece.to.y ) );

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

    const double weightCoefficientError =
      rational_ ?
        ( ( interpolation_.normInfinity() * weightError ) +
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
  double                          stripUMin_   = 0.0;
  double                          stripPeriod_ = 0.0;
  double                          tolerance_   = 0.0;
  mutable Counters                counters_;
};

}  // namespace conway::geometry
