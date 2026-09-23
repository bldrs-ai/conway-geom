#pragma once

#include <tinynurbs/tinynurbs.h>

#include <cmath>
#include <glm/glm.hpp>
#include <limits>
#include <tuple>
#include <vector>

namespace conway::geometry {

/**
 * Maximum B-spline degree supported by the stack-allocated basis buffers in
 * RationalSurfaceEvaluator. Degrees above this fall back to tinynurbs.
 */
constexpr uint32_t NURBS_MAX_STACK_DEGREE = 15;

/**
 * Allocation-free evaluator for rational B-spline surfaces.
 *
 * tinynurbs::surfacePoint / surfaceTangent on a RationalSurface re-convert
 * the ENTIRE control grid to homogeneous coordinates (heap allocating) on
 * every sample, which dominated STEP advanced-BREP tessellation cost. This
 * evaluator performs that conversion once, and each sample then costs
 * O((degree_u + 1) * (degree_v + 1)) with no heap allocation.
 *
 * The accumulation order of the basis/point/derivative loops deliberately
 * mirrors tinynurbs so results stay numerically identical to the previous
 * code path.
 */
struct RationalSurfaceEvaluator {

  explicit RationalSurfaceEvaluator( const tinynurbs::RationalSurface3d& srf )
      : surface_( srf ) {

    fastPath_ =
        srf.degree_u <= NURBS_MAX_STACK_DEGREE &&
        srf.degree_v <= NURBS_MAX_STACK_DEGREE &&
        srf.control_points.rows() > 0 &&
        srf.control_points.cols() > 0;

    if ( !fastPath_ ) {
      return;
    }

    size_t rows = srf.control_points.rows();

    cols_ = srf.control_points.cols();

    homogeneous_.resize( rows * cols_ );

    for ( size_t i = 0; i < rows; ++i ) {
      for ( size_t j = 0; j < cols_; ++j ) {

        const glm::dvec3& point  = srf.control_points( i, j );
        double            weight = srf.weights( i, j );

        homogeneous_[ i * cols_ + j ] =
            glm::dvec4( point * weight, weight );

        polynomial_ = polynomial_ && ( weight == 1.0 );

        // The magnitudes an evaluation's rounding error is actually governed
        // by - see maxHomogeneousNorm(). Taken over the WHOLE grid rather
        // than the local window: a window-wise maximum would be tighter, but
        // this is read once per candidate and the difference is nowhere near
        // the margin involved.
        maxHomogeneousNorm_ =
            std::max( maxHomogeneousNorm_,
                      glm::length( glm::dvec3( point * weight ) ) );

        maxWeight_ = std::max( maxWeight_, std::abs( weight ) );
      }
    }
  }

  /**
   * Knot-span lookup (The NURBS Book A2.1). tinynurbs::findSpan guards the
   * domain ends with an ABSOLUTE epsilon (numeric_limits::epsilon), which is
   * below one ULP for knot values > 2 - e.g. STEP surfaces whose parameter
   * range is a real length like [0, 200]. There `lastKnot - epsilon` rounds
   * back to lastKnot, the guard never fires for u == lastKnot, and its
   * binary search (`low = mid` with mid pinned by integer division) spins
   * forever. Exact >= / <= boundary guards make the search invariant
   * (knots[low] <= u < knots[high]) hold strictly, so it always terminates.
   */
  static int findSpan(
      uint32_t degree,
      const std::vector< double >& knots,
      double u ) {

    // Index of the last knot span start (n = knotCount - degree - 2).
    int n = static_cast< int >( knots.size() ) - static_cast< int >( degree ) - 2;

    // Degenerate knot vector (fewer than degree + 2 knots) - callers
    // validate with surfaceIsValid first so this shouldn't happen, but an
    // out-of-bounds read in wasm is silent garbage, so fail safe.
    if ( n < 0 ) {
      return static_cast< int >( degree );
    }

    if ( u >= knots[ n + 1 ] ) {
      return n;
    }

    if ( u <= knots[ degree ] ) {
      return degree;
    }

    int low  = degree;
    int high = n + 1;
    int mid  = ( low + high ) / 2;

    while ( u < knots[ mid ] || u >= knots[ mid + 1 ] ) {

      if ( u < knots[ mid ] ) {
        high = mid;
      } else {
        low = mid;
      }

      mid = ( low + high ) / 2;
    }

    return mid;
  }

  /**
   * Clamp a parameter strictly below the last knot so the tinynurbs
   * fallback (degree > NURBS_MAX_STACK_DEGREE) can never enter
   * tinynurbs::findSpan's non-terminating end-of-domain case (see
   * findSpan above).
   */
  static double clampBelowLastKnot(
      uint32_t degree, const std::vector< double >& knots, double u ) {

    if ( knots.size() < degree + 1 ) {
      return u;
    }

    double lastKnot = knots[ knots.size() - degree - 1 ];

    return u >= lastKnot ?
      std::nextafter( lastKnot, -std::numeric_limits< double >::infinity() ) : u;
  }

  /**
   * Outward surface normal at (u, v), from the first partial derivatives.
   *
   * No fast path: this is called once per triangle corner at emission time,
   * not inside the refinement loop, so it goes straight to tinynurbs rather
   * than duplicating the stack-allocated basis evaluation. Returns a zero
   * vector at a degenerate point — a collapsed pole, where the two partials
   * are parallel and the cross product vanishes — which the caller reads as
   * "no analytic normal here" (bldrs-ai/conway#667).
   *
   * @param u Surface parameter u.
   * @param v Surface parameter v.
   * @return The normal, not normalized, or zero where it is undefined.
   */
  glm::dvec3 normal( double u, double v ) const {

    glm::dvec3 result = tinynurbs::surfaceNormal(
      surface_,
      clampBelowLastKnot( surface_.degree_u, surface_.knots_u, u ),
      clampBelowLastKnot( surface_.degree_v, surface_.knots_v, v ) );

    if ( !std::isfinite( result.x ) ||
         !std::isfinite( result.y ) ||
         !std::isfinite( result.z ) ) {

      return glm::dvec3( 0.0 );
    }

    return result;
  }

  /**
   * The surface's DEGREES and KNOTS, and whether the weights are all one.
   *
   * Exposed for the deflection certificate in
   * `operations/deflection_certificate.h`, which needs the degrees to know
   * what polynomial degree the surface restricted to a line has, the knots to
   * know where that polynomial CHANGES, and the weights to know whether it is
   * a polynomial at all. Nothing else reads them; the evaluation paths below
   * use `surface_` directly.
   */
  uint32_t degreeU() const { return surface_.degree_u; }

  uint32_t degreeV() const { return surface_.degree_v; }

  const std::vector< double >& knotsU() const { return surface_.knots_u; }

  const std::vector< double >& knotsV() const { return surface_.knots_v; }

  /** False on the tinynurbs fallback, where `pointHomogeneous` has no path. */
  bool supportsFastPath() const { return fastPath_; }

  /**
   * Largest homogeneous control value, and largest weight, over the grid.
   *
   * These, not the magnitude of the surface POINTS, are what the rounding
   * error of an evaluation is governed by: `point` accumulates
   * `basis * controlPointW` over ( degreeU + 1 ) * ( degreeV + 1 ) terms, and
   * the basis functions are a partition of unity, so the absolute error of
   * the sum scales with the largest term in it. A patch whose surface values
   * are small can still have large control values that cancel - how much
   * larger is bounded by the basis's own condition number, but it is not
   * bounded by the surface values, which is why the certificate's error term
   * reads these rather than the chord endpoints.
   */
  double maxHomogeneousNorm() const { return maxHomogeneousNorm_; }

  double maxWeight() const { return maxWeight_; }

  /**
   * Is `span` a span the basis can be evaluated on - i.e. does it have a
   * non-empty knot interval?
   *
   * The Cox-de-Boor denominators are `knots[ span + a ] - knots[ span + b ]`
   * with `a >= 1 >= b`, so they are all at least `knots[ span + 1 ] -
   * knots[ span ]`; a non-empty span therefore divides by nothing near zero,
   * INCLUDING for a parameter outside the span, which is what
   * `pointHomogeneousAtSpan` relies on.
   */
  static bool spanIsEvaluable(
      const std::vector< double >& knots, int span ) {

    return span >= 0 &&
           static_cast< size_t >( span ) + 1 < knots.size() &&
           knots[ span ] < knots[ span + 1 ];
  }

  /**
   * True when every weight is exactly one, so `point` is a POLYNOMIAL map of
   * the parameters rather than a quotient of two.
   *
   * Exact comparison, not a tolerance: a weight an ulp off one leaves the
   * surface rational, and reporting it polynomial would put a rational
   * function inside a polynomial bound. STEP's
   * `B_SPLINE_SURFACE_WITH_KNOTS` carries no weights and is built with
   * literal 1.0, so the common case answers true on the nose.
   */
  bool isPolynomial() const { return polynomial_; }

  /**
   * The surface point in HOMOGENEOUS form - ( x w, y w, z w, w ) - before the
   * perspective divide `point` ends with.
   *
   * Both halves are polynomial in the parameters where `point` is not, which
   * is what lets the certificate bound a rational surface: see the RATIONAL
   * SURFACES note in `deflection_certificate.h`.
   */
  glm::dvec4 pointHomogeneous( double u, double v ) const {

    return pointHomogeneousAtSpan(
        findSpan( surface_.degree_u, surface_.knots_u, u ),
        findSpan( surface_.degree_v, surface_.knots_v, v ),
        u,
        v );
  }

  /**
   * The homogeneous point, evaluated as the polynomial of the NAMED spans
   * rather than of the spans `u` and `v` fall in.
   *
   * This is what makes a ONE-SIDED reading possible. `findSpan` resolves a
   * parameter sitting exactly on an interior knot to the span on its RIGHT,
   * and at a knot of multiplicity `degree + 1` the two sides are different
   * polynomials with a step between them - so a piece whose last node lands
   * on such a knot would otherwise be sampled from the polynomial it is not
   * certifying. Naming the span pins every node of a piece to that piece's
   * own polynomial; the basis extends outside its knot interval without any
   * division by zero, see spanIsEvaluable.
   *
   * Callers must have checked `spanIsEvaluable` for both spans.
   */
  glm::dvec4 pointHomogeneousAtSpan(
      int spanU, int spanV, double u, double v ) const {

    // The homogeneous control grid only exists on the fast path. Rather than
    // read an empty vector, hand back a value the caller's finiteness check
    // rejects - the certificate then declines instead of bounding garbage.
    if ( !fastPath_ ) {
      return glm::dvec4( std::numeric_limits< double >::quiet_NaN() );
    }

    uint32_t degreeU = surface_.degree_u;
    uint32_t degreeV = surface_.degree_v;

    double basisU[ NURBS_MAX_STACK_DEGREE + 1 ];
    double basisV[ NURBS_MAX_STACK_DEGREE + 1 ];

    basis( degreeU, spanU, surface_.knots_u, u, basisU );
    basis( degreeV, spanV, surface_.knots_v, v, basisV );

    glm::dvec4 pointw( 0.0 );

    for ( uint32_t l = 0; l <= degreeV; ++l ) {

      glm::dvec4 temp( 0.0 );

      for ( uint32_t k = 0; k <= degreeU; ++k ) {

        temp +=
            basisU[ k ] *
            controlPointW( spanU - degreeU + k, spanV - degreeV + l );
      }

      pointw += basisV[ l ] * temp;
    }

    return pointw;
  }

  /**
   * Largest step between ADJACENT homogeneous control points in the window
   * that supports span ( spanU, spanV ), along u or along v.
   *
   * This is what a B-spline's derivative is built from: the derivative is
   * itself a B-spline whose control points are
   * `degree * ( P[ i + 1 ] - P[ i ] ) / ( knot difference )`, and its basis
   * is a partition of unity, so the derivative over this span is bounded by
   * the largest of them. Taken over the LOCAL window rather than the whole
   * grid, which matters: on a surface whose coordinates span 1e9, a global
   * maximum would put a bound of 1e9 on a span whose own control points are
   * all of order one, and the certificate would decline chords it can
   * perfectly well bound.
   *
   * The window is `degree + 1` control points wide, so there are `degree`
   * steps across it and every index stays inside the span's own support.
   */
  void controlSpread(
      int     spanU,
      int     spanV,
      bool    alongU,
      double& pointSpread,
      double& weightSpread ) const {

    pointSpread  = 0.0;
    weightSpread = 0.0;

    if ( !fastPath_ ) {
      pointSpread  = std::numeric_limits< double >::infinity();
      weightSpread = std::numeric_limits< double >::infinity();
      return;
    }

    const int rowBase = spanU - static_cast< int >( surface_.degree_u );
    const int colBase = spanV - static_cast< int >( surface_.degree_v );

    const uint32_t lastRow =
        alongU ?
          ( surface_.degree_u > 0 ? surface_.degree_u - 1 : 0 ) :
          surface_.degree_u;

    const uint32_t lastCol =
        alongU ?
          surface_.degree_v :
          ( surface_.degree_v > 0 ? surface_.degree_v - 1 : 0 );

    if ( ( alongU && surface_.degree_u == 0 ) ||
         ( !alongU && surface_.degree_v == 0 ) ) {
      return;
    }

    for ( uint32_t a = 0; a <= lastRow; ++a ) {

      for ( uint32_t b = 0; b <= lastCol; ++b ) {

        const glm::dvec4& here =
            controlPointW( rowBase + a, colBase + b );

        const glm::dvec4& next =
            controlPointW( rowBase + a + ( alongU ? 1 : 0 ),
                           colBase + b + ( alongU ? 0 : 1 ) );

        pointSpread =
            std::max( pointSpread,
                      glm::length( glm::dvec3( next ) - glm::dvec3( here ) ) );

        weightSpread =
            std::max( weightSpread, std::abs( next.w - here.w ) );
      }
    }
  }

  /** `pointHomogeneousAtSpan` with the perspective divide applied. */
  glm::dvec3 pointAtSpan( int spanU, int spanV, double u, double v ) const {

    const glm::dvec4 homogeneous =
        pointHomogeneousAtSpan( spanU, spanV, u, v );

    return glm::dvec3( homogeneous ) / homogeneous.w;
  }

  /** Point on the surface, matching tinynurbs::surfacePoint( rational ). */
  glm::dvec3 point( double u, double v ) const {

    if ( !fastPath_ ) {
      return tinynurbs::surfacePoint(
        surface_,
        clampBelowLastKnot( surface_.degree_u, surface_.knots_u, u ),
        clampBelowLastKnot( surface_.degree_v, surface_.knots_v, v ) );
    }

    uint32_t degreeU = surface_.degree_u;
    uint32_t degreeV = surface_.degree_v;

    int spanU = findSpan( degreeU, surface_.knots_u, u );
    int spanV = findSpan( degreeV, surface_.knots_v, v );

    double basisU[ NURBS_MAX_STACK_DEGREE + 1 ];
    double basisV[ NURBS_MAX_STACK_DEGREE + 1 ];

    basis( degreeU, spanU, surface_.knots_u, u, basisU );
    basis( degreeV, spanV, surface_.knots_v, v, basisV );

    glm::dvec4 pointw( 0.0 );

    for ( uint32_t l = 0; l <= degreeV; ++l ) {

      glm::dvec4 temp( 0.0 );

      for ( uint32_t k = 0; k <= degreeU; ++k ) {

        temp +=
            basisU[ k ] *
            controlPointW( spanU - degreeU + k, spanV - degreeV + l );
      }

      pointw += basisV[ l ] * temp;
    }

    return glm::dvec3( pointw ) / pointw.w;
  }

  /**
   * Unit surface tangents along u and v, matching
   * tinynurbs::surfaceTangent( rational ).
   */
  std::tuple< glm::dvec3, glm::dvec3 > tangent( double u, double v ) const {

    if ( !fastPath_ ) {
      return tinynurbs::surfaceTangent(
        surface_,
        clampBelowLastKnot( surface_.degree_u, surface_.knots_u, u ),
        clampBelowLastKnot( surface_.degree_v, surface_.knots_v, v ) );
    }

    uint32_t degreeU = surface_.degree_u;
    uint32_t degreeV = surface_.degree_v;

    int spanU = findSpan( degreeU, surface_.knots_u, u );
    int spanV = findSpan( degreeV, surface_.knots_v, v );

    // ders[ 0 ] = basis values, ders[ 1 ] = first derivatives.
    double dersU[ 2 ][ NURBS_MAX_STACK_DEGREE + 1 ];
    double dersV[ 2 ][ NURBS_MAX_STACK_DEGREE + 1 ];

    derBasis( degreeU, spanU, surface_.knots_u, u, dersU );
    derBasis( degreeV, spanV, surface_.knots_v, v, dersV );

    // Homogeneous surface derivatives, mirroring
    // tinynurbs::internal::surfaceDerivatives with num_ders = 1.
    glm::dvec4 homoDers[ 2 ][ 2 ] = {
        { glm::dvec4( 0.0 ), glm::dvec4( 0.0 ) },
        { glm::dvec4( 0.0 ), glm::dvec4( 0.0 ) } };

    uint32_t du = std::min( 1u, degreeU );
    uint32_t dv = std::min( 1u, degreeV );

    glm::dvec4 temp[ NURBS_MAX_STACK_DEGREE + 1 ];

    for ( uint32_t k = 0; k <= du; ++k ) {

      for ( uint32_t s = 0; s <= degreeV; ++s ) {

        temp[ s ] = glm::dvec4( 0.0 );

        for ( uint32_t r = 0; r <= degreeU; ++r ) {

          temp[ s ] +=
              dersU[ k ][ r ] *
              controlPointW( spanU - degreeU + r, spanV - degreeV + s );
        }
      }

      uint32_t dd = std::min( 1u - k, dv );

      for ( uint32_t l = 0; l <= dd; ++l ) {

        for ( uint32_t s = 0; s <= degreeV; ++s ) {

          homoDers[ k ][ l ] += dersV[ l ][ s ] * temp[ s ];
        }
      }
    }

    // Rational correction (NURBS book eq. 4.20 truncated to first
    // derivatives), matching tinynurbs::surfaceDerivatives( rational ).
    // Note: tinynurbs multiplies by the reciprocal (der *= 1 / w) rather
    // than dividing; mirror that exactly so results stay bit-identical.
    double     wInv  = 1.0 / homoDers[ 0 ][ 0 ].w;
    glm::dvec3 s     = glm::dvec3( homoDers[ 0 ][ 0 ] ) * wInv;
    glm::dvec3 du3   = ( glm::dvec3( homoDers[ 1 ][ 0 ] ) -
                         homoDers[ 1 ][ 0 ].w * s ) * wInv;
    glm::dvec3 dv3   = ( glm::dvec3( homoDers[ 0 ][ 1 ] ) -
                         homoDers[ 0 ][ 1 ].w * s ) * wInv;

    double duLen = glm::length( du3 );
    double dvLen = glm::length( dv3 );

    if ( !tinynurbs::util::close( duLen, 0.0 ) ) {
      du3 /= duLen;
    }

    if ( !tinynurbs::util::close( dvLen, 0.0 ) ) {
      dv3 /= dvLen;
    }

    return std::make_tuple( du3, dv3 );
  }

 private:

  const glm::dvec4& controlPointW( size_t i, size_t j ) const {
    return homogeneous_[ i * cols_ + j ];
  }

  /** All weights exactly 1 - see isPolynomial(). */
  bool polynomial_ = true;

  /** Control-grid magnitudes - see maxHomogeneousNorm(). */
  double maxHomogeneousNorm_ = 0.0;
  double maxWeight_          = 0.0;

  /** Cox-de-Boor basis, mirroring tinynurbs::bsplineBasis. */
  static void basis(
      uint32_t degree,
      int span,
      const std::vector< double >& knots,
      double u,
      double* result ) {

    double left[ NURBS_MAX_STACK_DEGREE + 1 ]  = { 0.0 };
    double right[ NURBS_MAX_STACK_DEGREE + 1 ] = { 0.0 };

    double saved = 0.0;
    double temp  = 0.0;

    result[ 0 ] = 1.0;

    for ( int j = 1; j <= static_cast< int >( degree ); ++j ) {

      left[ j ]  = u - knots[ span + 1 - j ];
      right[ j ] = knots[ span + j ] - u;
      saved      = 0.0;

      for ( int r = 0; r < j; ++r ) {

        temp        = result[ r ] / ( right[ r + 1 ] + left[ j - r ] );
        result[ r ] = saved + right[ r + 1 ] * temp;
        saved       = left[ j - r ] * temp;
      }

      result[ j ] = saved;
    }
  }

  /**
   * Basis values + first derivatives, mirroring tinynurbs::bsplineDerBasis
   * with num_ders = 1.
   */
  static void derBasis(
      uint32_t degree,
      int span,
      const std::vector< double >& knots,
      double u,
      double ( &ders )[ 2 ][ NURBS_MAX_STACK_DEGREE + 1 ] ) {

    double left[ NURBS_MAX_STACK_DEGREE + 1 ]  = { 0.0 };
    double right[ NURBS_MAX_STACK_DEGREE + 1 ] = { 0.0 };

    double saved = 0.0;
    double temp  = 0.0;

    double ndu[ NURBS_MAX_STACK_DEGREE + 1 ][ NURBS_MAX_STACK_DEGREE + 1 ];

    ndu[ 0 ][ 0 ] = 1.0;

    for ( int j = 1; j <= static_cast< int >( degree ); ++j ) {

      left[ j ]  = u - knots[ span + 1 - j ];
      right[ j ] = knots[ span + j ] - u;
      saved      = 0.0;

      for ( int r = 0; r < j; ++r ) {

        ndu[ j ][ r ] = right[ r + 1 ] + left[ j - r ];
        temp          = ndu[ r ][ j - 1 ] / ndu[ j ][ r ];

        ndu[ r ][ j ] = saved + right[ r + 1 ] * temp;
        saved         = left[ j - r ] * temp;
      }

      ndu[ j ][ j ] = saved;
    }

    for ( int j = 0; j <= static_cast< int >( degree ); ++j ) {
      ders[ 0 ][ j ] = ndu[ j ][ degree ];
      ders[ 1 ][ j ] = 0.0;
    }

    // First derivative only (num_ders = 1 specialisation of the
    // triangular-table algorithm).
    double a[ 2 ][ NURBS_MAX_STACK_DEGREE + 1 ];

    for ( int r = 0; r <= static_cast< int >( degree ); ++r ) {

      int s1 = 0;
      int s2 = 1;

      a[ 0 ][ 0 ] = 1.0;

      constexpr int k  = 1;
      double        d  = 0.0;
      int           rk = r - k;
      int           pk = static_cast< int >( degree ) - k;
      int           j1 = 0;
      int           j2 = 0;

      if ( r >= k ) {
        a[ s2 ][ 0 ] = a[ s1 ][ 0 ] / ndu[ pk + 1 ][ rk ];
        d            = a[ s2 ][ 0 ] * ndu[ rk ][ pk ];
      }

      j1 = ( rk >= -1 ) ? 1 : -rk;
      j2 = ( r - 1 <= pk ) ? k - 1 : static_cast< int >( degree ) - r;

      for ( int j = j1; j <= j2; ++j ) {
        a[ s2 ][ j ] =
            ( a[ s1 ][ j ] - a[ s1 ][ j - 1 ] ) / ndu[ pk + 1 ][ rk + j ];
        d += a[ s2 ][ j ] * ndu[ rk + j ][ pk ];
      }

      if ( r <= pk ) {
        a[ s2 ][ k ] = -a[ s1 ][ k - 1 ] / ndu[ pk + 1 ][ r ];
        d += a[ s2 ][ k ] * ndu[ r ][ pk ];
      }

      ders[ 1 ][ r ] = d * static_cast< double >( degree );
    }
  }

  const tinynurbs::RationalSurface3d& surface_;

  std::vector< glm::dvec4 > homogeneous_;

  size_t cols_     = 0;
  bool   fastPath_ = false;
};

}  // namespace conway::geometry
