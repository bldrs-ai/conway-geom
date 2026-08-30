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
