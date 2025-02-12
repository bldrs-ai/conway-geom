#pragma once

#include <glm/glm.hpp>
#include <stdint.h>
#include <utility>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES 1
#endif

#include <math.h>

#if defined (_MSC_VER)

#pragma warning( push )
#pragma warning( disable : 26812 )

#endif

#include "predicates.h"

#if defined (_MSC_VER)

#pragma warning( pop )

#endif


namespace conway {
  
  constexpr size_t SECOND_AXIS_SHIFT = 2;
  constexpr size_t AXIS_MASK         = ( 1 << SECOND_AXIS_SHIFT ) - 1;
  constexpr size_t X_AXIS_INDEX      = 0;
  constexpr size_t Y_AXIS_INDEX      = 1;
  constexpr size_t Z_AXIS_INDEX      = 2;

  constexpr inline size_t make_axis_pair( size_t first, size_t second ) {

    return ( first ) | ( second << SECOND_AXIS_SHIFT );
  }

  enum class AxisPair : size_t {
    NONE = 0,
    X_Y = make_axis_pair( X_AXIS_INDEX, Y_AXIS_INDEX ),
    X_Z = make_axis_pair( X_AXIS_INDEX, Z_AXIS_INDEX ),
    Y_Z = make_axis_pair( Y_AXIS_INDEX, Z_AXIS_INDEX )

  };

  constexpr inline glm::length_t first_axis( AxisPair from ) {
  
    return static_cast< glm::length_t >( from ) & AXIS_MASK;
  }

  constexpr inline glm::length_t second_axis( AxisPair from ) {
  
    return ( static_cast< glm::length_t >( from ) >> SECOND_AXIS_SHIFT );
  }

  inline uint32_t largest_component( const glm::dvec3& of ) {

    uint32_t result = 0;

    double candidateValue = of.x;

    if (of.y > candidateValue) {

      result = 1;
      candidateValue = of.x;
    }

    if (of.z > candidateValue) {

      result = 2;
      candidateValue = of.z;
    }

    return result;
  }

  
  inline glm::dvec2 extract( const glm::dvec3& from, AxisPair axes ) {

    return glm::dvec2( from[ first_axis( axes ) ], from[ second_axis( axes ) ] );
  }

  inline std::pair< double, double > extract_pair( const glm::dvec3& from, AxisPair axes ) {

    return std::make_pair( from[ first_axis( axes ) ], from[ second_axis( axes ) ] );
  }

  inline double length2(const glm::dvec3& v) {

    return glm::dot(v, v);
  }

  inline bool same_point(const glm::dvec3& v0, const glm::dvec3& v1, double tolerance = 0) {

    glm::dvec3 comparison = glm::abs( v0 - v1 );

    return tolerance >= comparison.x && tolerance >= comparison.y && tolerance >= comparison.z;
  }
  
  namespace {

    // See here for this approximation https://mazzo.li/posts/vectorized-atan2.html#atan2-primer
    inline constexpr double fast_atan_scalar_approximation( double x ) {

      constexpr double A1  = 0.99997726;
      constexpr double A3  = -0.33262347;
      constexpr double A5  = 0.19354346;
      constexpr double A7  = -0.11643287;
      constexpr double A9  = 0.05265332;
      constexpr double A11 = -0.01172120;

      double x2 = x * x;

      return x * ( A1 + x2 * ( A3 + x2 * ( A5 + x2 * ( A7 + x2 * ( A9 + x2 * A11 ) ) ) ) );
    }
  }

  /** Fast approximation of atan2 */
  inline double fast_atan2( double y, double x ) {

    double absY = std::abs( y );
    double absX = std::abs( x );

    double result;

    if ( absX < absY ) {

      double atanInput = x / y;

      result = ( atanInput >= 0.0 ? M_PI_2 : -M_PI_2 ) -
        fast_atan_scalar_approximation( atanInput );
    }
    else {

      double atanInput = y / x;

      result = fast_atan_scalar_approximation( atanInput );
    }

    if ( x < 0 && y >= 0 ) {
      result += M_PI; 
    }
    else if ( x < 0 && y < 0 ) {
      result -= M_PI;
    }

    return result;
  }


  /**
   * Determinant for 3 vectors (matrix essentially) using a basic non-exact formula
   */
  inline double determinant3x3( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2 ) {

    return
      (v0.x * ((v1.y * v2.z) - (v1.z * v2.y))) -
      (v1.x * ((v0.y * v2.z) - (v0.z * v2.y))) +
      (v2.x * ((v0.y * v1.z) - (v0.z * v1.y)));
  }

  /** Get the value of the determinant for a pair of axes for a 3D triangle (i.e. a 2D triangle with an axis aligned ortho projection */
  inline double orient2D( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, AxisPair axes ) {

    glm::dvec2 v0t = extract( v0, axes );
    glm::dvec2 v1t = extract( v1, axes );
    glm::dvec2 v2t = extract( v2, axes );

    return predicates::adaptive::orient2d(
      &(v0t.x),
      &(v1t.x),
      &(v2t.x) );
  }

  /** Get the sign of the determinant for a pair of axes for a 3D triangle (i.e. a 2D triangle with an axis aligned ortho projection */
  inline int32_t orient2D(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, AxisPair axes, double tolerance ) {

    double orientation = orient2D( v0, v1, v2, axes );
    
    if ( orientation > tolerance ) {

      return 1;
    }

    if ( orientation < tolerance ) {

      return -1;
    }

    return 0;
  }

  /** Get the sign of the determinant for a pair of axes for a 3D triangle (i.e. a 2D triangle with an axis aligned ortho projection) */
  inline int32_t orient2D( const glm::dvec3(&vertices)[ 3 ], AxisPair axes, double tolerance = 0 ) {

    double determinant = orient2D( vertices[ 0 ], vertices[ 1 ], vertices[ 2 ], axes );

    if (determinant > tolerance) {

      return 1;
    }

    if (determinant < tolerance) {

      return -1;
    }

    return 0;
  }

  /** Is a 3D triangle zero area, computed with exact predicates, which means that for the tolerance == 0 case, it will give an exact answer */
  inline bool is_zero_area_triangle( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, double tolerance = 0 ) {

    double doubleTolerance = 2.0 * tolerance;

    return
      orient2D( v0, v1, v2, AxisPair::X_Y, doubleTolerance ) == 0 &&
      orient2D( v0, v1, v2, AxisPair::X_Z, doubleTolerance ) == 0 &&
      orient2D( v0, v1, v2, AxisPair::Y_Z, doubleTolerance ) == 0;
  } 

  /** Assuming 2 intersecting line segments, this gets the intersection point */
  inline glm::dvec3 line_segment_line_segment_intersection(
    const glm::dvec3& a0,
    const glm::dvec3& a1,
    const glm::dvec3& b0,
    const glm::dvec3& b1 ) {

    glm::dvec3 e0        = a1 - a0;
    glm::dvec3 direction = b1 - b0;

    // intersecting non-colinear line segments are 
    // coplanar, and have a well defined normal.
    glm::dvec3 e1     = glm::cross( e0, direction );
    glm::dvec3 origin = b0 - a0;
    glm::dvec3 normal = glm::cross( e0, e1 );

    double t = dot( origin, normal ) / -dot( direction, normal );

    // printf( "%f %f %f %f %d %d %d %d\n", t, predicates::adaptive::orient3d( &a0.x, &a1.x, &b0.x, &b1.x ), glm::length( normal ), glm::length( direction ), is_zero_area_triangle( a0, a1, b0 ) ? 1 : 0, is_zero_area_triangle( a0, a1, b1 ) ? 1 : 0, is_zero_area_triangle( a0, b0, b1 ) ? 1 : 0, is_zero_area_triangle( a1, b0, b1 ) ? 1 : 0 );

    return b0 + direction * t;
  }

  
  /** Calculate the signed solid angle of a triangle relative a point */
  inline double triangle_solid_angle( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& from ) {

    glm::dvec3 r0 = v0 - from;
    glm::dvec3 r1 = v1 - from;
    glm::dvec3 r2 = v2 - from;

    double r0l = glm::length( r0 );
    double r1l = glm::length( r1 );
    double r2l = glm::length( r2 );

    double x =
      r0l * r1l * r2l + 
      glm::dot( r0, r1 ) * r2l +
      glm::dot( r1, r2 ) * r0l +
      glm::dot( r2, r0 ) * r1l;

    // We use the exactly signed version of the determinant here because we really want the sign to be right,
    // as this gives us our accurate face winding up close for inside vs outside.
    // Because the dipoles in the BVH are far away, we don't need the precision, but here we do.
    double y = determinant3x3( r0, r1, r2 ); //predicates::adaptive::orient3d< double >( &v0.x, &v1.x, &v2.x, &from.x );

    return 2.0 * fast_atan2( y, x );
  }

  
  /**
   * Simple determininant based comparison given 4 points that will calculate 3 vectors relative to the first and then
   * work out the sign of the determinant within tolerance (-1 for negative, 0 for within tolerance of 0 and 1 for positive)
   */
  inline int32_t orient3D( const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3, double tolerance = 0 ) {

    double result = predicates::adaptive::orient3d< double >( &v0.x, &v1.x, &v2.x, &v3.x );

    if ( result > tolerance ) {

      return 1;
    }

    if ( result < tolerance ) {

      return -1;
    }

    return 0;
  }

  namespace {

    constexpr double ONE_THIRD = 1.0 / 3.0;
  }

}