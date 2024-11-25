#pragma once

#include <glm/glm.hpp>
#include <stdint.h>
#include <utility>

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

    glm::dvec3 comparison = glm::abs(v0 - v1);

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
   * Determinant for 3 vectors (matrix essentially)
   */
  inline double determinant3x3(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2) {

    return
      (v0.x * ((v1.y * v2.z) - (v1.z * v2.y))) -
      (v1.x * ((v0.y * v2.z) - (v0.z * v2.y))) +
      (v2.x * ((v0.y * v1.z) - (v0.z * v1.y)));
  }

}