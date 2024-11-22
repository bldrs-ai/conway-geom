#pragma once

#include <stdint.h>

#include <glm/glm.hpp>

namespace conway::geometry {

  constexpr uint32_t MAX_MORTON_COMPONENT = 1023;

  inline uint32_t pack3( uint32_t x ) {

    x &= 0x09249249;
    x  = ( x ^ ( x >>  2) ) & 0x030c30c3;
    x  = ( x ^ ( x >>  4) ) & 0x0300f00f;
    x  = ( x ^ ( x >>  8) ) & 0xff0000ff;
    x  = ( x ^ ( x >> 16) ) & 0x000003ff;
  
    return x;
  }

  inline uint32_t unpack3( uint32_t x ) {

    x = ( x ^ ( x << 16 ) ) & 0x030000ff;
    x = ( x ^ ( x << 8  ) ) & 0x0300f00f;
    x = ( x ^ ( x << 4  ) ) & 0x030c30c3;
    x = ( x ^ ( x << 2  ) ) & 0x09249249;

    return x;
  }

  inline uint32_t morton3( uint32_t x, uint32_t y, uint32_t z ) {

    return unpack3( x ) | ( unpack3( y ) << 1 ) | ( unpack3( z ) << 2 );
  }

  inline uint32_t morton3_x( uint32_t code ) {

    return pack3( code );
  }

  inline uint32_t morton3_y( uint32_t code ) {

    return pack3( code >> 1 );
  }
   
  inline uint32_t morton3_z( uint32_t code ) {

    return pack3( code >> 2 );
  }

  inline glm::uvec3 unpacked_coord3( const glm::dvec3& value, const glm::dvec3& origin, double inverseStep ) {

    return glm::uvec3( glm::clamp( glm::round( ( value - origin ) * inverseStep ), 0.0, 1023.0 ) );
  }

  inline uint32_t pack( const glm::dvec3& value, const glm::dvec3& origin, double inverseStep ) {

    glm::uvec3 coord = unpacked_coord3( value, origin, inverseStep );

    return morton3( coord.x, coord.y, coord.z );
  }

  inline uint32_t pack( const glm::uvec3& coord ) {

    return morton3( coord.x, coord.y, coord.z );
  }
}