#pragma once

#include <stdint.h>
#include <utility>
#include <bit>

#include <glm/glm.hpp>

namespace conway {

  inline uint32_t hash( uint32_t x ) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
  }

  
  inline size_t hash_mix( size_t left, size_t right ) {

    return std::rotl( left, 1 ) ^ std::rotr( right, 15 );
  } 

  inline size_t hash_mix( size_t a, size_t b, size_t c ) {

    return hash_mix( a, hash_mix( b, c ) );
  } 
}

namespace std {

  template <>
  struct hash< std::pair< uint32_t, uint32_t > > {

    size_t operator()( const std::pair< uint32_t, uint32_t >& value ) const {

      return conway::hash_mix( conway::hash( value.first ), conway::hash( value.second ) );
    }
  };

  template <>
  struct hash< glm::dvec3 > {

    size_t operator()( const glm::dvec3& value ) const {

      std::hash< double > doubleHasher;

      return conway::hash_mix( doubleHasher( value.x ), doubleHasher( value.y ), doubleHasher( value.z ) );
    }
  };

  template <>
  struct hash< glm::dvec2 > {

    size_t operator()( const glm::dvec2& value ) const {

      std::hash< double > doubleHasher;

      return conway::hash_mix( doubleHasher( value.x ), doubleHasher( value.y ) );
    }
  };

}
