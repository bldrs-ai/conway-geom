#include <stdint.h>
#include <utility>
#include <bit>

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
}

namespace std {

  template <>
  struct hash< std::pair< uint32_t, uint32_t > > {

    size_t operator()( const std::pair< uint32_t, uint32_t >& value ) const {

      return conway::hash_mix( conway::hash( value.first ), conway::hash( value.second ) );
    }
  };
}
