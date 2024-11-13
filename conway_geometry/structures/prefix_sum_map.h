#pragma once 

#include <vector>

namespace conway {

  struct PrefixSumMap {

    void reset() {

      counts.clear();
    }

    template < typename T, typename IdMappingFunction >
    void construct( const std::vector< T >& from, uint32_t idSize, IdMappingFunction idFunction ) {
      
      reset();

      counts.resize( idSize + 1, 0 );
      aggregate.resize( from.size(), 0 );

     static_assert(
        std::is_invocable_r_v< uint32_t, IdMappingFunction, T >,
        "Intersect requires a callable function that receives a uint32_t triangle indice." );

      for ( const T& item : from ) {

        ++counts[ idFunction( item ) ];
      }

      for ( uint32_t where = 1, end = idSize + 1; where < end; ++where ) {

        counts[ where ] += counts[ where - 1 ];
      }

      for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

        const T& item = from[ where ];
        uint32_t id   = idFunction( item );

        aggregate[ --counts[ id ] ] = where;
      }
    }

    std::span< const uint32_t > get( uint32_t id ) const {

      size_t offset = counts[ id ];
      size_t size   = counts[ id + 1 ] - offset;

      return std::span( aggregate.data() + offset, size );
    }

    std::vector< uint32_t > aggregate;
    std::vector< uint32_t > counts;
  };
}