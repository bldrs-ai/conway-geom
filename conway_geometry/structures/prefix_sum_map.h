#pragma once 

#include <vector>
#include <

namespace conway {

  struct PrefixSumMap {

    void reset() {

      counts.clear();
    }

    template < typename T, typename IdMappingFunction >
    void construct( const std::vector< T >& from, uint32_t idSize, IdMappingFunction idFunction ) {
      
      reset();

      counts.resize( idSize + 1, 0 );

     static_assert(
        std::is_invocable_r_v< uint32_t, IdMappingFunction, T >,
        "Intersect requires a callable function that receives a uint32_t triangle indice." );

      for ( const T& item : from ) {

        ++counts[ idFunction(  item ) ];
      }

      for ( uint32_t where = 1; where < idSize; ++where ) {

        counts[ where ] += counts[ where - 1 ];
      }

      for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

        const T& item = from[ where ];
        uint32_t id   = idFunction( item );

        aggregate[ --counts[ id ] ] = where;
      }
    }

    std::vector< uint32_t > aggregate;
    std::vector< uint32_t > counts;
  };
}