#pragma once 

#include <vector>

namespace conway {

  /** A dense index multi-map based on prefix sum bucketing,
   * where the index is relatively dense (keys -> secondary keys),
   * this is much faster and more cache efficient than a hash,
   * while maintaining ordering.
   */
  struct PrefixSumMap {

    void reset() {

      counts.clear();
      aggregate.clear();
    }

    /**
     * Sorts the values within each bucket.
     */
    template < typename ComparisonFunction >
    void bucketSort( ComparisonFunction comparison ) {

      static_assert(
        std::is_invocable_r_v< bool, ComparisonFunction, uint32_t, uint32_t >,
        "Bucket sort requires a comparison function that takes 2 indices from the constructed from vector" );

      for ( size_t where = 0, end = counts.size() - 1; where < end; ++where ) {

        auto bucketBegin = aggregate.begin() + counts[ where ];
        auto bucketEnd = aggregate.begin() + counts[ where + 1 ];
      
        std::sort( bucketBegin, bucketEnd, comparison );
      }
    }

    /** Construct a prefix sum map with enumeration, allowing a fixed number of id mappings per input T,
     *  this is good where you have say, multiple vertices in a triangle and want to create a map vertex->triangle.
     */
    template < typename T, typename IdMappingFunction >
    void construct( const std::vector< T >& from, uint32_t idSize, IdMappingFunction idFunction, uint32_t enumerationSize ) {

      reset();

      counts.resize( idSize + 1, 0 );

     static_assert(
        std::is_invocable_r_v< uint32_t, IdMappingFunction, T, uint32_t > ||
        std::is_invocable_r_v< std::pair< uint32_t, uint32_t >, IdMappingFunction, T, uint32_t >,
        "Construct requires a matching identity mapping function." );

      if constexpr ( std::is_invocable_r_v< uint32_t, IdMappingFunction, T, uint32_t > ) {
        
        aggregate.resize( from.size() * enumerationSize, 0 );

        for ( const T& item : from ) {

          for ( uint32_t keyIndex = 0; keyIndex < enumerationSize; ++keyIndex ) {

            uint32_t idValue = idFunction( item, keyIndex );

            assert( idValue < idSize + 1 );

            ++counts[ idValue ];
          }
        }
      } else {

         aggregate.resize( from.size() * 2 * enumerationSize, 0 );

         for ( const T& item : from ) {

          for ( uint32_t keyIndex = 0; keyIndex < enumerationSize; ++keyIndex ) {

            auto [idValue0, idValue1] = idFunction( item, keyIndex );

            assert( idValue0 < idSize + 1 );
            assert( idValue1 < idSize + 1 );

            ++counts[ idValue0 ];
            ++counts[ idValue1 ];
          }
        }
      }

      for ( uint32_t where = 1, end = idSize + 1; where < end; ++where ) {

        counts[ where ] += counts[ where - 1 ];
      }

      if constexpr ( std::is_invocable_r_v< uint32_t, IdMappingFunction, T, uint32_t > ) {

        for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

          for ( uint32_t keyIndex = 0; keyIndex < enumerationSize; ++keyIndex ) {
          
            const T& item = from[ where ];
            uint32_t id   = idFunction( item, keyIndex );

            aggregate[ --counts[ id ] ] = where;
          }
        }

      } else {

        for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

          for ( uint32_t keyIndex = 0; keyIndex < enumerationSize; ++keyIndex ) {

            const T& item     = from[ where ];
            auto [ id0, id1 ] = idFunction( item, keyIndex );

            aggregate[ --counts[ id0 ] ] = where;
            aggregate[ --counts[ id1 ] ] = where;
          }
        }
      }
    }

    template < typename T, typename IdMappingFunction >
    void construct( const std::vector< T >& from, uint32_t idSize, IdMappingFunction idFunction ) {
      
      reset();

      counts.resize( idSize + 1, 0 );

     static_assert(
        std::is_invocable_r_v< uint32_t, IdMappingFunction, T > || std::is_invocable_r_v< std::pair< uint32_t, uint32_t >, IdMappingFunction, T >,
        "Construct requires a matching identity mapping function." );

      if constexpr ( std::is_invocable_r_v< uint32_t, IdMappingFunction, T > ) {
        
        aggregate.resize( from.size(), 0 );

        for ( const T& item : from ) {

          uint32_t idValue = idFunction( item );

          assert( idValue < idSize + 1 );

          ++counts[ idValue ];
        }

      } else {

          aggregate.resize( from.size() * 2, 0 );

         for ( const T& item : from ) {

          auto [idValue0, idValue1] = idFunction( item );

          assert( idValue0 < idSize + 1 );
          assert( idValue1 < idSize + 1 );

          ++counts[ idValue0 ];
          ++counts[ idValue1 ];
        }

      }

      for ( uint32_t where = 1, end = idSize + 1; where < end; ++where ) {

        counts[ where ] += counts[ where - 1 ];
      }


      if constexpr ( std::is_invocable_r_v< uint32_t, IdMappingFunction, T > ) {

        for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

          const T& item = from[ where ];
          uint32_t id   = idFunction( item );

          aggregate[ --counts[ id ] ] = where;
        }

      } else {

        for ( uint32_t where = 0, end = static_cast< uint32_t >( from.size() ); where < end; ++where ) {

          const T& item   = from[ where ];
          auto [id0, id1] = idFunction( item );

          aggregate[ --counts[ id0 ] ] = where;
          aggregate[ --counts[ id1 ] ] = where;
        }
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