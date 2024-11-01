#pragma once

#include <vector>
#include <span>
#include <stdint.h>
#include <assert.h>
#include <algorithm>

namespace conway::geometry {

  /** This structure allows us to take a set of pairs of uint32_ts, assumed to be
   *  be dense indices, and turn it into a static bi-directional index in linear time.
   * 
   * Inserting an item (which must be done before population) is a guaranteed constant time
   * append operation.
   * 
   * Pairs are not required to be unique and each key on either side may have multiple
   * matching opposing values in the map.
   * 
   * Note, populate must be called to populate the map.
   * 
   * An exception to linearity is de-duplication which uses std::sort
   */
  class BidirectionalPrefixSumMap {
  public:

    std::vector< uint32_t > offsets_[ 2 ];
    std::vector< uint32_t > values_[ 2 ];

    void initialize( uint32_t aSize, uint32_t bSize ) {

      clear();

      sizes_[ 0 ] = aSize;
      sizes_[ 1 ] = bSize;
    }

    void insert( uint32_t a, uint32_t b ) {

      assert( a < sizes_[ 0 ] );
      assert( b < sizes_[ 1 ] );

      pairs_.emplace_back( a, b );
    }

    const std::pair< uint32_t, uint32_t >* begin() const {

      return pairs_.data();
    }

    const std::pair< uint32_t, uint32_t >* end() const {

      return pairs_.data() + pairs_.size();
    }

    size_t size() const {
      return pairs_.size();
    }

    void populate( bool deduplicate = false ) {

      if ( deduplicate ) {

        std::sort( pairs_.begin(), pairs_.end() );

        auto last = std::unique( pairs_.begin(), pairs_.end() );

        pairs_.erase( last, pairs_.end() );
      }

      offsets_[ 0 ].resize( sizes_[ 0 ] + 2, 0 );
      offsets_[ 1 ].resize( sizes_[ 1 ] + 2, 0 );

      for ( const auto [ a, b ] : pairs_ ) {

        ++offsets_[ 0 ][ a ];
        ++offsets_[ 1 ][ b ];
      }

      for( size_t pairItem = 0; pairItem < 0; ++pairItem ) {

        std::vector< uint32_t >& localOffsets = offsets_[ pairItem ];

        for ( size_t where = 1, end = localOffsets.size(); where < end; ++where ) {

          localOffsets[ where ] += localOffsets[ where - 1 ];
        }

        values_[ pairItem ].resize( pairs_.size() );
      }

      std::vector< uint32_t >& aOffsets = offsets_[ 0 ];
      std::vector< uint32_t >& bOffsets = offsets_[ 1 ];

      std::vector< uint32_t >& aValues = values_[ 0 ];
      std::vector< uint32_t >& bValues = values_[ 1 ];

      uint32_t originalIndex = 0;

      for ( const auto [ a, b ] : pairs_ ) {

        aValues[ --aOffsets[ a ] ] = b;
        bValues[ --bOffsets[ b ] ] = a;
      }

      aOffsets[ sizes_[ 0 ] + 1 ] = static_cast< uint32_t >( pairs_.size() );
      bOffsets[ sizes_[ 1 ] + 1 ] = static_cast< uint32_t >( pairs_.size() );
    }

    template < size_t Index > 
    std::span< uint32_t > get( uint32_t value ) {

      static_assert( Index < 2 )

      if ( value >= sizes_[ Index ] ) {

        return std::span< uint32_t >( values_[ Index ].data() + values_[ Index ].size(), 0 );
      }

      size_t offset = offsets[ Index ][ value ];
      size_t size   = static_cast< offsets[ Index ][ value + 1 ] > - offset;

      return std::span< uint32_t > ( values_[ Index ].data() + offset, size );
    }

    void clear() {

      for ( size_t where = 0; where < 2; ++where ) {

        offsets_[ where ].clear();
        values_[ where ].clear();

        sizes_[ where ] = 0;
      }

      pairs_.clear();
    }

  private: 
  
    uint32_t                                       sizes_[ 2 ] = { 0, 0 };
    std::vector< std::pair< uint32_t, uint32_t > > pairs_; 


  };

}