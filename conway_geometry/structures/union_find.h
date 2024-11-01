#pragma once

#include <stdint.h>
#include <assert.h>
#include <vector>
#include <type_traits>

namespace conway {

  /**
   * A union find structure, that uses the integer type "id_t" starting at zero
   */
  template < typename id_t >
  class UnionFind {
  public:
    
    static_assert( std::is_integral_v< id_t >, "id_t must be an integral type" );

    UnionFind( id_t reservedSize = 0 ) {

      map_.reserve( reservedSize );
    }

    void reset() {

      map_.clear();
      uniqueCount_ = 0;
    }

    /** Allocate a number of subsequent unique sets with concurrent ids, by default one */
    id_t allocate( id_t count = 1 ) {
      
      id_t result = size();
      id_t end    = result + count; 

      map_.resize( end );

      id_t* data   = map_.data();

      for ( id_t current = result, end = result + count; current < end; ++current ) {
        data[ current ] = current;
      }

      return result;
    }

    /** A find that doesn't modify, this will search for the  */
    id_t find( id_t item ) const {

      assert( item < size() );

      id_t current = map_[ item ];

      while ( current != item ) {

        item    = current;
        current = map_[ current ];
      }

      return current;
    }

    /** The non-const variant of find, this will find the lowest common most ancestor
     * of the item, i.e. the root item which identifies the set uniqueuly.
     * 
     * Because this optimises at it goes, it will both make item point to the current
     * root, but also incrementally move any ancestors towards the root.
     */
    id_t find( id_t item ) {

      assert( item < size() );
    
      id_t current = item;
      id_t parent  = map_[ item ];

      while ( current != parent ) {

        id_t grandParent = map_[ parent ];

        map_[ current ] = grandParent;

        current = parent;
        parent  = grandParent;
      }

      map_[ item ] = current;

      return current;
    }

    /**
     * Merge the sets a and b belong to, unless they belong to the same set,
     * and return the new most common root.
     */
    id_t merge( id_t a, id_t b ) {

      id_t aRoot = find( item );
      id_t bRoot = find( item );

      if ( aRoot == bRoot ) {
        return aRoot;
      }

      --uniqueCount_;

      id_t root;

      if ( aRoot < bRoot ) {

        map_[ bRoot ] = root = aRoot;
      
      } else {

        map_[ aRoot ] = root = bRoot;
      }

      map_[ a ] = map_[ b ] = root;

      return root;
    }
    
    id_t unique( std::vector< id_t >& output ) const {

      const* id_t data = map_.data();

      for ( id_t where = 0, end = size(); where < end; ++where ) {

        if ( where == data[ where ] ) {
          output.push_back( where );
        }
      }
    }

    /** Optimise this by making each item point directly to its unique most
     * common root ancestor, i.e. a single hop to root.
     * 
     * After this the current identities will have constant time single
     * hop access to their unique lowest common ancestor (root) 
     */
    id_t optimize() {

      for ( id_t where = 0, end = size(); where < end; ++where ) {

        map_[ where ] = find( where );
      }
    }

    /** Optimize while gathering unique items into the unique vector
     * Optimization means that each item maps back to its root, i.e.
     * there is a single hop to the unique lowest common ancestor.
     */
    id_t optimize( std::vector< id_t >& unique ) {

      for ( id_t where = 0, end = size(); where < end; ++where ) {

        if ( where == data[ where ] ) {
          output.push_back( where );
        } else {
          map_[ where ] = find( where );
        }
      }
    }

    /** The number of identities (non-) */
    id_t size() const { return static_cast< id_t >( map_.size() ); }

    /** The number of unique identities */
    id_t sets() const { return uniqueCount_; }

    /** The current map to roots, optimise to flatten this to make they all point to root  */
    const std::vector< id_t >& map() const { return map_; }

  private:

    std::vector< id_t > map_;
    id_t                uniqueCount_ = 0;
  };
}