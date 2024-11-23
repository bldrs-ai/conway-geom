#pragma once

#include <array>
#include <span>

namespace conway {

  template < typename T, size_t Size >
  class FixedStack {
  public:

    FixedStack() {}

    bool empty() const { return end_ == 0; }

    bool full() const { return end_ == Size; }

    static constexpr size_t max_size = Size;

    template < typename... Args >
    T& emplace( Args... args ) {
      
      assert( !full() );

      return ( *new (&values_[ end_++ ].t) T( args... ) );
    }

    T* begin() { return &values_[ 0 ].t; }

    T* end() { return &values_[ 0 ].t + end_; }
    
    const T* begin() const { return &values_[ 0 ].t; }

    const T* end() const { return &values_[ 0 ].t + end_; }

    const T& operator[]( size_t index ) const { return values_[ index ].t; }

    T& operator[]( size_t index ) { return values_[ index ].t; }

    std::span< T > values() { return std::span( &values_[ 0 ].t, end_); }

    std::span< const T > values() const { return std::span( &values_[ 0 ].t, end_ ); }

    T* data() { return values_; }

    const T* data() const { return values_; }

    T& top() {
      return values_[ end_ - 1 ].t;
    }

    const T& top() const {
      return values_[ end_ - 1 ].t;
    }

    void pop() {

      assert( !empty() );

      values_[ end_-- ].t.~T();
    }

    size_t size() const {

      return end_;
    }

    void push( const T& value ) {

      assert( !full() );

      (new (&values_[end_++]) T( value ) );
    }

    void clear() {

      while ( !empty() ) {

        // guarantees destruction.
        pop();
      }
    }

  private:

    union TorBytes {

      T t;
      uint8_t byte;

    };

    TorBytes values_[ Size ] = {};

    size_t end_ = 0;

  };
}