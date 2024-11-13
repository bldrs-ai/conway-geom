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

      return ( *new (&values_[ end_++ ]) T( args... ) );
    }

    T* begin() { return values_; }

    T* end() { return values_ + end_; }
    
    const T* begin() const { return values_; }

    const T* end() const { return values_ + end_; }

    const T& operator[]( size_t index ) const { return values_[ index ]; }

    T& operator[]( size_t index ) { return values_[ index ]; }

    std::span< T > values() { return std::span( values_, end_ ); }

    std::span< const T > values() const { return std::span( values_, end_ ); }

    T* data() { return values_; }

    const T* data() const { return values_; }

    T& top() {
      return values_[ end_ ];
    }

    const T& top() const {
      return values_[ end_ ];
    }

    void pop() {

      assert( !empty() );

      --end_;
    }

    size_t size() const {

      return end_;
    }

    void push( const T& value ) {

      assert( !full() );

      (new (&values_[end_++]) T( value ) );
    }

  private:

    T      values_[ Size ];
    size_t end_ = 0;

  };
}