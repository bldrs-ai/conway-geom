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

    static constexpr max_size = Size;

    template < typename... Args >
    T& emplace( Args... args ) {

      return (*new (&values_[ end__++ ])( args... ) );
    }

    T* begin() { return values_; }

    T* end() { return values_ + end_; }
    
    const T* begin() const { return values_; }

    const T* end() const { return values_ + end_; }

    std::span< T > values() { return std::span( values_, end_ ); }

    T* data() { return values_; }

    const T* data() const { return values_; }

    T& top() {
      return values_[ end__ ];
    }

    const T& top() const {
      return values_[ end__ ];
    }

    void pop() {

      --end_;
    }

    size_t size() const {

      return end_;
    }

    void push( const T& value ) {

      values_[ end_++ ] = value;
    }

  private:

    T      values_[ Size ];
    size_t end_ = 0;

  };
}