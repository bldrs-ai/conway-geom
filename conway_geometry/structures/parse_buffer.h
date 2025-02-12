#pragma once

#include <stdint.h>
#include <span>
#include <bit>

namespace conway {

  class ParseBuffer {
  public:

    static constexpr size_t DEFAULT_BASE_CAPACITY = 128 * 1024; 

    ParseBuffer() : capacity_( DEFAULT_BASE_CAPACITY ), size_( 0 ) {

      data_ = new uint8_t[ capacity_ ];
    }

    ~ParseBuffer() {

      delete [] data_;
    }

    size_t size() const { 
      return size_;
    }

    size_t capacity() const {
      return capacity_;
    }

    std::span< const uint8_t > range() const {

      return std::span< const uint8_t >( data_, size_ );
    }

    uintptr_t dataInteger() {

      return reinterpret_cast< uintptr_t >( data_ ); 
    }
    
    const uint8_t* data() const {

      return data_; 
    }
    
    uint8_t* data() {

      return data_; 
    }

    /** Resize the parse buffer, note this may actually destroy its current contents */
    uintptr_t resize( size_t size ) {

      if ( capacity_ < size ) {

        size_t newCapacity = std::bit_ceil( size );

        uint8_t* candidateData = new uint8_t[ newCapacity ];

        delete[] data_;

        data_     = candidateData;
        capacity_ = newCapacity;
      }

      size_ = size;

      return dataInteger();
    }

  private:

    uint8_t* data_;
    size_t capacity_;
    size_t size_;

  };

}