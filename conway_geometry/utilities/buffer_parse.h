#pragma once

#include <glm/glm.hpp>
#include <fast_float/fast_float.h>
#include <stdint.h>
#include <span>
#include <vector>
#include <optional>
#include "representation/Topology.h"


namespace conway {

  size_t skip_comment( size_t cursor, const std::span< const uint8_t >& input );

  std::optional< ptrdiff_t > step_begin_array( size_t cursor, const std::span< const uint8_t >& input );
  
  std::optional< ptrdiff_t > step_array_token( size_t cursor, const std::span< const uint8_t >& input );

  struct ParseResult {

    size_t cursor;

    std::errc error;
  };

  template < typename T >
  ParseResult parse_vector( const std::span< const uint8_t >& input, std::vector< T >& output );
  
  template < typename T, typename TParser >
  ParseResult parse_vector_custom( const std::span< const uint8_t >& input, std::vector< T >& output, TParser parser );
}


namespace fast_float {

  from_chars_result from_chars( const char* begin, const char* end, conway::geometry::Triangle& triangle );
  
  from_chars_result from_chars( const char* begin, const char* end, glm::dvec3& value );

  from_chars_result from_chars( const char* begin, const char* end, glm::dvec2& value );
}


inline size_t conway::skip_comment( size_t cursor, const std::span< const uint8_t >& input ) {

  size_t         size = input.size();
  const uint8_t* data = input.data();

  if ( cursor + 4 > size || data[ cursor ] != '/' || data[ cursor + 1 ] != '*' ) {

    return cursor;
  }

  cursor += 2;

  do {

    while ( cursor < size && data[ cursor ] != '*' ) {

      ++cursor;
    }

    if ( data[ cursor ] == '/' ) {

      break;
    }

  } while ( cursor < size );

  return cursor;
}


inline std::optional< ptrdiff_t > conway::step_begin_array( size_t cursor, const std::span< const uint8_t >& input ) {

  size_t         endCursor = input.size();
  const uint8_t* buffer    = input.data();

  if ( ( cursor + 1 ) >= endCursor || buffer[ cursor ] != '(' ) {

    return std::nullopt;
  }

  ++cursor;

  size_t previousCursor;

  do {

    previousCursor = cursor;

    while ( cursor < endCursor ) {

      uint8_t data = buffer[ cursor ];

      if ( data != ' ' && data != '\t' && data != '\r' && data != '\n' ) {
        break;
      }

      ++cursor;
    }

    cursor = skip_comment( cursor, input );
  }
  while ( previousCursor != cursor );

  if ( buffer[ cursor ] == ')' ) {
    return -( cursor + 1 );
  }

  return cursor;
}


inline std::optional< ptrdiff_t > conway::step_array_token( size_t cursor, const std::span< const uint8_t >& input ) {

  size_t         endCursor = input.size();
  const uint8_t* buffer    = input.data();

  size_t previousCursor;

  do {

    previousCursor = cursor;

    while ( cursor < endCursor ) {

      uint8_t data = buffer[ cursor ];

      if ( data != ' ' && data != '\t' && data != '\r' && data != '\n' ) {
        break;
      }

      ++cursor;
    }

    cursor = skip_comment( cursor, input );
  }
  while ( previousCursor != cursor );

  if ( cursor >= endCursor ) {

    return std::nullopt;
  }

  if ( buffer[ cursor ] == ')' ) {
    return -( cursor + 1 );
  }

  if ( buffer[ cursor ] != ',' ) {

    return std::nullopt;
  }

  ++cursor;

  do {

    previousCursor = cursor;

    while ( cursor < endCursor ) {

      uint8_t data = buffer[ cursor ];

      if ( data != ' ' && data != '\t' && data != '\r' && data != '\n' ) {
        break;
      }

      ++cursor;
    }

    cursor = skip_comment( cursor, input );
  }
  while ( previousCursor != cursor );

  return cursor;
}


template < typename T, typename TParser >
inline conway::ParseResult conway::parse_vector_custom( const std::span< const uint8_t >& input, std::vector< T >& output, TParser parser ) {

  std::optional< ptrdiff_t > arrayBegin = step_begin_array( 0, input );

  if ( !arrayBegin.has_value() ) {

    return { 0, std::errc::invalid_argument };
  }

  ptrdiff_t signedCursor = arrayBegin.value();

  const char* buffer    = reinterpret_cast< const char* >( input.data() );
  size_t         end    = input.size();
  const char* bufferEnd = buffer + end;

  while ( signedCursor >= 0 ) {

    size_t   cursor    = static_cast< size_t >( signedCursor );
    T        cellValue {};

    auto [ endPtr, errorCode ] = parser( buffer + cursor, bufferEnd, cellValue );

    if ( errorCode != std::errc() ) {

      return { cursor, errorCode };
    }

    cursor = static_cast< size_t >( endPtr - buffer );

    output.push_back( cellValue );

    std::optional< ptrdiff_t > nextCursor = step_array_token( cursor, input );

    if ( !nextCursor.has_value() ) {

      return { cursor, std::errc::invalid_argument };
    }

    signedCursor = nextCursor.value();
  }

  return { static_cast< size_t >( -signedCursor ), std::errc() };
}


template < typename T >
inline conway::ParseResult conway::parse_vector( const std::span< const uint8_t >& input, std::vector< T >& output ) {

  std::optional< ptrdiff_t > arrayBegin = step_begin_array( 0, input );

  if ( !arrayBegin.has_value() ) {

    return { 0, std::errc::invalid_argument };
  }

  ptrdiff_t signedCursor = arrayBegin.value();

  const char* buffer    = reinterpret_cast< const char* >( input.data() );
  size_t      end       = input.size();
  const char* bufferEnd = buffer + end;

  while ( signedCursor >= 0 ) {

    size_t   cursor    = static_cast< size_t >( signedCursor );
    T        cellValue {};

    auto [ endPtr, errorCode ] = fast_float::from_chars( buffer + cursor, bufferEnd, cellValue );

    if ( errorCode != std::errc() ) {

      return { cursor, errorCode };
    }

    cursor = static_cast< size_t >( endPtr - buffer );

    output.push_back( cellValue );

    std::optional< ptrdiff_t > nextCursor = step_array_token( cursor, input );

    if ( !nextCursor.has_value() ) {

      return { cursor, std::errc::invalid_argument };
    }

    signedCursor = nextCursor.value();
  }

  return { static_cast< size_t >( -signedCursor ), std::errc() };
}


inline fast_float::from_chars_result fast_float::from_chars(
  const char* begin,
  const char* end,
  conway::geometry::Triangle& triangle ) {

    std::span< const uint8_t > bufferSpan(
      reinterpret_cast< const uint8_t* >( begin ),
      static_cast< size_t >( end - begin ) );
  
    std::optional< ptrdiff_t > arrayBegin = conway::step_begin_array( 0, bufferSpan );
  
    if ( !arrayBegin.has_value() ) {
  
      return { begin, std::errc::invalid_argument };
    }
  
    ptrdiff_t signedCursor = arrayBegin.value();
  
    for ( size_t where = 0; where < 3; ++where )
    {
      if ( signedCursor < 0 ) {
  
        return { begin, std::errc::invalid_argument };
      }
      
      size_t cursor = abs( signedCursor );
  
      auto [ endPtr, errorCode ] = fast_float::from_chars( begin + cursor, end, triangle.vertices[ where ] );
  
      if ( errorCode != std::errc() ) {
  
        return { begin, std::errc::invalid_argument };
      }
  
      cursor = endPtr - begin;
  
      std::optional< ptrdiff_t > nextCursor = conway::step_array_token( cursor, bufferSpan );
  
      if ( !nextCursor.has_value() ) {
  
        return { begin, std::errc::invalid_argument };
      }
  
      signedCursor = nextCursor.value();
    }
  
    if ( signedCursor >= 0 ) {
  
      return { begin, std::errc::invalid_argument };
    }
  
    return { begin - signedCursor, std::errc() };
}


inline fast_float::from_chars_result fast_float::from_chars(
  const char* begin, const char* end, glm::dvec3& value ) {

  std::span< const uint8_t > bufferSpan(
    reinterpret_cast< const uint8_t* >( begin ),
    static_cast< size_t >( end - begin ) );

  std::optional< ptrdiff_t > arrayBegin = conway::step_begin_array( 0, bufferSpan );

  if ( !arrayBegin.has_value() ) {

    return { begin, std::errc::invalid_argument };
  }

  ptrdiff_t signedCursor = arrayBegin.value();

  for ( size_t where = 0; where < 3; ++where )
  {
    if ( signedCursor < 0 ) {

      return { begin, std::errc::invalid_argument };
    }
    
    size_t cursor = abs( signedCursor );

    auto [ endPtr, errorCode ] = fast_float::from_chars( begin + cursor, end, value[ where ] );

    if ( errorCode != std::errc() ) {

      return { begin, std::errc::invalid_argument };
    }

    cursor = endPtr - begin;

    std::optional< ptrdiff_t > nextCursor = conway::step_array_token( cursor, bufferSpan );

    if ( !nextCursor.has_value() ) {

      return { begin, std::errc::invalid_argument };
    }

    signedCursor = nextCursor.value();
  }

  if ( signedCursor >= 0 ) {

    return { begin, std::errc::invalid_argument };
  }

  return { begin - signedCursor, std::errc() };
}


inline fast_float::from_chars_result fast_float::from_chars(
  const char* begin, const char* end, glm::dvec2& value ) {

  std::span< const uint8_t > bufferSpan(
    reinterpret_cast< const uint8_t* >( begin ),
    static_cast< size_t >( end - begin ) );

  std::optional< ptrdiff_t > arrayBegin = conway::step_begin_array( 0, bufferSpan );

  if ( !arrayBegin.has_value() ) {

    return { begin, std::errc::invalid_argument };
  }

  ptrdiff_t signedCursor = arrayBegin.value();

  for ( size_t where = 0; where < 2; ++where )
  {
    if ( signedCursor < 0 ) {

      return { begin, std::errc::invalid_argument };
    }
    
    size_t cursor = abs( signedCursor );

    auto [ endPtr, errorCode ] = fast_float::from_chars( begin + cursor, end, value[ where ] );

    if ( errorCode != std::errc() ) {

      return { begin, std::errc::invalid_argument };
    }

    cursor = endPtr - begin;

    std::optional< ptrdiff_t > nextCursor = conway::step_array_token( cursor, bufferSpan );

    if ( !nextCursor.has_value() ) {

      return { begin, std::errc::invalid_argument };
    }

    signedCursor = nextCursor.value();
  }

  if ( signedCursor >= 0 ) {

    return { begin, std::errc::invalid_argument };
  }

  return { begin - signedCursor, std::errc() };
}