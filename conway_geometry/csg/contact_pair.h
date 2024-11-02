#include <stdint.h>

namespace conway::geometry {

  enum class ContactRegion : uint8_t {

    E0   = 0,
    E1   = 1,
    E2   = 2,
    FACE = 3,
    V0   = 4,
    V1   = 5,
    V2   = 6,

  };

  constexpr ContactRegion vertex0( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V0;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          static_cast< uint32_t >( a ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr ContactRegion vertex1( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V1;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          ( ( static_cast< uint32_t >( a ) + 1 ) % 3 ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr ContactRegion vertex1( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V1;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          ( ( static_cast< uint32_t >( a ) + 1 ) % 3 ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr ContactRegion vertex2( ContactRegion a ) {

    if ( a == ContactRegion::FACE ) {

      return ContactRegion::V1;
    }

    if ( a > ContactRegion::FACE ) {

      return a;
    }

    return
      static_cast< ContactRegion >(
          ( ( static_cast< uint32_t >( a ) + 2 ) % 3 ) +
          static_cast< uint32_t >( ContactRegion::V0 ) );
  }

  constexpr bool isSameEdge( ContactRegion a, ContactRegion b ) {

    if ( a == ContactRegion::FACE || b == ContactRegion::FACE ) {

      return false;
    }

    if ( a == b ) {
      return true;
    }

    return
      ( vertex0( a ) == vertex0( b ) ) &&
      ( vertex0( a ) == vertex0( b ) );
  }

  static constexpr ContactRegion vertex( uint32_t vertexInFace ) {

    return static_cast< ContactRegion >( vertexInFace + static_cast< uint32_t >( ContactRegion::V1 ) ) );
  }

  static constexpr ContactRegion edge( uint32_t edgeInFace ) {

    return static_cast< ContactRegion >( edgeInFace );
  }

  struct ContactPair {

    ContactRegion with : 3;
    ContactRegion against : 3;
    
    ContactPair( ContactRegion _with, ContactRegion _against ) : with( _with ), against( _against ) {}
  
    ContactPair() {}

    ContactPair( const ContactPair& ) = default;

    ContactPair& operator=( ContactPair& ) = default;
  };

}