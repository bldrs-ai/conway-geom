#pragma once

#include <glm/glm.hpp>

namespace conway::geometry {
  
  struct box3 {

    glm::dvec3 min = { DBL_MAX, DBL_MAX, DBL_MAX };
    glm::dvec3 max = { -DBL_MAX, -DBL_MAX, -DBL_MAX }; 

    void merge( const glm::dvec3& with ) {

      min = glm::min( with, min );
      max = glm::max( with, max );
    }

    void merge( const box3& with ) {

      min = glm::min( with.min, min );
      max = glm::max( with.max, max );
    }

    glm::dvec3 interval() const {

      return max - min;
    }
  };

  inline bool overlaps( const box3& left, const box3& right ) {

    for ( glm::length_t axis = 0; axis < 3; ++axis ) {

      if( left.max[ axis ] < right.min[ axis ] ) {
        return false;
      }
      
      if( right.max[ axis ] < left.min[ axis ] ) {
        return false;
      }
    }

    return true;
  }

  inline bool contains( const box3& left, const box3& right ) {

    for ( glm::length_t axis = 0; axis < 3; ++axis ) {

      if ( left.max[ axis ] < right.max[ axis ] ) {
        return false;
      }

      if ( left.min[ axis ] > right.min[ axis ] ) {
        return false;
      }
    }

    return true;
  }

  
  inline bool overlaps( const box3& left, const glm::dvec3& right ) {

    return
      left.max.x >= right.x &&
      left.min.x <= right.x &&
      left.max.y >= right.y &&
      left.min.y <= right.y &&
      left.max.z >= right.z &&
      left.min.z <= right.z;
  }

  inline box3 merge( const box3& left, const box3& right ) {

    return {
      glm::min( left.min, right.min ),
      glm::max( left.max, right.max )
    };
  }
}