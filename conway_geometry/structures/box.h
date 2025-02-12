#pragma once

#include <glm/glm.hpp>

#include "ray.h"

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

    void rescale( const glm::dvec3& scale, const glm::dvec3& origin ) {

      min = ( ( min - origin ) * scale ) + origin;
      max = ( ( max - origin ) * scale ) + origin;
    }

    glm::dvec3 interval() const {

      return max - min;
    }

    glm::dvec3 centre() const {

      return min + ( interval() * 0.5 );
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

  inline double overlaps( const box3& box, const ray3& ray, const glm::dvec3& inverseDirection, double& tMin, double& tMax ) {

    glm::dvec3 tMinExtent  = ( box.min - ray.origin ) * inverseDirection;
    glm::dvec3 tMaxTextent = ( box.max - ray.origin ) * inverseDirection;

    glm::dvec3 tMin3 = glm::min( tMinExtent, tMaxTextent );
    glm::dvec3 tMax3 = glm::max( tMinExtent, tMaxTextent );
   
    tMin = std::min( tMin3.x, std::min( tMin3.y, tMin3.z ) );
    tMax = std::max( tMax3.x, std::max( tMax3.y, tMax3.z ) );

    return ( tMin < tMax && tMax >= 0 );
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