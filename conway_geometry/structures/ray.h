#pragma once

#include <glm/glm.hpp>

namespace conway::geometry {
  
  struct ray3 {

    glm::dvec3 origin;
    glm::dvec3 direction;
  };

  inline glm::dvec3 inverse_direction( const ray3& from ) {

    return glm::dvec3( 1 ) / from.direction;
  }
}