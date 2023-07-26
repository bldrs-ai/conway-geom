#pragma once

#include <glm/glm.hpp>
#include <cmath>
#include <optional>

namespace conway::geometry {

  enum class BLEND_MODE {
    BLEND_OPAQUE = 0,
    BLEND  = 1,
    MASK   = 2
  };

  struct Material {

    // Note - base alpha is used for transparency 
    glm::dvec4 base        = glm::dvec4( 0.8, 0.8, 0.8, 1 );
    double     metallic    = 0.0f;
    double     roughness   = 1.0f;
    double     alphaCutoff = 0;
    double     ior         = 1.4;

    glm::dvec4 specular = glm::dvec4( 0, 0, 0, 1 );
    bool hasSpecular = false;

    BLEND_MODE alphaMode = BLEND_MODE::BLEND_OPAQUE;

    bool doubleSided = false;

  };
}