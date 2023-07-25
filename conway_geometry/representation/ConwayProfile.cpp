


#include "ConwayProfile.h"

#include <glm/glm.hpp>
#include <vector>

#include "../operations/geometryutils.h"

namespace conway::geometry {

std::string IfcProfile::getType() const {
  return type;
}

IfcCurve IfcProfile::getCurve() const {
  return curve;
}

std::vector<IfcCurve> IfcProfile::getHoles() const {
  return holes;
}

bool IfcProfile::getIsConvex() const {
  return isConvex;
}

bool IfcProfile::getIsComposite() const {
  return isComposite;
}

std::vector<IfcProfile> IfcProfile::getProfiles() const {
  return profiles;
}

} // namespace conway::geometry