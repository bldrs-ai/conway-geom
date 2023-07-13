#pragma once
#include <string>
#include <glm/glm.hpp>
#include <vector>
#include "ConwayCurve.h"


namespace conway::geometry {
struct IfcProfile {
  std::string type;
  IfcCurve curve;
  std::vector<IfcCurve> holes;
  bool isConvex;
  bool isComposite = false;
  std::vector<IfcProfile> profiles;

  std::string getType() const;
  IfcCurve getCurve() const;
  std::vector<IfcCurve> getHoles() const;
  bool getIsConvex() const;
  bool getIsComposite() const;
  std::vector<IfcProfile> getProfiles() const;

};
}