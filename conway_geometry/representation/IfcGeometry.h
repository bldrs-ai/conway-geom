/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.h
 * */

// Represents a single piece of IFC Geometry

#pragma once

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-function"

#include <fuzzy/geometry.h>

#pragma clang diagnostic pop

#include <glm/glm.hpp>
#include <optional>
#include <string>
#include <vector>

//#include "geometry.h"
#include "material.h"

/*namespace fuzzybools
{
    struct AABB; // Forward declaration for fuzzybools::AABB
}*/

namespace conway::geometry {

template <typename T>
constexpr size_t byteSize(const std::vector<T> &data) {
  return data.size() * sizeof(T);
}

inline bool computeSafeNormal(const glm::dvec3& v1, const glm::dvec3& v2,
                              const glm::dvec3& v3, glm::dvec3 &normal,
                              double eps = 0) {
  glm::dvec3 v12(v2 - v1);
  glm::dvec3 v13(v3 - v1);

  glm::dvec3 norm = glm::cross(v12, v13);

  double len = glm::length(norm);

  if (len <= eps) {
    return false;
  }

  normal = norm / len;

  return true;
}


struct IfcGeometry : fuzzybools::Geometry {
  //new additions 0.0.44
  bool halfSpace = false;
  std::vector<IfcGeometry> part;
  glm::dvec3 center = glm::dvec3(0,0,0);
  glm::dvec3 halfSpaceX = glm::dvec3(1, 0, 0);
  glm::dvec3 halfSpaceY = glm::dvec3(0, 1, 0);
  glm::dvec3 halfSpaceZ = glm::dvec3(0, 0, 1);
  glm::dvec3 halfSpaceOrigin = glm::dvec3(0, 0, 0);
  //end new additions 0.0.44

  std::vector<IfcGeometry> getParts();
  fuzzybools::AABB getAABB() const;
  glm::dvec3 GetAABBCenter() const;
  bool normalized = false;

  uint32_t GetAllocationSize() const;

  inline void AddFaceFloat(const float* af, const float* bf, const float* cf) {
    glm::dvec3 normal;

    glm::dvec3 a( af[ 0 ], af[ 1 ], af[ 2 ] );
    glm::dvec3 b( bf[ 0 ], bf[ 1 ], bf[ 2 ] );
    glm::dvec3 c( cf[ 0 ], cf[ 1 ], cf[ 2 ] );

    if (!conway::geometry::computeSafeNormal(a, b, c, normal, EPS_SMALL))
    {
      // bail out, zero area triangle
      //printf("zero tri");
      return;
    }

    AddPoint(a, normal);
    AddPoint(b, normal);
    AddPoint(c, normal);

    AddFace(numPoints - 3, numPoints - 2, numPoints - 1);
  }

  void ReverseFace(uint32_t index);
  void ReverseFaces();
  void AddPart(IfcGeometry geom);
  void AddPart(fuzzybools::Geometry geom);
  glm::dvec3 GetPoint(uint32_t index) const;
  glm::dvec3 GetNormal(uint32_t index) const;
  uint32_t GetVertexData();
  void AppendGeometry(IfcGeometry &geom);
  void AddGeometry(fuzzybools::Geometry geom, glm::dmat4 trans = glm::dmat4(1), double scx = 1, double scy = 1, double scz = 1, glm::dvec3 origin = glm::dvec3(0, 0, 0));
  void MergeGeometry(fuzzybools::Geometry geom);
  void AppendWithTransform(const IfcGeometry &geom, const glm::dmat4x4& transform);
  uint32_t GetVertexDataSize();
  uint32_t GetIndexData();
  uint32_t GetIndexDataSize();
  glm::dvec3 Normalize();
  void ApplyTransform(const glm::dmat4x4& transform);
  IfcGeometry Clone();
  std::string GeometryToObj( const std::string& preamble = "" ) const;

 private:
  bool computeSafeNormal(const glm::dvec3 v1, const glm::dvec3 v2,
                         const glm::dvec3 v3, glm::dvec3 &normal, double eps);
};

struct IfcGeometryCollection {
  std::vector<IfcGeometry *> components;
  std::vector<glm::dmat4x4> transforms;

  uint32_t materialIndex = 0;
  bool hasDefaultMaterial = true;

  void AddComponentWithTransform(IfcGeometry *geom,
                                 const glm::dmat4x4 &transform) {
    if (geom != nullptr) {
      components.push_back(geom);
      transforms.push_back(transform);

      currentSize += geom->GetAllocationSize();
    }
  }

  size_t currentSize = 0;
};

}  // namespace conway::geometry
