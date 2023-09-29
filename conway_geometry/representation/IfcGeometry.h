/* 
 * Decoupling: https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md 
 * Ref: https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.h
 * */

// Represents a single piece of IFC Geometry

#pragma once

#include <glm/glm.hpp>
#include <string>
#include <vector>
#include <optional>

#include "geometry.h"
#include "material.h"

namespace fuzzybools
{
    struct AABB; // Forward declaration for fuzzybools::AABB
}

namespace conway::geometry {

template< typename T >
constexpr size_t byteSize( const std::vector< T >& data ) {
  return data.size() * sizeof( T );
}

struct IfcGeometry {

  std::vector<float> fvertexData;
  std::vector<double> vertexData;
  std::vector<uint32_t> indexData;
  glm::dvec3 min = glm::dvec3(DBL_MAX, DBL_MAX, DBL_MAX);
  glm::dvec3 max = glm::dvec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  bool normalized = false;

  uint32_t numPoints = 0;
  uint32_t numFaces = 0;

  uint32_t GetAllocationSize() const;
  
  glm::dvec3 GetExtent() const;
  void Normalize();
  void NormalizeInPlace();
  void AddComponent(IfcGeometry &g);
  void AddPoint(const glm::dvec4 &pt, const glm::dvec3 &n);
  void AddPoint(const glm::dvec3 &pt, const glm::dvec3 &n);
  void AddFace(const glm::dvec3& a, const glm::dvec3& b, const glm::dvec3& c);
  void AddFace(uint32_t a, uint32_t b, uint32_t c);
  void ReverseFace(uint32_t index);
  void ReverseFaces();
  Face GetFace(uint32_t index) const;
  fuzzybools::AABB GetFaceBox(uint32_t index) const;
  glm::dvec3 GetPoint(uint32_t index) const;
  void GetCenterExtents(glm::dvec3 &center, glm::dvec3 &extents) const;
  IfcGeometry Normalize(glm::dvec3 center, glm::dvec3 extents) const;
  IfcGeometry DeNormalize(glm::dvec3 center, glm::dvec3 extents) const;
  uint32_t GetVertexData();
  void AppendGeometry(IfcGeometry &geom);
  void AppendWithTransform(IfcGeometry &geom, glm::dmat4x4 transform);
  uint32_t GetVertexDataSize();
  uint32_t GetIndexData();
  uint32_t GetIndexDataSize();
  void ApplyTransform(glm::dmat4x4 transform);
  IfcGeometry Clone();
  void ToObj(uint32_t expressID);

 private:
  bool computeSafeNormal(const glm::dvec3 v1, const glm::dvec3 v2,
                         const glm::dvec3 v3, glm::dvec3 &normal, double eps);
};


struct IfcGeometryCollection {

  std::vector<IfcGeometry*> components;
  std::vector<glm::dmat4x4> transforms;
  
  uint32_t materialIndex      = 0;
  bool     hasDefaultMaterial = true;

  void AddComponentWithTransform( IfcGeometry *geom, const glm::dmat4x4& transform) {

    if ( geom != nullptr ) {
      components.push_back( geom );
      transforms.push_back( transform );

      currentSize += geom->GetAllocationSize();
    }
  }
  
  size_t currentSize = 0;
};

}  // namespace conway::geometry
