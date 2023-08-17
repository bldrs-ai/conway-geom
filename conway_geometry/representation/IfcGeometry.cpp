/*
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.cpp
 * */

// Implementation for IfcGeometry

#include "IfcGeometry.h"

#include "fuzzy/aabb.h"

namespace conway::geometry {

glm::dvec3 IfcGeometry::GetExtent() const { return max - min; }

// set all vertices relative to min
void IfcGeometry::Normalize() {
  for (size_t i = 0; i < vertexData.size(); i += 6) {
    vertexData[i + 0] = vertexData[i + 0] - min.x;
    vertexData[i + 1] = vertexData[i + 1] - min.y;
    vertexData[i + 2] = vertexData[i + 2] - min.z;
  }

  normalized = true;
}

void IfcGeometry::AddComponent(IfcGeometry *g) { 
  components.push_back(g); 
}

void IfcGeometry::AddComponentTransform(glm::dmat4x4 transform) {
  componentTransforms.push_back(transform);
}

void IfcGeometry::AddComponentWithTransform(IfcGeometry *geom, glm::dmat4x4 transform) {
  for (uint32_t index = 0; index < numPoints; ++index) {
    glm::dvec4 t =
        transform *
        glm::dvec4(
            glm::dvec3(geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
                       geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
                       geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]),
            1);

    glm::dvec4 n =
        transform *
        glm::dvec4(
            glm::dvec3(geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 3],
                       geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 4],
                       geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 5]),
            0);

    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0] = t.x;
    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1] = t.y;
    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2] = t.z;
    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 3] = n.x;
    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 4] = n.y;
    geom->vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 5] = n.z;
  }

  components.push_back(geom);
}

void IfcGeometry::AddPoint(const glm::dvec4 &pt, const glm::dvec3 &n) {
  glm::dvec3 p = pt;
  AddPoint(p, n);
}

// just follow the ifc spec, damn
bool IfcGeometry::computeSafeNormal(const glm::dvec3 v1, const glm::dvec3 v2,
                                    const glm::dvec3 v3, glm::dvec3 &normal,
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

void IfcGeometry::AddPoint(const glm::dvec3 &pt, const glm::dvec3 &n) {
  // auto const source = std::vector<double>{pt.x, pt.y, pt.z, n.x, n.y, n.z};
  vertexData.insert(vertexData.end(), {pt.x, pt.y, pt.z, n.x, n.y, n.z});

  min = glm::min(min, pt);
  max = glm::max(max, pt);

  // if (std::isnan(pt.x) || std::isnan(pt.y) || std::isnan(pt.z)) {
  //   printf("NaN in geom!\n");
  // }

  // if (std::isnan(n.x) || std::isnan(n.y) || std::isnan(n.z)) {
  //   printf("NaN in geom!\n");
  // }

  numPoints += 1;
}

void IfcGeometry::AddFace(const glm::dvec3& a, const glm::dvec3& b, const glm::dvec3& c) {
  glm::dvec3 normal;
  if (!computeSafeNormal(a, b, c, normal)) {
    // bail out, zero area triangle
    printf("zero tri\n");
    return;
  }

  AddFace(numPoints + 0, numPoints + 1, numPoints + 2);

  AddPoint(a, normal);
  AddPoint(b, normal);
  AddPoint(c, normal);
}

void IfcGeometry::AddFace(uint32_t a, uint32_t b, uint32_t c) {
  indexData.insert(indexData.end(), {a, b, c});

  numFaces++;
}

void IfcGeometry::ReverseFace(uint32_t index) {
  Face f = GetFace(index);
  indexData[index * 3 + 0] = f.i2;
  indexData[index * 3 + 1] = f.i1;
  indexData[index * 3 + 2] = f.i0;
}

void IfcGeometry::ReverseFaces() {
  for (size_t i = 0; i < numFaces; i++) {
    ReverseFace(i);
  }
}

Face IfcGeometry::GetFace(uint32_t index) const {
  Face f;
  f.i0 = indexData[index * 3 + 0];
  f.i1 = indexData[index * 3 + 1];
  f.i2 = indexData[index * 3 + 2];
  return f;
}

glm::dvec3 IfcGeometry::GetPoint(uint32_t index) const {
  return glm::dvec3(vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
                    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
                    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]);
}

fuzzybools::AABB IfcGeometry::GetFaceBox(uint32_t index) const {
  fuzzybools::AABB aabb;
  aabb.index = index;
  glm::dvec3 a = GetPoint(indexData[index * 3 + 0]);
  glm::dvec3 b = GetPoint(indexData[index * 3 + 1]);
  glm::dvec3 c = GetPoint(indexData[index * 3 + 2]);
  aabb.min = glm::min(a, aabb.min);
  aabb.min = glm::min(b, aabb.min);
  aabb.min = glm::min(c, aabb.min);
  aabb.max = glm::max(a, aabb.max);
  aabb.max = glm::max(b, aabb.max);
  aabb.max = glm::max(c, aabb.max);
  aabb.center = (aabb.max + aabb.min) / 2.0;
  return aabb;
}

void IfcGeometry::GetCenterExtents(glm::dvec3 &center,
                                   glm::dvec3 &extents) const {
  glm::dvec3 min = glm::dvec3(DBL_MAX, DBL_MAX, DBL_MAX);
  glm::dvec3 max = glm::dvec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);

  for (size_t i = 0; i < numPoints; i++) {
    auto pt = GetPoint(i);
    min = glm::min(min, pt);
    max = glm::max(max, pt);
  }

  extents = (max - min);
  center = min + extents / 2.0;
}

IfcGeometry IfcGeometry::Normalize(glm::dvec3 center,
                                   glm::dvec3 extents) const {
  IfcGeometry newGeom;

  double scale = std::max(extents.x, std::max(extents.y, extents.z));

  for (size_t i = 0; i < numFaces; i++) {
    auto face = GetFace(i);
    auto a = (GetPoint(face.i0) - center) / scale;
    auto b = (GetPoint(face.i1) - center) / scale;
    auto c = (GetPoint(face.i2) - center) / scale;

    newGeom.AddFace(a, b, c);
  }

  return newGeom;
}

IfcGeometry IfcGeometry::DeNormalize(glm::dvec3 center,
                                     glm::dvec3 extents) const {
  IfcGeometry newGeom;

  double scale = std::max(extents.x, std::max(extents.y, extents.z));

  for (size_t i = 0; i < numFaces; i++) {
    auto face = GetFace(i);
    auto a = GetPoint(face.i0) * scale + center;
    auto b = GetPoint(face.i1) * scale + center;
    auto c = GetPoint(face.i2) * scale + center;

    newGeom.AddFace(a, b, c);
  }

  return newGeom;
}

uint32_t IfcGeometry::GetVertexData() {
  // unfortunately webgl can't do doubles
  if (fvertexData.size() != vertexData.size()) {
    fvertexData.resize(vertexData.size());
    for (size_t i = 0; i < vertexData.size(); i += 6) {
      fvertexData[i + 0] = vertexData[i + 0];
      fvertexData[i + 1] = vertexData[i + 1];
      fvertexData[i + 2] = vertexData[i + 2];

      fvertexData[i + 3] = vertexData[i + 3];
      fvertexData[i + 4] = vertexData[i + 4];
      fvertexData[i + 5] = vertexData[i + 5];
    }

    // cleanup
    // vertexData = {};
  }

  if (fvertexData.empty()) {
    return 0;
  }

  return (uint32_t)(size_t)&fvertexData[0];
}

IfcGeometry IfcGeometry::Clone() { return *this; }

void IfcGeometry::AppendWithTransform(IfcGeometry &geom,
                                      glm::dmat4x4 transform) {
  uint32_t maxIndex = numPoints;
  numPoints += geom.numPoints;
  min = glm::min(min, geom.min);
  max = glm::max(max, geom.max);

  size_t vertexDataStart = vertexData.size();
  vertexData.insert(vertexData.end(), geom.vertexData.begin(),
                    geom.vertexData.end());

  double *vertexDataCursor = vertexData.data() + vertexDataStart;
  double *vertexDataEnd = vertexData.data() + vertexData.size();

  for (; vertexDataCursor < vertexDataEnd;
       vertexDataCursor += VERTEX_FORMAT_SIZE_FLOATS) {
    glm::dvec4 t = transform * glm::dvec4(
      glm::dvec3(
        vertexDataCursor[0],
        vertexDataCursor[1],
        vertexDataCursor[2]),
      1);

    glm::dvec4 n = transform * glm::dvec4(
      glm::dvec3(
        vertexDataCursor[3],
        vertexDataCursor[4],
        vertexDataCursor[5]),
      0);

    vertexDataCursor[0] = t.x;
    vertexDataCursor[1] = t.y;
    vertexDataCursor[2] = t.z;
    vertexDataCursor[3] = n.x;
    vertexDataCursor[4] = n.y;
    vertexDataCursor[5] = n.z;
  }

  size_t indexDataStart = indexData.size();

  indexData.insert(indexData.end(), geom.indexData.begin(),
                   geom.indexData.end());

  uint32_t *indexDataCursor = indexData.data() + indexDataStart;
  uint32_t *indexDataEnd = indexData.data() + indexData.size();

  for (; indexDataCursor < indexDataEnd; indexDataCursor += 3) {
    indexDataCursor[0] += maxIndex;
    indexDataCursor[1] += maxIndex;
    indexDataCursor[2] += maxIndex;
  }
}

void IfcGeometry::AppendGeometry(IfcGeometry &geom) {
  uint32_t maxIndex = numPoints;
  numPoints += geom.numPoints;
  min = glm::min(min, geom.min);
  max = glm::max(max, geom.max);
  vertexData.insert(vertexData.end(), geom.vertexData.begin(),
                    geom.vertexData.end());
  for (uint32_t k = 0; k < geom.numFaces; k++) {
    AddFace(maxIndex + geom.indexData[k * 3 + 0],
            maxIndex + geom.indexData[k * 3 + 1],
            maxIndex + geom.indexData[k * 3 + 2]);
  }
}

uint32_t IfcGeometry::GetVertexDataSize() {
  return (uint32_t)fvertexData.size();
}

uint32_t IfcGeometry::GetIndexData() { return (uint32_t)(size_t)&indexData[0]; }

uint32_t IfcGeometry::GetIndexDataSize() { return (uint32_t)indexData.size(); }

void IfcGeometry::ApplyTransform(glm::dmat4 transform) {
  for (uint32_t index = 0; index < numPoints; ++index) {
    glm::dvec4 t =
        transform *
        glm::dvec4(
            glm::dvec3(vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
                       vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
                       vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]),
            1);

    glm::dvec4 n =
        transform *
        glm::dvec4(
            glm::dvec3(vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 3],
                       vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 4],
                       vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 5]),
            0);

    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0] = t.x;
    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1] = t.y;
    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2] = t.z;
    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 3] = n.x;
    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 4] = n.y;
    vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 5] = n.z;
  }
}

}  // namespace conway::geometry
