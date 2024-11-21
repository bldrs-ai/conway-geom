/*
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.cpp
 * */

// Implementation for IfcGeometry

#include "IfcGeometry.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "structures/hash_functions.h"
#include "csg/csg_utils.h"
#include "csg/csg.h"
#include "structures/vertex_welder.h"

namespace conway::geometry {

VertexWelder welder;
  
IfcGeometry::IfcGeometry( WingedEdgeDV3&& from ) : fuzzybools::Geometry(), wingedEdgeMesh() {

  const WingedEdgeDV3&             mesh     = wingedEdgeMesh.emplace( from );
  const std::vector< glm::dvec3 >& vertices = mesh.vertices;

  for ( const ConnectedTriangle& triangle : from.triangles ) {

    if (
      vertices[ triangle.vertices[ 0 ] ] == vertices[ triangle.vertices[ 1 ] ] ||
      vertices[ triangle.vertices[ 1 ] ] == vertices[ triangle.vertices[ 2 ] ] ||
      vertices[ triangle.vertices[ 2 ] ] == vertices[ triangle.vertices[ 0 ] ] ) {

      continue;
    }

    AddFace( vertices[ triangle.vertices[ 0 ] ], vertices[ triangle.vertices[ 1 ] ], vertices[ triangle.vertices[ 2 ] ] );
  }
}

void IfcGeometry::ReverseFace(uint32_t index) {
  fuzzybools::Face f = GetFace(index);
  indexData[index * 3 + 0] = f.i2;
  indexData[index * 3 + 1] = f.i1;
  indexData[index * 3 + 2] = f.i0;
}

void IfcGeometry::ReverseFaces() {
  for (size_t i = 0; i < numFaces; i++) {
    ReverseFace(i);
  }
}

glm::dvec3 IfcGeometry::Normalize()
{
  if (!normalized)
  {
    glm::dvec3 extents;
    GetCenterExtents(center,extents);

    for (size_t i = 0; i < vertexData.size(); i += 6)
    {
      vertexData[i + 0] = vertexData[i + 0] - center.x;
      vertexData[i + 1] = vertexData[i + 1] - center.y;
      vertexData[i + 2] = vertexData[i + 2] - center.z;
    }
    normalized = true;
  }
  return center;
}

glm::dvec3 IfcGeometry::GetAABBCenter() const {
   glm::dvec3 _center(0,0,0);
   glm::dvec3 extents;

   GetCenterExtents(_center, extents);

   return _center;
}

glm::dvec3 IfcGeometry::GetPoint(uint32_t index) const {
  return glm::dvec3(
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 0],
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 1],
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 2]);
}

glm::dvec3 IfcGeometry::GetNormal(uint32_t index) const {
  return glm::dvec3(
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 3],
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 4],
      vertexData[index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS + 5]);
}

uint32_t IfcGeometry::GetVertexData() {
  // unfortunately webgl can't do doubles
  if (fvertexData.size() != vertexData.size()) {
    fvertexData.resize(vertexData.size());
    for (size_t i = 0; i < vertexData.size(); i++) {
      // The vector was previously copied in batches of 6, but
      // copying single entry at a time is more resilient if the
      // underlying geometry lib changes the treatment of normals
      fvertexData[i] = static_cast< float >( vertexData[i] );
    }
  }
  if (fvertexData.empty()) {
    return 0;
  }
  return (uint32_t)(size_t)&fvertexData[0];
}

std::string IfcGeometry::GeometryToObj(
  const std::string& preamble) const {
  
  std::string obj;
  obj.reserve(numPoints * 32 + numFaces * 32);  // preallocate memory

  const char *vFormat = "v %.6f %.6f %.6f\nvn %.6f %.6f %.6f\n";
  const char *fFormat = "f %zu//%zu %zu//%zu %zu//%zu\n";

  obj.append( preamble );

  for (uint32_t i = 0; i < numPoints; ++i) {
    glm::dvec3 t = GetPoint(i);
    glm::dvec3 n = GetNormal(i);
    char vBuffer[128];
    snprintf(vBuffer, sizeof(vBuffer), vFormat, t.x, t.y, t.z, n.x, n.y, n.z);
    obj.append(vBuffer);    
  }

  for (uint32_t i = 0; i < numFaces; ++i) {
    size_t f1 = indexData[i * 3 + 0] + 1;
    size_t f2 = indexData[i * 3 + 1] + 1;
    size_t f3 = indexData[i * 3 + 2] + 1;

    char fBuffer[128];
    snprintf(fBuffer, sizeof(fBuffer), fFormat, f1, f1, f2, f2, f3, f3);
    obj.append(fBuffer);
  }

  return obj;
}

IfcGeometry IfcGeometry::Clone() { return *this; }

void IfcGeometry::AppendWithTransform(const IfcGeometry &geom,
                                      const glm::dmat4x4& transform) {
  uint32_t maxIndex = numPoints;
  numPoints += geom.numPoints;

  size_t vertexDataStart = vertexData.size();
  vertexData.insert(vertexData.end(), geom.vertexData.begin(),
                    geom.vertexData.end());

  double *vertexDataCursor = vertexData.data() + vertexDataStart;
  double *vertexDataEnd = vertexData.data() + vertexData.size();

  for (; vertexDataCursor < vertexDataEnd;
       vertexDataCursor += fuzzybools::VERTEX_FORMAT_SIZE_FLOATS) {
    glm::dvec4 t = transform * glm::dvec4(glm::dvec3(vertexDataCursor[0],
                                                     vertexDataCursor[1],
                                                     vertexDataCursor[2]),
                                          1);

    glm::dvec4 n = transform * glm::dvec4(glm::dvec3(vertexDataCursor[3],
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

  for ( ; indexDataCursor < indexDataEnd; indexDataCursor += 3 ) {
    indexDataCursor[0] += maxIndex;
    indexDataCursor[1] += maxIndex;
    indexDataCursor[2] += maxIndex;
  }

  wingedEdgeMesh.reset();
  //TODO: see if this is needed 
  //AddPart(geom);
}

conway::geometry::WingedEdgeDV3& IfcGeometry::GetWingedEdgeMesh() {

  if ( !wingedEdgeMesh.has_value() ) {

    wingedEdgeMesh.emplace();
 
    WingedEdgeDV3& mesh = *wingedEdgeMesh;

    for ( uint32_t where = 0; where < numPoints; ++where ) {
    
      glm::dvec3 from = GetPoint( where );

      from *= exp2( 24 );

      from = glm::round( from );

      from *= exp2( -24 );

      mesh.vertices.push_back( from );
    }

    for ( uint32_t where = 0, end = numFaces * 3; where < end; where += 3 ) {

      size_t f1 = indexData[ where + 0 ];
      size_t f2 = indexData[ where + 1 ];
      size_t f3 = indexData[ where + 2 ];

      mesh.makeTriangle( f1, f2, f3 );
    }

    welder.weld( *wingedEdgeMesh, exp2( -23 ) );

    {
      CSG cleaner;
 
      cleaner.clean( *wingedEdgeMesh );
    }
  }

  return *wingedEdgeMesh;
}

void IfcGeometry::AppendGeometry(IfcGeometry &geom) {
  uint32_t maxIndex = numPoints;
  numPoints += geom.numPoints;
  vertexData.insert(vertexData.end(), geom.vertexData.begin(),
                    geom.vertexData.end());
  for (uint32_t k = 0; k < geom.numFaces; k++) {
    AddFace(maxIndex + geom.indexData[k * 3 + 0],
            maxIndex + geom.indexData[k * 3 + 1],
            maxIndex + geom.indexData[k * 3 + 2]);
  }

  if ( wingedEdgeMesh.has_value() ) {

    wingedEdgeMesh.reset();
  }

  //TODO: see if this is needed 
  //AddPart(geom);
}

std::vector<IfcGeometry> IfcGeometry::getParts() { return part; }

void IfcGeometry::AddPart( const IfcGeometry& geom ) { part.push_back(geom); }

void IfcGeometry::AddPart( const fuzzybools::Geometry& geom)
{
  part.emplace_back(geom);
}

void IfcGeometry::AddPart( fuzzybools::Geometry&& geom)
{
  part.emplace_back(geom);
}

fuzzybools::AABB IfcGeometry::getAABB() const
{
  return GetAABB();
}

void IfcGeometry::AddGeometry(const fuzzybools::Geometry& geom, const glm::dmat4& trans,
                              double scx, double scy, double scz,
                              const glm::dvec3& origin) {
  
  if ( wingedEdgeMesh.has_value() ) {

    wingedEdgeMesh.reset();
  }
  
  for (uint32_t i = 0; i < geom.numFaces; i++) {
    fuzzybools::Face f = geom.GetFace(i);
    glm::dvec3 a = geom.GetPoint(f.i0);
    glm::dvec3 b = geom.GetPoint(f.i1);
    glm::dvec3 c = geom.GetPoint(f.i2);
    if (scx != 1 || scy != 1 || scz != 1) {
      double aax = glm::dot(trans[0], glm::dvec4(a - origin, 1)) * scx;
      double aay = glm::dot(trans[1], glm::dvec4(a - origin, 1)) * scy;
      double aaz = glm::dot(trans[2], glm::dvec4(a - origin, 1)) * scz;
      a = origin + glm::dvec3(aax * trans[0]) + glm::dvec3(aay * trans[1]) +
          glm::dvec3(aaz * trans[2]);
      double bbx = glm::dot(trans[0], glm::dvec4(b - origin, 1)) * scx;
      double bby = glm::dot(trans[1], glm::dvec4(b - origin, 1)) * scy;
      double bbz = glm::dot(trans[2], glm::dvec4(b - origin, 1)) * scz;
      b = origin + glm::dvec3(bbx * trans[0]) + glm::dvec3(bby * trans[1]) +
          glm::dvec3(bbz * trans[2]);
      double ccx = glm::dot(trans[0], glm::dvec4(c - origin, 1)) * scx;
      double ccy = glm::dot(trans[1], glm::dvec4(c - origin, 1)) * scy;
      double ccz = glm::dot(trans[2], glm::dvec4(c - origin, 1)) * scz;
      c = origin + glm::dvec3(ccx * trans[0]) + glm::dvec3(ccy * trans[1]) +
          glm::dvec3(ccz * trans[2]);
    }
    AddFace(a, b, c);
  }

  AddPart(geom);
}

/*
TODO: change over to copy indices with a base vertex position added and append points including normals
      so normals don't get invalidated 
*/
void IfcGeometry::MergeGeometry(const fuzzybools::Geometry& geom)
{
  for (uint32_t i = 0; i < geom.numFaces; i++)
  {
    fuzzybools::Face f = geom.GetFace(i);
    glm::dvec3 a = geom.GetPoint(f.i0);
    glm::dvec3 b = geom.GetPoint(f.i1);
    glm::dvec3 c = geom.GetPoint(f.i2);
    AddFace(a, b, c);
  }
  
  if ( wingedEdgeMesh.has_value() ) {

    wingedEdgeMesh.reset();
  }
}

uint32_t IfcGeometry::GetAllocationSize() const {
  return byteSize(fvertexData) + byteSize(vertexData) + byteSize(indexData);
}

uint32_t IfcGeometry::GetVertexDataSize() {
  return (uint32_t)vertexData.size();
}

uint32_t IfcGeometry::GetIndexData() { return (uint32_t)(size_t)&indexData[0]; }

uint32_t IfcGeometry::GetIndexDataSize() { return (uint32_t)indexData.size(); }

void IfcGeometry::ApplyTransform(const glm::dmat4& transform) {

  if (halfSpace) {
    halfSpaceOrigin = transform * glm::dvec4(halfSpaceOrigin, 1);
    halfSpaceX = transform * glm::dvec4(halfSpaceX, 1);
    halfSpaceY = transform * glm::dvec4(halfSpaceY, 1);
    halfSpaceZ = transform * glm::dvec4(halfSpaceZ, 1);
  }

  for (uint32_t index = 0; index < numPoints; ++index) {

    size_t floatIndex = index * fuzzybools::VERTEX_FORMAT_SIZE_FLOATS;

    glm::dvec4 t =
        transform *
        glm::dvec4(
            glm::dvec3(
                vertexData[floatIndex + 0],
                vertexData[floatIndex + 1],
                vertexData[floatIndex + 2]),
            1);

    glm::dvec4 n =
        transform *
        glm::dvec4(
            glm::dvec3(
                vertexData[floatIndex + 3],
                vertexData[floatIndex + 4],
                vertexData[floatIndex + 5]),
            0);

    vertexData[floatIndex + 0] = t.x;
    vertexData[floatIndex + 1] = t.y;
    vertexData[floatIndex + 2] = t.z;
    vertexData[floatIndex + 3] = n.x;
    vertexData[floatIndex + 4] = n.y;
    vertexData[floatIndex + 5] = n.z;
  }

  if ( glm::determinant( transform ) < 0 ) {

    for ( uint32_t* indexPtr = indexData.data(), 
                  * end = indexData.data() + indexData.size(); 
          indexPtr < end; 
          indexPtr += 3 ) {

      std::swap( *indexPtr, indexPtr[ 2 ] );
    }
  }

  for ( auto& localPart : part ) {
    localPart.ApplyTransform( transform );
  }

  if ( wingedEdgeMesh.has_value() ) {

    wingedEdgeMesh.reset();
  }
}

}  // namespace conway::geometry
