/*
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.cpp
 * */

// Implementation for IfcGeometry

#include "Geometry.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "structures/hash_functions.h"
#include "csg/csg_utils.h"
#include "csg/csg.h"
#include "structures/vertex_welder.h"

namespace conway::geometry {

constexpr double SMOOTHING_GROUP_ANGLE = glm::radians( 40.0 );

VertexWelder welder;

void Geometry::ReverseFace( uint32_t index ) {

  Triangle& triangle = triangles[ index ];

  std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );

  if ( bvh.has_value() ) {

    bvh->clearDipoles();
  }
}

void Geometry::ReverseFaces() {

  for ( Triangle& triangle : triangles ) {

    std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
  }
  
  if ( bvh.has_value() ) {

    bvh->clearDipoles();
  }
}

void Geometry::EnableConnectivity() {

  if ( !hasConnectivity_ ) {

    hasConnectivity_ = true;

    edges.clear();
    edge_map.clear();
    triangle_edges.clear();

    // Guess using Euler's formula approximation.
    triangle_edges.reserve( 3 * vertices.size() );

    uint32_t triangleIndex = 0;

    for ( const Triangle& triangle : triangles ) {

      triangle_edges.emplace_back();

      TriangleEdges& triangleEdges = triangle_edges.back();

      triangleEdges.edges[ 0 ] = MakeEdge( triangle.vertices[ 0 ], triangle.vertices[ 1 ], triangleIndex );
      triangleEdges.edges[ 1 ] = MakeEdge( triangle.vertices[ 1 ], triangle.vertices[ 2 ], triangleIndex );
      triangleEdges.edges[ 2 ] = MakeEdge( triangle.vertices[ 2 ], triangle.vertices[ 0 ], triangleIndex );

      ++triangleIndex;
    }

    assert( triangles.size() == triangle_edges.size() );
  }
}

void Geometry::Reify() {

  if ( isReified_ ) {

    return;
  }

  if ( !cleanedUp_ ) {

    welder.weld( *this, 0 );
  }

  isReified_ = true;

  PrefixSumMap vertexTriangles;

  uint32_t vertexCount   = static_cast< uint32_t >( vertices.size() );
  uint32_t triangleCount = static_cast< uint32_t >( triangles.size() );

  std::vector< glm::dvec3 > faceNormals;

  faceNormals.resize( triangles.size() );

  for ( uint32_t triangleIndex = 0; triangleIndex < triangleCount; ++triangleIndex ) {
  
    const Triangle& triangle = triangles[ triangleIndex ];

    const glm::dvec3& v0 = vertices[ triangle.vertices[ 0 ] ];
    const glm::dvec3& v1 = vertices[ triangle.vertices[ 1 ] ];
    const glm::dvec3& v2 = vertices[ triangle.vertices[ 2 ] ];

    if ( is_zero_area_triangle( v0, v1, v2 ) ) {

      faceNormals[ triangleIndex ] = glm::dvec3( 0 );
      continue;
    }
    
    faceNormals[ triangleIndex ] = glm::cross( v1 - v0, v2 - v0 );
  }

  // Create a prefix sum map that maps vertices to their respective triangles.
  vertexTriangles.construct(
    triangles,
    static_cast< uint32_t >( vertices.size() ),
    []( const Triangle& triangle, uint32_t vertexIndex ) {

      return triangle.vertices[ vertexIndex ];
    },
    3 );

  indexData_.clear();
  floatVertexData_.clear();

  indexData_.resize( 3 * triangles.size(), EMPTY_INDEX );
  floatVertexData_.reserve( triangles.size() * 3 );

  uint32_t outputVertexIndex = 0;
  
  double cosineCutoff = cos( SMOOTHING_GROUP_ANGLE );

  // Greedy vertex smoothing.
  for ( uint32_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {

    std::span< const uint32_t > trianglesPerVertex = vertexTriangles.get( vertexIndex );

    for ( size_t triangleInSpan = 0, end = trianglesPerVertex.size(); triangleInSpan < end; ++triangleInSpan ) {

      uint32_t triangleIndex = trianglesPerVertex[ triangleInSpan ];

      const Triangle& triangle         = triangles[ triangleIndex ];
      uint32_t        vertexInTriangle = triangle.vertexInTriangle( vertexIndex );      
      uint32_t        indexDataOffset  = triangleIndex * 3 + vertexInTriangle;

      if ( indexData_[ indexDataOffset ] != EMPTY_INDEX ) {

        continue;
      }

      glm::fvec3 vertex( vertices[ vertexIndex ] );

      floatVertexData_.push_back( vertex.x );
      floatVertexData_.push_back( vertex.y );
      floatVertexData_.push_back( vertex.z );
    
      const glm::dvec3& normal     = faceNormals[ triangleIndex ];
      double            doubleArea = glm::length( normal );

      indexData_[ indexDataOffset ] = outputVertexIndex;

      if ( doubleArea == 0 || isnan( doubleArea ) ) {

        // Push back a fake normal for this degenerate case.
        floatVertexData_.push_back( 0 );
        floatVertexData_.push_back( 0 );
        floatVertexData_.push_back( 1 );

        ++outputVertexIndex;
        continue;
      }

      glm::dvec3 accumulator = normal;

      // Probe forwards looking for matches within the cutoff angle of this normal greedily.
      for ( size_t nextTriangleInSpan = triangleInSpan + 1; nextTriangleInSpan < end; ++nextTriangleInSpan ) {

        uint32_t        nextTriangleIndex     = trianglesPerVertex[ nextTriangleInSpan ];
        const Triangle& nextTriangle          = triangles[ nextTriangleIndex ];
        uint32_t        vertexInNextTriangle  = nextTriangle.vertexInTriangle( vertexIndex );
        uint32_t        nextindexDataOffset   = nextTriangleIndex * 3 + vertexInNextTriangle;

        // This has already been merged with another triangle's normal.
        if ( indexData_[ nextindexDataOffset ] != EMPTY_INDEX ) {

          continue;
        }

        const glm::dvec3& opposingNormal     = faceNormals[ nextTriangleIndex ];
        double            doubleOpposingArea = glm::length( opposingNormal );

        if ( doubleOpposingArea < DBL_EPSILON || ( doubleOpposingArea * doubleArea ) < DBL_EPSILON ) {

          continue;
        }

        double cosBetweenNormals = glm::dot( normal, opposingNormal ) / ( doubleOpposingArea * doubleArea );

        if ( cosBetweenNormals < cosineCutoff ) {

          continue;
        }

        indexData_[ nextindexDataOffset ] = outputVertexIndex;

        accumulator += opposingNormal;
      }

      // Normal weighted by triangle area.
      glm::fvec3 outputNormal = glm::normalize( accumulator );

      floatVertexData_.push_back( outputNormal.x );
      floatVertexData_.push_back( outputNormal.y );
      floatVertexData_.push_back( outputNormal.z );

      ++outputVertexIndex;
    }
  }
}

uint32_t Geometry::GetVertexData() {
  
  Reify();

  if ( floatVertexData_.empty() ) {
    return 0;
  }

  return (uint32_t)(size_t)floatVertexData_.data();
}

std::string Geometry::GeometryToObj(
  const std::string& preamble) {
  
  bool isReified = isReified_;

  if ( !isReified ) {

    Reify();
  }

  std::string obj;
  
  obj.reserve( floatVertexData_.size() * 6 + triangles.size() * 32 );  // preallocate memory

  const char *vFormat = "v %.6f %.6f %.6f\nvn %.6f %.6f %.6f\n";
  const char *fFormat = "f %zu//%zu %zu//%zu %zu//%zu\n";

  obj.append( preamble );

  for (size_t i = 0, end = floatVertexData_.size(); i < end; i += 6 ) {

    char vBuffer[128];
    snprintf(
      vBuffer,
      sizeof( vBuffer ),
      vFormat,
      floatVertexData_[ i + 0 ],
      floatVertexData_[ i + 1 ],
      floatVertexData_[ i + 2 ],
      floatVertexData_[ i + 3 ],
      floatVertexData_[ i + 4 ],
      floatVertexData_[ i + 5 ] );
    obj.append(vBuffer);    
  }

  // for (uint32_t i = 0; i < numPoints; ++i) {
  //   glm::dvec3 t = GetPoint(i);
  //   glm::dvec3 n = GetNormal(i);
  //   char vBuffer[128];
  //   snprintf(vBuffer, sizeof(vBuffer), vFormat, t.x, t.y, t.z, n.x, n.y, n.z);
  //   obj.append(vBuffer);    
  // }

  for (size_t i = 0, end = indexData_.size(); i < end; i += 3 ) {

    uint32_t f1 = indexData_[ i + 0 ] + 1;
    uint32_t f2 = indexData_[ i + 1 ] + 1;
    uint32_t f3 = indexData_[ i + 2 ] + 1;

    char fBuffer[ 128 ];
    snprintf( fBuffer, sizeof(fBuffer), fFormat, f1, f1, f2, f2, f3, f3 );
    obj.append(fBuffer);
  }

  if ( !isReified ) {

    ClearReification();
  }

  return obj;
}

Geometry Geometry::Clone() { return *this; }

void Geometry::AppendWithTransform(
  const Geometry &geom,
  const glm::dmat4x4& transform ) {

  size_t currentVertexCount = vertices.size();

  AppendGeometry( geom );

  for ( glm::dvec3* where = vertices.data() + currentVertexCount, *end = vertices.data() + vertices.size(); where < end; ++where ) {

    *where = transform * glm::dvec4( *where, 1 );
  }
}

void Geometry::AppendWithScalingTransform(
    const Geometry& geom,
    const glm::dmat4& trans,
    double scx,
    double scy,
    double scz,
    const glm::dvec3& origin ) {

  ClearReification();

  bvh.reset();

  cleanedUp_ = false;

  size_t currentVertexCount = vertices.size();
  
  AppendGeometry( geom );

  if (scx != 1 || scy != 1 || scz != 1) {

    for (
      glm::dvec3
        *where = vertices.data() + currentVertexCount,
        *end = vertices.data() + vertices.size();
        where < end;
        ++where ) {
    
      double x = glm::dot( trans[0], glm::dvec4( *where - origin, 1 ) ) * scx;
      double y = glm::dot( trans[1], glm::dvec4( *where - origin, 1 ) ) * scy;
      double z = glm::dot( trans[2], glm::dvec4( *where - origin, 1 ) ) * scz;
     
      *where =
        origin +
        glm::dvec3( x * trans[ 0 ] ) +
        glm::dvec3( y * trans[ 1 ] ) +
        glm::dvec3( z * trans[ 2 ] );
    }
  }
}

void Geometry::Cleanup() {

  if ( !cleanedUp_ ) {

    if ( !halfSpace ) {

      for ( glm::dvec3& from : vertices ) {

        from *= exp2( 24 );
        from = glm::round( from );
        from *= exp2( -24 );
      }
    }

    welder.weld( *this, exp2( -23 ) );

    EnableConnectivity();

    assert( triangles.size() == triangle_edges.size() );

    if ( !halfSpace ) {
      CSG cleaner;

      cleaner.clean( *this );
      EnableConnectivity();
    }
    
    cleanedUp_ = true;

  } else {

    // for ( glm::dvec3& from : vertices ) {

    //   from *= exp2( 24 );
    //   from = glm::round( from );
    //   from *= exp2( -24 );
    // }

    // welder.weld( *this, exp2( -23 ) );

    EnableConnectivity();
  }
}

glm::dvec3 Geometry::Normalize() {

  glm::dvec3 centre = GetAABB().centre();
 
  if ( !normalized_ ) {

    for ( glm::dvec3& vertex : vertices ) {

      vertex -= centre;
    }

    normalized_ = true;
  }

  return center;
}

void Geometry::AppendGeometry( const Geometry &geom ) {

  bool hasConnectivity = hasConnectivity_;

  if ( hasConnectivity ) {

    ClearConnectivity();
  }

  uint32_t maxIndex = static_cast< uint32_t >( vertices.size() );

  vertices.insert( vertices.end(), geom.vertices.begin(), geom.vertices.end() );

  size_t maxTriangle = triangles.size();

  triangles.insert( triangles.end(), geom.triangles.begin(), geom.triangles.end() );

  for ( auto where = triangles.begin() + maxTriangle, end = triangles.end(); where < end; ++where ) {

    uint32_t (&triangleVertices)[ 3 ] = where->vertices;

    triangleVertices[ 0 ] += maxIndex;
    triangleVertices[ 1 ] += maxIndex;
    triangleVertices[ 2 ] += maxIndex;
  }

  if ( hasConnectivity ) {

    EnableConnectivity();
  }

  ClearReification();
  bvh.reset();
  cleanedUp_ = false;
  normalized_ = false;

  //TODO: see if this is needed 
  //AddPart(geom);
}


// void Geometry::AddGeometry(
//   const fuzzybools::Geometry& geom, const glm::dmat4& trans,
//                               double scx, double scy, double scz,
//                               const glm::dvec3& origin) {
  
//   if ( wingedEdgeMesh.has_value() ) {

//     wingedEdgeMesh.reset();
//   }
  
//   for (uint32_t i = 0; i < geom.numFaces; i++) {
//     fuzzybools::Face f = geom.GetFace(i);
//     glm::dvec3 a = geom.GetPoint(f.i0);
//     glm::dvec3 b = geom.GetPoint(f.i1);
//     glm::dvec3 c = geom.GetPoint(f.i2);
//     if (scx != 1 || scy != 1 || scz != 1) {
//       double aax = glm::dot(trans[0], glm::dvec4(a - origin, 1)) * scx;
//       double aay = glm::dot(trans[1], glm::dvec4(a - origin, 1)) * scy;
//       double aaz = glm::dot(trans[2], glm::dvec4(a - origin, 1)) * scz;
//       a = origin + glm::dvec3(aax * trans[0]) + glm::dvec3(aay * trans[1]) +
//           glm::dvec3(aaz * trans[2]);
//       double bbx = glm::dot(trans[0], glm::dvec4(b - origin, 1)) * scx;
//       double bby = glm::dot(trans[1], glm::dvec4(b - origin, 1)) * scy;
//       double bbz = glm::dot(trans[2], glm::dvec4(b - origin, 1)) * scz;
//       b = origin + glm::dvec3(bbx * trans[0]) + glm::dvec3(bby * trans[1]) +
//           glm::dvec3(bbz * trans[2]);
//       double ccx = glm::dot(trans[0], glm::dvec4(c - origin, 1)) * scx;
//       double ccy = glm::dot(trans[1], glm::dvec4(c - origin, 1)) * scy;
//       double ccz = glm::dot(trans[2], glm::dvec4(c - origin, 1)) * scz;
//       c = origin + glm::dvec3(ccx * trans[0]) + glm::dvec3(ccy * trans[1]) +
//           glm::dvec3(ccz * trans[2]);
//     }
//     AddFace(a, b, c);
//   }

//   //AddPart(geom);
// }

/*
TODO: change over to copy indices with a base vertex position added and append points including normals
      so normals don't get invalidated 
*/

uint32_t Geometry::GetAllocationSize() const {
  return
    static_cast< uint32_t >(
      byteSize( floatVertexData_ ) +
      byteSize( vertices ) +
      byteSize( triangles ) +
      byteSize( edges ) +
      byteSize( triangle_edges ) );
}

uint32_t Geometry::GetVertexDataSize() {
  
  Reify();

  return (uint32_t)( floatVertexData_.size() );
}

box3 Geometry::GetAABB() const {

  if ( bvh.has_value() ) {

    return bvh->bounds();
  }

  box3 result;

  for ( const glm::dvec3& vertex : vertices ) {

    result.merge( vertex );
  }

  return result;
}

uint32_t Geometry::GetIndexData() { Reify(); return (uint32_t)(size_t)indexData_.data(); }

uint32_t Geometry::GetIndexDataSize() { Reify(); return static_cast< uint32_t >( indexData_.size() ); }

void Geometry::ApplyTransform( const glm::dmat4& transform ) {

  if (halfSpace) {
    halfSpaceOrigin = transform * glm::dvec4(halfSpaceOrigin, 1);
    halfSpaceX = transform * glm::dvec4(halfSpaceX, 1);
    halfSpaceY = transform * glm::dvec4(halfSpaceY, 1);
    halfSpaceZ = transform * glm::dvec4(halfSpaceZ, 1);
  }

  for ( glm::dvec3& vertex : vertices ) {

    glm::dvec4 t = transform * glm::dvec4( vertex, 1 );

    vertex = t;
  }

  if ( glm::determinant( transform ) < 0 ) {

    for ( Triangle& triangle : triangles ) {

      std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
    }

    if ( hasConnectivity_ ) {

      // to flip in the corresponding way, we reverse the first 2 edges ((0,1),(1,2),(2,0)) becomes ((2,1),(1,0),(0,2)) which given the previous 
      // flip matches the original identity. 
      for ( TriangleEdges& triangleEdges : triangle_edges ) {

        std::swap( triangleEdges.edges[ 0 ], triangleEdges.edges[ 1 ] );
      }
    }
  }

  ClearReification();
  bvh.reset();
  cleanedUp_ = false;
}

}  // namespace conway::geometry
