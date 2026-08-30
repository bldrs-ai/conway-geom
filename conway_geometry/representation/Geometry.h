/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/representation/IfcGeometry.h
 * */

// Represents a single piece of IFC Geometry

#pragma once

#include <glm/glm.hpp>
#include <optional>
#include <string>
#include <vector>
#include <unordered_map>

#include "structures/aabb_tree.h"
#include "material.h"
#include "Topology.h"

namespace conway {

  class ParseBuffer;
}

namespace conway::geometry {

template <typename T>
constexpr size_t byteSize(const std::vector<T> &data) {
  return data.size() * sizeof(T);
}


struct Geometry {

  std::vector< glm::dvec3 >                vertices;
  std::vector< Triangle >                  triangles;
  std::vector< TriangleEdges >             triangle_edges;
  std::vector< Edge >                      edges;
  std::unordered_map< uint64_t, uint32_t > edge_map;

  /**
   * Per-triangle-corner analytic surface normals, three per triangle in the
   * corner order of `triangles`, in this geometry's own space.
   *
   * INVARIANT: either empty, or exactly `3 * triangles.size()` long. Every
   * mutator that changes the triangle list either maintains that (MakeTriangle
   * appends zeros, DeleteTriangle mirrors the swap-and-pop, ReverseFace mirrors
   * the corner swap) or clears the whole array to fall back — see the
   * per-function notes. A ZERO vector at a corner means "no analytic normal
   * here"; Reify() falls back to the triangle's face normal for that corner, so
   * a partially-populated geometry is legal and degrades locally.
   *
   * Why corners rather than vertices: the weld merges the coincident vertices
   * of two faces that meet at an edge, so a per-vertex array could only hold
   * one of the two faces' normals. Corners survive the weld with their face
   * identity intact, which is what makes the crease split fall out for free
   * (bldrs-ai/conway#667).
   */
  std::vector< glm::vec3 >                 corner_normals;

  std::optional< AABBTree > bvh;
  //new additions 0.0.44
  bool halfSpace = false;

  glm::dvec3 center          = glm::dvec3(0,0,0);
  glm::dvec3 halfSpaceX      = glm::dvec3(1, 0, 0);
  glm::dvec3 halfSpaceY      = glm::dvec3(0, 1, 0);
  glm::dvec3 halfSpaceZ      = glm::dvec3(0, 0, 1);
  glm::dvec3 halfSpaceOrigin = glm::dvec3(0, 0, 0);
  //end new additions 0.0.44

  //std::vector<Geometry> getParts();
 // fuzzybools::AABB getAABB() const;
 // glm::dvec3 GetAABBCenter() const;
  //bool normalized = false;

  uint32_t GetAllocationSize() const;

  bool IsEmpty() const { return triangles.empty() || vertices.empty(); } 

  Geometry() {}

  Geometry( Geometry&& ) = default;
  Geometry( const Geometry& ) = default;

  Geometry& operator=( const Geometry& ) = default;
  Geometry& operator=( Geometry&& ) = default;

  box3 GetAABB() const;

  const glm::dvec3& GetPoint( uint32_t index ) const { return vertices[ index ]; }

  void ReverseFace( uint32_t index );
  void ReverseFaces();
  uint32_t GetVertexData();

  void AppendGeometry( const Geometry &geom );
 
  void AppendWithScalingTransform(
    const Geometry& geom,
    const glm::dmat4& trans = glm::dmat4(1),
    double scx = 1,
    double scy = 1,
    double scz = 1,
    const glm::dvec3& origin = glm::dvec3(0, 0, 0) );

  void AppendWithTransform( const Geometry &geom, const glm::dmat4x4& transform );
  uint32_t GetVertexDataSize();
  uint32_t GetIndexData();
  uint32_t GetIndexDataSize();
  glm::dvec3 Normalize();
  uint32_t GetVertexCount() const { return static_cast< uint32_t >( vertices.size() ); }
  uint32_t GetTriangleCount() const { return static_cast< uint32_t >( triangles.size() ); }
  void ApplyTransform( const glm::dmat4x4& transform );

  void ApplyRescale( const glm::dvec3& scale, const glm::dvec3& origin = glm::dvec3( 0 ) );

  void ExtractVertices( const ParseBuffer& buffer );

  void ExtractTriangles( const ParseBuffer& buffer );

  void ExtractVerticesAndTriangles( const ParseBuffer& verticesBuffer, const ParseBuffer& triangleBuffer );

  Geometry Clone();

  std::string GeometryToObj( const std::string& preamble = "" );

  bool HasConnectivity() const {

    return hasConnectivity_;
  }

  void ClearConnectivity( bool releaseMemory = false ) {

    if ( hasConnectivity_ ) {

      edge_map.clear();
      edges.clear();
      triangle_edges.clear();

      hasConnectivity_ = false;
      cleanedUp_       = false;
    }

    if ( releaseMemory ) {

      if ( edge_map.bucket_count() > 0 ) {

        std::unordered_map< uint64_t, uint32_t > empty;

        std::swap( edge_map, empty );
      }

      edges.shrink_to_fit();
      triangle_edges.shrink_to_fit();
    }
  }

  void EnableConnectivity();

  void MakeTriangle( uint32_t a, uint32_t b, uint32_t c, uint32_t index );

  /** Cleanup this mesh for CSG */
  void Cleanup( bool forSubtract = false );

  void MarkedCleanedup() { cleanedUp_ = true; }

  void Clear();

  /** Construct a BVH for this. */
  AABBTree& MakeBVH() {

    if ( !bvh.has_value() ) {
      bvh.emplace( *this );
    }

    return *bvh;
  }

  std::optional< uint32_t > GetEdge( uint32_t v0, uint32_t v1 ) const;

  void DeleteTriangle( uint32_t index );

  uint32_t MakeVertex( const glm::dvec3& value ) ;

  uint32_t MakeVertex( double x, double y, double z );

  uint32_t MakeTriangle( uint32_t a, uint32_t b, uint32_t c );

  /**
   * Record the analytic surface normals of one triangle's three corners,
   * activating `corner_normals` on first use (every other corner in the
   * geometry then reads as zero, i.e. "fall back to the face normal").
   *
   * Normals need not be unit length, only correctly directed and non-zero;
   * Reify() normalizes. Passing a zero vector for a corner leaves that corner
   * on the face-normal fallback.
   *
   * @param triangleIndex Index into `triangles`, as returned by MakeTriangle.
   * @param n0 Analytic normal at corner 0.
   * @param n1 Analytic normal at corner 1.
   * @param n2 Analytic normal at corner 2.
   */
  void SetCornerNormals(
    uint32_t triangleIndex,
    const glm::vec3& n0,
    const glm::vec3& n1,
    const glm::vec3& n2 ) {

    if ( corner_normals.empty() ) {

      // First use: size to the whole triangle list so MakeTriangle's cheap
      // append keeps the invariant from here on. O(n) once, not per call.
      corner_normals.resize( triangles.size() * 3, glm::vec3( 0.0f ) );
    }

    uint32_t base = triangleIndex * 3;

    corner_normals[ base     ] = n0;
    corner_normals[ base + 1 ] = n1;
    corner_normals[ base + 2 ] = n2;
  }

  /**
   * Drop the analytic normals, returning this geometry to Reify()'s
   * dihedral-guessing fallback. Called by every operation that invalidates
   * them — CSG, non-uniform rescale, scaling appends.
   */
  void ClearCornerNormals() {

    corner_normals.clear();
    corner_normals.shrink_to_fit();
  }

  /**
   * Mirror a triangle-list swap-and-pop in `corner_normals`. No-op when the
   * analytic normals are not in use.
   *
   * @param from Triangle index being moved (the back of the list).
   * @param to Triangle index being overwritten.
   */
  /**
   * Rotate the corner normals of triangles in `[fromTriangle, toTriangle)` by a
   * transform's normal matrix (inverse transpose of the upper-left 3x3), so a
   * baked placement leaves them pointing outward.
   *
   * A singular linear part has no normal matrix; the range is zeroed rather
   * than filled with NaN, which puts those corners back on the face-normal
   * fallback.
   *
   * @param transform The transform applied to the positions.
   * @param fromTriangle First triangle in the range.
   * @param toTriangle One past the last triangle in the range.
   */
  void TransformCornerNormalRange(
    const glm::dmat4& transform,
    size_t fromTriangle,
    size_t toTriangle ) {

    if ( corner_normals.empty() || fromTriangle >= toTriangle ) {

      return;
    }

    glm::dmat3 linear( transform );

    bool singular = std::abs( glm::determinant( linear ) ) < DBL_EPSILON;

    glm::dmat3 normalMatrix =
      singular ? glm::dmat3( 1.0 ) : glm::transpose( glm::inverse( linear ) );

    for ( size_t corner = fromTriangle * 3, end = toTriangle * 3; corner < end; ++corner ) {

      if ( singular ) {

        corner_normals[ corner ] = glm::vec3( 0.0f );
        continue;
      }

      glm::dvec3 rotated = normalMatrix * glm::dvec3( corner_normals[ corner ] );

      corner_normals[ corner ] = glm::vec3( rotated );
    }
  }

  void MoveCornerNormals( uint32_t from, uint32_t to ) {

    if ( corner_normals.empty() ) {

      return;
    }

    for ( uint32_t corner = 0; corner < 3; ++corner ) {

      corner_normals[ ( to * 3 ) + corner ] = corner_normals[ ( from * 3 ) + corner ];
    }
  }

  uint32_t MakeEdge( uint32_t v1, uint32_t v2, uint32_t triangleIndex );

  void Reify( const glm::dvec3& offset = glm::dvec3( 0 ) );

  void ClearReification() {

    floatVertexData_.clear();
    floatVertexData_.shrink_to_fit();
    indexData_.clear();
    indexData_.shrink_to_fit();

    isReified_ = false;
  }

  bool isScaled = false;

 private:

  bool hasConnectivity_ = false;
  bool isReified_       = false;
  bool cleanedUp_       = false;
  bool normalized_      = false;

  glm::dvec3 previousReificationOffset_ = glm::dvec3( 0 );

  std::vector< float >    floatVertexData_;
  std::vector< uint32_t > indexData_;
};


inline void Geometry::MakeTriangle( uint32_t a, uint32_t b, uint32_t c, uint32_t index ) {

  Triangle& triangle = triangles[ index ];
  
  assert( a < vertices.size() );
  assert( b < vertices.size() );
  assert( c < vertices.size() );

  triangle.vertices[ 0 ] = a;
  triangle.vertices[ 1 ] = b;
  triangle.vertices[ 2 ] = c;

  if ( hasConnectivity_ ) {

    TriangleEdges& triangleEdges = triangle_edges[ index ];

    triangleEdges.edges[ 0 ] = MakeEdge( a, b, index );
    triangleEdges.edges[ 1 ] = MakeEdge( b, c, index );
    triangleEdges.edges[ 2 ] = MakeEdge( c, a, index );
  }
}

inline void Geometry::Clear() {

  bvh.reset();

  ClearConnectivity();

  vertices.clear();
  triangles.clear();

  ClearCornerNormals();

  halfSpace       = false;
  center          = glm::dvec3( 0, 0, 0 );
  halfSpaceX      = glm::dvec3( 1, 0, 0 );
  halfSpaceY      = glm::dvec3( 0, 1, 0 );
  halfSpaceZ      = glm::dvec3( 0, 0, 1 );
  halfSpaceOrigin = glm::dvec3( 0, 0, 0 );
}

inline std::optional< uint32_t > Geometry::GetEdge( uint32_t v0, uint32_t v1 ) const {

  auto mapIterator =
    edge_map.find( edgeCompoundID( v0, v1 ) );

  if ( mapIterator == edge_map.end() ) {
    return std::nullopt;
  }

  return std::optional< uint32_t >( mapIterator->second );
}

inline void Geometry::DeleteTriangle( uint32_t index ) {

  if ( hasConnectivity_ ) {

    const TriangleEdges& toDeleteEdges = triangle_edges[ index ];

    for ( uint32_t localEdge = 0; localEdge < 3; ++localEdge ) {

      Edge& edge = edges[ toDeleteEdges.edges[ localEdge ] ];

      if ( edge.triangles[ 0 ] == index ) {

        edge.triangles[ 0 ] = edge.triangles[ 1 ];

      }

      edge.triangles[ 1 ] = EMPTY_INDEX;
    }

    uint32_t backIndex = static_cast< uint32_t >( triangles.size() - 1 );

    if ( index != backIndex ) {

      const TriangleEdges& back = triangle_edges.back();

      for ( uint32_t localEdge = 0; localEdge < 3; ++localEdge ) {

        Edge& edge = edges[ back.edges[ localEdge ] ];

        for ( uint32_t onEdge = 0; onEdge < 2; ++onEdge ) {

          if ( edge.triangles[ onEdge ] == backIndex ) {
            edge.triangles[ onEdge ] = index;
            break;
          }
        }
      }

      triangles[ index ]      = triangles.back();
      triangle_edges[ index ] = triangle_edges.back();

      MoveCornerNormals( backIndex, index );
    }

    triangle_edges.pop_back();

  } else if (
    uint32_t backIndex = static_cast< uint32_t >( triangles.size() - 1 );
    index != backIndex ) {

    triangles[ index ] = triangles.back();

    MoveCornerNormals( backIndex, index );
  }

  triangles.pop_back();

  if ( !corner_normals.empty() ) {

    corner_normals.resize( triangles.size() * 3 );
  }

  if ( bvh.has_value() ) {

    bvh.reset();
  }
}

inline uint32_t Geometry::MakeVertex( const glm::dvec3& value ) {

  uint32_t index = static_cast<uint32_t>( vertices.size() );

  vertices.push_back( value );

  return index;
}

inline uint32_t Geometry::MakeVertex( double x, double y, double z ) {

  uint32_t index = static_cast<uint32_t>( vertices.size() );

  vertices.emplace_back( x, y, z );

  return index;
}

inline uint32_t Geometry::MakeTriangle( uint32_t a, uint32_t b, uint32_t c ) {

  assert( a < vertices.size() );
  assert( b < vertices.size() );
  assert( c < vertices.size() );

  uint32_t index = static_cast< uint32_t >( triangles.size() );

  Triangle& triangle = triangles.emplace_back();

  triangle.vertices[ 0 ] = a;
  triangle.vertices[ 1 ] = b;
  triangle.vertices[ 2 ] = c;

  // Keep the corner_normals invariant (3 per triangle) without requiring every
  // caller to know about it: a triangle appended with no analytic normal gets
  // three zeros, which Reify() reads as "use the face normal here".
  if ( !corner_normals.empty() ) {

    corner_normals.insert( corner_normals.end(), 3, glm::vec3( 0.0f ) );
  }

  if ( hasConnectivity_ ) {

    TriangleEdges& triangleEdges = triangle_edges.emplace_back();

    triangleEdges.edges[ 0 ] = MakeEdge( a, b, index );
    triangleEdges.edges[ 1 ] = MakeEdge( b, c, index );
    triangleEdges.edges[ 2 ] = MakeEdge( c, a, index );
  }

  return index;
}

inline uint32_t Geometry::MakeEdge( uint32_t v1, uint32_t v2, uint32_t triangleIndex ) {

  uint64_t edgeIdentifier = edgeCompoundID( v1, v2 );

  while ( true ) {

    auto [ mapIterator, emplaced ] = edge_map.try_emplace( edgeIdentifier, static_cast<uint32_t>( edges.size() ) );

    if ( emplaced ) {
      
      edges.push_back( Edge{ { triangleIndex, EMPTY_INDEX }, { v1, v2 } });

    }
    else {

      Edge& currentEdge = edges[ mapIterator->second ];

      // Note, this allows non manifold edges to exist,
      // but edges will not point back to triangles
      if ( currentEdge.triangles[ 0 ] == EMPTY_INDEX || currentEdge.triangles[ 1 ] == EMPTY_INDEX ) {

        currentEdge.triangles[ currentEdge.triangles[ 0 ] == EMPTY_INDEX ? 0 : 1 ] = triangleIndex;
      }
      else {

        currentEdge.non_manifold_linked_edge = static_cast< uint32_t >( edges.size() );

        // This is a work around for non-manifold cases.
        edge_map.erase( mapIterator );
        continue;

      }
    }

    return mapIterator->second;
  }
}


struct GeometryCollection {
  std::vector<Geometry *> components;
  std::vector<glm::dmat4x4> transforms;

  uint32_t materialIndex = 0;
  bool hasDefaultMaterial = true;

  void AddComponentWithTransform(Geometry *geom,
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
