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
#include "structures/alloc_telemetry.h"
#include "structures/parse_buffer.h"
#include "utilities/buffer_parse.h"

#include <unordered_map>

namespace conway::geometry {

constexpr int32_t TOLERANCE = -23;

constexpr double SMOOTHING_GROUP_ANGLE = glm::radians( 40.0 );

VertexWelder welder;

void Geometry::ReverseFace( uint32_t index ) {

  Triangle& triangle = triangles[ index ];

  std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );

  // Corner normals are indexed by corner position, so reversing the winding has
  // to carry them along or corner 0 would read the normal of the old corner 2 —
  // and they have to be NEGATED, because this API exists to invert the face's
  // outward direction, not to renumber its corners. Leaving them alone would
  // shade a reversed analytic face against its own winding, and Reify() prefers
  // the analytic vector, so the flipped face normal would not save it.
  ReverseCornerNormals( index );

  if ( bvh.has_value() ) {

    bvh->clearDipoles();
  }
}

void Geometry::ReverseFaces() {

  for ( Triangle& triangle : triangles ) {

    std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
  }

  // Same reordering AND negation as ReverseFace, for the same reason.
  if ( !corner_normals.empty() ) {

    for ( uint32_t index = 0, end = static_cast< uint32_t >( triangles.size() ); index < end; ++index ) {

      ReverseCornerNormals( index );
    }
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

void Geometry::Reify( const glm::dvec3& offset ) {

  if ( isReified_ && offset == previousReificationOffset_ ) {

    return;
  }

  previousReificationOffset_ = offset;

  floatVertexData_.clear();
  indexData_.clear();

  if ( !cleanedUp_ ) {

    conway::AllocTagScope weldTag( conway::AllocSite::VertexWeld );

    welder.weld( *this, DBL_EPSILON );
  
  } else {

    size_t triangleCursor = 0;

    while ( triangleCursor < triangles.size() ) {

      const Triangle& triangle = triangles[ triangleCursor ];

      uint32_t i0 = triangle.vertices[ 0 ];
      uint32_t i1 = triangle.vertices[ 1 ];
      uint32_t i2 = triangle.vertices[ 2 ];

      if ( i0 == i1 || i0 == i2 || i1 == i2 ) {

        DeleteTriangle( triangleCursor );
        continue;
      }
      
      const glm::dvec3& v0 = vertices[ triangle.vertices[ 0 ] ];
      const glm::dvec3& v1 = vertices[ triangle.vertices[ 1 ] ];
      const glm::dvec3& v2 = vertices[ triangle.vertices[ 2 ] ];

      if ( is_zero_area_triangle( v0, v1, v2 ) ) {

        DeleteTriangle( triangleCursor );
        continue;
      }

      ++triangleCursor;
    }
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

  bool hasAnalyticNormals = !corner_normals.empty();

  /*
   * The shading normal for one triangle corner: the analytic surface normal the
   * tessellator recorded, when there is one, and otherwise the triangle's own
   * face normal exactly as before.
   *
   * This one substitution is the whole of the bldrs-ai/conway#667 fix. The
   * grouping below is greedy — every candidate is compared against the FIRST
   * triangle's normal, never against the running average — so its result
   * depends on triangle order, and a single near-degenerate sliver whose face
   * normal points anywhere can capture or exclude a whole fan. Along a trimmed
   * face's boundary row, where the constrained edge forces exactly such
   * slivers, that produced a normal error that cycled 0..64 degrees once per
   * trim-polyline segment: the reported "regular, broken dark line". Feeding
   * the analytic normal in instead removes the guess at the source — a sliver's
   * SURFACE normal is still correct even when its face normal is garbage, and
   * every normal within one smooth face agrees to within a rounding error, so
   * the grouping stops depending on order at all. Creases still split, because
   * two faces meeting at an edge genuinely disagree by the dihedral angle.
   *
   * Note the two return paths differ in magnitude on purpose: the fallback
   * returns the UNNORMALIZED cross product, so the accumulation below stays
   * area-weighted and byte-identical to the previous behaviour for every
   * geometry with no analytic normals. The analytic path returns a unit vector,
   * giving each contributing corner equal say — which is what we want when a
   * sliver's area would otherwise erase a correct normal.
   *
   * Because the two magnitudes mean different things, they must never be summed
   * into one accumulator: a unit vector added to an area-weighted one is
   * neither, and the mix is scale-dependent, since only one of the two terms
   * carries area units. So the choice below is made once for a whole local
   * smoothing group, not per corner — see the two-attempt loop.
   */
  auto hasAnalyticCorner =
    [&]( uint32_t triangleIndex, uint32_t vertexInTriangle ) -> bool {

      if ( !hasAnalyticNormals ) {

        return false;
      }

      const glm::vec3& analytic =
        corner_normals[ ( triangleIndex * 3 ) + vertexInTriangle ];

      return analytic.x != 0.0f || analytic.y != 0.0f || analytic.z != 0.0f;
    };

  auto cornerShadingNormal =
    [&]( uint32_t triangleIndex, uint32_t vertexInTriangle, bool analyticMode ) -> glm::dvec3 {

      if ( analyticMode && hasAnalyticCorner( triangleIndex, vertexInTriangle ) ) {

        return glm::normalize(
          glm::dvec3( corner_normals[ ( triangleIndex * 3 ) + vertexInTriangle ] ) );
      }

      return faceNormals[ triangleIndex ];
    };

  // Corners merged into the group being built, as offsets into indexData_.
  // Held rather than written straight through, because a group that turns out
  // to be mixed is discarded and rebuilt on the face-normal path; declared out
  // here so the allocation is made once for the whole reification.
  std::vector< uint32_t > groupMembers;

  // Greedy vertex smoothing.
  for ( uint32_t vertexIndex = 0; vertexIndex < vertexCount; ++vertexIndex ) {

    std::span< const uint32_t > trianglesPerVertex =
      vertexTriangles.get( vertexIndex );

    for ( size_t triangleInSpan = 0, end = trianglesPerVertex.size(); triangleInSpan < end; ++triangleInSpan ) {

      uint32_t triangleIndex = trianglesPerVertex[ triangleInSpan ];

      const Triangle& triangle         = triangles[ triangleIndex ];
      uint32_t        vertexInTriangle = triangle.vertexInTriangle( vertexIndex );      
      uint32_t        indexDataOffset  = triangleIndex * 3 + vertexInTriangle;

      if ( indexData_[ indexDataOffset ] != EMPTY_INDEX ) {

        continue;
      }

      const glm::dvec3& vertex = vertices[ vertexIndex ];

      floatVertexData_.push_back( static_cast< float >( vertex.x - offset.x ) );
      floatVertexData_.push_back( static_cast< float >( vertex.y - offset.y ) );
      floatVertexData_.push_back( static_cast< float >( vertex.z - offset.z ) );

      indexData_[ indexDataOffset ] = outputVertexIndex;

      /*
       * Build the local smoothing group, in a mode chosen for the WHOLE group.
       *
       * At a welded vertex where an analytic face meets an uncovered one — a
       * cylinder running into a swept or boolean surface, say — a per-corner
       * choice would sum unit analytic normals with unnormalized area-weighted
       * face normals in a single accumulator. That result is neither semantics
       * and, because only the fallback term carries area units, it changes with
       * the model's scale. So when a group in analytic mode turns out to
       * include even one corner without an analytic normal, the whole group is
       * discarded and rebuilt on the face-normal path: that corner's answer is
       * then exactly the pre-#667 one, which is the right thing to degrade to.
       *
       * At most two attempts: the second runs with analyticMode false, which
       * cannot report a mix.
       */
      bool       analyticMode = hasAnalyticCorner( triangleIndex, vertexInTriangle );
      bool       degenerate   = false;
      glm::dvec3 accumulator( 0.0 );

      for ( int attempt = 0; attempt < 2; ++attempt ) {

        groupMembers.clear();

        glm::dvec3 normal     = cornerShadingNormal( triangleIndex, vertexInTriangle, analyticMode );
        double     doubleArea = glm::length( normal );

        if ( doubleArea == 0 || isnan( doubleArea ) ) {

          degenerate = true;
          break;
        }

        accumulator = normal;

        bool mixedGroup = false;

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

          glm::dvec3 opposingNormal =
            cornerShadingNormal( nextTriangleIndex, vertexInNextTriangle, analyticMode );
          double     doubleOpposingArea = glm::length( opposingNormal );

          if ( doubleOpposingArea < DBL_EPSILON || ( doubleOpposingArea * doubleArea ) < DBL_EPSILON ) {

            continue;
          }

          double cosBetweenNormals = glm::dot( normal, opposingNormal ) / ( doubleOpposingArea * doubleArea );

          if ( cosBetweenNormals < cosineCutoff ) {

            continue;
          }

          // A corner that JOINS the group and has no analytic normal is what
          // makes the group mixed — one that fails the cutoff is in a different
          // group and does not constrain this one's semantics.
          if ( analyticMode && !hasAnalyticCorner( nextTriangleIndex, vertexInNextTriangle ) ) {

            mixedGroup = true;
            break;
          }

          groupMembers.push_back( nextindexDataOffset );

          accumulator += opposingNormal;
        }

        if ( !mixedGroup ) {

          break;
        }

        analyticMode = false;
      }

      if ( degenerate ) {

        // Push back a fake normal for this degenerate case.
        floatVertexData_.push_back( 0 );
        floatVertexData_.push_back( 0 );
        floatVertexData_.push_back( 1 );

        ++outputVertexIndex;
        continue;
      }

      for ( uint32_t memberOffset : groupMembers ) {

        indexData_[ memberOffset ] = outputVertexIndex;
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

  Reify( previousReificationOffset_ );

  if ( floatVertexData_.empty() ) {
    return 0;
  }

  return (uint32_t)(size_t)floatVertexData_.data();
}

std::string Geometry::GeometryToObj(
  const std::string& preamble) {
  
  bool isReified = isReified_;

  if ( !isReified ) {

    Reify( previousReificationOffset_ );
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

  size_t currentVertexCount   = vertices.size();
  size_t currentTriangleCount = triangles.size();

  AppendGeometry( geom );

  for ( glm::dvec3* where = vertices.data() + currentVertexCount, *end = vertices.data() + vertices.size(); where < end; ++where ) {

    *where = transform * glm::dvec4( *where, 1 );
  }

  TransformCornerNormalRange( transform, currentTriangleCount, triangles.size() );
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

  // The per-axis scale below is applied in `trans`'s own frame rather than as a
  // single matrix, so there is no normal matrix to rotate the analytic normals
  // by. Rather than reconstruct one for a path this rare, drop them and let
  // Reify() fall back — correctness over coverage (bldrs-ai/conway#667).
  if ( ( scx != 1 || scy != 1 || scz != 1 ) && !corner_normals.empty() ) {

    ClearCornerNormals();
  }

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

void Geometry::Cleanup( bool forSubtract ) {

  // CSG preparation. The boolean that follows clips triangles and synthesises
  // new ones along the intersection curve, and those have no face to inherit an
  // analytic normal from — so a CSG result is a fallback geometry by
  // construction. Dropping them here rather than after the boolean keeps the
  // invariant simple: nothing downstream of Cleanup has to reason about
  // partially-valid normals (bldrs-ai/conway#667).
  ClearCornerNormals();

  {
    conway::AllocTagScope weldTag( conway::AllocSite::VertexWeld );

    welder.weld( *this, exp2( TOLERANCE ), forSubtract );
  }

  EnableConnectivity();

  // if ( !cleanedUp_ ) {

  //   assert( triangles.size() == triangle_edges.size() );

  //   if ( !halfSpace ) {
  //     CSG cleaner;

  //     cleaner.clean( *this );
  //     EnableConnectivity();
  //   }
    
     cleanedUp_ = true;
  // }
  // } else {

  //   // welder.weld( *this, exp2( TOLERANCE ), forSubtract );//exp2( -23 ) );

  //   // EnableConnectivity();

  // }
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

  // Corner normals have to be merged BEFORE the triangle insert below, because
  // both sides are sized off their own triangle counts. Either side may be
  // empty (no analytic normals); the missing side contributes zeros so the
  // result keeps the 3-per-triangle invariant and the un-normalled half simply
  // falls back to face normals in Reify().
  if ( !corner_normals.empty() || !geom.corner_normals.empty() ) {

    corner_normals.resize( maxTriangle * 3, glm::vec3( 0.0f ) );

    if ( geom.corner_normals.empty() ) {

      corner_normals.insert(
        corner_normals.end(), geom.triangles.size() * 3, glm::vec3( 0.0f ) );

    } else {

      corner_normals.insert(
        corner_normals.end(),
        geom.corner_normals.begin(),
        geom.corner_normals.end() );
    }
  }

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


void Geometry::ExtractVertices( const ParseBuffer& buffer ) {

  parse_vector( buffer.range(), vertices );

}

void Geometry::ExtractTriangles( const ParseBuffer& buffer ) {

  bool hasConnectivity = hasConnectivity_;

  if ( hasConnectivity ) {

    ClearConnectivity();
  }

  parse_vector( buffer.range(), triangles );

  // The triangle list is rewritten from the buffer — parse_vector appends, and
  // the decrement below then renumbers every triangle in it, so this expects to
  // own the list — while the parsed indices carry no surface information at
  // all. Either way `corner_normals` stops matching the triangle count, and
  // this is the other path (with the per-face rollback) that changes
  // `triangles` outside MakeTriangle/DeleteTriangle. Drop rather than pad:
  // "no analytic normals" is always a legal answer, a misaligned one is not
  // (bldrs-ai/conway#667).
  ClearCornerNormals();

  for ( Triangle& triangle : triangles ) {

    --triangle.vertices[ 0 ];
    --triangle.vertices[ 1 ];
    --triangle.vertices[ 2 ];
  }

  if ( hasConnectivity ) {

    EnableConnectivity();
  }
}

void Geometry::ExtractVerticesAndTriangles( const ParseBuffer& verticesBuffer, const ParseBuffer& triangleBuffer ) {

  ExtractVertices( verticesBuffer );
  ExtractTriangles( triangleBuffer );
}

uint32_t Geometry::GetVertexDataSize() {
  
  Reify( previousReificationOffset_ );

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

uint32_t Geometry::GetIndexData() { Reify( previousReificationOffset_ ); return (uint32_t)(size_t)indexData_.data(); }

uint32_t Geometry::GetIndexDataSize() { Reify( previousReificationOffset_ ); return static_cast< uint32_t >( indexData_.size() ); }


void Geometry::ApplyRescale( const glm::dvec3& scale, const glm::dvec3& origin ) {

  if (halfSpace) {
    halfSpaceOrigin = ( ( halfSpaceOrigin - origin ) * scale ) + origin;
  }

  for ( glm::dvec3& vertex : vertices ) {

    vertex = ( ( vertex - origin ) * scale ) + origin;
  }

  // A per-axis rescale shears normals unless they are pushed through the
  // inverse scale; a zero component makes that undefined outright. Not worth
  // reconstructing here — fall back.
  ClearCornerNormals();

  ClearReification();

  if ( bvh.has_value() ) {

    bvh->applyRescale( scale, origin );
  }

  bvh.reset();
  cleanedUp_ = false;
}

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

  TransformCornerNormalRange( transform, 0, triangles.size() );

  // The mirroring invariant, stated once on TransformMirrors: a mirroring
  // transform leaves the winding normal opposed to the transformed surface
  // normal, and reversing the winding — not negating the normals — is what puts
  // them back in agreement.
  if ( TransformMirrors( transform ) ) {

    for ( Triangle& triangle : triangles ) {

      std::swap( triangle.vertices[ 0 ], triangle.vertices[ 2 ] );
    }

    // Same corner-order argument as ReverseFace: the winding flip above
    // renumbers corners 0 and 2, so their normals have to move with them.
    if ( !corner_normals.empty() ) {

      for ( size_t index = 0, end = triangles.size(); index < end; ++index ) {

        std::swap( corner_normals[ index * 3 ], corner_normals[ ( index * 3 ) + 2 ] );
      }
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
