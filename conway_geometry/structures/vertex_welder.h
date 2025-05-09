#pragma once

#include "hash_functions.h"
#include "representation/Geometry.h"
#include "csg/csg_utils.h"
#include "structures/box.h"
#include "operations/morton_code_utils.h"
#include <utility>

namespace conway::geometry {

  struct VertexWelder {

    UnionFind< uint32_t >                         unified;
    std::unordered_multimap< uint32_t, uint32_t > spatial_hash;
    std::vector< uint32_t >                       remap;
    std::vector< Triangle >                       old_triangles;
    std::vector< glm::dvec3 >                     converged;
    std::vector< uint32_t >                       unique_mapping;
    std::vector< glm::dvec3 >                     area_weighted_normals;

    void weld( Geometry& toWeld, double tolerance, bool outset = false ) {

      std::vector< glm::dvec3 >& vertices = toWeld.vertices;

      unified.reset();
      spatial_hash.clear();
      remap.resize( toWeld.vertices.size() );
      old_triangles.clear();
      spatial_hash.clear();
      converged.clear();
      unique_mapping.clear();

      box3 boundingBox;

      for ( glm::dvec3 vertex : vertices ) {

        boundingBox.merge( vertex );
      }

      glm::dvec3 interval = boundingBox.interval();

      uint32_t largestComponent = largest_component( interval );

      if ( interval[ largestComponent ] == 0 ) {

        interval += DBL_EPSILON;
      }

      double inverseStep = 1.0 / std::max( interval[ largestComponent ] / double( MAX_MORTON_COMPONENT ), tolerance );

      uint32_t verticesSize = static_cast< uint32_t >( vertices.size() );

      unified.allocate( verticesSize );
      converged.reserve( verticesSize );

      for (
        uint32_t where = 0;
        where < verticesSize;
        ++where) {

        glm::dvec3& vertex = vertices[ where ];

        converged.push_back( vertex );

        glm::uvec3 coord = unpacked_coord3( vertex, boundingBox.min, inverseStep );

        for ( int32_t z = -1; z <= 1; ++z ) {

          int32_t probeZ = static_cast< int32_t >( coord.z ) + z;

          if ( probeZ < 0 || probeZ > static_cast< int32_t >( MAX_MORTON_COMPONENT ) ) {
         
            continue;
          }

          for ( int32_t y = -1; y <= 1; ++y ) {

            int32_t probeY = static_cast<int32_t>( coord.y ) + y;

            if ( probeY < 0 || probeY > static_cast<int32_t>( MAX_MORTON_COMPONENT ) ) {

              continue;
            }

            for ( int32_t x = -1; x <= 1; ++x ) {

              int32_t probeX = static_cast<int32_t>( coord.x ) + x;

              if ( probeX < 0 || probeX > static_cast<int32_t>( MAX_MORTON_COMPONENT ) ) {

                continue;
              }

              glm::uvec3 probeCoord =
                glm::uvec3( 
                  static_cast< uint32_t >( probeX ),
                  static_cast< uint32_t >( probeY ), 
                  static_cast< uint32_t >( probeZ ) );

              // Using morton makes up for "dumb" integer hash code implementations.
              uint32_t mortonProbeCoord = pack( probeCoord );

              auto [ coordMatch, endCoordMatch ] = spatial_hash.equal_range( mortonProbeCoord );

              for ( ; coordMatch != endCoordMatch; ++coordMatch) {

                const glm::dvec3& candidateCoord = vertices[ coordMatch->second ];

                if ( same_point( candidateCoord, vertex, tolerance ) ) {

                  uint32_t foundCandidate = unified.find( coordMatch->second );
                  uint32_t foundThis      = unified.find( where );

                  if ( foundCandidate != foundThis ) {

                    uint32_t mergeIndex = unified.merge( foundCandidate, foundThis );

                    converged[ mergeIndex ] = vertices[ mergeIndex ];
                  }
                }
              }
            }
          }
        }

        spatial_hash.emplace( pack( coord ), where );
      }

      unique_mapping.reserve( unified.sets() );

      unified.optimize( unique_mapping );
      remap.resize( vertices.size() );

      vertices.resize( unique_mapping.size() );

      uint32_t vertexIndex = 0;

      for ( uint32_t uniqueItem : unique_mapping ) {

        remap[ uniqueItem ]       = vertexIndex;
        vertices[ vertexIndex++ ] = converged[ uniqueItem ];
      }

      std::swap( old_triangles, toWeld.triangles );

      toWeld.triangle_edges.clear();
      toWeld.edges.clear();
      toWeld.edge_map.clear();

      for ( const Triangle& triangle : old_triangles ) {

        uint32_t i0 = remap[ unified.find( triangle.vertices[ 0 ] ) ];
        uint32_t i1 = remap[ unified.find( triangle.vertices[ 1 ] ) ];
        uint32_t i2 = remap[ unified.find( triangle.vertices[ 2 ] ) ];

        if ( i0 == i1 || i1 == i2 || i2 == i0 ) {

          continue;
        }

        if ( is_zero_area_triangle(
          vertices[ i0 ],
          vertices[ i1 ],
          vertices[ i2 ] ) ) {

          continue;
        }

        toWeld.MakeTriangle( i0, i1, i2 );
      }

      if ( outset ) {

        area_weighted_normals.clear();
        area_weighted_normals.resize( vertices.size() );

        for ( const Triangle& triangle : toWeld.triangles ) {
              
          uint32_t i0 = triangle.vertices[ 0 ];
          uint32_t i1 = triangle.vertices[ 1 ];
          uint32_t i2 = triangle.vertices[ 2 ];

          const glm::dvec3& v0 = vertices[ i0 ];
          const glm::dvec3& v1 = vertices[ i1 ];
          const glm::dvec3& v2 = vertices[ i2 ];

          glm::dvec3 normal = glm::cross( v1 - v0, v2 - v0 );

          area_weighted_normals[ i0 ] += normal;
          area_weighted_normals[ i1 ] += normal;
          area_weighted_normals[ i2 ] += normal;
        }

        for ( uint32_t vertexIndex = 0, end = vertices.size(); vertexIndex < end; ++vertexIndex ) {

          glm::dvec3 vertexOffset = glm::sign( area_weighted_normals[ vertexIndex ] );

          //printf( "%f %f %f\n", vertexOffset.x, vertexOffset.y, vertexOffset.z );

          glm::dvec3& vertex = vertices[ vertexIndex ];

          vertex += 2.0 * ( tolerance  ) * vertexOffset;

          vertex /= 0.5 * tolerance;
          vertex = glm::round( vertex );
          vertex *= 0.5 * tolerance;

        }

      } else {

        for ( uint32_t vertexIndex = 0, end = vertices.size(); vertexIndex < end; ++vertexIndex ) {

          glm::dvec3& vertex = vertices[ vertexIndex ];

          vertex /= 0.5 * tolerance;
          vertex = glm::round( vertex );
          vertex *= 0.5 * tolerance;
        }
      }
        
      size_t triangleIndex = 0;

      while ( triangleIndex < toWeld.triangles.size() ) {
              
        const Triangle& triangle = toWeld.triangles[ triangleIndex ];

        uint32_t i0 = triangle.vertices[ 0 ];
        uint32_t i1 = triangle.vertices[ 1 ];
        uint32_t i2 = triangle.vertices[ 2 ];

        if( i0 == i1 || i0 == i2 || i1 == i2 ) {

          toWeld.DeleteTriangle( static_cast< uint32_t >( triangleIndex ) );
          continue;
        }

        const glm::dvec3& v0 = vertices[ i0 ];
        const glm::dvec3& v1 = vertices[ i1 ];
        const glm::dvec3& v2 = vertices[ i2 ];

        if ( is_zero_area_triangle(
          v0,
          v1,
          v2 ) ) {

          toWeld.DeleteTriangle( static_cast< uint32_t >( triangleIndex ) );
          continue;
        }

        ++triangleIndex;
      }
    }
  };
}