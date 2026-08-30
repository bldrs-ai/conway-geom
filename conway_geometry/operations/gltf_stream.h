#pragma once

/*
 * The one place a reified component becomes glTF vertex attributes.
 *
 * This lives in a header, taking plain typed vectors, for two reasons that came
 * straight out of review (bldrs-ai/conway#667):
 *
 *   1. It is testable without linking the whole processor, so the rules below
 *      — bias applied exactly once, normals by the normal matrix, winding
 *      reversed for a mirroring placement — are pinned by a native test rather
 *      than by reading the writer.
 *   2. Typed `std::vector` in the signature makes the 32-bit-address bug
 *      unrepresentable here: there is no pointer to truncate.
 */

#include "representation/Geometry.h"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

namespace conway::geometry {

/**
 * Floats per vertex in the reified vertex stream: interleaved position then
 * normal, [px, py, pz, nx, ny, nz]. This is a CONTRACT, not an implementation
 * detail — Share's conway-direct loader reads the same buffer with the same
 * stride (src/viewer/ifc/flatMeshToBufferGeometry.js), so it cannot change
 * without changing the viewer in lockstep.
 */
constexpr uint32_t GLTF_STREAM_STRIDE = 6;

/**
 * Transform one component's reified stream into the collection's glTF buffers.
 *
 * Positions are placed by `transform` and then rebased by `positionBias`
 * EXACTLY ONCE — the writer hands the same bias back as the node's translation,
 * so subtracting it a second time anywhere (the Draco encoder, say) displaces
 * the whole primitive by -bias relative to everything else in the file.
 *
 * Normals go through the normal matrix (inverse transpose), not the placement,
 * or a non-uniform scale shears them off the surface. A singular linear part has
 * no normal matrix; those normals are emitted as +Z rather than as NaN, which
 * would fail glTF validation.
 *
 * A MIRRORING placement (negative determinant) gets its index triples reversed
 * — see TransformMirrors for why that, and not negated normals, is the correct
 * correction. Baking the transform into the positions is precisely what makes
 * this the writer's job: a consumer reverses winding itself only when the NODE
 * transform is mirroring, and by then this one is gone.
 *
 * @param vertexStream Interleaved [px, py, pz, nx, ny, nz] per vertex.
 * @param indexStream Three corner indices per triangle, into vertexStream.
 * @param transform Placement for this component, already including the
 *   collection transform.
 * @param positionBias Rebasing offset, subtracted once from every position.
 * @param pointOffset Index of this component's first vertex in `positions`.
 * @param positions Output positions, appended to.
 * @param normals Output normals, appended to.
 * @param indexData Output indices, appended to, rebased by pointOffset.
 */
inline void appendComponentToGltfStream(
  const std::vector< float >&    vertexStream,
  const std::vector< uint32_t >& indexStream,
  const glm::dmat4&              transform,
  const glm::dvec3&              positionBias,
  uint32_t                       pointOffset,
  std::vector< float >&          positions,
  std::vector< float >&          normals,
  std::vector< uint32_t >&       indexData ) {

  glm::dmat3 linear( transform );

  bool singularPlacement = std::abs( glm::determinant( linear ) ) < DBL_EPSILON;

  glm::dmat3 normalMatrix =
    singularPlacement ? glm::dmat3( 1.0 ) : glm::transpose( glm::inverse( linear ) );

  size_t vertexCount = vertexStream.size() / GLTF_STREAM_STRIDE;

  for ( size_t vertex = 0; vertex < vertexCount; ++vertex ) {

    const float* source = vertexStream.data() + ( vertex * GLTF_STREAM_STRIDE );

    glm::dvec3 placed =
      glm::dvec3( transform * glm::dvec4( source[ 0 ], source[ 1 ], source[ 2 ], 1 ) ) -
      positionBias;

    positions.push_back( static_cast< float >( placed.x ) );
    positions.push_back( static_cast< float >( placed.y ) );
    positions.push_back( static_cast< float >( placed.z ) );

    glm::dvec3 normal =
      singularPlacement ?
        glm::dvec3( 0.0 ) :
        normalMatrix * glm::dvec3( source[ 3 ], source[ 4 ], source[ 5 ] );

    double length = glm::length( normal );

    if ( !std::isfinite( length ) || length < DBL_EPSILON ) {

      normals.push_back( 0.0f );
      normals.push_back( 0.0f );
      normals.push_back( 1.0f );

    } else {

      normal /= length;

      normals.push_back( static_cast< float >( normal.x ) );
      normals.push_back( static_cast< float >( normal.y ) );
      normals.push_back( static_cast< float >( normal.z ) );
    }
  }

  bool mirrored = TransformMirrors( transform );

  size_t triangleCount = indexStream.size() / 3;

  for ( size_t triangle = 0; triangle < triangleCount; ++triangle ) {

    const uint32_t* corners = indexStream.data() + ( triangle * 3 );

    if ( mirrored ) {

      indexData.push_back( corners[ 2 ] + pointOffset );
      indexData.push_back( corners[ 1 ] + pointOffset );
      indexData.push_back( corners[ 0 ] + pointOffset );

    } else {

      indexData.push_back( corners[ 0 ] + pointOffset );
      indexData.push_back( corners[ 1 ] + pointOffset );
      indexData.push_back( corners[ 2 ] + pointOffset );
    }
  }

  // Any trailing partial triangle is dropped rather than emitted: Reify()
  // always writes three indices per triangle, so a remainder would mean a
  // corrupt stream, and emitting it would make the primitive's index count
  // indivisible by three — invalid glTF for a consumer to even read.
}

}  // namespace conway::geometry
