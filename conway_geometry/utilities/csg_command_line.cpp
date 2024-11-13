// Change this to 1 to randomize voxel sampling - CS
#define RANDOMIZED_VOXELS 0

#include "structures/winged_edge.h"
#include "csg/csg.h"

#include <vector>

#define TINYOBJLOADER_USE_DOUBLE
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include <stdio.h>
#include <iostream>
#include <format>

#define GLM_EXT_INCLUDED 1
#include <glm/gtx/transform.hpp>
#include <glm/gtx/euler_angles.hpp>

using namespace conway;
using namespace conway::geometry;

bool loadOBJ( const char* name, WingedEdgeMesh< glm::dvec3 >& output ) {

  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  std::string warn;
  std::string err;

  bool result = tinyobj::LoadObj( &attrib, &shapes, &materials, &warn, &err, name );

  if (!err.empty()) {

    printf( "obj error: %s\n", err.c_str() );
  }

  if (!warn.empty()) {

      printf("obj warning: %s\n", err.c_str());
  }

  if (!result) {

    return false;
  }

  if (shapes.size() != 1) {

    printf( "warning: obj doesn't have single shape (this may cause strange behaviour)\n" );
  }
  
  std::vector< glm::dvec3 >& outVertices = output.vertices;
  std::vector< double >&     inVertices  = attrib.vertices;

  size_t vertexCount = inVertices.size() / 3;

  outVertices.reserve(vertexCount);

  const double* inVertData = inVertices.data();

  for (size_t where = 0; where < vertexCount; ++where, inVertData += 3) {
          
    outVertices.emplace_back( inVertData[ 0 ], inVertData[ 1 ], inVertData[ 2 ] );
  }

  for ( const tinyobj::shape_t& shape : shapes ) {

    const std::vector< tinyobj::index_t >& indices = shape.mesh.indices;

    size_t triangleCount = indices.size() / 3;

    const tinyobj::index_t* inIndexData = indices.data();

    for ( size_t where = 0; where < triangleCount; ++where, inIndexData += 3 ) {

      output.makeTriangle( inIndexData[ 0 ].vertex_index, inIndexData[ 1 ].vertex_index, inIndexData[ 2 ].vertex_index );
    }
  }

  return true;
}

struct PointCloudPoint {

  glm::dvec3 position;
  glm::fvec4 color;
};

void dumpPointCloud(const char* outputName, std::vector< PointCloudPoint >& data) {

   std::ofstream file( outputName );

   file <<
    std::format(
    "ply\n"
    "format ascii 1.0\n"
    "element vertex {}\n"
    "property float x\n"
    "property float y\n"
    "property float z\n"
    "property uchar red\n"
    "property uchar green\n"
    "property uchar blue\n"
    "property uchar alpha\n"
    "end_header\n", data.size());

  for ( const PointCloudPoint& point : data ) {

    glm::uvec4 color( glm::clamp( glm::round( glm::pow( point.color, glm::fvec4( 1.0f / 2.2f ) ) * 255.0f ), 0.0f, 255.0f ) );

    const glm::dvec3& position = point.position;

    file <<
      std::format(
        "{} {} {} {} {} {} {}\n",
        position.x,
        position.y,
        position.z,
        color.r,
        color.g,
        color.b,
        color.a );
  }

  file.close();
}

void usage() {

  printf("usage: csg_command_line\n\t"
         "( \"v\" <obj to voxelise> <resolution> <output name> ) |\n\t"
         "( <\"u\" | \"i\" | \"d\"> <obj operand a> <obj operand b> <output name(no extension)> [ x y z (of b) ] [ x-angle y-angle z-angle (of b)] [ scale (of b)] )\n" );
}

int voxelise(int argc, const char* argv[]) {

  if (argc < 5) {
  
    usage();
    return -1;
  }

  WingedEdgeMesh< glm::dvec3 > mesh;

  size_t resolution = atoll(argv[2]);

  if (resolution < 4) {

      printf("Resolution %s is not a number or too small (< 4)\n", argv[2]);
      return -1;
  }

  if (!loadOBJ(argv[3], mesh)) {

      return -1;
  }

  mesh.makeBVH();

  AABBTree& bvh = *mesh.bvh;

  bvh.dipoles( mesh );

  box3 meshBox = bvh.bounds();

  if (meshBox.min.x > meshBox.max.x) {

      printf("Mesh has no valid triangles\n");
      return -1;
  }

  glm::dvec3 boxInterval = meshBox.interval();

  double longestEdge = std::max({ boxInterval.x, boxInterval.y, boxInterval.z });

  if (longestEdge <= 0) {

      printf("Mesh has zero extents\n");
  }

  double step = longestEdge / (resolution - 3);

  glm::uvec3 boxResolution =
    glm::clamp(
      glm::uvec3( glm::ceil( boxInterval / step ) ),
      3u,
      static_cast<unsigned int>(resolution) - 3u) +
    glm::uvec3(3);

  glm::dvec3 start = meshBox.min - glm::dvec3(step);

  std::ofstream file( std::string(argv[4]) + ".ply" );

  file <<
    std::format(
      "ply\n"
      "format ascii 1.0\n"
      "element vertex {}\n"
      "property float x\n"
      "property float y\n"
      "property float z\n"
      "property uchar red\n"
      "property uchar green\n"
      "property uchar blue\n"
      "property uchar alpha\n"
      "end_header\n", (boxResolution.z * boxResolution.y * boxResolution.x));

  for (uint32_t z = 0; z < boxResolution.z; ++z) {
    for (uint32_t y = 0; y < boxResolution.y; ++y) {
      for (uint32_t x = 0; x < boxResolution.x; ++x) {

#if ( RANDOMIZED_VOXELS == 0 )
        glm::dvec3 coord = start + glm::dvec3(x, y, z) * step;
#else
        glm::dvec3 coord =
          start + 
          ( boxInterval + 2 * step ) * 
          glm::dvec3(
              double( rand() ) / double( RAND_MAX ),
              double( rand() ) / double( RAND_MAX ),
              double( rand() ) / double( RAND_MAX ) );
#endif
        double gwn = bvh.gwn( mesh, coord );

        if (gwn < -0.5) {

          file <<
            std::format(
              "{} {} {} 0 255 0 255\n",
              coord.x,
              coord.y,
              coord.z);
        }
        else if (gwn > 0.5) {

          file <<
            std::format(
              "{} {} {} 0 0 255 255\n",
              coord.x,
              coord.y,
              coord.z);
        }
        else {
          file <<
            std::format(
              "{} {} {} 255 128 0 0\n",
              coord.x,
              coord.y,
              coord.z);
        }
      }
    }
  }

  file.close();

  return 0;
}

int csg( int argc, const char* argv[] ) {

  if (argc < 5) {

    usage();
    return -1;
  }

  CSG::Operation operation{};

  bool append = false;

  switch ( argv[1][0] ) {

  case 'u':

    printf( "union\n" );
    operation = CSG::Operation::UNION;
    break;

  case 'i':

    printf( "intersection\n" );
    operation = CSG::Operation::INTERSECTION;
    break;

  case 'd':

    printf( "difference\n" );
    operation = CSG::Operation::DIFFERENCE;
    break;

  case 'a':

    printf( "append\n" );
    append = true;
    break;

  }

  bool performTransform = false;

  glm::dmat4x4 transform( 1 );

  if ( argc >= 8 ) {

    performTransform = true;
  
    glm::dvec3 translation;

    translation.x = atof( argv[ 5 ] );
    translation.y = atof( argv[ 6 ] );
    translation.z = atof( argv[ 7 ] );

    transform = glm::translate( transform, translation );
  }

  if ( argc >= 11 ) {

    glm::dvec3 euler;

    euler.x = glm::degrees( atof( argv[ 8 ] ) );
    euler.y = glm::degrees( atof( argv[ 9 ] ) );
    euler.z = glm::degrees( atof( argv[ 10 ] ) );

    transform *= glm::eulerAngleYXZ( euler.y, euler.x, euler.z );
  }
    
  if ( argc > 11 ) {

    double scale = atof( argv[ 12 ] );

    transform = glm::scale( transform, glm::dvec3( scale ) );
  }

  WingedEdgeMesh< glm::dvec3 > a;
  WingedEdgeMesh< glm::dvec3 > b;
  WingedEdgeMesh< glm::dvec3 > output;

  if ( !loadOBJ( argv[ 2 ], a ) ) {
    return -1;
  }

  if ( !loadOBJ( argv[ 3 ], b ) ) {
    return -1;
  }

  if ( performTransform ) {

    for ( glm::dvec3& vertex : b.vertices ) {

      vertex = glm::dvec3( transform * glm::dvec4( vertex, 1 ) );
    }
  }

  if ( !append ) {
  
    CSG runner;

    runner.run( operation, a, b, output );

    std::ofstream outputFile( ( std::string( argv[ 4 ] ) + ".ply" ), std::ofstream::out | std::ofstream::binary );

    std::string novelVertices = runner.dumpNovelVertices();

    outputFile.write( novelVertices.data(), novelVertices.size() );

    outputFile.close();
  }
  else {

    output = std::move( a );

    output.append( b );
  }

  std::string result = output.dumpToOBJ();

  std::ofstream outputFile( ( std::string( argv[ 4 ] ) + ".obj" ), std::ofstream::out | std::ofstream::binary);

  outputFile.write( result.data(), result.size() );
  outputFile.close();

  return 0;
}

int main( int argc, const char* argv[] ) {

  if ( argc >= 2 ) {
  
    std::string mode = argv[ 1 ];
  
    if ( mode == "v" || mode == "voxelize" || mode == "voxelise" ) {

      return voxelise( argc, argv );
    }
    else if (mode == "u" || mode == "union" || mode == "i" || mode == "intersection" || mode == "d" || mode == "difference" || mode == "a" || mode == "append") {

      return csg( argc, argv );
    }
  }
  else {

    usage();
  }

  return -1;
}