#include <TinyCppTest.hpp>
#include <glm/glm.hpp>

#include "../conway_geometry/ConwayGeometryProcessor.h"

std::vector<conway::geometry::Geometry> genIndexIfcTest() {
  std::vector<conway::geometry::Geometry> geometryVec;

  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
      parametersGetPolygonalFaceSetGeometry;
  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 76.0000;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 76.0000;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 76.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 76.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 86.0000;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 86.0000;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 86.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 86.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry);

  // clear
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 48.0000;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 48.0000;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 48.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 48.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 58.0000;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 58.0000;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 58.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 58.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry2 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry2);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 10.0000;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 10.0000;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 10.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 10.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry3 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry3);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 0.00232305;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -12.647637;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 0.00232305;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -24.098042;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 0.00232305;
  parametersGetPolygonalFaceSetGeometry.points[2].y = -24.098042;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 0.00232305;
  parametersGetPolygonalFaceSetGeometry.points[3].y = -12.647637;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 10.002323;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -12.647637;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 10.002323;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -24.098042;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 10.002323;
  parametersGetPolygonalFaceSetGeometry.points[6].y = -24.098042;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 10.002323;
  parametersGetPolygonalFaceSetGeometry.points[7].y = -12.647637;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(5);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry4 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry4);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 24.0000;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 24.0000;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 24.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 24.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 34.0000;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 34.0000;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 34.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 34.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry5 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry5);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 47.859639;
  parametersGetPolygonalFaceSetGeometry.points[0].y = 0.973380;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 47.859639;
  parametersGetPolygonalFaceSetGeometry.points[1].y = 0.973380;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 47.859639;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 12.423785;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 47.859639;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 12.423785;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 57.859639;
  parametersGetPolygonalFaceSetGeometry.points[4].y = 0.973380;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 57.859639;
  parametersGetPolygonalFaceSetGeometry.points[5].y = 0.973380;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 57.859639;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 12.423785;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 15.000000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 57.859639;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 12.423785;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.000000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry6 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry6);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  parametersGetPolygonalFaceSetGeometry.points.resize(8);
  parametersGetPolygonalFaceSetGeometry.points[0].x = 62.0000;
  parametersGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[1].x = 62.0000;
  parametersGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[1].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[2].x = 62.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[2].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[3].x = 62.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[4].x = 72.0000;
  parametersGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[4].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[5].x = 72.0000;
  parametersGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  parametersGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.points[6].x = 72.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[6].z = 15.0000;

  parametersGetPolygonalFaceSetGeometry.points[7].x = 72.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  parametersGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  parametersGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  parametersGetPolygonalFaceSetGeometry.faces.resize(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[0].indices.push_back(4);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[1].indices.push_back(6);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(2);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[2].indices.push_back(3);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(3);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(7);
  parametersGetPolygonalFaceSetGeometry.faces[3].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(1);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(4);
  parametersGetPolygonalFaceSetGeometry.faces[4].indices.push_back(8);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(5);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(6);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(8);
  parametersGetPolygonalFaceSetGeometry.faces[5].indices.push_back(7);

  parametersGetPolygonalFaceSetGeometry.faces[0].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[1].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[2].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[3].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[4].face_starts.push_back(0);
  parametersGetPolygonalFaceSetGeometry.faces[5].face_starts.push_back(0);

  conway::geometry::Geometry geometry7 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry7);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.faces.clear();

  return geometryVec;
}

TEST(GeometryVectorSizeTest) {
  std::vector<conway::geometry::Geometry> geometryVec = genIndexIfcTest();

  ASSERT_EQ(geometryVec.size(), 7);
}
