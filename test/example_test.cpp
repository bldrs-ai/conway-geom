#include <TinyCppTest.hpp>
#include <glm/glm.hpp>

#include "../conway_geometry/ConwayGeometryProcessor.h"

std::vector<conway::geometry::IfcGeometry> genIndexIfcTest() {
  std::vector<conway::geometry::IfcGeometry> geometryVec;

  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
      parametersGetPolygonalFaceSetGeometry;
  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 7;

  conway::geometry::IfcGeometry geometry =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry);

  // clear
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 6;

  conway::geometry::IfcGeometry geometry2 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry2);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 5;

  conway::geometry::IfcGeometry geometry3 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry3);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 5;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 8;

  conway::geometry::IfcGeometry geometry4 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry4);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 7;

  conway::geometry::IfcGeometry geometry5 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry5);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 5;

  conway::geometry::IfcGeometry geometry6 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry6);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  parametersGetPolygonalFaceSetGeometry.numPoints = 8;
  parametersGetPolygonalFaceSetGeometry.points.resize(
      parametersGetPolygonalFaceSetGeometry.numPoints);
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
  parametersGetPolygonalFaceSetGeometry.numIndices =
      6 * parametersGetPolygonalFaceSetGeometry.indicesPerFace;
  parametersGetPolygonalFaceSetGeometry.indices.resize(
      parametersGetPolygonalFaceSetGeometry.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[0] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[1] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[2] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[4] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[5] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[6] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[8] = 2;
  parametersGetPolygonalFaceSetGeometry.indices[9] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[10] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[12] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[13] = 3;
  parametersGetPolygonalFaceSetGeometry.indices[14] = 7;
  parametersGetPolygonalFaceSetGeometry.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[16] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[17] = 1;
  parametersGetPolygonalFaceSetGeometry.indices[18] = 4;
  parametersGetPolygonalFaceSetGeometry.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersGetPolygonalFaceSetGeometry.indices[20] = 5;
  parametersGetPolygonalFaceSetGeometry.indices[21] = 6;
  parametersGetPolygonalFaceSetGeometry.indices[22] = 8;
  parametersGetPolygonalFaceSetGeometry.indices[23] = 7;

  conway::geometry::IfcGeometry geometry7 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry7);

  // free memory
  parametersGetPolygonalFaceSetGeometry.points.clear();
  parametersGetPolygonalFaceSetGeometry.indices.clear();

  return geometryVec;
}

TEST(GeometryVectorSizeTest) {
  std::vector<conway::geometry::IfcGeometry> geometryVec = genIndexIfcTest();

  ASSERT_EQ(geometryVec.size(), 7);
}
