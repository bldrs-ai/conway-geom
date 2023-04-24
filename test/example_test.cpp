#include <TinyCppTest.hpp>
#include <glm/glm.hpp>

#include "../conway_geometry/ConwayGeometryProcessor.h"

std::vector<conway::geometry::IfcGeometry> genIndexIfcTest() {
  std::vector<conway::geometry::IfcGeometry> geometryVec;

  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet
      parametersPolygonalFaceset;
  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 76.0000;
  parametersPolygonalFaceset.points[0].y = -11.4504;
  parametersPolygonalFaceset.points[0].z = 0.0000;

  parametersPolygonalFaceset.points[1].x = 76.0000;
  parametersPolygonalFaceset.points[1].y = -11.4504;
  parametersPolygonalFaceset.points[1].z = 15.0000;

  parametersPolygonalFaceset.points[2].x = 76.0000;
  parametersPolygonalFaceset.points[2].y = 0.0000;
  parametersPolygonalFaceset.points[2].z = 15.0000;

  parametersPolygonalFaceset.points[3].x = 76.0000;
  parametersPolygonalFaceset.points[3].y = 0.0000;
  parametersPolygonalFaceset.points[3].z = 0.0000;

  parametersPolygonalFaceset.points[4].x = 86.0000;
  parametersPolygonalFaceset.points[4].y = -11.4504;
  parametersPolygonalFaceset.points[4].z = 15.0000;

  parametersPolygonalFaceset.points[5].x = 86.0000;
  parametersPolygonalFaceset.points[5].y = -11.4504;
  parametersPolygonalFaceset.points[5].z = 0.0000;

  parametersPolygonalFaceset.points[6].x = 86.0000;
  parametersPolygonalFaceset.points[6].y = 0.0000;
  parametersPolygonalFaceset.points[6].z = 15.0000;

  parametersPolygonalFaceset.points[7].x = 86.0000;
  parametersPolygonalFaceset.points[7].y = 0.0000;
  parametersPolygonalFaceset.points[7].z = 0.0000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 5;
  parametersPolygonalFaceset.indices[5] = 2;
  parametersPolygonalFaceset.indices[6] = 1;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 2;
  parametersPolygonalFaceset.indices[9] = 5;
  parametersPolygonalFaceset.indices[10] = 7;
  parametersPolygonalFaceset.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 6;
  parametersPolygonalFaceset.indices[17] = 1;
  parametersPolygonalFaceset.indices[18] = 4;
  parametersPolygonalFaceset.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 5;
  parametersPolygonalFaceset.indices[21] = 6;
  parametersPolygonalFaceset.indices[22] = 8;
  parametersPolygonalFaceset.indices[23] = 7;

  conway::geometry::IfcGeometry geometry =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry);

  // clear
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 48.0000;
  parametersPolygonalFaceset.points[0].y = -11.4504;
  parametersPolygonalFaceset.points[0].z = 0.0000;

  parametersPolygonalFaceset.points[1].x = 48.0000;
  parametersPolygonalFaceset.points[1].y = -11.4504;
  parametersPolygonalFaceset.points[1].z = 30.0000;

  parametersPolygonalFaceset.points[2].x = 48.0000;
  parametersPolygonalFaceset.points[2].y = 0.0000;
  parametersPolygonalFaceset.points[2].z = 30.0000;

  parametersPolygonalFaceset.points[3].x = 48.0000;
  parametersPolygonalFaceset.points[3].y = 0.0000;
  parametersPolygonalFaceset.points[3].z = 0.0000;

  parametersPolygonalFaceset.points[4].x = 58.0000;
  parametersPolygonalFaceset.points[4].y = -11.4504;
  parametersPolygonalFaceset.points[4].z = 0.0000;

  parametersPolygonalFaceset.points[5].x = 58.0000;
  parametersPolygonalFaceset.points[5].y = -11.4504;
  parametersPolygonalFaceset.points[5].z = 30.0000;

  parametersPolygonalFaceset.points[6].x = 58.0000;
  parametersPolygonalFaceset.points[6].y = 0.0000;
  parametersPolygonalFaceset.points[6].z = 30.0000;

  parametersPolygonalFaceset.points[7].x = 58.0000;
  parametersPolygonalFaceset.points[7].y = 0.0000;
  parametersPolygonalFaceset.points[7].z = 0.0000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 2;
  parametersPolygonalFaceset.indices[5] = 1;
  parametersPolygonalFaceset.indices[6] = 5;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 7;
  parametersPolygonalFaceset.indices[9] = 3;
  parametersPolygonalFaceset.indices[10] = 2;
  parametersPolygonalFaceset.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 4;
  parametersPolygonalFaceset.indices[17] = 8;
  parametersPolygonalFaceset.indices[18] = 5;
  parametersPolygonalFaceset.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 5;
  parametersPolygonalFaceset.indices[21] = 8;
  parametersPolygonalFaceset.indices[22] = 7;
  parametersPolygonalFaceset.indices[23] = 6;

  conway::geometry::IfcGeometry geometry2 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry2);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 0.0000;
  parametersPolygonalFaceset.points[0].y = -11.4504;
  parametersPolygonalFaceset.points[0].z = 0.0000;

  parametersPolygonalFaceset.points[1].x = 0.0000;
  parametersPolygonalFaceset.points[1].y = -11.4504;
  parametersPolygonalFaceset.points[1].z = 30.0000;

  parametersPolygonalFaceset.points[2].x = 0.0000;
  parametersPolygonalFaceset.points[2].y = 0.0000;
  parametersPolygonalFaceset.points[2].z = 30.0000;

  parametersPolygonalFaceset.points[3].x = 0.0000;
  parametersPolygonalFaceset.points[3].y = 0.0000;
  parametersPolygonalFaceset.points[3].z = 0.0000;

  parametersPolygonalFaceset.points[4].x = 10.0000;
  parametersPolygonalFaceset.points[4].y = -11.4504;
  parametersPolygonalFaceset.points[4].z = 0.0000;

  parametersPolygonalFaceset.points[5].x = 10.0000;
  parametersPolygonalFaceset.points[5].y = -11.4504;
  parametersPolygonalFaceset.points[5].z = 30.0000;

  parametersPolygonalFaceset.points[6].x = 10.0000;
  parametersPolygonalFaceset.points[6].y = 0.0000;
  parametersPolygonalFaceset.points[6].z = 30.0000;

  parametersPolygonalFaceset.points[7].x = 10.0000;
  parametersPolygonalFaceset.points[7].y = 0.0000;
  parametersPolygonalFaceset.points[7].z = 0.0000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 2;
  parametersPolygonalFaceset.indices[5] = 1;
  parametersPolygonalFaceset.indices[6] = 5;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 7;
  parametersPolygonalFaceset.indices[9] = 3;
  parametersPolygonalFaceset.indices[10] = 2;
  parametersPolygonalFaceset.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 4;
  parametersPolygonalFaceset.indices[17] = 8;
  parametersPolygonalFaceset.indices[18] = 5;
  parametersPolygonalFaceset.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 8;
  parametersPolygonalFaceset.indices[21] = 7;
  parametersPolygonalFaceset.indices[22] = 6;
  parametersPolygonalFaceset.indices[23] = 5;

  conway::geometry::IfcGeometry geometry3 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry3);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 0.00232305;
  parametersPolygonalFaceset.points[0].y = -12.647637;
  parametersPolygonalFaceset.points[0].z = 0.000000;

  parametersPolygonalFaceset.points[1].x = 0.00232305;
  parametersPolygonalFaceset.points[1].y = -24.098042;
  parametersPolygonalFaceset.points[1].z = 0.000000;

  parametersPolygonalFaceset.points[2].x = 0.00232305;
  parametersPolygonalFaceset.points[2].y = -24.098042;
  parametersPolygonalFaceset.points[2].z = 15.000000;

  parametersPolygonalFaceset.points[3].x = 0.00232305;
  parametersPolygonalFaceset.points[3].y = -12.647637;
  parametersPolygonalFaceset.points[3].z = 15.000000;

  parametersPolygonalFaceset.points[4].x = 10.002323;
  parametersPolygonalFaceset.points[4].y = -12.647637;
  parametersPolygonalFaceset.points[4].z = 0.000000;

  parametersPolygonalFaceset.points[5].x = 10.002323;
  parametersPolygonalFaceset.points[5].y = -24.098042;
  parametersPolygonalFaceset.points[5].z = 0.000000;

  parametersPolygonalFaceset.points[6].x = 10.002323;
  parametersPolygonalFaceset.points[6].y = -24.098042;
  parametersPolygonalFaceset.points[6].z = 15.000000;

  parametersPolygonalFaceset.points[7].x = 10.002323;
  parametersPolygonalFaceset.points[7].y = -12.647637;
  parametersPolygonalFaceset.points[7].z = 15.000000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 2;
  parametersPolygonalFaceset.indices[5] = 1;
  parametersPolygonalFaceset.indices[6] = 5;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 2;
  parametersPolygonalFaceset.indices[9] = 6;
  parametersPolygonalFaceset.indices[10] = 7;
  parametersPolygonalFaceset.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 1;
  parametersPolygonalFaceset.indices[17] = 4;
  parametersPolygonalFaceset.indices[18] = 8;
  parametersPolygonalFaceset.indices[19] = 5;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 7;
  parametersPolygonalFaceset.indices[21] = 6;
  parametersPolygonalFaceset.indices[22] = 5;
  parametersPolygonalFaceset.indices[23] = 8;

  conway::geometry::IfcGeometry geometry4 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry4);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 24.0000;
  parametersPolygonalFaceset.points[0].y = -11.4504;
  parametersPolygonalFaceset.points[0].z = 0.0000;

  parametersPolygonalFaceset.points[1].x = 24.0000;
  parametersPolygonalFaceset.points[1].y = -11.4504;
  parametersPolygonalFaceset.points[1].z = 30.0000;

  parametersPolygonalFaceset.points[2].x = 24.0000;
  parametersPolygonalFaceset.points[2].y = 0.0000;
  parametersPolygonalFaceset.points[2].z = 30.0000;

  parametersPolygonalFaceset.points[3].x = 24.0000;
  parametersPolygonalFaceset.points[3].y = 0.0000;
  parametersPolygonalFaceset.points[3].z = 0.0000;

  parametersPolygonalFaceset.points[4].x = 34.0000;
  parametersPolygonalFaceset.points[4].y = -11.4504;
  parametersPolygonalFaceset.points[4].z = 30.0000;

  parametersPolygonalFaceset.points[5].x = 34.0000;
  parametersPolygonalFaceset.points[5].y = -11.4504;
  parametersPolygonalFaceset.points[5].z = 0.0000;

  parametersPolygonalFaceset.points[6].x = 34.0000;
  parametersPolygonalFaceset.points[6].y = 0.0000;
  parametersPolygonalFaceset.points[6].z = 30.0000;

  parametersPolygonalFaceset.points[7].x = 34.0000;
  parametersPolygonalFaceset.points[7].y = 0.0000;
  parametersPolygonalFaceset.points[7].z = 0.0000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 5;
  parametersPolygonalFaceset.indices[5] = 2;
  parametersPolygonalFaceset.indices[6] = 1;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 2;
  parametersPolygonalFaceset.indices[9] = 5;
  parametersPolygonalFaceset.indices[10] = 7;
  parametersPolygonalFaceset.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 6;
  parametersPolygonalFaceset.indices[17] = 1;
  parametersPolygonalFaceset.indices[18] = 4;
  parametersPolygonalFaceset.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 5;
  parametersPolygonalFaceset.indices[21] = 6;
  parametersPolygonalFaceset.indices[22] = 8;
  parametersPolygonalFaceset.indices[23] = 7;

  conway::geometry::IfcGeometry geometry5 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry5);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 47.859639;
  parametersPolygonalFaceset.points[0].y = 0.973380;
  parametersPolygonalFaceset.points[0].z = 0.000000;

  parametersPolygonalFaceset.points[1].x = 47.859639;
  parametersPolygonalFaceset.points[1].y = 0.973380;
  parametersPolygonalFaceset.points[1].z = 15.000000;

  parametersPolygonalFaceset.points[2].x = 47.859639;
  parametersPolygonalFaceset.points[2].y = 12.423785;
  parametersPolygonalFaceset.points[2].z = 15.000000;

  parametersPolygonalFaceset.points[3].x = 47.859639;
  parametersPolygonalFaceset.points[3].y = 12.423785;
  parametersPolygonalFaceset.points[3].z = 0.000000;

  parametersPolygonalFaceset.points[4].x = 57.859639;
  parametersPolygonalFaceset.points[4].y = 0.973380;
  parametersPolygonalFaceset.points[4].z = 0.000000;

  parametersPolygonalFaceset.points[5].x = 57.859639;
  parametersPolygonalFaceset.points[5].y = 0.973380;
  parametersPolygonalFaceset.points[5].z = 15.000000;

  parametersPolygonalFaceset.points[6].x = 57.859639;
  parametersPolygonalFaceset.points[6].y = 12.423785;
  parametersPolygonalFaceset.points[6].z = 15.000000;

  parametersPolygonalFaceset.points[7].x = 57.859639;
  parametersPolygonalFaceset.points[7].y = 12.423785;
  parametersPolygonalFaceset.points[7].z = 0.000000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 2;
  parametersPolygonalFaceset.indices[5] = 1;
  parametersPolygonalFaceset.indices[6] = 5;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 7;
  parametersPolygonalFaceset.indices[9] = 3;
  parametersPolygonalFaceset.indices[10] = 2;
  parametersPolygonalFaceset.indices[11] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 4;
  parametersPolygonalFaceset.indices[17] = 8;
  parametersPolygonalFaceset.indices[18] = 5;
  parametersPolygonalFaceset.indices[19] = 1;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 8;
  parametersPolygonalFaceset.indices[21] = 7;
  parametersPolygonalFaceset.indices[22] = 6;
  parametersPolygonalFaceset.indices[23] = 5;

  conway::geometry::IfcGeometry geometry6 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry6);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  parametersPolygonalFaceset.numPoints = 8;
  parametersPolygonalFaceset.points.resize(
      parametersPolygonalFaceset.numPoints);
  parametersPolygonalFaceset.points[0].x = 62.0000;
  parametersPolygonalFaceset.points[0].y = -11.4504;
  parametersPolygonalFaceset.points[0].z = 0.0000;

  parametersPolygonalFaceset.points[1].x = 62.0000;
  parametersPolygonalFaceset.points[1].y = -11.4504;
  parametersPolygonalFaceset.points[1].z = 15.0000;

  parametersPolygonalFaceset.points[2].x = 62.0000;
  parametersPolygonalFaceset.points[2].y = 0.0000;
  parametersPolygonalFaceset.points[2].z = 15.0000;

  parametersPolygonalFaceset.points[3].x = 62.0000;
  parametersPolygonalFaceset.points[3].y = 0.0000;
  parametersPolygonalFaceset.points[3].z = 0.0000;

  parametersPolygonalFaceset.points[4].x = 72.0000;
  parametersPolygonalFaceset.points[4].y = -11.4504;
  parametersPolygonalFaceset.points[4].z = 15.0000;

  parametersPolygonalFaceset.points[5].x = 72.0000;
  parametersPolygonalFaceset.points[5].y = -11.4504;
  parametersPolygonalFaceset.points[5].z = 0.0000;

  parametersPolygonalFaceset.points[6].x = 72.0000;
  parametersPolygonalFaceset.points[6].y = 0.0000;
  parametersPolygonalFaceset.points[6].z = 15.0000;

  parametersPolygonalFaceset.points[7].x = 72.0000;
  parametersPolygonalFaceset.points[7].y = 0.0000;
  parametersPolygonalFaceset.points[7].z = 0.0000;

  parametersPolygonalFaceset.indicesPerFace = 4;
  parametersPolygonalFaceset.numIndices =
      6 * parametersPolygonalFaceset.indicesPerFace;
  parametersPolygonalFaceset.indices.resize(
      parametersPolygonalFaceset.numIndices);
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[0] = 1;
  parametersPolygonalFaceset.indices[1] = 2;
  parametersPolygonalFaceset.indices[2] = 3;
  parametersPolygonalFaceset.indices[3] = 4;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[4] = 5;
  parametersPolygonalFaceset.indices[5] = 2;
  parametersPolygonalFaceset.indices[6] = 1;
  parametersPolygonalFaceset.indices[7] = 6;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[8] = 2;
  parametersPolygonalFaceset.indices[9] = 5;
  parametersPolygonalFaceset.indices[10] = 7;
  parametersPolygonalFaceset.indices[11] = 3;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[12] = 4;
  parametersPolygonalFaceset.indices[13] = 3;
  parametersPolygonalFaceset.indices[14] = 7;
  parametersPolygonalFaceset.indices[15] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[16] = 6;
  parametersPolygonalFaceset.indices[17] = 1;
  parametersPolygonalFaceset.indices[18] = 4;
  parametersPolygonalFaceset.indices[19] = 8;
  // IFCINDEXEDPOLYGONALFACE
  parametersPolygonalFaceset.indices[20] = 5;
  parametersPolygonalFaceset.indices[21] = 6;
  parametersPolygonalFaceset.indices[22] = 8;
  parametersPolygonalFaceset.indices[23] = 7;

  conway::geometry::IfcGeometry geometry7 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          parametersPolygonalFaceset);

  geometryVec.push_back(geometry7);

  // free memory
  parametersPolygonalFaceset.points.clear();
  parametersPolygonalFaceset.indices.clear();

  return geometryVec;
}

TEST(GeometryVectorSizeTest) {
  std::vector<conway::geometry::IfcGeometry> geometryVec = genIndexIfcTest();

  ASSERT_EQ(geometryVec.size(), 7);
}
