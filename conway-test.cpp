/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <filesystem>
#include <fstream>
#include <iostream>

/*#include "include/conway-geometry.h"
#include "include/web-ifc.h"
// #include "include/web-ifc-geometry.h"*/
#include "conway_geometry/ConwayGeometryProcessor.h"
//#include "include/ifc-schema.h"
//#include "include/math/triangulate-with-boundaries.h"

namespace conway::statistics {
bool exportGltfs = false;
bool exportGlbs = false;
bool exportObjs = false;
bool exportSingleGeometry = false;
bool exportIndividualGeometryFiles = false;
bool exportDraco = false;
}  // namespace conway::statistics

std::string ReadFile(std::string filename) {
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

void writeFile(std::wstring filename, std::string data) {
  std::ofstream out(filename.c_str());
  out << data;
  out.close();
}

void DumpRefs(std::unordered_map<uint32_t, std::vector<uint32_t>> &refs) {
  std::ofstream of("refs.txt");

  int32_t prev = 0;
  for (auto &it : refs) {
    if (!it.second.empty()) {
      for (auto &i : it.second) {
        of << (((int32_t)i) - (prev));
        prev = i;
      }
    }
  }
}

struct BenchMarkResult {
  std::string file;
  long long timeMS;
  long long sizeBytes;
};

long long ms() {
  using namespace std::chrono;
  milliseconds millis =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch());

  return millis.count();
}

void genIndexIfc() {
  std::vector<conway::geometry::IfcGeometry> geometryVec;

  // taken from web ifc obj dump code
  glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                          glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));

  // from outside debugger
  // std::string content = ReadFile("../../../index.ifc");

  auto start = ms();

  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  // generated from webifc_geom output
  conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
      parametersAxis2Placement3D;
  parametersAxis2Placement3D.zAxisRef.x = 0.000;
  parametersAxis2Placement3D.zAxisRef.y = 0.000;
  parametersAxis2Placement3D.zAxisRef.z = 1.000;
  parametersAxis2Placement3D.normalizeZ = true;
  parametersAxis2Placement3D.xAxisRef.x = 1.000;
  parametersAxis2Placement3D.xAxisRef.y = 0.000;
  parametersAxis2Placement3D.xAxisRef.z = 0.000;
  parametersAxis2Placement3D.normalizeX = true;
  parametersAxis2Placement3D.position.x = 0.000;
  parametersAxis2Placement3D.position.y = 0.000;
  parametersAxis2Placement3D.position.z = 0.000;

  glm::dmat4 localPlacementMatrix =
      conwayGeometryProcessor.GetAxis2Placement3D(parametersAxis2Placement3D);

  /*printf("localPlacementMatrix[0][1] = %.3f, localPlacementMatrix[0][2] =
  %.3f, localPlacementMatrix[0][3] = %.3f, localPlacementMatrix[0][4] =
  %.3f,\n", localPlacementMatrix[0][0], localPlacementMatrix[0][1],
  localPlacementMatrix[0][2], localPlacementMatrix[0][3]);
  printf("localPlacementMatrix[1][1] = %.3f, localPlacementMatrix[1][2] = %.3f,
  localPlacementMatrix[1][3] = %.3f, localPlacementMatrix[1][4] = %.3f,\n",
  localPlacementMatrix[1][0], localPlacementMatrix[1][1],
  localPlacementMatrix[1][2], localPlacementMatrix[1][3]);
  printf("localPlacementMatrix[2][1] = %.3f, localPlacementMatrix[2][2] = %.3f,
  localPlacementMatrix[2][3] = %.3f, localPlacementMatrix[2][4] = %.3f,\n",
  localPlacementMatrix[2][0], localPlacementMatrix[2][1],
  localPlacementMatrix[2][2], localPlacementMatrix[2][3]);
  printf("localPlacementMatrix[3][1] = %.3f, localPlacementMatrix[3][2] = %.3f,
  localPlacementMatrix[3][3] = %.3f, localPlacementMatrix[3][4] = %.3f,\n",
  localPlacementMatrix[3][0], localPlacementMatrix[3][1],
  localPlacementMatrix[3][2], localPlacementMatrix[3][3]);
  */

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

  auto time = ms() - start;

  std::cout << "Processing geometry took " << time << "ms" << std::endl;

  if (conway::statistics::exportObjs) {
    std::cout << "Testing OBJ export..." << std::endl;
  }

  if (conway::statistics::exportGltfs) {
    if (conway::statistics::exportDraco) {
      std::cout << "Testing GLTF export (Draco)..." << std::endl;
    } else {
      std::cout << "Testing GLTF export..." << std::endl;
    }
  }

  if (conway::statistics::exportGlbs) {
    if (conway::statistics::exportDraco) {
      std::cout << "Testing GLB export (Draco)..." << std::endl;
    } else {
      std::cout << "Testing GLB export..." << std::endl;
    }
  }

  for (int geometryIndex = 0; geometryIndex < geometryVec.size();
       geometryIndex++) {
    size_t offset = 0;
    std::string singleObj = conwayGeometryProcessor.GeometryToObj(
        geometryVec[geometryIndex], offset, NormalizeMat);

    // filthy but just testing quick
    std::string fileNameGltf = "./";
    fileNameGltf += std::to_string(geometryIndex);
    fileNameGltf += "_conway";

    if (conway::statistics::exportGltfs &&
        conway::statistics::exportIndividualGeometryFiles) {
      if (conway::statistics::exportDraco) {
        printf("Writing GLTF (Draco)...\n");
      } else {
        printf("Writing GLTF...\n");
      }

      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometryVec[geometryIndex], false,
              conway::statistics::exportDraco, fileNameGltf, true,
              NormalizeMat);

      if (!results.success) {
        printf("Error writing GLTF...");
      }
    }

    if (conway::statistics::exportGlbs &&
        conway::statistics::exportIndividualGeometryFiles) {
      if (conway::statistics::exportDraco) {
        printf("Writing GLB (Draco)\n");
      } else {
        printf("Writing GLB...\n");
      }

      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometryVec[geometryIndex], true, conway::statistics::exportDraco,
              fileNameGltf, true, NormalizeMat);

      if (!results.success) {
        printf("Error writing GLTF...");
      }
    }

    if (conway::statistics::exportObjs &&
        conway::statistics::exportIndividualGeometryFiles) {
      // filthy but just testing quick
      std::string fileNameObj = "./";
      fileNameObj += std::to_string(geometryIndex);
      fileNameObj += "_conway.obj";

      std::wstring wsTmp(fileNameObj.begin(), fileNameObj.end());

      printf("Writing OBJ...\n");
      writeFile(wsTmp, singleObj);
    }
  }

  std::string completeObj = "";
  size_t offset = 0;

  if (conway::statistics::exportSingleGeometry) {
    conway::geometry::IfcGeometry fullGeometry;

    for (int geometryIndex = 0; geometryIndex < geometryVec.size();
         geometryIndex++) {
      fullGeometry.AddGeometry(geometryVec[geometryIndex]);
    }

    std::string fileNameGltf = "./index_ifc_full_conway";
    if (conway::statistics::exportDraco) {
      fileNameGltf += "_draco";
    }

    if (conway::statistics::exportGltfs) {
      if (conway::statistics::exportDraco) {
        printf("Writing Complete GLTF (Draco)...\n");
      } else {
        printf("Writing Complete GLTF...\n");
      }

      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              fullGeometry, false, conway::statistics::exportDraco,
              fileNameGltf, true, NormalizeMat);

      if (!results.success) {
        printf("Error writing GLTF...");
      }
    }

    if (conway::statistics::exportGlbs) {
      if (conway::statistics::exportDraco) {
        printf("Writing Complete GLB (Draco)...\n");
      } else {
        printf("Writing Complete GLB...\n");
      }

      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              fullGeometry, true, conway::statistics::exportDraco, fileNameGltf,
              true, NormalizeMat);

      if (!results.success) {
        printf("Error writing GLB.");
      }
    }

    if (conway::statistics::exportObjs) {
      std::cout << "Writing Complete OBJ..." << std::endl;
      completeObj = conwayGeometryProcessor.GeometryToObj(fullGeometry, offset,
                                                          NormalizeMat);

      std::string fileName = "./index_ifc_full_conway.obj";

      std::wstring wsTmp(fileName.begin(), fileName.end());

      writeFile(wsTmp, completeObj);
    }
  }

  geometryVec.clear();
}

int main(int argc, char *argv[]) {
  std::cout << "Hello Conway-Geom test!\n";

  if (argc < 2) {
    std::cout << "conway_geom_native Usage\n\tLaunching executable with no "
                 "arguments displays help"
              << "\n\t-obj   - Outputs geometry from index.ifc to obj file(s)"
              << "\n\t-gltf  - Outputs geometry from index.ifc to gltf file(s)"
              << "\n\t-glb   - Outputs geometry from index.ifc to glb file(s)"
              << "\n\t-draco - Applies default Draco compression"
              << "\n\t-full  - Outputs geometry from index.ifc to a single "
                 "geometry file"
              << "\n\t(-full is the default if -full and -split not specified)"
              << "\n\t-split - Outputs geometry from index.ifc to individual "
                 "geometry files"
              << "\n\t-h     - Displays help\n";
    return 0;
  }

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "-gltf") {
      conway::statistics::exportGltfs = true;
    } else if (arg == "-obj") {
      conway::statistics::exportObjs = true;
    } else if (arg == "-glb") {
      conway::statistics::exportGlbs = true;
    } else if (arg == "-draco") {
      conway::statistics::exportDraco = true;
    } else if (arg == "-full") {
      conway::statistics::exportSingleGeometry = true;
    } else if (arg == "-split") {
      conway::statistics::exportIndividualGeometryFiles = true;
    } else if (arg == "-h") {
      std::cout
          << "conway_geom_native Usage\n\tLaunching executable with no "
             "arguments displays help"
          << "\n\t-obj   - Outputs geometry from index.ifc to obj file(s)"
          << "\n\t-gltf  - Outputs geometry from index.ifc to gltf file(s)"
          << "\n\t-glb   - Outputs geometry from index.ifc to glb file(s)"
          << "\n\t-draco - Applies default Draco compression"
          << "\n\t-full  - Outputs geometry from index.ifc to a single "
             "geometry file"
          << "\n\t(-full is the default if -full and -split not specified)"
          << "\n\t-split - Outputs geometry from index.ifc to individual "
             "geometry files"
          << "\n\t-h     - Displays help\n";
      return 0;
    } else {
      std::cerr << "Error: invalid argument " << arg << std::endl;
      return 1;
    }
  }

  if (conway::statistics::exportDraco &&
      (!conway::statistics::exportGltfs && !conway::statistics::exportGlbs)) {
    std::cout << "Must choose -gltf or -glb with -draco switch." << std::endl;
  }

  if (!conway::statistics::exportIndividualGeometryFiles &&
      !conway::statistics::exportSingleGeometry) {
    conway::statistics::exportSingleGeometry = true;
  }

  // generate simple index.ifc geometry
  genIndexIfc();

  std::cout << "Done" << std::endl;
  return 0;
}
