#include <filesystem>
#include <fstream>
#include <iostream>

#include "conway_geometry/ConwayGeometryProcessor.h"

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

#ifdef _WIN32
#define STRING_TYPE std::wstring
#else
#define STRING_TYPE std::string
#endif

void writeFile(STRING_TYPE filename, std::string data) {
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

  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
      paramsGetPolygonalFaceSetGeometry;

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 76.0000;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 76.0000;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 76.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 76.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 86.0000;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 86.0000;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 86.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 86.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test1\n");

  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  printf("test2\n");

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 3;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 7;

  conway::geometry::IfcGeometry geometry =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry);

  // clear
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 48.0000;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 48.0000;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 48.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 48.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 58.0000;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 58.0000;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 58.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 58.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test2\n");

  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 1;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 6;

  conway::geometry::IfcGeometry geometry2 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry2);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 10.0000;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 10.0000;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 10.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 10.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test3\n");
  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 1;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 5;

  conway::geometry::IfcGeometry geometry3 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry3);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 0.00232305;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -12.647637;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 0.00232305;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -24.098042;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 0.00232305;
  paramsGetPolygonalFaceSetGeometry.points[2].y = -24.098042;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 0.00232305;
  paramsGetPolygonalFaceSetGeometry.points[3].y = -12.647637;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 10.002323;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -12.647637;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 10.002323;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -24.098042;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 10.002323;
  paramsGetPolygonalFaceSetGeometry.points[6].y = -24.098042;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 10.002323;
  paramsGetPolygonalFaceSetGeometry.points[7].y = -12.647637;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test4\n");
  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 3;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 5;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 8;

  conway::geometry::IfcGeometry geometry4 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry4);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 24.0000;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 24.0000;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 24.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 24.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 34.0000;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 34.0000;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 34.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 30.0000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 34.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test5\n");
  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 3;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 7;

  conway::geometry::IfcGeometry geometry5 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry5);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 47.859639;
  paramsGetPolygonalFaceSetGeometry.points[0].y = 0.973380;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 47.859639;
  paramsGetPolygonalFaceSetGeometry.points[1].y = 0.973380;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 47.859639;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 12.423785;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 47.859639;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 12.423785;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 57.859639;
  paramsGetPolygonalFaceSetGeometry.points[4].y = 0.973380;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 57.859639;
  paramsGetPolygonalFaceSetGeometry.points[5].y = 0.973380;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 57.859639;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 12.423785;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 15.000000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 57.859639;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 12.423785;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.000000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test6\n");
  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }
  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 1;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 5;

  conway::geometry::IfcGeometry geometry6 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry6);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

  paramsGetPolygonalFaceSetGeometry.points.resize(8);
  paramsGetPolygonalFaceSetGeometry.points[0].x = 62.0000;
  paramsGetPolygonalFaceSetGeometry.points[0].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[0].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[1].x = 62.0000;
  paramsGetPolygonalFaceSetGeometry.points[1].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[1].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[2].x = 62.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[2].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[3].x = 62.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[3].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[4].x = 72.0000;
  paramsGetPolygonalFaceSetGeometry.points[4].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[4].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[5].x = 72.0000;
  paramsGetPolygonalFaceSetGeometry.points[5].y = -11.4504;
  paramsGetPolygonalFaceSetGeometry.points[5].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.points[6].x = 72.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[6].z = 15.0000;

  paramsGetPolygonalFaceSetGeometry.points[7].x = 72.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].y = 0.0000;
  paramsGetPolygonalFaceSetGeometry.points[7].z = 0.0000;

  paramsGetPolygonalFaceSetGeometry.indicesPerFace = 4;
  paramsGetPolygonalFaceSetGeometry.faces.resize(6);

  printf("test7\n");
  for (size_t faceIndex = 0;
       faceIndex < paramsGetPolygonalFaceSetGeometry.faces.size();
       faceIndex++) {
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].face_starts.push_back(0);
    paramsGetPolygonalFaceSetGeometry.faces[faceIndex].indices.resize(4);
  }

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[0] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[2] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[0].indices[3] = 4;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[1] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[2] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[1].indices[3] = 6;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[0] = 2;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[1] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[2].indices[3] = 3;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[0] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[1] = 3;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[2] = 7;
  paramsGetPolygonalFaceSetGeometry.faces[3].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[0] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[1] = 1;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[2] = 4;
  paramsGetPolygonalFaceSetGeometry.faces[4].indices[3] = 8;

  // IFCINDEXEDPOLYGONALFACE
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[0] = 5;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[1] = 6;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[2] = 8;
  paramsGetPolygonalFaceSetGeometry.faces[5].indices[3] = 7;

  conway::geometry::IfcGeometry geometry7 =
      conwayGeometryProcessor.getPolygonalFaceSetGeometry(
          paramsGetPolygonalFaceSetGeometry);

  geometryVec.push_back(geometry7);

  // free memory
  paramsGetPolygonalFaceSetGeometry.points.clear();
  paramsGetPolygonalFaceSetGeometry.faces.clear();

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

      std::vector<conway::geometry::IfcGeometry> geometrySingle;
      geometrySingle.push_back(geometryVec[0]);
      std::vector<conway::geometry::Material> materials;
      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometrySingle, materials, false, conway::statistics::exportDraco,
              fileNameGltf, true, NormalizeMat);

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

      std::vector<conway::geometry::IfcGeometry> geometrySingle;
      geometrySingle.push_back(geometryVec[0]);
      std::vector<conway::geometry::Material> materials;
      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometrySingle, materials, true, conway::statistics::exportDraco,
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

      printf("Writing OBJ...\n");
#ifdef _WIN32
      std::wstring wsTmp(fileNameObj.begin(), fileNameObj.end());
      writeFile(wsTmp, singleObj);
#else
      writeFile(fileNameObj, singleObj);
#endif
    }
  }

  std::string completeObj = "";
  size_t offset = 0;

  if (conway::statistics::exportSingleGeometry) {
    conway::geometry::IfcGeometry fullGeometry;

    /*for (int geometryIndex = 0; geometryIndex < geometryVec.size();
         geometryIndex++) {
      fullGeometry.AppendGeometry(geometryVec[geometryIndex]);
    }*/

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

      std::vector<conway::geometry::Material> materials;
      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometryVec, materials, false, conway::statistics::exportDraco,
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

      std::vector<conway::geometry::Material> materials;
      conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
          conwayGeometryProcessor.GeometryToGltf(
              geometryVec, materials, true, conway::statistics::exportDraco,
              fileNameGltf, true, NormalizeMat);

      if (!results.success) {
        printf("Error writing GLB.");
      }
    }

    if (conway::statistics::exportObjs) {
      std::cout << "Writing Complete OBJ..." << std::endl;
      completeObj = conwayGeometryProcessor.GeometryToObj(fullGeometry, offset,
                                                          NormalizeMat);

      std::string fileName = "./index_ifc_full_conway.obj";

#ifdef _WIN32
      std::wstring wsTmp(fileName.begin(), fileName.end());
      writeFile(wsTmp, singleObj);
#else
      writeFile(fileName, completeObj);
#endif
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
