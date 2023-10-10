/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>

// #include "geometry/IfcGeometryProcessor.h"
#include "parsing/IfcLoader.h"
#include "schema/IfcSchemaManager.h"
#include "schema/ifc-schema.h"
#include "test/io_helpers.h"
#include "utility/LoaderError.h"
#include "utility/LoaderSettings.h"
#include "utility/ifcStatistics.h"

using namespace webifc::io;

long long ms() {
  using namespace std::chrono;
  milliseconds millis =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch());

  return millis.count();
}

double RandomDouble(double lo, double hi) {
  return lo + static_cast<double>(rand()) /
                  (static_cast<double>(RAND_MAX / (hi - lo)));
}

std::string ReadFile(std::string filename) {
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

void SpecificLoadTest(webifc::parsing::IfcLoader &loader,
                      webifc::geometry::IfcGeometryProcessor &geometryLoader,
                      uint64_t num) {
  auto walls = loader.GetExpressIDsWithType(webifc::schema::IFCSLAB);

  bool writeFiles = true;

  auto mesh = geometryLoader.GetMesh(num);

  if (writeFiles) {
    DumpMesh(mesh, geometryLoader, "TEST.obj");
  }
}

std::vector<webifc::geometry::IfcAlignment> GetAlignments(
    webifc::parsing::IfcLoader &loader,
    webifc::geometry::IfcGeometryProcessor &geometryLoader) {
  std::vector<webifc::geometry::IfcAlignment> alignments;

  auto type = webifc::schema::IFCALIGNMENT;

  auto elements = loader.GetExpressIDsWithType(type);

  for (unsigned int i = 0; i < elements.size(); i++) {
    auto alignment = geometryLoader.GetLoader().GetAlignment(elements[i]);
    alignment.transform(geometryLoader.GetCoordinationMatrix());
    alignments.push_back(alignment);
  }

  bool writeFiles = true;

  if (writeFiles) {
    DumpAlignment(alignments, "V_ALIGN.obj", "H_ALIGN.obj");
  }

  return alignments;
}

std::vector<webifc::geometry::IfcCrossSections> GetCrossSections3D(
    webifc::parsing::IfcLoader &loader,
    webifc::geometry::IfcGeometryProcessor &geometryLoader) {
  std::vector<webifc::geometry::IfcCrossSections> crossSections;

  std::vector<uint32_t> typeList;
  typeList.push_back(webifc::schema::IFCSECTIONEDSOLID);
  typeList.push_back(webifc::schema::IFCSECTIONEDSURFACE);
  typeList.push_back(webifc::schema::IFCSECTIONEDSOLIDHORIZONTAL);

  for (auto &type : typeList) {
    auto elements = loader.GetExpressIDsWithType(type);

    for (unsigned int i = 0; i < elements.size(); i++) {
      auto crossSection =
          geometryLoader.GetLoader().GetCrossSections3D(elements[i]);
      crossSections.push_back(crossSection);
    }
  }

  bool writeFiles = true;

  if (writeFiles) {
    DumpCrossSections(crossSections, "CrossSection.obj");
  }

  return crossSections;
}

std::string GeometryToObj(const webifc::geometry::IfcGeometry &geom,
                          size_t &offset, std::string materialName,
                          glm::dmat4 transform = glm::dmat4(1)) {
  std::string obj;
  obj.reserve(geom.numPoints * 32 + geom.numFaces * 32);  // preallocate memory

  obj.append("mtllib " + materialName + ".mtl\n");
  obj.append("usemtl " + materialName + "\n");

  const char *vFormat = "v %.6f %.6f %.6f\n";
  const char *fFormat = "f %zu// %zu// %zu//\n";

  for (uint32_t i = 0; i < geom.numPoints; ++i) {
    glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
    char vBuffer[64];
    snprintf(vBuffer, sizeof(vBuffer), vFormat, t.x, t.y, t.z);
    obj.append(vBuffer);
  }

  for (uint32_t i = 0; i < geom.numFaces; ++i) {
    size_t f1 = geom.indexData[i * 3 + 0] + 1 + offset;
    size_t f2 = geom.indexData[i * 3 + 1] + 1 + offset;
    size_t f3 = geom.indexData[i * 3 + 2] + 1 + offset;
    char fBuffer[64];
    snprintf(fBuffer, sizeof(fBuffer), fFormat, f1, f2, f3);
    obj.append(fBuffer);
  }

  offset += geom.numPoints;

  return obj;
}

/*void writeFile(std::string filename, std::string data) {
  std::ofstream out(filename.c_str());
  out << data;
  out.close();
}*/

std::string ColorToMtl(const glm::dvec3 &color,
                       const std::string &materialName) {
  std::string mtl;
  const char *mtlFormat =
      "newmtl %s\n"
      "Ka %f %f %f\n"            // Ambient color
      "Kd %f %f %f\n"            // Diffuse color
      "Ks 0.000 0.000 0.000\n";  // Specular color (set to black for simplicity)

  char mtlBuffer[128];
  snprintf(mtlBuffer, sizeof(mtlBuffer), mtlFormat, materialName.c_str(),
           color.r, color.g, color.b, color.r, color.g, color.b);
  mtl.append(mtlBuffer);

  return mtl;
}

std::string DumpMeshToStr(webifc::geometry::IfcComposedMesh &mesh,
              webifc::geometry::IfcGeometryProcessor &processor,
              size_t &offset) {
 // writeFile(filename,
      //      ToObj(mesh, processor, offset, webifc::geometry::NormalizeIFC));
   return ToObj(mesh, processor, offset, webifc::geometry::NormalizeIFC);
}

glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                        glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));
void DumpGeom(const webifc::geometry::IfcGeometry &geom, glm::dvec3 &color,
              std::string name) {
  size_t offset = 0;
  std::string fileName = "./";
  std::string objFileName = fileName + name + ".obj";
  std::string mtlFileName = fileName + name + ".mtl";
  writeFile(objFileName, GeometryToObj(geom, offset, name, NormalizeMat));

  printf("name: %s\n", name.c_str());
  writeFile(mtlFileName, ColorToMtl(color, name));
}

std::vector<webifc::geometry::IfcFlatMesh> LoadAllTest(
    webifc::parsing::IfcLoader &loader,
    webifc::geometry::IfcGeometryProcessor &geometryLoader) {
  std::vector<webifc::geometry::IfcFlatMesh> meshes;
  webifc::schema::IfcSchemaManager schema;

  webifc::geometry::IfcGeometry fullGeometry;

  for (auto type : schema.GetIfcElementList()) {
    auto elements = loader.GetExpressIDsWithType(type);

    for (unsigned int i = 0; i < elements.size(); i++) {
      auto mesh = geometryLoader.GetFlatMesh(elements[i]);

      for (auto &geom : mesh.geometries) {
        auto flatGeom = geometryLoader.GetGeometry(geom.geometryExpressID);
        glm::dvec3 color(geom.color.r, geom.color.g, geom.color.b);
        if (webifc::statistics::exportObjs) {
          webifc::geometry::IfcGeometry tmpGeometry;
          std::string fileName = std::to_string(geom.geometryExpressID);
          fileName += "_webifc";

          DumpGeom(flatGeom, color, fileName);

          std::cout << "Dumped mesh to file: " << fileName.c_str() << "\n";
        }
      }

      meshes.push_back(mesh);
    }
  }

  /*if (webifc::statistics::exportObjs) {
    webifc::geometry::IfcGeometry tmpGeometry;
    size_t offset = 0;
    std::string completeObj = "";
    for (auto type : schema.GetIfcElementList()) {
      if (type == webifc::schema::IFCOPENINGELEMENT ||
          type == webifc::schema::IFCSPACE ||
          type == webifc::schema::IFCOPENINGSTANDARDCASE) {
        continue;
      }
      auto elements = loader.GetExpressIDsWithType(type);

      for (unsigned int i = 0; i < elements.size(); i++) {
        auto mesh = geometryLoader.GetMesh(elements[i]);
        completeObj += DumpMeshToStr(mesh, geometryLoader, offset );
      }
    }

    glm::dvec3 color(0, 0, 0);
    std::string fullFileName = "fullOBJ_webifc_test.obj";
    writeFile(fullFileName, completeObj);
  }*/

  return meshes;
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

void Benchmark() {
  std::vector<BenchMarkResult> results;
  std::string path = "../../../benchmark/ifcfiles";
  for (const auto &entry : std::filesystem::directory_iterator(path)) {
    if (entry.path().extension().string() != ".ifc") {
      continue;
    }

    std::string filePath = entry.path().string();
    std::string filename = entry.path().filename().string();

    std::string content = ReadFile(filePath);

    auto start = ms();
    {
      // loader.LoadFile(content);
    }
    auto time = ms() - start;

    BenchMarkResult result;
    result.file = filename;
    result.timeMS = time;
    result.sizeBytes = entry.file_size();
    results.push_back(result);

    std::cout << "Reading " << result.file << " took " << time << "ms"
              << std::endl;
  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Results:" << std::endl;

  double avgMBsec = 0;
  for (auto &result : results) {
    double MBsec = result.sizeBytes / 1000.0 / result.timeMS;
    avgMBsec += MBsec;
    std::cout << result.file << ": " << MBsec << " MB/sec" << std::endl;
  }

  avgMBsec /= results.size();

  std::cout << std::endl;
  std::cout << "Average: " << avgMBsec << " MB/sec" << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
}

void TestTriangleDecompose() {
  const int NUM_TESTS = 100;
  const int PTS_PER_TEST = 100;
  const int EDGE_PTS_PER_TEST = 10;

  const double scaleX = 650;
  const double scaleY = 1;

  glm::dvec2 a(0, 0);
  glm::dvec2 b(scaleX, 0);
  glm::dvec2 c(0, scaleY);

  for (int i = 0; i < NUM_TESTS; i++) {
    srand(i);

    std::vector<glm::dvec2> points;

    // random points
    for (unsigned int j = 0; j < PTS_PER_TEST; j++) {
      points.push_back({RandomDouble(0, scaleX), RandomDouble(0, scaleY)});
    }

    // points along the edges
    for (unsigned int j = 0; j < EDGE_PTS_PER_TEST; j++) {
      glm::dvec2 e1 = b - a;
      glm::dvec2 e2 = c - a;
      glm::dvec2 e3 = b - c;

      points.push_back(a + e1 * RandomDouble(0, 1));
      points.push_back(a + e2 * RandomDouble(0, 1));
      points.push_back(c + e3 * RandomDouble(0, 1));
    }

    std::cout << "Start test " << i << std::endl;

    bool swapped = false;

    // webifc::IsValidTriangulation(triangles, points);

    std::vector<webifc::io::Point> pts;

    for (auto &pt : points) {
      webifc::io::Point p;
      p.x = pt.x;
      p.y = pt.y;
      pts.push_back(p);
    }
  }
}

int main(int argc, char *argv[]) {
  std::cout << "Hello web IFC test 0.44!\n";

  if (argc < 2) {
    std::cout
        << "webifc_native Usage\n\tLaunching executable with no arguments "
           "loads index.ifc from repository root and parses"
        << " ifc type information.\n\t-i /Path/To/IFC.ifc - Loads external "
           "ifc files and parses ifc type information"
        << "\n\t-objs - Outputs geometry from ifc to individudal obj files"
        << "\n\t-stats - Outputs currently supported ifc type coverage for a "
           "given ifc file"
        << "\n\t-statsv - Outputs currently supported ifc type coverage + "
           "prints type name info for a given ifc file"
        << "\n\t-t - Traces and prints IFC type information gathered from "
           "parsing a given ifc file"
        << "\n\t-h - Displays help\n";
  }

  std::string file_path;

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "-i") {
      // The next argument should be the file path
      if (i + 1 < argc) {
        file_path = argv[i + 1];
        i++;  // Skip the next argument
      } else {
        std::cerr << "Error: file path not provided" << std::endl;
        return 1;
      }
    } else if (arg == "-objs") {
      webifc::statistics::exportObjs = true;
    } else if (arg == "-stats") {
      webifc::statistics::collectStats = true;
      webifc::statistics::uniqueTypeDefs.clear();
    } else if (arg == "-statsv") {
      webifc::statistics::collectStats = true;
      webifc::statistics::verboseStats = true;
      webifc::statistics::uniqueTypeDefs.clear();
    } else if (arg == "-t") {
      webifc::statistics::shouldPrintTypeInfo = true;
    } else if (arg == "-h") {
      std::cout
          << "webifc_native Usage\n\tLaunching executable with no arguments "
             "loads index.ifc from repository root and parses"
          << " ifc type information.\n\t-i /Path/To/IFC.ifc - Loads external "
             "ifc files and parses ifc type information"
          << "\n\t-objs - Outputs geometry from ifc to individudal obj files"
          << "\n\t-stats - Outputs currently supported ifc type coverage for a "
             "given ifc file"
          << "\n\t-statsv - Outputs currently supported ifc type coverage + "
             "prints type name info for a given ifc file"
          << "\n\t-t - Traces and prints IFC type information gathered from "
             "parsing a given ifc file"
          << "\n\t-h - Displays help\n";
      return 0;
    } else {
      std::cerr << "Error: invalid argument " << arg << std::endl;
      return 1;
    }
  }

  // from outside debugger
  std::string content;
  if (file_path.empty()) {
    file_path = "../../../index.ifc";
    content = ReadFile("../../../index.ifc");
  } else {
    content = ReadFile(file_path);
  }

  if (content.empty()) {
    std::cerr << "Error: invalid filepath." << std::endl;
    return 1;
  }

  // from inside debugger
  // std::string content = ReadFile("index.ifc");

  if (webifc::statistics::shouldPrintTypeInfo) {
    std::cout << "Tracing IFC Schema for input file...\n";
  }

  webifc::utility::LoaderSettings set;
  set.COORDINATE_TO_ORIGIN = true;
  set.OPTIMIZE_PROFILES = true;

  webifc::utility::LoaderErrorHandler errorHandler;
  webifc::schema::IfcSchemaManager schemaManager;
  webifc::parsing::IfcLoader loader(set.TAPE_SIZE, set.MEMORY_LIMIT,
                                    errorHandler, schemaManager);

  auto start = ms();
  loader.LoadFile([&](char *dest, size_t sourceOffset, size_t destSize) {
    uint32_t length = std::min(content.size() - sourceOffset, destSize);
    memcpy(dest, &content[sourceOffset], length);

    return length;
  });
  auto time = ms() - start;

  std::cout << "Reading took " << time << "ms" << std::endl;

  webifc::geometry::IfcGeometryProcessor geometryLoader(
      loader, errorHandler, schemaManager, set.CIRCLE_SEGMENTS,
      set.COORDINATE_TO_ORIGIN, set.OPTIMIZE_PROFILES);

  start = ms();
  auto meshes = LoadAllTest(loader, geometryLoader);
  auto errors = errorHandler.GetErrors();
  errorHandler.ClearErrors();

  for (auto error : errors) {
    std::cout << error.expressID << " " << error.ifcType << " "
              << std::to_string((int)error.type) << " " << error.message
              << std::endl;
  }

  time = ms() - start;

  std::cout << "Generating geometry took " << time << "ms" << std::endl;

  if (webifc::statistics::collectStats) {
    std::vector<uint32_t> currentlySupportedTypes = {
        webifc::schema::IFCPLANE,
        webifc::schema::IFCAXIS2PLACEMENT2D,
        webifc::schema::IFCCARTESIANTRANSFORMATIONOPERATOR2D,
        webifc::schema::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM,
        webifc::schema::IFCPOLYLINE,
        webifc::schema::IFCLINE,
        webifc::schema::IFCINDEXEDPOLYCURVE,
        webifc::schema::IFCCIRCLE,
        webifc::schema::IFCELLIPSE,
        webifc::schema::IFCBSPLINECURVE,
        webifc::schema::IFCBSPLINECURVEWITHKNOTS,
        webifc::schema::IFCRATIONALBSPLINECURVEWITHKNOTS,
        webifc::schema::IFCAXIS1PLACEMENT,
        webifc::schema::IFCAXIS2PLACEMENT3D,
        webifc::schema::IFCLOCALPLACEMENT,
        webifc::schema::IFCCARTESIANTRANSFORMATIONOPERATOR3D,
        webifc::schema::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM,
        webifc::schema::IFCCONNECTEDFACESET,
        webifc::schema::IFCCLOSEDSHELL,
        webifc::schema::IFCOPENSHELL,
        webifc::schema::IFCFACE,
        webifc::schema::IFCPOLYLOOP,
        webifc::schema::IFCINDEXEDPOLYGONALFACE,
        webifc::schema::IFCPOLYGONALFACESET,
        webifc::schema::IFCMAPPEDITEM,
        webifc::schema::IFCBOOLEANCLIPPINGRESULT,
        webifc::schema::IFCBOOLEANRESULT,
        webifc::schema::IFCHALFSPACESOLID,
        webifc::schema::IFCPOLYGONALBOUNDEDHALFSPACE,
        webifc::schema::IFCFACEBASEDSURFACEMODEL,
        webifc::schema::IFCSHELLBASEDSURFACEMODEL,
        webifc::schema::IFCSURFACESTYLE,
        webifc::schema::IFCSURFACESTYLERENDERING,
        webifc::schema::IFCCOLOURRGB,
        webifc::schema::IFCSHAPEREPRESENTATION,
        webifc::schema::IFCINDEXEDPOLYGONALFACEWITHVOIDS,
        webifc::schema::IFCPRODUCTDEFINITIONSHAPE,
        webifc::schema::IFCEXTRUDEDAREASOLID,
        webifc::schema::IFCBUILDINGELEMENTPROXY,
        webifc::schema::IFCARBITRARYCLOSEDPROFILEDEF,
        webifc::schema::IFCSTYLEDITEM,
        webifc::schema::IFCMATERIAL,
        webifc::schema::IFCOPENINGELEMENT,
        webifc::schema::IFCMATERIALCONSTITUENT,
        webifc::schema::IFCMATERIALCONSTITUENTSET,
        webifc::schema::IFCMATERIALPROFILESET,
        webifc::schema::IFCMATERIALPROFILE,
        webifc::schema::IFCCOMPOSITEPROFILEDEF,
        webifc::schema::IFCREPRESENTATIONMAP,
        webifc::schema::IFCFACEOUTERBOUND,
        webifc::schema::IFCPRESENTATIONSTYLEASSIGNMENT,
        webifc::schema::IFCFACETEDBREP,
        webifc::schema::IFCMATERIALLAYERSET,
        webifc::schema::IFCMATERIALLAYER,
        webifc::schema::IFCMATERIALLAYERSETUSAGE,
        webifc::schema::IFCMATERIALLIST};

    uint32_t uniqueTypeDefsSize = webifc::statistics::uniqueTypeDefs.size();
    std::vector<std::pair<unsigned int, unsigned int>> supportedTypes;
    std::vector<std::pair<unsigned int, unsigned int>> unsupportedTypes;

    std::cout << "\n\n********** IFC Statistics **********" << std::endl;
    std::cout << "Input File: " << file_path.c_str() << std::endl;
    std::cout << "Unique IFC Types: " << uniqueTypeDefsSize << std::endl;

    // Copy the elements of the map into a vector
    std::vector<std::pair<unsigned int, unsigned int>> vectorUniqueTypeDefs(
        webifc::statistics::uniqueTypeDefs.begin(),
        webifc::statistics::uniqueTypeDefs.end());

    // Sort the vector in descending order by the second value
    std::sort(vectorUniqueTypeDefs.begin(), vectorUniqueTypeDefs.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    if (webifc::statistics::verboseStats) {
      std::cout << "Ifc Type\tFrequency\n" << std::endl;
    }

    for (const auto &ifcType : vectorUniqueTypeDefs) {
      if (webifc::statistics::verboseStats) {
        std::cout << schemaManager.IfcTypeCodeToType(ifcType.first) << "\t"
                  << ifcType.second << std::endl;
        // std::cout << ifcType.second << "\t"
        //       << schema.IfcTypeCodeToType(ifcType.first) << std::endl;
      }

      if (std::find(currentlySupportedTypes.begin(),
                    currentlySupportedTypes.end(),
                    ifcType.first) != currentlySupportedTypes.end()) {
        supportedTypes.push_back(ifcType);
      } else {
        unsupportedTypes.push_back(ifcType);
      }
    }

    uint32_t supportedSize = supportedTypes.size();
    uint32_t unsupportedSize = unsupportedTypes.size();
    std::cout << "\nSupported Types: " << supportedSize << std::endl;
    if (webifc::statistics::verboseStats) {
      std::cout << "Ifc Type\tFrequency\n" << std::endl;

      for (const auto &ifcType : supportedTypes) {
        std::cout << schemaManager.IfcTypeCodeToType(ifcType.first) << "\t"
                  << ifcType.second << std::endl;
        // std::cout << ifcType.second << "\t"
        //        << schema.IfcTypeCodeToType(ifcType.first) << std::endl;
      }
    }

    std::cout << "\nUnsupported Types: " << unsupportedSize << std::endl;
    if (webifc::statistics::verboseStats) {
      std::cout << "Ifc Type\tFrequency\n" << std::endl;
      for (const auto &ifcType : unsupportedTypes) {
        std::cout << schemaManager.IfcTypeCodeToType(ifcType.first) << "\t"
                  << ifcType.second << std::endl;
        // std::cout << ifcType.second << "\t"
        //         << schema.IfcTypeCodeToType(ifcType.first) << std::endl;
      }
    }

    std::cout << "\nPercentage Supported: "
              << ((double)supportedSize / (double)uniqueTypeDefsSize) * 100
              << "%" << std::endl;
    std::cout << "************************************" << std::endl;

    std::cout << "Done" << std::endl;
  }
}
