/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

// #include "include/web-ifc.h"
// #include "include/conway-geometry.h"
// #include "include/ifc-schema.h"
// #include "include/math/triangulate-with-boundaries.h"
// #include "include/web-ifc-geometry.h"
#include "fuzzy/obj-exporter.h"
#include "geometry/IfcGeometryProcessor.h"
#include "parsing/IfcLoader.h"
#include "schema/IfcSchemaManager.h"
#include "utility/LoaderSettings.h"
#include "utility/Logging.h"
#include "utility/ifcStatistics.h"
// using namespace webifc::io;

webifc::schema::IfcSchemaManager schema;

long long ms() {
  using namespace std::chrono;
  milliseconds millis =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch());

  return millis.count();
}

std::string ReadFile(std::string filename) {
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

std::string GeometryToObj(const fuzzybools::Geometry &geom, size_t &offset,
                          glm::dmat4 transform = glm::dmat4(1)) {
  std::string obj;
  obj.reserve(geom.numPoints * 32 + geom.numFaces * 32);  // preallocate memory

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

void writeFile(std::wstring filename, std::string data) {
  std::ofstream out(filename.c_str());
  out << data;
  out.close();
}

glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                          glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));
void DumpGeom(const fuzzybools::Geometry &geom,
                         std::wstring filename) {
  size_t offset = 0;
  writeFile(filename, GeometryToObj(geom, offset, NormalizeMat));
}

std::vector<webifc::geometry::IfcFlatMesh> LoadAllTest(
    webifc::parsing::IfcLoader &loader,
    webifc::geometry::IfcGeometryProcessor &geometryLoader) {
  std::vector<webifc::geometry::IfcFlatMesh> meshes;

  for (auto type : schema.GetIfcElementList()) {
    auto elements = loader.GetExpressIDsWithType(type);

    for (unsigned int i = 0; i < elements.size(); i++) {
      auto mesh = geometryLoader.GetFlatMesh(elements[i]);

      if (webifc::statistics::exportObjs) {

        for (auto& geom : mesh.geometries)
        {
            auto ifc_geom = geometryLoader.GetGeometry(geom.geometryExpressID);
            fuzzybools::Geometry fbGeom;
            std::string fileName = "./";
            fileName += std::to_string(i);
            fileName += "_webifc.obj";

             std::wstring wsTmp(fileName.begin(), fileName.end());

            for (size_t j = 0; j < ifc_geom.numFaces; j++) {
            const webifc::geometry::Face &f = ifc_geom.GetFace(j);

            auto a = ifc_geom.GetPoint(f.i0);
            auto b = ifc_geom.GetPoint(f.i1);
            auto c = ifc_geom.GetPoint(f.i2);

            fbGeom.AddFace(a, b, c);
            }

            DumpGeom(fbGeom, wsTmp);

            std::cout << "Dumped mesh to file: " << fileName.c_str() << "\n";
        }
      }

      /*for (auto &geom : mesh.geometries) {
        auto flatGeom = geometryLoader.GetGeometry(geom.geometryExpressID);
      }*/

      meshes.push_back(mesh);
    }
  }

  return meshes;
}

int main(int argc, char *argv[]) {
  std::cout << "Hello web IFC test!\n";

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

  // TestTriangleDecompose();

  // return 0;

  // Benchmark();

  // return 0;

  // std::string content =
  // ReadFile(L"C:/Users/qmoya/Desktop/PROGRAMES/VSCODE/IFC.JS/issues/#83
  // processing/05111002_IFCR2_Geo_Columns_1.ifc");

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
  set.USE_FAST_BOOLS = true;

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
      loader, errorHandler, schemaManager, set.CIRCLE_SEGMENTS_HIGH,
      set.COORDINATE_TO_ORIGIN);

  start = ms();
  // SpecificLoadTest(loader, geometryLoader, 2591);
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
        webifc::schema::IFCBSPLINESURFACE,
        webifc::schema::IFCBSPLINESURFACEWITHKNOTS,
        webifc::schema::IFCRATIONALBSPLINESURFACEWITHKNOTS,
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
        webifc::schema::IFCADVANCEDFACE,
        webifc::schema::IFCPOLYLOOP,
        webifc::schema::IFCINDEXEDPOLYGONALFACE,
        webifc::schema::IFCPOLYGONALFACESET,
        webifc::schema::IFCMAPPEDITEM,
        webifc::schema::IFCBOOLEANCLIPPINGRESULT,
        webifc::schema::IFCBOOLEANRESULT,
        webifc::schema::IFCHALFSPACESOLID,
        webifc::schema::IFCPOLYGONALBOUNDEDHALFSPACE,
        webifc::schema::IFCFACEBASEDSURFACEMODEL,
        webifc::schema::IFCSHELLBASEDSURFACEMODEL};

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
        std::cout << schema.IfcTypeCodeToType(ifcType.first) << "\t" << ifcType.second << std::endl;
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
         std::cout << schema.IfcTypeCodeToType(ifcType.first) << "\t" << ifcType.second << std::endl;
       // std::cout << ifcType.second << "\t"
          //        << schema.IfcTypeCodeToType(ifcType.first) << std::endl;
      }
    }

    std::cout << "\nUnsupported Types: " << unsupportedSize << std::endl;
    if (webifc::statistics::verboseStats) {
      std::cout << "Ifc Type\tFrequency\n" << std::endl;
      for (const auto &ifcType : unsupportedTypes) {
         std::cout << schema.IfcTypeCodeToType(ifcType.first) << "\t" << ifcType.second << std::endl;
       // std::cout << ifcType.second << "\t"
         //         << schema.IfcTypeCodeToType(ifcType.first) << std::endl;
      }
    }

    std::cout << "\nPercentage Supported: "
              << ((double)supportedSize / (double)uniqueTypeDefsSize) * 100
              << "%" << std::endl;
    std::cout << "************************************" << std::endl;
  }

  std::cout << "Done." << std::endl;
}
