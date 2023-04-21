/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <emscripten/bind.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stack>
#include <string>

#include "conway_geometry/ConwayGeometryProcessor.h"

std::map<uint32_t, std::unique_ptr<conway::geometry::ConwayGeometryProcessor>> processors;

uint32_t GLOBAL_MODEL_ID_COUNTER = 0;

#ifdef __EMSCRIPTEN_PTHREADS__
constexpr bool MT_ENABLED = true;
#else
constexpr bool MT_ENABLED = false;
#endif

// use to construct API placeholders
int main() {
  processors.emplace();

  return 0;
}

// taken from web ifc obj dump code
glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                        glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));

conway::geometry::ConwayGeometryProcessor::ResultsGltf GeometryToGltf(
    uint32_t modelID, conway::geometry::IfcGeometry geom, bool isGlb, bool outputDraco,
    std::string filePath) {
  auto& conwayProcessor = processors[modelID];

  conway::geometry::ConwayGeometryProcessor::ResultsGltf results =
      conwayProcessor->GeometryToGltf(geom, isGlb, outputDraco, filePath, false,
                                      NormalizeMat);

  return results;
}

std::string GeometryToObj(uint32_t modelID, conway::geometry::IfcGeometry geom,
                          size_t offset) {
  auto& conwayProcessor = processors[modelID];
  return conwayProcessor->GeometryToObj(geom, offset, NormalizeMat);
}

conway::geometry::IfcGeometry GetGeometry(
    uint32_t modelID,
    conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet parameters) {
  auto& conwayProcessor = processors[modelID];
  return conwayProcessor->getPolygonalFaceSetGeometry(parameters);
}

uint32_t InitializeGeometryProcessor() {
  uint32_t modelID = GLOBAL_MODEL_ID_COUNTER++;
  auto conwayProcessor = std::make_unique<conway::geometry::ConwayGeometryProcessor>();
  processors.emplace(modelID, std::move(conwayProcessor));

  return modelID;
}

bool FreeGeometryProcessor(uint32_t modelID) {
  processors.erase(modelID);
  return true;
}

emscripten::val GetUint8Array(std::vector<uint8_t>& buffer) {
  return emscripten::val(
      emscripten::typed_memory_view(buffer.size(), buffer.data()));
}

typedef glm::vec3 glmVec3;
typedef std::vector<glm::vec3> glmVec3Array;

EMSCRIPTEN_BINDINGS(my_module) {
  emscripten::class_<conway::geometry::IfcGeometry>("IfcGeometry")
      .constructor<>()
      .function("GetVertexData", &conway::geometry::IfcGeometry::GetVertexData)
      .function("GetVertexDataSize", &conway::geometry::IfcGeometry::GetVertexDataSize)
      .function("GetIndexData", &conway::geometry::IfcGeometry::GetIndexData)
      .function("GetIndexDataSize", &conway::geometry::IfcGeometry::GetIndexDataSize)
      .function("AddGeometry", &conway::geometry::IfcGeometry::AddGeometry);

  emscripten::value_object<glm::dvec4>("dvec4")
      .field("x", &glm::dvec4::x)
      .field("y", &glm::dvec4::y)
      .field("z", &glm::dvec4::z)
      .field("w", &glm::dvec4::w);

  emscripten::value_object<glm::vec3>("glmVec3")
      .field("x", &glm::vec3::x)
      .field("y", &glm::vec3::y)
      .field("z", &glm::vec3::z);

  emscripten::register_vector<glm::vec3>("glmVec3Array");

  // conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet>(
      "ParamsPolygonalFaceSet")
      .field(
          "numPoints",
          &conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet::numPoints)
      .field(
          "numIndices",
          &conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet::numIndices)
      .field("indicesPerFace", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsPolygonalFaceSet::indicesPerFace)
      .field("indexedPolygonalFaceWithVoids",
             &conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet::
                 indexedPolygonalFaceWithVoids)
      .field("points",
             &conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet::points)
      .field("indices",
             &conway::geometry::ConwayGeometryProcessor::ParamsPolygonalFaceSet::indices);

  // Define the ResultsGltf object
  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::ResultsGltf>(
      "ResultsGltf")
      .field("success", &conway::geometry::ConwayGeometryProcessor::ResultsGltf::success)
      .field("bufferUris",
             &conway::geometry::ConwayGeometryProcessor::ResultsGltf::bufferUris)
      .field("buffers", &conway::geometry::ConwayGeometryProcessor::ResultsGltf::buffers);

  emscripten::value_array<std::array<double, 16>>("array_double_16")
      .element(emscripten::index<0>())
      .element(emscripten::index<1>())
      .element(emscripten::index<2>())
      .element(emscripten::index<3>())
      .element(emscripten::index<4>())
      .element(emscripten::index<5>())
      .element(emscripten::index<6>())
      .element(emscripten::index<7>())
      .element(emscripten::index<8>())
      .element(emscripten::index<9>())
      .element(emscripten::index<10>())
      .element(emscripten::index<11>())
      .element(emscripten::index<12>())
      .element(emscripten::index<13>())
      .element(emscripten::index<14>())
      .element(emscripten::index<15>());

  emscripten::register_vector<std::string>("stringVector");
  emscripten::register_vector<uint32_t>("UintVector");
  emscripten::register_vector<uint8_t>("VectorUint8");
  emscripten::register_vector<std::vector<uint8_t>>("VectorVectorUint8");
  emscripten::function("GetGeometry", &GetGeometry);
  emscripten::function("InitializeGeometryProcessor",
                       &InitializeGeometryProcessor);
  emscripten::function("FreeGeometryProcessor", &FreeGeometryProcessor);
  emscripten::function("GeometryToObj", &GeometryToObj);
  emscripten::function("GeometryToGltf", &GeometryToGltf);
  emscripten::function("GetUint8Array", &GetUint8Array,
                       emscripten::allow_raw_pointers());
}
