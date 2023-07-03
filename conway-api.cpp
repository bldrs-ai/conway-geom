#include <emscripten/bind.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stack>
#include <string>

#include "conway_geometry/ConwayGeometryProcessor.h"

std::unique_ptr<conway::geometry::ConwayGeometryProcessor> processor;

#ifdef __EMSCRIPTEN_PTHREADS__
constexpr bool MT_ENABLED = true;
#else
constexpr bool MT_ENABLED = false;
#endif

// use to construct API placeholders
int main() { return 0; }

// taken from web ifc obj dump code
glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                        glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));

conway::geometry::ConwayGeometryProcessor::ResultsGltf GeometryToGltf(
    conway::geometry::IfcGeometry geom, bool isGlb, bool outputDraco,
    std::string filePath) {
  conway::geometry::ConwayGeometryProcessor::ResultsGltf results;
  if (processor) {
    results = processor->GeometryToGltf(geom, isGlb, outputDraco, filePath,
                                        false, NormalizeMat);
  }

  return results;
}

std::string GeometryToObj(conway::geometry::IfcGeometry geom, size_t offset) {
  std::string result;
  if (processor) {
    return processor->GeometryToObj(geom, offset, NormalizeMat);
  } else {
    return result;
  }
}

conway::geometry::IfcGeometry GetGeometry(
    conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
        parameters) {
  conway::geometry::IfcGeometry geom;

  if (processor) {
    return processor->getPolygonalFaceSetGeometry(parameters);
  } else {
    return geom;
  }
}

conway::geometry::IfcCurve GetIndexedPolyCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetIfcIndexedPolyCurve
        parameters) {
  conway::geometry::IfcCurve curve;
  if (processor) {
    return processor->getIndexedPolyCurve(parameters);
  } else {
    return curve;
  }
}

glm::dmat4 GetLocalPlacement(
    conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement
        parameters) {
  glm::dmat4 resultMat;
  if (processor) {
    resultMat = processor->GetLocalPlacement(parameters);
  }

  return resultMat;
}

glm::dmat3 GetAxis2Placement2D(
    conway::geometry::ConwayGeometryProcessor::ParamsGetAxis2Placement2D
        parameters) {
  glm::dmat3 resultMat;
  if (processor) {
    resultMat = processor->GetAxis2Placement2D(parameters);
  }

  return resultMat;
}

glm::dmat4 GetAxis2Placement3D(
    conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
        parameters) {
  glm::dmat4 resultMat;
  if (processor) {
    resultMat = processor->GetAxis2Placement3D(parameters);
  }

  return resultMat;
}

bool InitializeGeometryProcessor() {
  processor = std::make_unique<conway::geometry::ConwayGeometryProcessor>();

  return true;
}

bool FreeGeometryProcessor() {
  processor.release();
  return true;
}

emscripten::val GetUint8Array(std::vector<uint8_t>& buffer) {
  return emscripten::val(
      emscripten::typed_memory_view(buffer.size(), buffer.data()));
}

// Helper function to convert glm::dmat4x4 to a linear array.
emscripten::val getMatrixValues(const glm::dmat4& mat) {
  emscripten::val array = emscripten::val::array();
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      array.set(i * 4 + j, mat[i][j]);
    }
  }
  return array;
}

// Helper function to set values of a glm::dmat4x4 from a linear array.
void setMatrixValues(glm::dmat4& mat, emscripten::val array) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat[i][j] = array[i * 4 + j].as<double>();
    }
  }
}

typedef glm::vec3 glmVec3;
typedef std::vector<glm::vec3> glmVec3Array;

EMSCRIPTEN_BINDINGS(my_module) {
  emscripten::class_<conway::geometry::IfcGeometry>("IfcGeometry")
      .constructor<>()
      .function("getVertexData", &conway::geometry::IfcGeometry::GetVertexData)
      .function("getVertexDataSize",
                &conway::geometry::IfcGeometry::GetVertexDataSize)
      .function("getIndexData", &conway::geometry::IfcGeometry::GetIndexData)
      .function("getIndexDataSize",
                &conway::geometry::IfcGeometry::GetIndexDataSize)
      .function("appendGeometry",
                &conway::geometry::IfcGeometry::AppendGeometry)
      .function("applyTransform",
                &conway::geometry::IfcGeometry::ApplyTransform)
      .function("clone", &conway::geometry::IfcGeometry::Clone);

  emscripten::class_<conway::geometry::IfcCurve>("IfcCurve")
      .constructor<>()
      .function("add2d", &conway::geometry::IfcCurve::Add2d)
      .function("add3d", &conway::geometry::IfcCurve::Add3d)
      .function("get2d", &conway::geometry::IfcCurve::Get2d)
      .function("get3d", &conway::geometry::IfcCurve::Get3d)
      .function("invert", &conway::geometry::IfcCurve::Invert)
      .function("isCCW", &conway::geometry::IfcCurve::IsCCW);

  /**
   *  void Add(glm::dvec3 pt);
void Add(glm::dvec2 pt);
glm::dvec2 Get2d(size_t i) const;
glm::dvec3 Get3d(size_t i) const;
void Invert();
bool IsCCW() const;

  */

  emscripten::value_object<glm::dvec4>("dvec4")
      .field("x", &glm::dvec4::x)
      .field("y", &glm::dvec4::y)
      .field("z", &glm::dvec4::z)
      .field("w", &glm::dvec4::w);

  emscripten::value_object<glm::vec3>("glmVec3")
      .field("x", &glm::vec3::x)
      .field("y", &glm::vec3::y)
      .field("z", &glm::vec3::z);

  emscripten::value_object<glm::dvec3>("glmdVec3")
      .field("x", &glm::dvec3::x)
      .field("y", &glm::dvec3::y)
      .field("z", &glm::dvec3::z);

  emscripten::value_object<glm::vec2>("vec2")
      .field("x", &glm::vec2::x)
      .field("y", &glm::vec2::y);

  emscripten::register_vector<glm::vec2>("vec2Array");

  emscripten::class_<glm::dmat4>("glmdmat4")
      .constructor<>()
      .function("getValues", &getMatrixValues)
      .function("setValues", &setMatrixValues);

  emscripten::register_vector<glm::vec3>("glmVec3Array");

  // conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>(
      "IndexedPolygonalFace")
      .field("indices", &conway::geometry::ConwayGeometryProcessor::
                            IndexedPolygonalFace::indices)
      .field("face_starts", &conway::geometry::ConwayGeometryProcessor::
                                IndexedPolygonalFace::face_starts);

  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::Segment>(
      "Segment")
      .field("isArcType",
             &conway::geometry::ConwayGeometryProcessor::Segment::isArcType)
      .field("indices",
             &conway::geometry::ConwayGeometryProcessor::Segment::indices);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetIfcIndexedPolyCurve
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetIfcIndexedPolyCurve>(
      "ParamsGetIfcIndexedPolyCurve")
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetIfcIndexedPolyCurve::dimensions)
      .field("segments", &conway::geometry::ConwayGeometryProcessor::
                             ParamsGetIfcIndexedPolyCurve::segments)
      .field("points", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetIfcIndexedPolyCurve::points);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::
                               ParamsGetPolygonalFaceSetGeometry>(
      "ParamsGetPolygonalFaceSetGeometry")
      .field("indicesPerFace",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetPolygonalFaceSetGeometry::indicesPerFace)
      .field("points", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetPolygonalFaceSetGeometry::points)
      .field("faces", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetPolygonalFaceSetGeometry::faces);

  // conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement2D
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetAxis2Placement2D>(
      "ParamsAxis2Placement2D")
      .field("isAxis2Placement2D",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetAxis2Placement2D::isAxis2Placement2D)
      .field("isCartesianTransformationOperator2D",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetAxis2Placement2D::isCartesianTransformationOperator2D)
      .field("isCartesianTransformationOperator2DNonUniform",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetAxis2Placement2D::
                     isCartesianTransformationOperator2DNonUniform)
      .field("position2D", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetAxis2Placement2D::position2D)
      .field("customAxis1Ref", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetAxis2Placement2D::customAxis1Ref)
      .field("axis1Ref", &conway::geometry::ConwayGeometryProcessor::
                             ParamsGetAxis2Placement2D::axis1Ref)
      .field("customAxis2Ref", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetAxis2Placement2D::customAxis2Ref)
      .field("axis2Ref", &conway::geometry::ConwayGeometryProcessor::
                             ParamsGetAxis2Placement2D::axis2Ref)
      .field("customScale", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetAxis2Placement2D::customScale)
      .field("scale1", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetAxis2Placement2D::scale1)
      .field("customScale2", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetAxis2Placement2D::customScale2)
      .field("scale2", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetAxis2Placement2D::scale2);

  // conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D>(
      "ParamsAxis2Placement3D")
      .field("position", &conway::geometry::ConwayGeometryProcessor::
                             ParamsAxis2Placement3D::position)
      .field("zAxisRef", &conway::geometry::ConwayGeometryProcessor::
                             ParamsAxis2Placement3D::zAxisRef)
      .field("xAxisRef", &conway::geometry::ConwayGeometryProcessor::
                             ParamsAxis2Placement3D::xAxisRef)
      .field("normalizeZ", &conway::geometry::ConwayGeometryProcessor::
                               ParamsAxis2Placement3D::normalizeZ)
      .field("normalizeX", &conway::geometry::ConwayGeometryProcessor::
                               ParamsAxis2Placement3D::normalizeX);

  // conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement>(
      "ParamsLocalPlacement")
      .field("useRelPlacement", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsLocalPlacement::useRelPlacement)
      .field("axis2Placement", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsLocalPlacement::axis2Placement)
      .field("relPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsLocalPlacement::relPlacement);

  // Define the ResultsGltf object
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ResultsGltf>("ResultsGltf")
      .field("success",
             &conway::geometry::ConwayGeometryProcessor::ResultsGltf::success)
      .field(
          "bufferUris",
          &conway::geometry::ConwayGeometryProcessor::ResultsGltf::bufferUris)
      .field("buffers",
             &conway::geometry::ConwayGeometryProcessor::ResultsGltf::buffers);

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
  emscripten::register_vector<size_t>("ULongVector");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>(
      "VectorIndexedPolygonalFace");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::Segment>("VectorSegment");
  emscripten::function("getGeometry", &GetGeometry);
  emscripten::function("getIndexedPolyCurve", &GetIndexedPolyCurve);
  emscripten::function("initializeGeometryProcessor",
                       &InitializeGeometryProcessor);
  emscripten::function("freeGeometryProcessor", &FreeGeometryProcessor);
  emscripten::function("geometryToObj", &GeometryToObj);
  emscripten::function("geometryToGltf", &GeometryToGltf);
  emscripten::function("getAxis2Placement2D", &GetAxis2Placement2D);
  emscripten::function("getAxis2Placement3D", &GetAxis2Placement3D);
  emscripten::function("getLocalPlacement", &GetLocalPlacement);
  emscripten::function("getUint8Array", &GetUint8Array,
                       emscripten::allow_raw_pointers());
}
