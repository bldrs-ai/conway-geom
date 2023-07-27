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
    std::vector<conway::geometry::IfcGeometry>& geoms,
    std::vector<conway::geometry::Material>& materials, bool isGlb,
    bool outputDraco, std::string filePath) {
  conway::geometry::ConwayGeometryProcessor::ResultsGltf results;
  if (processor) {
    results = processor->GeometryToGltf(geoms, materials, isGlb, outputDraco,
                                        filePath, false, NormalizeMat);
  }

  return results;
}

std::string GeometryToObj(conway::geometry::IfcGeometry geom, size_t offset) {
  if (processor) {
    return processor->GeometryToObj(geom, offset, NormalizeMat);
  }

  std::string result;
  return result;
}

conway::geometry::IfcGeometry GetPolygonalFaceSetGeometry(
    conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
        parameters) {
  if (processor) {
    return processor->getPolygonalFaceSetGeometry(parameters);
  }

  conway::geometry::IfcGeometry geom;
  return geom;
}

conway::geometry::IfcCurve GetIndexedPolyCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetIfcIndexedPolyCurve
        parameters) {
  if (processor) {
    return processor->getIndexedPolyCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetCircleCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetCircleCurve
        parameters) {
  if (processor) {
    return processor->getCircleCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetTrimmedCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetIfcTrimmedCurve
        parameters) {
  if (processor) {
    return processor->getTrimmedCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcGeometry GetExtrudedAreaSolid(
    conway::geometry::ConwayGeometryProcessor::ParamsGetExtrudedAreaSolid
        parameters) {
  if (processor) {
    return processor->getExtrudedAreaSolid(parameters);
  }

  conway::geometry::IfcGeometry geom;
  return geom;
}

glm::dmat4 GetLocalPlacement(
    conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement
        parameters) {
  if (processor) {
    return processor->GetLocalPlacement(parameters);
  }

  glm::dmat4 resultMat;

  return resultMat;
}

glm::dmat3 GetAxis2Placement2D(
    conway::geometry::ConwayGeometryProcessor::ParamsGetAxis2Placement2D
        parameters) {
  if (processor) {
    return processor->GetAxis2Placement2D(parameters);
  }

  glm::dmat3 resultMat;
  return resultMat;
}

glm::dmat4 GetAxis2Placement3D(
    conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
        parameters) {
  if (processor) {
    return processor->GetAxis2Placement3D(parameters);
  }

  glm::dmat4 resultMat;
  return resultMat;
}

conway::geometry::IfcGeometry GetBooleanResult(
    conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult
        parameters) {
  if (processor) {
    return processor->GetBooleanResult(parameters);
  }
  conway::geometry::IfcGeometry geometry;
  return geometry;
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
emscripten::val getMatrixValues4x4(const glm::dmat4& mat) {
  emscripten::val array = emscripten::val::array();
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      array.set(i * 4 + j, mat[i][j]);
    }
  }
  return array;
}

// Helper function to set values of a glm::dmat4x4 from a linear array.
void setMatrixValues4x4(glm::dmat4& mat, emscripten::val array) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mat[i][j] = array[i * 4 + j].as<double>();
    }
  }
}

// Helper function to convert glm::dmat3x3 to a linear array.
emscripten::val getMatrixValues3x3(const glm::dmat3& mat) {
  emscripten::val array = emscripten::val::array();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      array.set(i * 3 + j, mat[i][j]);
    }
  }
  return array;
}

// Helper function to set values of a glm::dmat3x3 from a linear array.
void setMatrixValues3x3(glm::dmat3& mat, emscripten::val array) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      mat[i][j] = array[i * 3 + j].as<double>();
    }
  }
}

struct ParamsCreateNativeIfcProfile {
  conway::geometry::IfcCurve curve;
  std::vector<conway::geometry::IfcCurve> holes;
  bool isConvex;
  bool isComposite;
  std::vector<conway::geometry::IfcProfile> profiles;
};
conway::geometry::IfcProfile createNativeIfcProfile(
    ParamsCreateNativeIfcProfile parameters) {
  conway::geometry::IfcProfile profile;

  profile.type = "testType";
  profile.curve = parameters.curve;
  profile.holes = parameters.holes;
  profile.isConvex = parameters.isConvex;
  profile.isComposite = parameters.isComposite;
  profile.profiles = parameters.profiles;

  return profile;
}

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
      .function("clone", &conway::geometry::IfcGeometry::Clone)
      .property("materialIndex", &conway::geometry::IfcGeometry::materialIndex)
      .property("hasDefaultMaterial",
                &conway::geometry::IfcGeometry::hasDefaultMaterial);

  emscripten::class_<conway::geometry::IfcCurve>("IfcCurve")
      .constructor<>()
      .function("add2d", &conway::geometry::IfcCurve::Add2d)
      .function("add3d", &conway::geometry::IfcCurve::Add3d)
      .function("get2d", &conway::geometry::IfcCurve::Get2d)
      .function("get3d", &conway::geometry::IfcCurve::Get3d)
      .function("invert", &conway::geometry::IfcCurve::Invert)
      .function("isCCW", &conway::geometry::IfcCurve::IsCCW);

  emscripten::class_<conway::geometry::IfcProfile>("IfcProfile")
      .constructor<>()
      .function("getType", &conway::geometry::IfcProfile::getType)
      .function("getCurve", &conway::geometry::IfcProfile::getCurve)
      .function("getHoles", &conway::geometry::IfcProfile::getHoles)
      .function("isConvex", &conway::geometry::IfcProfile::getIsConvex)
      .function("isComposite", &conway::geometry::IfcProfile::getIsComposite)
      .function("getProfiles", &conway::geometry::IfcProfile::getProfiles);

  emscripten::class_<glm::dmat4>("glmdmat4")
      .constructor<>()
      .function("getValues", &getMatrixValues4x4)
      .function("setValues", &setMatrixValues4x4);

  emscripten::class_<glm::dmat3>("glmdmat3")
      .constructor<>()
      .function("getValues", &getMatrixValues3x3)
      .function("setValues", &setMatrixValues3x3);

  emscripten::enum_<conway::geometry::BLEND_MODE>("BlendMode")
      .value("OPAQUE", conway::geometry::BLEND_MODE::BLEND_OPAQUE)
      .value("BLEND", conway::geometry::BLEND_MODE::BLEND)
      .value("MASK", conway::geometry::BLEND_MODE::MASK);

  emscripten::value_object<conway::geometry::Material>("MaterialObject")
      .field("base", &conway::geometry::Material::base)
      .field("metallic", &conway::geometry::Material::metallic)
      .field("roughness", &conway::geometry::Material::roughness)
      .field("alphaCutoff", &conway::geometry::Material::alphaCutoff)
      .field("ior", &conway::geometry::Material::ior)
      .field("specular", &conway::geometry::Material::specular)
      .field("alphaMode", &conway::geometry::Material::alphaMode)
      .field("doubleSided", &conway::geometry::Material::doubleSided);

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

  emscripten::value_object<glm::dvec2>("glmdVec2")
      .field("x", &glm::dvec2::x)
      .field("y", &glm::dvec2::y);

  emscripten::value_object<glm::vec2>("glmVec2")
      .field("x", &glm::vec2::x)
      .field("y", &glm::vec2::y);

  emscripten::register_vector<conway::geometry::Material>("materialArray");

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

  // conway::geometry::ConwayGeometryProcessor::ParamsGetCircleCurve {
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetCircleCurve>(
      "ParamsGetCircleCurve")
      .field("radius", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetCircleCurve::radius)
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetCircleCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetCircleCurve::placement);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetExtrudedAreaSolid
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetExtrudedAreaSolid>(
      "ParamsGetExtrudedAreaSolid")
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetExtrudedAreaSolid::depth)
      .field("dir", &conway::geometry::ConwayGeometryProcessor::
                        ParamsGetExtrudedAreaSolid::dir)
      .field("profile", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetExtrudedAreaSolid::profile);

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
  /*    glm::dvec3 position;
      glm::dvec3 axis1Ref;
      glm::dvec3 axis2Ref;
      glm::dvec3 axis3Ref;
      bool normalizeAxis1 = false;
      bool normalizeAxis2 = false;
      bool normalizeAxis3 = false;
      bool nonUniform = false;
      bool realScale = false;
      double scale1_ = 0;
      double scale2_ = 0;
      double scale3_ = 0;*/
  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::
                               ParamsCartesianTransformationOperator3D>(
      "ParamsCartesianTransformationOperator3D")
      .field("position", &conway::geometry::ConwayGeometryProcessor::
                             ParamsCartesianTransformationOperator3D::position)
      .field("axis1Ref", &conway::geometry::ConwayGeometryProcessor::
                             ParamsCartesianTransformationOperator3D::axis1Ref)
      .field("axis2Ref", &conway::geometry::ConwayGeometryProcessor::
                             ParamsCartesianTransformationOperator3D::axis2Ref)
      .field("axis3Ref", &conway::geometry::ConwayGeometryProcessor::
                             ParamsCartesianTransformationOperator3D::axis3Ref)
      .field("normalizeAxis1",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsCartesianTransformationOperator3D::normalizeAxis1)
      .field("normalizeAxis2",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsCartesianTransformationOperator3D::normalizeAxis2)
      .field("normalizeAxis3",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsCartesianTransformationOperator3D::normalizeAxis3)
      .field("nonUniform",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsCartesianTransformationOperator3D::nonUniform)
      .field("realScale",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsCartesianTransformationOperator3D::realScale)
      .field("scale1_", &conway::geometry::ConwayGeometryProcessor::
                            ParamsCartesianTransformationOperator3D::scale1_)
      .field("scale2_", &conway::geometry::ConwayGeometryProcessor::
                            ParamsCartesianTransformationOperator3D::scale2_)
      .field("scale3_", &conway::geometry::ConwayGeometryProcessor::
                            ParamsCartesianTransformationOperator3D::scale3_);

  // ParamsCreateNativeIfcProfile
  emscripten::value_object<ParamsCreateNativeIfcProfile>(
      "ParamsCreateNativeIfcProfile")
      .field("curve", &ParamsCreateNativeIfcProfile::curve)
      .field("holes", &ParamsCreateNativeIfcProfile::holes)
      .field("isConvex", &ParamsCreateNativeIfcProfile::isConvex)
      .field("isComposite", &ParamsCreateNativeIfcProfile::isComposite)
      .field("profiles", &ParamsCreateNativeIfcProfile::profiles);

  // ParamsGetIfcTrimmedCurve
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetIfcTrimmedCurve>(
      "ParamsCreateNativeIfcProfile")
      .field("basisCurve", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetIfcTrimmedCurve::basisCurve)
      .field("masterRepresentation",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetIfcTrimmedCurve::masterRepresentation)
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetIfcTrimmedCurve::dimensions)
      .field("senseAgreement", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetIfcTrimmedCurve::senseAgreement)
      .field("trim1Vec3", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetIfcTrimmedCurve::trim1Vec3)
      .field("trim1VecDouble", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetIfcTrimmedCurve::trim1VecDouble)
      .field("trim2Vec3", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetIfcTrimmedCurve::trim2Vec3)
      .field("trim2VecDouble", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetIfcTrimmedCurve::trim2VecDouble);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult>(
      "ParamsGetBooleanResult")
      .field("flatFirstMesh", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsGetBooleanResult::flatFirstMesh)
      .field("flatSecondMesh", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetBooleanResult::flatSecondMesh)
      .field("operatorType", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetBooleanResult::operatorType);

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

  emscripten::register_vector<glm::vec2>("vec2Array");
  emscripten::register_vector<glm::vec3>("glmVec3Array");
  emscripten::register_vector<std::string>("stringVector");
  emscripten::register_vector<uint32_t>("UintVector");
  emscripten::register_vector<uint8_t>("VectorUint8");
  emscripten::register_vector<std::vector<uint8_t>>("VectorVectorUint8");
  emscripten::register_vector<size_t>("ULongVector");
  emscripten::register_vector<conway::geometry::IfcCurve>("curveArray");
  emscripten::register_vector<conway::geometry::IfcProfile>("profileArray");
  emscripten::register_vector<conway::geometry::IfcGeometry>("geometryArray");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>(
      "VectorIndexedPolygonalFace");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::Segment>("VectorSegment");

  emscripten::function("getPolygonalFaceSetGeometry",
                       &GetPolygonalFaceSetGeometry);
  emscripten::function("getIndexedPolyCurve", &GetIndexedPolyCurve);
  emscripten::function("getCircleCurve", &GetCircleCurve);
  emscripten::function("initializeGeometryProcessor",
                       &InitializeGeometryProcessor);
  emscripten::function("freeGeometryProcessor", &FreeGeometryProcessor);
  emscripten::function("geometryToObj", &GeometryToObj);
  emscripten::function("geometryToGltf", &GeometryToGltf);
  emscripten::function("getAxis2Placement2D", &GetAxis2Placement2D);
  emscripten::function("getAxis2Placement3D", &GetAxis2Placement3D);
  emscripten::function("getLocalPlacement", &GetLocalPlacement);
  emscripten::function("getCartesianTransformationOperator3D",
                       &conway::geometry::ConwayGeometryProcessor::
                           GetCartesianTransformationOperator3D);
  emscripten::function("getUint8Array", &GetUint8Array,
                       emscripten::allow_raw_pointers());
  emscripten::function("createNativeIfcProfile", &createNativeIfcProfile);
  emscripten::function("getExtrudedAreaSolid", &GetExtrudedAreaSolid);
  emscripten::function("getBooleanResult", &GetBooleanResult);
}
