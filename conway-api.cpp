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
    std::vector< conway::geometry::IfcGeometry >& geoms,
    std::vector< conway::geometry::Material >& materials,
    bool isGlb,
    bool outputDraco, std::string filePath) {
  conway::geometry::ConwayGeometryProcessor::ResultsGltf results;
  if (processor) {
    results = processor->GeometryToGltf(geoms, materials, isGlb, outputDraco, filePath,
                                        false, NormalizeMat);
  }

  return results;
}

std::string GeometryToObj(conway::geometry::IfcGeometry geom,
                          size_t offset) {
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

glm::dmat4 GetLocalPlacement(
    conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement
        parameters) {
  glm::dmat4 resultMat;
  if (processor) {
    resultMat = processor->GetLocalPlacement(parameters);
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
                .function("clone", 
                &conway::geometry::IfcGeometry::Clone)
      .property("materialIndex", &conway::geometry::IfcGeometry::materialIndex )
      .property("hasDefaultMaterial", &conway::geometry::IfcGeometry::hasDefaultMaterial );

/*    // Note - base alpha is used for transparency 
    glm::dvec4 base        = glm::dvec4( 0.8, 0.8, 0.8, 1 );
    double     metallic    = 0.0f;
    double     roughness   = 1.0f;
    double     alphaCutoff = 0;
    double     ior         = 1.4;

    std::optional< glm::dvec4 > specular;

    BLEND_MODE alphaMode = BLEND_MODE::BLEND_OPAQUE;

    bool doubleSized = false;*/

    emscripten::enum_< conway::geometry::BLEND_MODE >( "BlendMode" )
      .value( "OPAQUE", conway::geometry::BLEND_MODE::BLEND_OPAQUE )
      .value( "BLEND", conway::geometry::BLEND_MODE::BLEND )
      .value( "MASK", conway::geometry::BLEND_MODE::MASK );

    emscripten::value_object<conway::geometry::Material>("MaterialObject")
      .field("base", &conway::geometry::Material::base )
      .field("metallic", &conway::geometry::Material::metallic )
      .field("roughness", &conway::geometry::Material::roughness )
      .field("alphaCutoff", &conway::geometry::Material::alphaCutoff )
      .field("ior", &conway::geometry::Material::ior )
      .field("specular", &conway::geometry::Material::specular )
      .field("alphaMode", &conway::geometry::Material::alphaMode )
      .field("doubleSided", &conway::geometry::Material::doubleSided );

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

  emscripten::class_<glm::dmat4>("glmdmat4")
      .constructor<>()
      .function("getValues", &getMatrixValues)
      .function("setValues", &setMatrixValues);

  emscripten::register_vector<glm::vec3>("glmVec3Array");
  
  emscripten::register_vector<conway::geometry::Material>("MaterialVector");

  // conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>(
      "IndexedPolygonalFace")
      .field("indices", &conway::geometry::ConwayGeometryProcessor::
                            IndexedPolygonalFace::indices)
      .field("face_starts", &conway::geometry::ConwayGeometryProcessor::
                                IndexedPolygonalFace::face_starts);

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
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D >( "ParamsCartesianTransformationOperator3D" )
      .field("position", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::position )
      .field("axis1Ref", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::axis1Ref )
      .field("axis2Ref", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::axis2Ref )
      .field("axis3Ref", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::axis3Ref )
      .field("normalizeAxis1", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::normalizeAxis1 )
      .field("normalizeAxis2", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::normalizeAxis2 )
      .field("normalizeAxis3", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::normalizeAxis3 )
      .field("nonUniform", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::nonUniform )
      .field("realScale", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::realScale )
      .field("scale1_", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::scale1_ )
      .field("scale2_", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::scale2_ )
      .field("scale3_", &conway::geometry::ConwayGeometryProcessor::ParamsCartesianTransformationOperator3D::scale3_ );

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
  emscripten::register_vector<conway::geometry::IfcGeometry>("GeometryVector");
  emscripten::function("getGeometry", &GetGeometry);
  emscripten::function("initializeGeometryProcessor",
    &InitializeGeometryProcessor);
  emscripten::function("freeGeometryProcessor", &FreeGeometryProcessor);
  emscripten::function("geometryToObj", &GeometryToObj);
  emscripten::function("geometryToGltf", &GeometryToGltf);
  emscripten::function("getAxis2Placement3D", &GetAxis2Placement3D);
  emscripten::function("getLocalPlacement", &GetLocalPlacement);
  emscripten::function("getCartesianTransformationOperator3D", 
    &conway::geometry::ConwayGeometryProcessor::GetCartesianTransformationOperator3D );
  emscripten::function("getUint8Array", &GetUint8Array,
    emscripten::allow_raw_pointers());
}
