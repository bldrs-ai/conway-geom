#include <emscripten/bind.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stack>
#include <string>

#include "utilities/buffer_parse.h"
#include "structures/parse_buffer.h"
#include "conway_geometry/ConwayGeometryProcessor.h"
#include "logging/Logger.h"


std::unique_ptr<conway::geometry::ConwayGeometryProcessor> processor;

// use to construct API placeholders
int main() { return 0; }

// taken from web ifc obj dump code
glm::dmat4 NormalizeMat(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 0, -1, 0),
                        glm::dvec4(0, 1, 0, 0), glm::dvec4(0, 0, 0, 1));

conway::geometry::ConwayGeometryProcessor::ResultsGltf GeometryToGltf(
    std::vector<conway::geometry::GeometryCollection>& geoms,
    std::vector<conway::geometry::Material>& materials, bool isGlb,
    bool outputDraco, std::string filePath, size_t offset, size_t count) {
  conway::geometry::ConwayGeometryProcessor::ResultsGltf results;

  if ( processor ) {

    offset = std::min( offset, geoms.size() );
    count  = std::min( count, geoms.size() - offset );

    results = processor->GeometryToGltf(
        std::span{geoms.data() + offset, count},
        std::span{materials.data(), materials.size()}, isGlb, outputDraco,
        filePath, false, NormalizeMat);
  }

  return results;
}

std::string GeometryToObj(conway::geometry::Geometry geom, size_t offset) {
  if (processor) {
    return processor->GeometryToObj(geom, offset, NormalizeMat);
  }

  std::string result;
  return result;
}

glm::dmat4 multiplyNativeMatrices(glm::dmat4 mat1, glm::dmat4 mat2) {
  return mat1 * mat2;
}

conway::geometry::Geometry GetSweptDiskSolid(
  conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid& parameters) {
    if (processor) {
      return processor->getSweptDiskSolid(parameters);
    }

    conway::geometry::Geometry geom;
    return geom;
  }

conway::geometry::Geometry GetPolygonalFaceSetGeometry(
    conway::geometry::ConwayGeometryProcessor::
        ParamsGetPolygonalFaceSetGeometry& parameters) {
  if (processor) {
    return processor->getPolygonalFaceSetGeometry(parameters);
  }

  conway::geometry::Geometry geom;
  return geom;
}

conway::geometry::IfcCurve GetIndexedPolyCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetIfcIndexedPolyCurve&
        parameters) {
  if (processor) {
    return processor->getIndexedPolyCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetCircleCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetCircleCurve&
        parameters) {
  if (processor) {
    return processor->getCircleCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetCircleHoleCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetCircleCurve&
        parameters) {
  if (processor) {
    return processor->getCircleHoleCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetIfcCircle(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetIfcCircle& parameters) {
  if (processor) {
    return processor->getIfcCircle(parameters);
  }

  conway::geometry::IfcCurve curve;

  return curve;
}

conway::geometry::IfcCurve GetIfcLine(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetIfcLine& parameters) {
  if (processor) {
    return processor->getIfcLine(parameters);
  }

  conway::geometry::IfcCurve curve;

  return curve;
}

conway::geometry::IfcCurve GetBSplineCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetBSplineCurve&
        parameters) {
  if (processor) {
    return processor->getBSplineCurve(parameters);
  }

  conway::geometry::IfcCurve curve;

  return curve;
}

conway::geometry::Geometry GetPolygonalBoundedHalfspace(
  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace
  parameters) {
    if (processor) {
      return processor->GetPolygonalBoundedHalfspace(parameters);
    }

    return conway::geometry::Geometry();
  }

conway::geometry::Geometry GetExtrudedAreaSolid(
    conway::geometry::ConwayGeometryProcessor::ParamsGetExtrudedAreaSolid
        parameters) {
  if (processor) {
    return processor->getExtrudedAreaSolid(parameters);
  }

  conway::geometry::Geometry geom;
  return geom;
}

conway::geometry::IfcCurve GetLoop(
    conway::geometry::ConwayGeometryProcessor::ParamsGetLoop parameters) {
  if (processor) {
    return processor->GetLoop(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

void AddFaceToGeometry(
    conway::geometry::ConwayGeometryProcessor::ParamsAddFaceToGeometry&
        parameters,
    conway::geometry::Geometry& geometry) {
  if (processor) {
    return processor->AddFaceToGeometry(parameters, geometry);
  }
}

struct ParamsCreateBound3D {
  conway::geometry::IfcCurve curve;
  bool orientation;
  uint32_t type;
};

conway::geometry::IfcBound3D createBound3D(const ParamsCreateBound3D& parameters) {
  conway::geometry::IfcBound3D bounds3D;
  bounds3D.curve = parameters.curve;
  bounds3D.orientation = parameters.orientation;

  if (!bounds3D.orientation) {
    std::reverse(bounds3D.curve.points.begin(), bounds3D.curve.points.end());
  }

  switch (parameters.type) {
    case 0:
      bounds3D.type = conway::geometry::IfcBoundType::OUTERBOUND;
      break;
    case 1:
      bounds3D.type = conway::geometry::IfcBoundType::BOUND;
      break;
    default:
      Logger::logWarning("Invalid value for IfcBoundType enum!\n");
      // Handle the case when the uint32_t value doesn't correspond to any enum
      // value. You might want to provide a default or throw an exception here.
      break;
  }

  return bounds3D;
}

conway::geometry::Geometry GetHalfSpaceSolid(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetHalfspaceSolid&
        parameters) {
  if (processor) {
    return processor->GetHalfSpaceSolid(parameters);
  }

  conway::geometry::Geometry geom;
  return geom;
}

glm::dmat4 GetLocalPlacement(
    const conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement&
        parameters) {
  if (processor) {
    return processor->GetLocalPlacement(parameters);
  }

  glm::dmat4 resultMat;

  return resultMat;
}

glm::dmat3 GetAxis2Placement2D(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetAxis2Placement2D&
        parameters) {
  if (processor) {
    return processor->GetAxis2Placement2D(parameters);
  }

  glm::dmat3 resultMat;
  return resultMat;
}

glm::dmat4 GetAxis1Placement(
    const conway::geometry::ConwayGeometryProcessor::ParamsAxis1Placement3D&
        parameters) {
  if (processor) {
    return processor->GetAxis1Placement(parameters);
  }

  glm::dmat4 resultMat;
  return resultMat;
}

glm::dmat4 GetAxis2Placement3D(
    const conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D&
        parameters) {
  if (processor) {
    return processor->GetAxis2Placement3D(parameters);
  }

  glm::dmat4 resultMat;
  return resultMat;
}

conway::geometry::Geometry GetBooleanResult(
    conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult*
        parameters) {
  if (processor) {
    return processor->GetBooleanResult(parameters);
  }
  conway::geometry::Geometry geometry;
  return geometry;
}

conway::geometry::Geometry RelVoidSubtract(
    conway::geometry::ConwayGeometryProcessor::ParamsRelVoidSubtract&
        parameters) {
  if (processor) {
    return processor->RelVoidSubtract(parameters);
  }
  conway::geometry::Geometry geometry;
  return geometry;
}

conway::geometry::IfcCurve GetCShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetCShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetCShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetIShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetIShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetIShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetLShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetLShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetLShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetTShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetTShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetTShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetUShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetUShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetUShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetZShapeCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetZShapeCurve&
        parameters) {
  if (processor) {
    return processor->GetZShapeCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetEllipseCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetEllipseCurve&
        parameters) {
  if (processor) {
    return processor->getEllipseCurve(parameters);
  }
  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetRectangleProfileCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetRectangleProfileCurve&
        parameters) {
  if (processor) {
    return processor->GetRectangleProfileCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

conway::geometry::IfcCurve GetRectangleHollowProfileHole(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetRectangleProfileCurve&
        parameters) {
  if (processor) {
    return processor->GetRectangleHollowProfileHole(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
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

conway::geometry::IfcCurve GetPolyCurve(
    const conway::geometry::ConwayGeometryProcessor::ParamsGetPolyCurve&
        parameters) {
  if (processor) {
    return processor->getPolyCurve(parameters);
  }

  conway::geometry::IfcCurve curve;
  return curve;
}

struct ParamsCreateNativeIfcProfile {
  conway::geometry::IfcCurve curve;
  std::vector<conway::geometry::IfcCurve> holes;
  bool isConvex;
  bool isComposite;
  std::vector<conway::geometry::IfcProfile> profiles;
};
conway::geometry::IfcProfile createNativeIfcProfile(
    const ParamsCreateNativeIfcProfile& parameters) {
  conway::geometry::IfcProfile profile;

  profile.type = "testType";
  profile.curve = parameters.curve;
  profile.holes = parameters.holes;
  profile.isConvex = parameters.isConvex;
  profile.isComposite = parameters.isComposite;
  profile.profiles = parameters.profiles;

  return profile;
}

conway::geometry::IfcProfile TransformProfile(
    conway::geometry::ConwayGeometryProcessor::ParamsTransformProfile*
        parameters) {
  if (processor) {
    return processor->transformProfile(parameters);
  }

  conway::geometry::IfcProfile transformProfile;
  return transformProfile;
}

glm::dmat4 getIdentityTransform() {
  return glm::dmat4(glm::dvec4(1, 0, 0, 0), glm::dvec4(0, 1, 0, 0),
                    glm::dvec4(0, 0, 1, 0), glm::dvec4(0, 0, 0, 1));
}

void deleteParamsGetPolyCurve(
    conway::geometry::ConwayGeometryProcessor::ParamsGetPolyCurve* params) {
  delete params;
}

void deleteParamsGetBooleanResult(
    conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult* params) {
  delete params;
}

void deleteParamsTransformProfile(
    conway::geometry::ConwayGeometryProcessor::ParamsTransformProfile* params) {
  delete params;
}

void deleteParamsGetTriangulatedFaceSetGeometry(
    conway::geometry::ConwayGeometryProcessor::
        ParamsGetTriangulatedFaceSetGeometry* params) {
  delete params;
}

glm::dmat3 placementIdentity2D(1);

glm::dmat3 GetIdentity2DMatrix() { return placementIdentity2D; }

glm::dmat4 placementIdentity3D(1);
glm::dmat4 GetIdentity3DMatrix() { return placementIdentity3D; }

std::vector<glm::dvec3> createVertexVector(uintptr_t verticesArray_,
                                          size_t length) {
  const double* verticesArray = reinterpret_cast<double*>(verticesArray_);
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  std::vector<glm::dvec3> vec;
  vec.resize(length / 3);
  for (size_t i = 0; i < length / 3; ++i) {
    glm::dvec3 vertex = glm::dvec3(verticesArray[i * 3], verticesArray[i * 3 + 1],
                                 verticesArray[i * 3 + 2]);

    vec[i] = vertex;
  }

  return vec;
}

std::vector<glm::dvec3> parseVertexVector( const conway::ParseBuffer& data ) {
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  std::vector<glm::dvec3> vec;

  if ( auto err = conway::parse_vector( data.range(), vec ).error; err != std::errc() ) {

    Logger::logError( "Error encountered parsing glm::dvec3 vector data %s", std::error_condition( err ).message().c_str() );
  }

  return vec;
}

double parseDouble( const conway::ParseBuffer& data ) {
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  const char* buffer    = reinterpret_cast< const char* >( data.data() );
  size_t         end    = data.size();
  const char* bufferEnd = buffer + end;
  double result;

  if ( auto err = fast_float::from_chars( buffer, bufferEnd, result ).ec; err != std::errc() ) {

    return NAN;
  }

  return result;
}

std::vector<glm::dvec2> parsePoint2DVector( const conway::ParseBuffer& data ) {
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  std::vector<glm::dvec2> vec;

  if ( auto err = conway::parse_vector( data.range(), vec ).error; err != std::errc() ) {

    Logger::logError( "Error encountered parsing glm::dvec2 vector data %s", std::error_condition( err ).message().c_str() );
  }

  return vec;
}

std::vector<glm::dvec2> parsePoint3Dto2DVector( const conway::ParseBuffer& data ) {
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  std::vector<glm::dvec2> vec;

  if ( auto err = conway::parse_vector_custom(
    data.range(),
    vec,
    []( const char* begin, const char* end, glm::dvec2& value ) {

      glm::dvec3 intermediate;

      auto result = fast_float::from_chars( begin, end, intermediate );

      value.x = intermediate.x;
      value.y = intermediate.y;

      return result;

    } ).error; err != std::errc() ) {

    Logger::logError( "Error encountered parsing vector of glm::dvec3 truncated to glm::dvec2 %s", std::error_condition( err ).message().c_str() );
  }

  return vec;
}

glm::dvec2 parsePoint2D( const conway::ParseBuffer& data ) {

  const char* buffer    = reinterpret_cast< const char* >( data.data() );
  size_t         end    = data.size();
  const char* bufferEnd = buffer + end;
  glm::dvec2 result;

  if ( auto err = fast_float::from_chars( buffer, bufferEnd, result ).ec; err != std::errc() ) {

    Logger::logError( "Error encountered parsing glm::dvec2 %s", std::error_condition( err ).message().c_str() );
  }

  return result;
}

glm::dvec2 parsePoint3D( const conway::ParseBuffer& data ) {

  const char* buffer    = reinterpret_cast< const char* >( data.data() );
  size_t         end    = data.size();
  const char* bufferEnd = buffer + end;
  glm::dvec2 result;

  if ( auto err = fast_float::from_chars( buffer, bufferEnd, result ).ec; err != std::errc() ) {

    Logger::logError( "Error encountered parsing glm::dvec2 %s", std::error_condition( err ).message().c_str() );
  }

  return result;
}


std::vector<uint32_t> parseUInt32Vector( const conway::ParseBuffer& data ) {
  // Assume 'vec' is a std::vector<glm::vec3> that's part of your class or
  // accessible here
  std::vector<uint32_t> vec;

  if ( auto err = conway::parse_vector( data.range(), vec ).error; err != std::errc() ) {

    Logger::logError( "Error encountered parsing glm::dvec3 data %s", std::error_condition( err ).message().c_str() );
  }

  return vec;
}

conway::geometry::Geometry GetTriangulatedFaceSetGeometry(
    const conway::geometry::ConwayGeometryProcessor::
        ParamsGetTriangulatedFaceSetGeometry& parameters) {
  if (processor) {
    return processor->getTriangulatedFaceSetGeometry(parameters);
  }

  conway::geometry::Geometry geom;
  return geom;
}

std::vector<conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>
buildIndexedPolygonalFaceVector(uintptr_t indicesArray_,
                                uint32_t indicesArrayLength,
                                uintptr_t startIndicesArray_,
                                uintptr_t polygonalFaceBufferOffsetsArray_,
                                uint32_t polygonalFaceBufferOffsetsLength,
                                uintptr_t startIndicesBufferOffsets_,
                                uint32_t startIndicesBufferOffsetsLength) {

  
  const uint32_t* indicesArray = reinterpret_cast<uint32_t*>(indicesArray_);
  const uint32_t* startIndicesArray =
      reinterpret_cast<uint32_t*>(startIndicesArray_);
  const uint32_t* polygonalFaceBufferOffsetsArray =
      reinterpret_cast<uint32_t*>(polygonalFaceBufferOffsetsArray_);
  const uint32_t* startIndicesBufferOffsets =
      reinterpret_cast<uint32_t*>(startIndicesBufferOffsets_);
  std::vector<conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>
      result;

  // Loop through each polygonal face buffer offset
  for (uint32_t index = 0; index < polygonalFaceBufferOffsetsLength - 1; ++index) {
    // Create a new IndexedPolygonalFace
    conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace
        indexedPolygonalFace;

    // The starting point for this face in the startIndicesArray
    uint32_t startOffset = startIndicesBufferOffsets[index];
    // The ending point for this face in the startIndicesArray
    uint32_t endOffset = startIndicesBufferOffsets[index + 1];

    //  Populate the face_starts vector with indices from the startIndicesArray
    for (uint32_t j = startOffset; j < endOffset; ++j) {
      indexedPolygonalFace.face_starts.push_back(startIndicesArray[j]);
    }

    // Now, populate the indices vector
    // If this is not the last face, the end index is the start index of the
    // next face If this is the last face, the end index is the total length of
    // indicesArray
    uint32_t indicesStart = polygonalFaceBufferOffsetsArray[index];
    uint32_t indicesEnd = polygonalFaceBufferOffsetsArray[index + 1];

    for (uint32_t k = indicesStart; k < indicesEnd; ++k) {
      indexedPolygonalFace.indices.push_back(indicesArray[k]);
    }

    // Add the constructed face to the result vector
    result.push_back(indexedPolygonalFace);
  }

  return result;
}

void resizeVectorVectorDouble(std::vector<std::vector<double>>& vec, size_t newSize) {
    vec.resize(newSize);
}

EMSCRIPTEN_BINDINGS(my_module) {
  /*
    active: boolean
    direction: NativeTransform
    profile: IfcProfile3D*/

  emscripten::value_object<conway::geometry::Extrusion>("ExtrusionSurface")
      .field("active", &conway::geometry::Extrusion::Active)
      .field("direction", &conway::geometry::Extrusion::Direction)
      .field("profile", &conway::geometry::Extrusion::Profile)
      .field("length", &conway::geometry::Extrusion::Length);

  emscripten::value_object<conway::geometry::Revolution>("RevolutionSurface")
      .field("active", &conway::geometry::Revolution::Active)
      .field("direction", &conway::geometry::Revolution::Direction)
      .field("profile", &conway::geometry::Revolution::Profile);

  emscripten::value_object<conway::geometry::Cylinder>("CylinderSurface")
      .field("active", &conway::geometry::Cylinder::Active)
      .field("radius", &conway::geometry::Cylinder::Radius);

  emscripten::value_object<conway::geometry::BSpline>("BSplineSurface")
      .field("active", &conway::geometry::BSpline::Active)
      .field("uDegree", &conway::geometry::BSpline::UDegree)
      .field("vDegree", &conway::geometry::BSpline::VDegree)
      .field("closedU", &conway::geometry::BSpline::ClosedU)
      .field("closedV", &conway::geometry::BSpline::ClosedV)
      .field("controlPoints", &conway::geometry::BSpline::ControlPoints)
      .field("uMultiplicity", &conway::geometry::BSpline::UMultiplicity)
      .field("vMultiplicity", &conway::geometry::BSpline::VMultiplicity)
      .field("uKnots", &conway::geometry::BSpline::UKnots)
      .field("vKnots", &conway::geometry::BSpline::VKnots)
      .field("weightPoints", &conway::geometry::BSpline::WeightPoints);

  emscripten::class_<conway::geometry::IfcSurface>("IfcSurface")
      .constructor<>()
      .property("transformation", &conway::geometry::IfcSurface::transformation)
      .property("bspline", &conway::geometry::IfcSurface::BSplineSurface)
      .property("cylinder", &conway::geometry::IfcSurface::CylinderSurface)
      .property("revolution", &conway::geometry::IfcSurface::RevolutionSurface)
      .property("extrusion", &conway::geometry::IfcSurface::ExtrusionSurface);

  emscripten::class_<conway::geometry::IfcBound3D>("IfcBound3D")
      .constructor<>();

  emscripten::class_<conway::geometry::Geometry>("IfcGeometry")
      .constructor<>()
      .function("GetVertexData", &conway::geometry::Geometry::GetVertexData)
      .function("GetVertexDataSize",
                &conway::geometry::Geometry::GetVertexDataSize)
      .function("getPoint", &conway::geometry::Geometry::GetPoint)
      .function("GetIndexData", &conway::geometry::Geometry::GetIndexData)
      .function("GetIndexDataSize",
                &conway::geometry::Geometry::GetIndexDataSize)
      .function("normalize",
                &conway::geometry::Geometry::Normalize)
      .function("appendGeometry",
                &conway::geometry::Geometry::AppendGeometry)
      .function("applyTransform",
                &conway::geometry::Geometry::ApplyTransform)
      .function("appendWithTransform",
                &conway::geometry::Geometry::AppendWithTransform)
      .function("getAllocationSize",
                &conway::geometry::Geometry::GetAllocationSize)
      .function("extractVertices",
                &conway::geometry::Geometry::ExtractVertices)
      .function("extractTriangles",
                &conway::geometry::Geometry::ExtractTriangles)
      .function("extractVerticesAndTriangles",
                &conway::geometry::Geometry::ExtractVerticesAndTriangles)
      // .function("getAABBCenter", &conway::geometry::Geometry::GetAABBCenter)
      .function("clone", &conway::geometry::Geometry::Clone)
      .function("dumpToOBJ", &conway::geometry::Geometry::GeometryToObj );

  emscripten::class_<conway::geometry::GeometryCollection>(
      "IfcGeometryCollection")
      .constructor<>()
      .function(
          "addComponentWithTransform",
          &conway::geometry::GeometryCollection::AddComponentWithTransform,
          emscripten::allow_raw_pointers())
      .property("materialIndex",
                &conway::geometry::GeometryCollection::materialIndex)
      .property("hasDefaultMaterial",
                &conway::geometry::GeometryCollection::hasDefaultMaterial)
      .property("currentSize",
                &conway::geometry::GeometryCollection::currentSize);

  emscripten::class_<conway::geometry::IfcCurve>("IfcCurve")
      .constructor<>()
      .function("add2d", &conway::geometry::IfcCurve::Add2d)
      .function("add3d", &conway::geometry::IfcCurve::Add3d)
      .function("getPointsSize", &conway::geometry::IfcCurve::GetPointsSize)
      .function("get2d", &conway::geometry::IfcCurve::Get2d)
      .function("get3d", &conway::geometry::IfcCurve::Get3d)
      .function("invert", &conway::geometry::IfcCurve::Invert)
      .function("isCCW", &conway::geometry::IfcCurve::IsCCW)
      .function("dumpToOBJ", &conway::geometry::IfcCurve::DumpToOBJ)
      .function("dumpToSVG", &conway::geometry::IfcCurve::DumpToSVG)
      .property("indices", &conway::geometry::IfcCurve::indices)
      .function("clone", &conway::geometry::IfcCurve::Clone);

  emscripten::class_<conway::geometry::IfcProfile>("IfcProfile")
      .constructor<>()
      .function("getType", &conway::geometry::IfcProfile::getType)
      .function("getCurve", &conway::geometry::IfcProfile::getCurve)
      .function("getHoles", &conway::geometry::IfcProfile::getHoles)
      .function("isConvex", &conway::geometry::IfcProfile::getIsConvex)
      .function("isComposite", &conway::geometry::IfcProfile::getIsComposite)
      .function("getProfiles", &conway::geometry::IfcProfile::getProfiles)
      .function("dumpToOBJ", &conway::geometry::IfcProfile::DumpToOBJ)
      .function("dumpToSVG", &conway::geometry::IfcProfile::DumpToSVG);

  emscripten::class_<glm::dmat4>("Glmdmat4")
      .constructor<>()
      .function("getValues", &getMatrixValues4x4)
      .function("setValues", &setMatrixValues4x4);

  emscripten::class_<glm::dmat3>("Glmdmat3")
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
  emscripten::register_vector<conway::geometry::Geometry*>(
      "geometryPointerArray");
  emscripten::register_vector<conway::geometry::GeometryCollection>(
      "geometryCollectionArray");

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsAxis1Placement3D>(
      "ParamsAxis1Placement3D")
      .field("position", &conway::geometry::ConwayGeometryProcessor::
                             ParamsAxis1Placement3D::position)
      .field("zAxisRef", &conway::geometry::ConwayGeometryProcessor::
                             ParamsAxis1Placement3D::zAxisRef)
      .field("normalizeZ", &conway::geometry::ConwayGeometryProcessor::
                               ParamsAxis1Placement3D::normalizeZ);

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
                              ParamsGetCircleCurve::placement)
      .field("thickness", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetCircleCurve::thickness);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetEllipseCurve>(
      "ParamsGetEllipseCurve")
      .field("radiusX", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetEllipseCurve::radiusX)
      .field("radiusY", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetEllipseCurve::radiusY)
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetEllipseCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetEllipseCurve::placement)
      .field("circleSegments", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetEllipseCurve::circleSegments);

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

  // conway::geometry::ConwayGeometryProcessor::ParamsGetHalfspaceSolid
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetHalfspaceSolid>(
      "ParamsGetHalfspaceSolid")
      .field("flipWinding", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetHalfspaceSolid::flipWinding)
      .field("optionalLinearScalingFactor",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetHalfspaceSolid::optionalLinearScalingFactor);

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

  // conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
  emscripten::class_<conway::geometry::ConwayGeometryProcessor::
                         ParamsGetTriangulatedFaceSetGeometry>(
      "ParamsGetTriangulatedFaceSetGeometry")
      .constructor<>()
      .property("points",
                &conway::geometry::ConwayGeometryProcessor::
                    ParamsGetTriangulatedFaceSetGeometry::pointsArray_)
      .property("pointsArrayLength",
                &conway::geometry::ConwayGeometryProcessor::
                    ParamsGetTriangulatedFaceSetGeometry::pointsArrayLength)
      .property("indices",
                &conway::geometry::ConwayGeometryProcessor::
                    ParamsGetTriangulatedFaceSetGeometry::indicesArray_)
      .property("indicesArrayLength",
                &conway::geometry::ConwayGeometryProcessor::
                    ParamsGetTriangulatedFaceSetGeometry::indicesArrayLength);

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
      "ParamsGetIfcTrimmedCurve")
      .field("masterRepresentation",
             &conway::geometry::ConwayGeometryProcessor::
                 ParamsGetIfcTrimmedCurve::masterRepresentation)
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetIfcTrimmedCurve::dimensions)
      .field("senseAgreement", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetIfcTrimmedCurve::senseAgreement)
      .field("trim1Cartesian2D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcTrimmedCurve::trim1Cartesian2D)
      .field("trim1Cartesian3D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcTrimmedCurve::trim1Cartesian3D)
      .field("trim1Double", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetIfcTrimmedCurve::trim1Double)
      .field("trim2Cartesian2D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcTrimmedCurve::trim2Cartesian2D)
      .field("trim2Cartesian3D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcTrimmedCurve::trim2Cartesian3D)
      .field("trim2Double", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetIfcTrimmedCurve::trim2Double)
      .field("trimExists", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetIfcTrimmedCurve::trimExists);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetIfcCircle
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetIfcCircle>(
      "ParamsGetIfcCircle")
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetIfcCircle::dimensions)
      .field("axis2Placement2D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcCircle::axis2Placement2D)
      .field("axis2Placement3D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcCircle::axis2Placement3D)
      .field("radius", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetIfcCircle::radius)
      .field("radius2", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetIfcCircle::radius2)
      .field("paramsGetIfcTrimmedCurve",
             &conway::geometry::ConwayGeometryProcessor::ParamsGetIfcCircle::
                 paramsGetIfcTrimmedCurve);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetIfcLine
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetIfcLine>(
      "ParamsGetIfcCircle")
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetIfcLine::dimensions)
      .field("cartesianPoint2D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcLine::cartesianPoint2D)
      .field("cartesianPoint3D", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetIfcLine::cartesianPoint3D)
      .field("vectorOrientation", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetIfcLine::vectorOrientation)
      .field("vectorMagnitude", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetIfcLine::vectorMagnitude)
      .field("isEdge", &conway::geometry::ConwayGeometryProcessor::ParamsGetIfcLine::isEdge)
      .field("paramsGetIfcTrimmedCurve",
             &conway::geometry::ConwayGeometryProcessor::ParamsGetIfcLine::
                 paramsGetIfcTrimmedCurve);

  emscripten::class_<
      conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult>(
      "ParamsGetBooleanResult")
      .constructor<>()
      .property("flatFirstMesh", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetBooleanResult::flatFirstMesh)
      .property("flatSecondMesh", &conway::geometry::ConwayGeometryProcessor::
                                      ParamsGetBooleanResult::flatSecondMesh)
      .property("operatorType", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsGetBooleanResult::operatorType)
      .property("isSubtractOperand", &conway::geometry::ConwayGeometryProcessor::
                                          ParamsGetBooleanResult::isSubtractOperand);

  emscripten::class_<
      conway::geometry::ConwayGeometryProcessor::ParamsTransformProfile>(
      "ParamsTransformProfile")
      .constructor<>()
      .property("transformation", &conway::geometry::ConwayGeometryProcessor::
                                      ParamsTransformProfile::transformation)
      .property("profile", &conway::geometry::ConwayGeometryProcessor::
                               ParamsTransformProfile::profile);

  emscripten::class_<
      conway::geometry::ConwayGeometryProcessor::ParamsGetPolyCurve>(
      "ParamsGetPolyCurve")
      .constructor<>()
      .property("points", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetPolyCurve::points_)
      .property("pointsLength", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsGetPolyCurve::pointsLength)
      .property("dimensions", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsGetPolyCurve::dimensions)
      .property("senseAgreement", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsGetPolyCurve::senseAgreement)
      .property("isEdge", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsGetPolyCurve::isEdge);

  /*emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult>(
      "ParamsGetBooleanResult")
      .field("flatFirstMesh", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsGetBooleanResult::flatFirstMesh)
      .field("flatSecondMesh", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsGetBooleanResult::flatSecondMesh)
      .field("operatorType", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetBooleanResult::operatorType);*/

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsRelVoidSubtract>(
      "ParamsRelVoidSubtract")
      .field("flatFirstMesh", &conway::geometry::ConwayGeometryProcessor::
                                  ParamsRelVoidSubtract::flatFirstMesh)
      .field("flatSecondMesh", &conway::geometry::ConwayGeometryProcessor::
                                   ParamsRelVoidSubtract::flatSecondMesh)
      .field("operatorType", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsRelVoidSubtract::operatorType)
      .field("parentMatrix", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsRelVoidSubtract::parentMatrix);

  // conway::geometry::ConwayGeometryProcessor::ParamsGetLoop
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetLoop>("ParamsGetLoop")
      .field("points",
             &conway::geometry::ConwayGeometryProcessor::ParamsGetLoop::points)
      .field("edges",
             &conway::geometry::ConwayGeometryProcessor::ParamsGetLoop::edges);

  // ParamsCreateBound3D
  emscripten::value_object<ParamsCreateBound3D>("ParamsCreateBound3D")
      .field("curve", &ParamsCreateBound3D::curve)
      .field("orientation", &ParamsCreateBound3D::orientation)
      .field("type", &ParamsCreateBound3D::type);

  // conway::geometry::ConwayGeometryProcessor::ParamsAddFaceToGeometry
  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsAddFaceToGeometry>(
      "ParamsAddFaceToGeometry")
      .field("boundsArray", &conway::geometry::ConwayGeometryProcessor::
                                ParamsAddFaceToGeometry::boundsArray)
      .field("advancedBrep", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsAddFaceToGeometry::advancedBrep)
      .field("surface", &conway::geometry::ConwayGeometryProcessor::
                            ParamsAddFaceToGeometry::surface)
      .field("scaling", &conway::geometry::ConwayGeometryProcessor::
                            ParamsAddFaceToGeometry::scaling);

  // ifc::IFCRECTANGLEPROFILEDEF
  // ifc::IFCROUNDEDRECTANGLEPROFILEDEF
  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::
                               ParamsGetRectangleProfileCurve>(
      "ParamsGetRectangleProfileCurve")
      .field("xDim", &conway::geometry::ConwayGeometryProcessor::
                         ParamsGetRectangleProfileCurve::xDim)
      .field("yDim", &conway::geometry::ConwayGeometryProcessor::
                         ParamsGetRectangleProfileCurve::yDim)
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetRectangleProfileCurve::hasPlacement)
      .field("matrix", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetRectangleProfileCurve::matrix)
      .field("thickness", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetRectangleProfileCurve::thickness);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetCShapeCurve>(
      "ParamsGetCShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetCShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetCShapeCurve::placement)
      .field("hasFillet", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetCShapeCurve::hasFillet)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetCShapeCurve::depth)
      .field("width", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetCShapeCurve::width)
      .field("thickness", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetCShapeCurve::thickness)
      .field("girth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetCShapeCurve::girth)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetCShapeCurve::filletRadius);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetIShapeCurve>(
      "ParamsGetIShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetIShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetIShapeCurve::placement)
      .field("hasFillet", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetIShapeCurve::hasFillet)
      .field("width", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetIShapeCurve::width)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetIShapeCurve::depth)
      .field("webThickness", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetIShapeCurve::webThickness)
      .field("flangeThickness", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsGetIShapeCurve::flangeThickness)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetIShapeCurve::filletRadius);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetLShapeCurve>(
      "ParamsGetLShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetLShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetLShapeCurve::placement)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetLShapeCurve::filletRadius)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetLShapeCurve::depth)
      .field("width", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetLShapeCurve::width)
      .field("thickness", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetLShapeCurve::thickness)
      .field("edgeRadius", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetLShapeCurve::edgeRadius)
      .field("legSlope", &conway::geometry::ConwayGeometryProcessor::
                             ParamsGetLShapeCurve::legSlope);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetTShapeCurve>(
      "ParamsGetTShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetTShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetTShapeCurve::placement)
      .field("hasFillet", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetTShapeCurve::hasFillet)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetTShapeCurve::depth)
      .field("width", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetTShapeCurve::width)
      .field("webThickness", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetTShapeCurve::webThickness)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetTShapeCurve::filletRadius)
      .field("flangeEdgeRadius", &conway::geometry::ConwayGeometryProcessor::
                                     ParamsGetTShapeCurve::flangeEdgeRadius)
      .field("flangeScope", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetTShapeCurve::flangeScope);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetUShapeCurve>(
      "ParamsGetUShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetUShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetUShapeCurve::placement)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetUShapeCurve::depth)
      .field("flangeWidth", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetUShapeCurve::flangeWidth)
      .field("webThickness", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetUShapeCurve::webThickness)
      .field("flangeThickness", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsGetUShapeCurve::flangeThickness)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetUShapeCurve::filletRadius)
      .field("edgeRadius", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetUShapeCurve::edgeRadius)
      .field("flangeScope", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetUShapeCurve::flangeScope);

  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetZShapeCurve>(
      "ParamsGetZShapeCurve")
      .field("hasPlacement", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetZShapeCurve::hasPlacement)
      .field("placement", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetZShapeCurve::placement)
      .field("hasFillet", &conway::geometry::ConwayGeometryProcessor::
                              ParamsGetZShapeCurve::hasFillet)
      .field("depth", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetZShapeCurve::depth)
      .field("flangeWidth", &conway::geometry::ConwayGeometryProcessor::
                                ParamsGetZShapeCurve::flangeWidth)
      .field("webThickness", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetZShapeCurve::webThickness)
      .field("flangeThickness", &conway::geometry::ConwayGeometryProcessor::
                                    ParamsGetZShapeCurve::flangeThickness)
      .field("filletRadius", &conway::geometry::ConwayGeometryProcessor::
                                 ParamsGetZShapeCurve::filletRadius)
      .field("edgeRadius", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetZShapeCurve::edgeRadius);

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

  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid>(
    "ParamsGetSweptDiskSolid")
    .field("directrix", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::directrix)
    .field("radius", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::radius)
    .field("innerRadius", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::innerRadius)
    .field("startParam", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::startParam)
    .field("endParam", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::endParam)
    .field("closed", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::closed)
    .field("circleSegments", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::circleSegments)
    .field("scalingFactor", &conway::geometry::ConwayGeometryProcessor::ParamsGetSweptDiskSolid::scalingFactor);


  emscripten::value_object<
      conway::geometry::ConwayGeometryProcessor::ParamsGetBSplineCurve>(
      "ParamsGetBSplineCurve")
      .field("dimensions", &conway::geometry::ConwayGeometryProcessor::
                               ParamsGetBSplineCurve::dimensions)
      .field("degree", &conway::geometry::ConwayGeometryProcessor::
                           ParamsGetBSplineCurve::degree)
      .field("points2", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetBSplineCurve::points2)
      .field("points3", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetBSplineCurve::points3)
      .field("knots", &conway::geometry::ConwayGeometryProcessor::
                          ParamsGetBSplineCurve::knots)
      .field("weights", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetBSplineCurve::weights)
      .field("senseAgreement", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetBSplineCurve::senseAgreement)
      .field("isEdge", &conway::geometry::ConwayGeometryProcessor::
                            ParamsGetBSplineCurve::isEdge);

  emscripten::value_object<conway::geometry::IfcTrimmingSelect>(
      "TrimmingSelect")
      .field("hasParam", &conway::geometry::IfcTrimmingSelect::hasParam)
      .field("hasPos", &conway::geometry::IfcTrimmingSelect::hasPos)
      .field("hasLength", &conway::geometry::IfcTrimmingSelect::hasLength)
      .field("param", &conway::geometry::IfcTrimmingSelect::param)
      .field("pos", &conway::geometry::IfcTrimmingSelect::pos)
      .field("pos3D", &conway::geometry::IfcTrimmingSelect::pos3D);

  emscripten::value_object<conway::geometry::IfcTrimmingArguments>(
      "TrimmingArguments")
      .field("exist", &conway::geometry::IfcTrimmingArguments::exist)
      .field("start", &conway::geometry::IfcTrimmingArguments::start)
      .field("end", &conway::geometry::IfcTrimmingArguments::end);

  /**
   * bool agreement = false;
    double scaleFactor = 1.0;
    conway::geometry::IfcCurve curve;
    IfcSurface surface;
    glm::dmat4 position;
  */

  emscripten::value_object<conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace>(
    "ParamsGetPolygonalBoundedHalfspace")
    .field("agreement", &conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace::agreement)
    .field("scaleFactor", &conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace::scaleFactor)
    .field("surface", &conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace::surface)
    .field("curve", &conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace::curve)
    .field("position", &conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace::position);

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

  emscripten::class_< conway::ParseBuffer >( "ParseBuffer" )
    .constructor<>()
    .function( "resize", &conway::ParseBuffer::resize )
    .function( "data", &conway::ParseBuffer::dataInteger )
    .function( "size", &conway::ParseBuffer::size )
    .function( "capacity", &conway::ParseBuffer::capacity );

  emscripten::register_vector<double>("vectorDouble");
  emscripten::register_vector<std::vector<double>>("vectorVectorDouble");
  emscripten::register_vector<glm::vec2>("vec2Array");
  emscripten::register_vector<glm::vec3>("glmVec3Array");
  emscripten::register_vector<glm::dvec3>("glmdVec3Array");
  emscripten::register_vector<std::vector<glm::dvec3>>("glmdVec3ArrayArray");
  emscripten::register_vector<glm::dvec2>("glmdVec2Array");
  emscripten::register_vector<std::string>("stringVector");
  emscripten::register_vector<uint32_t>("UintVector");
  emscripten::register_vector<uint8_t>("VectorUint8");
  emscripten::register_vector<std::vector<uint8_t>>("VectorVectorUint8");
  emscripten::register_vector<size_t>("ULongVector");
  emscripten::register_vector<conway::geometry::IfcCurve>("curveArray");
  emscripten::register_vector<conway::geometry::IfcProfile>("profileArray");
  emscripten::register_vector<conway::geometry::Geometry>("geometryArray");
  emscripten::register_vector<conway::geometry::IfcBound3D>("Bound3DArray");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace>(
      "VectorIndexedPolygonalFace");
  emscripten::register_vector<
      conway::geometry::ConwayGeometryProcessor::Segment>("VectorSegment");

  emscripten::function("createVertexVector", &createVertexVector,
                       emscripten::allow_raw_pointers());
  emscripten::function("parseVertexVector", &parseVertexVector,
                       emscripten::allow_raw_pointers());
  emscripten::function("parseDouble", &parseDouble,
                       emscripten::allow_raw_pointers());
  emscripten::function("parsePoint2DVector", &parsePoint2DVector,
                       emscripten::allow_raw_pointers());
  emscripten::function("parsePoint3DTo2DVector", &parsePoint3Dto2DVector,
                       emscripten::allow_raw_pointers());                       
  emscripten::function("parseUint32Vector", &parseUInt32Vector,
                       emscripten::allow_raw_pointers());                          
  emscripten::function("getPolygonalFaceSetGeometry",
                       &GetPolygonalFaceSetGeometry);
  emscripten::function("getTriangulatedFaceSetGeometry",
                       &GetTriangulatedFaceSetGeometry,
                       emscripten::allow_raw_pointers());
  emscripten::function("getIndexedPolyCurve", &GetIndexedPolyCurve);
  emscripten::function("getCircleCurve", &GetCircleCurve);
  emscripten::function("initializeGeometryProcessor",
                       &InitializeGeometryProcessor);
  emscripten::function("freeGeometryProcessor", &FreeGeometryProcessor);
  emscripten::function("geometryToObj", &GeometryToObj);
  emscripten::function("geometryToGltf", &GeometryToGltf);
  emscripten::function("getAxis1Placement", &GetAxis1Placement);
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
  emscripten::function("getPolygonalBoundedHalfspace", &GetPolygonalBoundedHalfspace);
  emscripten::function("getHalfSpaceSolid", &GetHalfSpaceSolid);
  emscripten::function("getBooleanResult", &GetBooleanResult,
                       emscripten::allow_raw_pointers());
  emscripten::function("deleteParamsGetBooleanResult",
                       &deleteParamsGetBooleanResult,
                       emscripten::allow_raw_pointers());
  emscripten::function("deleteParamsGetPolyCurve", &deleteParamsGetPolyCurve,
                       emscripten::allow_raw_pointers());
  emscripten::function("transformProfile", &TransformProfile,
                       emscripten::allow_raw_pointers());
  emscripten::function("deleteParamsTransformProfile",
                       &deleteParamsTransformProfile,
                       emscripten::allow_raw_pointers());
  emscripten::function("deleteParamsGetTriangulatedFaceSetGeometry",
                       &deleteParamsGetTriangulatedFaceSetGeometry,
                       emscripten::allow_raw_pointers());
  emscripten::function("relVoidSubtract", &RelVoidSubtract);
  emscripten::function("getIfcCircle", &GetIfcCircle);
  emscripten::function("getIfcLine", &GetIfcLine);
  emscripten::function("getBSplineCurve", &GetBSplineCurve);
  emscripten::function("getLoop", &GetLoop);
  emscripten::function("createBound3D", &createBound3D);
  emscripten::function("addFaceToGeometry", &AddFaceToGeometry);
  emscripten::function("getRectangleProfileCurve", &GetRectangleProfileCurve);
  emscripten::function("getIdentityTransform", &getIdentityTransform);
  emscripten::function("multiplyNativeMatrices", &multiplyNativeMatrices);
  emscripten::function("getRectangleHollowProfileHole",
                       &GetRectangleHollowProfileHole);
  emscripten::function("getCircleHoleCurve", &GetCircleHoleCurve);
  emscripten::function("getEllipseCurve", &GetEllipseCurve);
  emscripten::function("getCShapeCurve", &GetCShapeCurve);
  emscripten::function("getIShapeCurve", &GetIShapeCurve);
  emscripten::function("getTShapeCurve", &GetTShapeCurve);
  emscripten::function("getLShapeCurve", &GetLShapeCurve);
  emscripten::function("getUShapeCurve", &GetUShapeCurve);
  emscripten::function("getZShapeCurve", &GetZShapeCurve);
  emscripten::function("getIdentity2DMatrix", &GetIdentity2DMatrix);
  emscripten::function("getIdentity3DMatrix", &GetIdentity3DMatrix);
  emscripten::function("buildIndexedPolygonalFaceVector",
                       &buildIndexedPolygonalFaceVector,
                       emscripten::allow_raw_pointers());
  emscripten::function("getPolyCurve", &GetPolyCurve,
                       emscripten::allow_raw_pointers());

  emscripten::function("resizeVectorVectorDouble", 
  &resizeVectorVectorDouble, emscripten::allow_raw_pointers());
  emscripten::function("getSweptDiskSolid", 
  &GetSweptDiskSolid, emscripten::allow_raw_pointers());
}
