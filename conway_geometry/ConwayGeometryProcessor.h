/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/IfcGeometryProcessor.cpp
 */

#pragma once

#include <tinynurbs/tinynurbs.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <filesystem>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <iostream>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#pragma clang diagnostic ignored "-Wunused-function"
#include <mapbox/earcut.hpp>
#pragma clang diagnostic push

#include <span>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

// draco
#include <draco/compression/config/compression_shared.h>
#include <draco/compression/encode.h>
#include <draco/compression/expert_encode.h>
#include <draco/core/cycle_timer.h>
#include <draco/io/file_utils.h>
#include <draco/io/mesh_io.h>
#include <draco/io/point_cloud_io.h>

// GLTFSDK
#include <GLTFSDK/BufferBuilder.h>
#include <GLTFSDK/ExtensionsKHR.h>
#include <GLTFSDK/GLBResourceWriter.h>
#include <GLTFSDK/GLTF.h>
#include <GLTFSDK/GLTFResourceWriter.h>
#include <GLTFSDK/IStreamWriter.h>
#include <GLTFSDK/Serialize.h>

#include "operations/geometry_utils.h"
// #include "representation/IfcGeometry.h"
#include "representation/Geometry.h"
#include "structures/parse_buffer.h"

namespace fuzzybools {
struct Geometry;
}

namespace {
// The glTF SDK is decoupled from all file I/O by the IStreamWriter (and
// IStreamReader) interface(s) and the C++ stream-based I/O library. This allows
// the glTF SDK to be used in sandboxed environments, such as WebAssembly
// modules and UWP apps, where any file I/O code must be platform or use-case
// specific.
class StreamWriter : public Microsoft::glTF::IStreamWriter {
 public:
  StreamWriter(const std::filesystem::path& pathBase) : m_pathBase(pathBase) {
    // printf( "%s\n", m_pathBase.c_str() );
    // assert(m_pathBase.has_root_path());
  }

  const std::vector<std::string>& getUris() const { return uris; }

  // Resolves the relative URIs of any external resources declared in the glTF
  // manifest
  std::shared_ptr<std::ostream> GetOutputStream(
      const std::string& filename) const override {
    auto stream = std::make_shared<std::ostringstream>();
    uris.push_back(filename);
    return stream;
  }

 private:
  std::filesystem::path m_pathBase;
  mutable std::vector<std::string> uris;
};
}  // namespace

namespace conway::geometry {
// TODO(nickcastel50): Pass these into Geometry to GLTF + GLB as a parameter
struct DracoOptions {
  bool isPointCloud = false;
  int posQuantizationBits = 11;
  int texCoordsQuantizationBits = 8;
  bool texCoordsDeleted = false;
  int normalsQuantizationBits = 8;
  bool normalsDeleted = false;
  int genericQuantizationBits = 8;
  bool genericDeleted = false;
  int compressionLevel = 7;
  bool preservePolygons = false;
  bool useMetadata = false;
  bool deduplicateInputValues = true;
};

/*
 * The purpose of this class is to process all of the geometry extraction /
 * transform request queries as required by the IFC file format. This will later
 * be extended to handle other 3D file formats.
 */
class ConwayGeometryProcessor {
 public:
  ConwayGeometryProcessor() {}

  // case ifc::IFCMAPPEDITEM:
  struct ParamsGetMappedItem {
    glm::dmat4 transformation = glm::dmat4();
    IfcComposedMesh ifcPresentationMesh;
  };

  Geometry BoolSubtract(std::vector<Geometry>& firstGroups,
                        std::vector<Geometry>& secondGroups,
                        bool isSubtractOperand);
  Geometry BoolSubtractLegacy(const std::vector<Geometry>& firstGeoms,
                              std::vector<Geometry>& secondGeoms);

  struct ParamsGetBlock {
    uint32_t xLength = 0;
    uint32_t yLength = 0;
    uint32_t zLength = 0;
    glm::dmat4 placement = glm::dmat4(1.0);
  };

  Geometry GetBlock(const ParamsGetBlock& parameters);

  // case ifc::IFCBOOLEANCLIPPINGRESULT:
  // case ifc::IFCBOOLEANRESULT:
  struct ParamsGetBooleanResult {
    std::vector<Geometry> flatFirstMesh;
    std::vector<Geometry> flatSecondMesh;
    int operatorType = 2;
    bool isSubtractOperand = false;
  };
  Geometry GetBooleanResult(ParamsGetBooleanResult* parameters);

  struct ParamsRelVoidSubtract {
    std::vector<Geometry> flatFirstMesh;
    std::vector<Geometry> flatSecondMesh;
    int operatorType = 2;
    glm::dmat4 parentMatrix;
  };

  Geometry RelVoidSubtract(ParamsRelVoidSubtract& parameters);

  // case ifc::IFCHALFSPACESOLID:
  struct ParamsGetHalfspaceSolid {
    bool flipWinding = false;
    double optionalLinearScalingFactor = 1.0;
  };

  Geometry GetHalfSpaceSolid(const ParamsGetHalfspaceSolid& parameters);

  // case ifc::IFCPOLYGONALBOUNDEDHALFSPACE
  struct ParamsGetPolygonalBoundedHalfspace {
    bool agreement = false;
    double scaleFactor = 1.0;
    conway::geometry::IfcCurve curve;
    IfcSurface surface;
    glm::dmat4 position;
  };

  /**
   * scaleFactor:number
  agreement:boolean
  curve:CurveObject | undefined
  surface:SurfaceObject | undefined
  matrix: any // glm::dmat4
  */

  Geometry GetPolygonalBoundedHalfspace(
      const ParamsGetPolygonalBoundedHalfspace& parameters);

  // case ifc::IFCREPRESENTATIONMAP
  // TODO(nickcastel50) : see if this is needed

  // case ifc::IFCCONNECTEDFACESET:
  // case ifc::IFCCLOSEDSHELL:
  // case ifc::IFCOPENSHELL:
  // These cases are handled by getBrep()
  struct ParamsAddFaceToGeometry {
    std::vector<IfcBound3D> boundsArray;
    bool advancedBrep = false;
    IfcSurface surface;
    double scaling;
  };

  struct ParamsAddFaceToGeometrySimple {
    std::vector<IfcBound3D> boundsArray;
    double scaling;
  };

  void AddFaceToGeometrySimple(ParamsAddFaceToGeometrySimple& parameters,
                               Geometry& geometry);

  struct ParamsGetBrep {
    uint32_t boundsSize = 0;
    uint32_t indicesPerFace = 0;
    size_t numIndices = 0;
    std::vector<uint32_t> indices;
    std::vector<IfcBound3D> boundsArray;
    bool advancedBrep = false;
    IfcSurface surface;
  };

  // IfcGeometry getBrep(ParamsGetBrep parameters);

  // case ifc::IFCFACE:
  // case ifc::IFCADVANCEDFACE:
  void AddFaceToGeometry(ParamsAddFaceToGeometry& parameters,
                         Geometry& geometry);

  // case ifc::IFCRECTANGLEPROFILEDEF:
  // case ifc::IFCROUNDEDRECTANGLEPROFILEDEF:
  struct ParamsGetRectangleProfileCurve {
    double xDim = 0.0f;
    double yDim = 0.0f;
    bool hasPlacement = false;
    bool hasRoundingRadius = false;
    double roundingRadius = 0.0f;
    glm::dmat3 matrix;
    double thickness = -1.0f;
    int circleSegments = 12;
  };

  struct ParamsGetCShapeCurve {
    bool hasPlacement = false;
    glm::dmat3 placement;
    bool hasFillet = false;
    double depth = 0.0f;
    double width = 0.0f;
    double thickness = 0.0f;
    double girth = 0.0f;
    double filletRadius = 0.0f;
  };

  IfcCurve GetCShapeCurve(const ParamsGetCShapeCurve& parameters);

  struct ParamsGetIShapeCurve {
    bool hasPlacement = false;
    // Assuming placement is some sort of matrix or object
    glm::dmat3 placement;
    bool hasFillet = false;
    double width = 0.0f;
    double depth = 0.0f;
    double webThickness = 0.0f;
    double flangeThickness = 0.0f;
    double filletRadius = 0.0f;
  };

  IfcCurve GetIShapeCurve(const ParamsGetIShapeCurve& parameters);

  struct ParamsGetLShapeCurve {
    bool hasPlacement = false;
    glm::dmat3 placement;
    bool hasFillet = false;
    double filletRadius = 0.0f;
    double depth = 0.0f;
    double width = 0.0f;
    double thickness = 0.0f;
    double edgeRadius = 0.0f;
    double legSlope = 0.0f;
  };

  IfcCurve GetLShapeCurve(const ParamsGetLShapeCurve& parameters);

  struct ParamsGetTShapeCurve {
    bool hasPlacement = false;
    glm::dmat3 placement;
    bool hasFillet = false;
    double depth = 0.0f;
    double width = 0.0f;
    double webThickness = 0.0f;
    double filletRadius = 0.0f;
    double flangeEdgeRadius = 0.0f;
    double flangeScope = 0.0f;
  };

  IfcCurve GetTShapeCurve(const ParamsGetTShapeCurve& parameters);

  struct ParamsGetUShapeCurve {
    bool hasPlacement = false;
    glm::dmat3 placement;
    double depth = 0.0f;
    double flangeWidth = 0.0f;
    double webThickness = 0.0f;
    double flangeThickness = 0.0f;
    double filletRadius = 0.0f;
    double edgeRadius = 0.0f;
    double flangeScope = 0.0f;
  };

  IfcCurve GetUShapeCurve(const ParamsGetUShapeCurve& parameters);

  struct ParamsGetZShapeCurve {
    bool hasPlacement = false;
    glm::dmat3 placement;
    bool hasFillet = false;
    double depth = 0.0f;
    double flangeWidth = 0.0f;
    double webThickness = 0.0f;
    double flangeThickness = 0.0f;
    double filletRadius = 0.0f;
    double edgeRadius = 0.0f;
  };

  IfcCurve GetZShapeCurve(const ParamsGetZShapeCurve& parameters);

  IfcCurve GetRectangleProfileCurve(
      const ParamsGetRectangleProfileCurve& parameters);

  IfcCurve GetRectangleHollowProfileHole(
      const ParamsGetRectangleProfileCurve& parameters);

  // case ifc::IFCFACEBASEDSURFACEMODEL:
  // case ifc::IFCSHELLBASEDSURFACEMODEL:
  struct ParamsGetSurfaceModel {
    uint32_t numShellRefs = 0;
    std::vector<ParamsGetBrep> shells;
  };

  // std::vector<IfcGeometry> GetSurfaceModel(ParamsGetSurfaceModel parameters);

  // case ifc::IFCPLANE:
  // case ifc::IFCBSPLINESURFACE:
  // case ifc::IFCBSPLINESURFACEWITHKNOTS:
  // case ifc::IFCRATIONALBSPLINESURFACEWITHKNOTS:
  struct ParamsGetSurface {
    bool isPlane = false;
    glm::dmat4 transformation;
    bool isBsplineSurface = false;
    double Udegree = 0;
    double Vdegree = 0;
    // TODO(nickcastel50): How do we pass these across?
    std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;
    std::string curveType;
    bool closedU;
    bool closedV;
    std::string selfIntersect;
    bool isBsplineSurfaceWithKnots = false;
    std::vector<glm::f64> UMultiplicity;
    std::vector<glm::f64> VMultiplicity;
    std::vector<glm::f64> UKnots;
    std::vector<glm::f64> VKnots;
    bool isRationalBsplineSurfaceWithKnots = false;
    std::vector<std::vector<glm::f64>> weightPts;
    bool isCylindricalSurface = false;
    double radius = 0;
    bool isSurfaceOfRevolution = false;
    glm::dmat4 revolutionDirection;
    IfcProfile revolutionProfile;
    bool includeTransformation = false;
    bool isSurfaceOfLinearExtrusion = false;
    glm::dvec3 extrusionDirection;
    IfcProfile extrusionProfile;
    bool customLength = false;
    double length = 0;
  };
  IfcSurface GetSurface(const ParamsGetSurface& parameters);

  // case ifc::IFCAXIS2PLACEMENT2D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM:
  struct ParamsGetAxis2Placement2D {
    bool isAxis2Placement2D = false;
    bool isCartesianTransformationOperator2D = false;
    bool isCartesianTransformationOperator2DNonUniform = false;
    glm::dvec2 position2D;
    bool customAxis1Ref = false;
    glm::dvec2 axis1Ref;
    bool customAxis2Ref = false;
    glm::dvec2 axis2Ref;
    bool customScale = false;
    double scale1 = 0;
    bool customScale2 = false;
    double scale2 = 0;
  };

  glm::dmat3 GetAxis2Placement2D(const ParamsGetAxis2Placement2D& parameters);

  // case ifc::IFCAXIS1PLACEMENT:
  struct ParamsAxis1Placement3D {
    glm::dvec3 position;
    glm::dvec3 zAxisRef;
    bool normalizeZ = false;
  };
  glm::dmat4 GetAxis1Placement(const ParamsAxis1Placement3D& parameters);

  // case ifc::IFCAXIS2PLACEMENT3D:
  struct ParamsAxis2Placement3D {
    glm::dvec3 position;
    glm::dvec3 zAxisRef;
    glm::dvec3 xAxisRef;
    bool normalizeZ = false;
    bool normalizeX = false;
  };
  glm::dmat4 GetAxis2Placement3D(const ParamsAxis2Placement3D& parameters);

  // case ifc::IFCLOCALPLACEMENT:
  // This case just recursively calls GetLocalPlacement, not sure if needed. See
  // GetLocalPlacement
  struct ParamsLocalPlacement {
    bool useRelPlacement = false;
    glm::dmat4 axis2Placement;
    glm::dmat4 relPlacement;
  };

  glm::dmat4 GetLocalPlacement(const ParamsLocalPlacement& parameters);

  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM:
  struct ParamsCartesianTransformationOperator3D {
    glm::dvec3 position;
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
    double scale3_ = 0;
  };

  static glm::dmat4 GetCartesianTransformationOperator3D(
      const ParamsCartesianTransformationOperator3D& parameters);

  // case ifc::IFCPOLYLOOP:
  // case ifc::IFCEDGELOOP:
  struct ParamsGetLoop {
    bool isEdgeLoop = false;
    std::vector<glm::dvec3> points;
    std::vector<conway::geometry::IfcCurve> edges;
  };

  IfcCurve GetLoop(const ParamsGetLoop& parameters);

  struct ParamsGetBound {
    bool isFaceOuterBound = false;
    bool orient = false;
    ParamsGetLoop parametersGetLoop;
  };
  IfcBound3D GetBound(const ParamsGetBound& parameters);

  struct IndexedPolygonalFace {
    std::vector<uint32_t> indices;
    std::vector<size_t> face_starts;
  };

  // case ifc::IFCINDEXEDPOLYGONALFACEWITHVOIDS:
  // case ifc::IFCINDEXEDPOLYGONALFACE:
  struct ParamsReadIndexedPolygonalFace {
    const std::vector<glm::dvec3>& points;
    const IndexedPolygonalFace& face;

    ParamsReadIndexedPolygonalFace(const std::vector<glm::dvec3>& points_ref,
                                   const IndexedPolygonalFace& face_ref)
        : points(points_ref), face(face_ref) {}
  };
  std::vector<IfcBound3D> ReadIndexedPolygonalFace(
      const ParamsReadIndexedPolygonalFace& parameters);

  struct ResultsGltf {
    bool success = false;
    std::vector<std::string> bufferUris;
    std::vector<std::vector<uint8_t>> buffers;
  };
  ResultsGltf GeometryToGltf(
      std::span<conway::geometry::GeometryCollection> geom,
      std::span<conway::geometry::Material> materials, bool isGlb,
      bool outputDraco, std::string filePath, bool outputFile,
      glm::dmat4 transform = glm::dmat4(1));

  std::string GeometryToObj(const conway::geometry::Geometry& geom,
                            size_t& offset,
                            glm::dmat4 transform = glm::dmat4(1));

  struct ParamsTransformProfile {
    glm::dmat3 transformation;
    IfcProfile profile;
  };

  IfcProfile transformProfile(ParamsTransformProfile* parameters);

  struct ParamsGetPolyCurve {
    uintptr_t points_ = 0;
    uint32_t pointsLength = 0;
    uint32_t dimensions = 0;
    bool senseAgreement = true;
    bool isEdge = false;
  };

  IfcCurve getPolyCurve(const ParamsGetPolyCurve& parameters);

  // casae ifc::IFCTRIANGULATEDFACESET
  struct ParamsGetTriangulatedFaceSetGeometry {
    uintptr_t indicesArray_;
    uint32_t indicesArrayLength = 0;
    uintptr_t pointsArray_;
    uint32_t pointsArrayLength = 0;
  };

  Geometry getTriangulatedFaceSetGeometry(
      const ParamsGetTriangulatedFaceSetGeometry& parameters);

  // case ifc::IFCPOLYGONALFACESET:
  struct ParamsGetPolygonalFaceSetGeometry {
    uint32_t indicesPerFace = 0;
    std::vector<glm::dvec3> points;
    std::vector<IndexedPolygonalFace> faces;
  };
  Geometry getPolygonalFaceSetGeometry(
      const ParamsGetPolygonalFaceSetGeometry& parameters);

  struct Segment {
    bool isArcType = false;
    std::vector<uint32_t> indices;
  };

  // case ifc::IFCSWEPTDISKSOLID
  struct ParamsGetSweptDiskSolid {
    IfcCurve directrix;
    double radius;
    double innerRadius;
    double startParam;
    double endParam;
    bool closed;
    uint32_t circleSegments;
    double scalingFactor;
  };

  conway::geometry::Geometry getSweptDiskSolid(
      const ParamsGetSweptDiskSolid& parameters);

  // case ifc::IFCINDEXEDPOLYCURVE
  struct ParamsGetIfcIndexedPolyCurve {
    uint32_t dimensions = 2;
    std::vector<Segment> segments;
    std::vector<glm::dvec2> points;
  };

  struct ParamsGetIfcIndexedPolyCurve3D {
    uint32_t dimensions = 2;
    std::vector<Segment> segments;
    std::vector<glm::dvec3> points;
  };

  conway::geometry::IfcCurve getIndexedPolyCurve(
      const ParamsGetIfcIndexedPolyCurve& parameters);

  conway::geometry::IfcCurve getIndexedPolyCurve3D(
      const ParamsGetIfcIndexedPolyCurve3D& parameters);

  // case ifc::CIRCLEPROFILEDEF
  struct ParamsGetCircleCurve {
    double radius;
    bool hasPlacement = true;
    glm::dmat3 placement;
    double thickness = -1.0f;
  };

  conway::geometry::IfcCurve getCircleCurve(
      const ParamsGetCircleCurve& parameters);

  conway::geometry::IfcCurve getCircleHoleCurve(
      const ParamsGetCircleCurve& parameters);

  // case ifc::EllipseProfileDef
  struct ParamsGetEllipseCurve {
    double radiusX;
    double radiusY;
    bool hasPlacement = true;
    glm::dmat3 placement;
    int circleSegments = 12;
  };

  conway::geometry::IfcCurve getEllipseCurve(
      const ParamsGetEllipseCurve& parameters);

  // case ifc::IFCTRIMMEDCURVE
  struct ParamsGetIfcTrimmedCurve {
    uint32_t masterRepresentation;
    uint32_t dimensions;
    bool senseAgreement;
    glm::dvec2 trim1Cartesian2D;
    glm::dvec3 trim1Cartesian3D;
    double trim1Double;
    glm::dvec2 trim2Cartesian2D;
    glm::dvec3 trim2Cartesian3D;
    double trim2Double;
    bool trimExists;
  };
  conway::geometry::IfcCurve getTrimmedCurve(
      const ParamsGetIfcTrimmedCurve& parameters);

  // case ifc::IFCCIRCLE
  struct ParamsGetIfcCircle {
    uint32_t dimensions;
    glm::dmat3 axis2Placement2D;
    glm::dmat4 axis2Placement3D;
    double radius;
    double radius2;
    ParamsGetIfcTrimmedCurve paramsGetIfcTrimmedCurve;
  };

  struct ParamsGetIfcLine {
    uint32_t dimensions;
    glm::dvec2 cartesianPoint2D;
    glm::dvec3 cartesianPoint3D;
    glm::dvec3 vectorOrientation;
    double vectorMagnitude;
    bool isEdge;
    ParamsGetIfcTrimmedCurve paramsGetIfcTrimmedCurve;
  };

  conway::geometry::IfcCurve getIfcLine(const ParamsGetIfcLine& parameters);

  conway::geometry::IfcCurve getIfcCircle(const ParamsGetIfcCircle& parameters);

  struct ParamsGetBSplineCurve {
    uint32_t dimensions;
    uint32_t degree;
    std::vector<glm::dvec2> points2;
    std::vector<glm::dvec3> points3;
    std::vector<double> knots;
    std::vector<double> weights;
    bool senseAgreement;
    bool isEdge;
  };

  conway::geometry::IfcCurve getBSplineCurve(
      const ParamsGetBSplineCurve& parameters);

  // case ifc::IFCREVOLVEDAREASOLID:
  struct ParamsGetRevolvedAreaSolid {
    glm::dmat4 placement = glm::dmat4(1.0);
    glm::dvec3 axisPosition = glm::dvec3(0, 0, 1);
    glm::dvec3 axis = glm::dvec3(0, 0, 1);
    double angle = 0;
    IfcProfile profile;
    float scalingFactor = 1.0f;
    uint32_t circleSegments = 12;
  };

  conway::geometry::Geometry getRevolvedAreaSolid(
      const ParamsGetRevolvedAreaSolid& parameters);

  // case ifc::IFCEXTRUDEDAREASOLID:
  struct ParamsGetExtrudedAreaSolid {
    double depth = 0.0f;
    glm::dvec3 dir;
    IfcProfile profile;
  };

  conway::geometry::Geometry getExtrudedAreaSolid(
      const ParamsGetExtrudedAreaSolid& parameters);

 private:
  fuzzybools::Geometry GeomToFBGeom(const Geometry& geom);
  Geometry FBGeomToGeom(const fuzzybools::Geometry& fbGeom);

  bool COORDINATE_TO_ORIGIN = false;
  bool USE_FAST_BOOLS = true;

  bool DUMP_CSG_MESHES = false;

 public:
  static constexpr int CIRCLE_SEGMENTS_LOW = 6;
  static constexpr int CIRCLE_SEGMENTS_MEDIUM = 24;
  static constexpr int CIRCLE_SEGMENTS_HIGH = 96;

  static constexpr double CIRCLE_SEGMENT_TO_DIAMETER = 0.05;
  static constexpr double CIRCLE_SEGMENT_TO_RADIUS_RATIO = 0.05 / (2.0 * M_PI);

 private:
  bool MESH_CACHE = false;
  int BOOL_ABORT_THRESHOLD = 10000;  // 10k verts
};
}  // namespace conway::geometry
