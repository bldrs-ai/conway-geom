#pragma once

#include <tinynurbs/tinynurbs.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <iostream>
#include <mapbox/earcut.hpp>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "operations/geometryutils.h"
#include "representation/IfcGeometry.h"
#include "representation/geometry.h"

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

#include <filesystem>

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
  StreamWriter(std::filesystem::path pathBase)
      : m_pathBase(std::move(pathBase)) {
    assert(m_pathBase.has_root_path());
  }

  const std::vector<std::string> getUris() const { return uris; }

  // Resolves the relative URIs of any external resources declared in the glTF
  // manifest
  std::shared_ptr<std::ostream> GetOutputStream(
      const std::string &filename) const override {
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
    glm::dmat4 transformation;
    IfcComposedMesh ifcPresentationMesh;
  };

  IfcComposedMesh getMappedItem(ParamsGetMappedItem parameters);
  IfcGeometry BoolSubtract(const std::vector<IfcGeometry> &firstGroups,
                           std::vector<IfcGeometry> &secondGroups);

  // case ifc::IFCBOOLEANCLIPPINGRESULT:
  // case ifc::IFCBOOLEANRESULT:
  struct ParamsGetBooleanResult {
    std::vector<IfcGeometry> &flatFirstMesh;
    std::vector<IfcGeometry> &flatSecondMesh;
    bool clippingResult;
  };
  IfcGeometry GetBooleanResult(ParamsGetBooleanResult parameters);

  // case ifc::IFCHALFSPACESOLID:
  struct ParamsGetHalfspaceSolid {
    IfcSurface &surface;
    bool flipWinding;
    double optionalLinearScalingFactor = 1.0;
  };

  IfcGeometry GetHalfSpaceSolid(ParamsGetHalfspaceSolid parameters);

  // case ifc::IFCPOLYGONALBOUNDEDHALFSPACE
  struct ParamsGetPolygonalBoundedHalfspace {
    IfcSurface &surface;
    glm::dmat4 position;
    conway::geometry::IfcCurve curve;
    bool halfSpaceInPlaneDirection;  // agreement != "T"
    double optionalLinearScalingFactor = 1.0;
  };

  IfcGeometry GetPolygonalBoundedHalfspace(
      ParamsGetPolygonalBoundedHalfspace parameters);

  // case ifc::IFCREPRESENTATIONMAP
  // TODO(nickcastel50) : see if this is needed

  // case ifc::IFCCONNECTEDFACESET:
  // case ifc::IFCCLOSEDSHELL:
  // case ifc::IFCOPENSHELL:
  // These cases are handled by getBrep()
  struct ParamsAddFaceToGeometry {
    uint32_t boundsSize;
    uint32_t *indices;
    uint32_t indicesPerFace;
    IfcBound3D *boundsArray;
    bool advancedBrep;
    IfcSurface *surface;
  };

  struct ParamsGetBrep {
    uint32_t boundsSize;
    uint32_t indicesPerFace;
    size_t numIndices;
    uint32_t *indices;
    IfcBound3D *boundsArray;
    bool advancedBrep;
    IfcSurface *surface;
  };

  IfcGeometry getBrep(ParamsGetBrep parameters);

  // case ifc::IFCFACE:
  // case ifc::IFCADVANCEDFACE:
  void AddFaceToGeometry(ParamsAddFaceToGeometry parameters,
                         IfcGeometry &geometry);

  // case ifc::IFCFACEBASEDSURFACEMODEL:
  // case ifc::IFCSHELLBASEDSURFACEMODEL:
  struct ParamsGetSurfaceModel {
    uint32_t numShellRefs;
    ParamsGetBrep *shells;
  };

  struct ResponseGeometryArray {
    uint32_t numRepresentations;
    IfcGeometry *geometryArray;
  };

  ResponseGeometryArray GetSurfaceModel(ParamsGetSurfaceModel parameters);

  // case ifc::IFCPLANE:
  // case ifc::IFCBSPLINESURFACE:
  // case ifc::IFCBSPLINESURFACEWITHKNOTS:
  // case ifc::IFCRATIONALBSPLINESURFACEWITHKNOTS:
  struct ParamsGetSurface {
    bool isPlane;
    glm::dmat4 transformation;
    bool isBsplineSurface;
    double Udegree;
    double Vdegree;
    // TODO(nickcastel50): How do we pass these across?
    std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;
    std::string curveType;
    std::string closedU;
    std::string closedV;
    std::string selfIntersect;
    bool isBsplineSurfaceWithKnots;
    std::vector<glm::f64> UMultiplicity;
    std::vector<glm::f64> VMultiplicity;
    std::vector<glm::f64> UKnots;
    std::vector<glm::f64> VKnots;
    bool isRationalBsplineSurfaceWithKnots;
    std::vector<std::vector<glm::f64>> weightPts;
    bool isCylindricalSurface;
    double radius;
    bool isSurfaceOfRevolution;
    glm::dmat4 revolutionDirection;
    IfcProfile3D revolutionProfile;
    bool includeTransformation;
    bool isSurfaceOfLinearExtrusion;
    glm::dvec3 extrusionDirection;
    IfcProfile extrusionProfile;
    bool customLength;
    double length;
  };
  IfcSurface GetSurface(ParamsGetSurface parameters);

  // case ifc::IFCAXIS2PLACEMENT2D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM:
  struct ParamsGetAxis2Placement2D {
    bool isAxis2Placement2D;
    bool isCartesianTransformationOperator2D;
    bool isCartesianTransformationOperator2DNonUniform;
    glm::dvec2 position2D;
    bool customAxis1Ref;
    glm::dvec2 axis1Ref;
    bool customAxis2Ref;
    glm::dvec2 axis2Ref;
    bool customScale;
    double scale1;
    bool customScale2;
    double scale2;
  };

  glm::dmat3 GetAxis2Placement2D(ParamsGetAxis2Placement2D parameters);

  // case ifc::IFCAXIS1PLACEMENT:
  struct ParamsAxis1Placement3D {
    glm::dvec3 position;
    glm::dvec3 xAxisRef;
    glm::dvec3 zAxisRef;
    bool normalizeZ;
  };
  glm::dmat4 GetAxis1Placement(ParamsAxis1Placement3D parameters);

  // case ifc::IFCAXIS2PLACEMENT3D:
  struct ParamsAxis2Placement3D {
    glm::dvec3 position;
    glm::dvec3 xAxisRef;
    glm::dvec3 zAxisRef;
    bool normalizeZ;
    bool normalizeX;
  };
  glm::dmat4 GetAxis2Placement3D(ParamsAxis2Placement3D parameters);

  // case ifc::IFCLOCALPLACEMENT:
  // This case just recursively calls GetLocalPlacement, not sure if needed. See
  // GetLocalPlacement
  struct ParamsLocalPlacement {
    bool useRelPlacement;
    glm::dmat4 axis2Placement;
    glm::dmat4 relPlacement;
  };

  glm::dmat4 GetLocalPlacement(ParamsLocalPlacement parameters);

  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3D:
  // case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM:
  struct ParamsCartesianTransformationOperator3D {
    glm::dvec3 position;
    glm::dvec3 axis1Ref;
    glm::dvec3 axis2Ref;
    glm::dvec3 axis3Ref;
    bool normalizeAxis1;
    bool normalizeAxis2;
    bool normalizeAxis3;
    bool nonUniform;
    bool realScale;
    double scale1_;
    double scale2_;
    double scale3_;
  };

  glm::dmat4 GetCartesianTransformationOperator3D(
      ParamsCartesianTransformationOperator3D parameters);

  // case ifc::IFCPOLYLOOP:
  // case ifc::IFCEDGELOOP:
  struct ParamsGetLoop {
    bool isEdgeLoop;
    size_t numPoints;
    std::vector<glm::dvec3> points;
  };

  IfcCurve GetLoop(ParamsGetLoop parameters);

  struct ParamsGetBound {
    bool isFaceOuterBound;
    bool orient;
    ParamsGetLoop parametersGetLoop;
  };
  IfcBound3D GetBound(ParamsGetBound parameters);

  // case ifc::IFCINDEXEDPOLYGONALFACEWITHVOIDS:
  // case ifc::IFCINDEXEDPOLYGONALFACE:
  struct ParamsReadIndexedPolygonalFace {
    size_t numPoints;
    size_t numIndices;
    bool indexedPolygonalFaceWithVoids;
    glm::vec3 *points;
    uint32_t *indices;
  };
  std::vector<IfcBound3D> ReadIndexedPolygonalFace(
      ParamsReadIndexedPolygonalFace parameters);

  struct ResultsGltf {
    bool success = false;
    std::vector<std::string> bufferUris;
    std::vector<std::vector<uint8_t>> buffers;
  };
  ResultsGltf GeometryToGltf(conway::geometry::IfcGeometry geom, bool isGlb,
                             bool outputDraco, std::string filePath,
                             bool outputFile,
                             glm::dmat4 transform = glm::dmat4(1));

  std::string GeometryToObj(const conway::geometry::IfcGeometry &geom,
                            size_t &offset,
                            glm::dmat4 transform = glm::dmat4(1));

  // case ifc::IFCPOLYGONALFACESET:
  struct ParamsPolygonalFaceSet {
    size_t numPoints;
    size_t numIndices;
    uint32_t indicesPerFace;
    bool indexedPolygonalFaceWithVoids;
    std::vector<glm::vec3> points;
    std::vector<uint32_t> indices;
  };
  IfcGeometry getPolygonalFaceSetGeometry(ParamsPolygonalFaceSet parameters);

 private:
  fuzzybools::Geometry GeomToFBGeom(const IfcGeometry &geom);
  IfcGeometry FBGeomToGeom(const fuzzybools::Geometry &fbGeom);

  bool COORDINATE_TO_ORIGIN = false;
  bool USE_FAST_BOOLS = true;

  bool DUMP_CSG_MESHES = false;
  int CIRCLE_SEGMENTS_LOW = 5;
  int CIRCLE_SEGMENTS_MEDIUM = 8;
  int CIRCLE_SEGMENTS_HIGH = 12;
  bool MESH_CACHE = false;
  int BOOL_ABORT_THRESHOLD = 10000;  // 10k verts
};
}  // namespace conway::geometry