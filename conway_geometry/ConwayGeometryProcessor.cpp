#include "ConwayGeometryProcessor.h"

#include <glm/glm.hpp>

#include "fuzzy/fuzzy-bools.h"
#include "operations/curve-utils.h"
#include "operations/geometryutils.h"
#include "operations/mesh_utils.h"
#include "representation/geometry.h"

namespace conway::geometry {
IfcComposedMesh ConwayGeometryProcessor::getMappedItem(
    ParamsGetMappedItem parameters) {
  IfcComposedMesh mesh;
  mesh.transformation = parameters.transformation;
  mesh.children.push_back(parameters.ifcPresentationMesh);

  return mesh;
}

fuzzybools::Geometry ConwayGeometryProcessor::GeomToFBGeom(
    const IfcGeometry &geom) {
  fuzzybools::Geometry fbGeom;

  for (size_t i = 0; i < geom.numFaces; i++) {
    const Face &f = geom.GetFace(i);

    auto a = geom.GetPoint(f.i0);
    auto b = geom.GetPoint(f.i1);
    auto c = geom.GetPoint(f.i2);

    fbGeom.AddFace(a, b, c);
  }

  return fbGeom;
}

IfcGeometry ConwayGeometryProcessor::FBGeomToGeom(
    const fuzzybools::Geometry &fbGeom) {
  IfcGeometry geom;

  for (size_t i = 0; i < fbGeom.numFaces; i++) {
    const fuzzybools::Face &f = fbGeom.GetFace(i);

    auto a = fbGeom.GetPoint(f.i0);
    auto b = fbGeom.GetPoint(f.i1);
    auto c = fbGeom.GetPoint(f.i2);

    geom.AddFace(a, b, c);
  }

  return geom;
}

IfcGeometry ConwayGeometryProcessor::BoolSubtract(
    const std::vector<IfcGeometry> &firstGroups,
    std::vector<IfcGeometry> &secondGroups) {
  std::vector<IfcGeometry> results;

  std::vector<IfcGeometry> firstGeoms;
  for (auto &firsts : firstGroups) {
    if (firsts.components.size() < 2) {
      firstGeoms.push_back(firsts);
    } else {
      firstGeoms.insert(firstGeoms.end(), firsts.components.begin(),
                        firsts.components.end());
    }
  }

  std::vector<IfcGeometry> secondGeoms;
  for (auto &seconds : secondGroups) {
    if (seconds.components.size() < 2) {
      secondGeoms.push_back(seconds);
    } else {
      secondGeoms.insert(secondGeoms.end(), seconds.components.begin(),
                         seconds.components.end());
    }
  }

  for (auto &firstGeom : firstGeoms) {
    IfcGeometry result = firstGeom;
    for (auto &secondGeom : secondGeoms) {
      bool doit = true;
      glm::dvec3 center;
      glm::dvec3 extents;
      result.GetCenterExtents(center, extents);

      glm::dvec3 s_center;
      glm::dvec3 s_extents;
      secondGeom.GetCenterExtents(s_center, s_extents);

      if (secondGeom.numFaces == 0) {
        printf("bool aborted due to empty source or target");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        doit = false;
      }

      if (result.numFaces == 0) {
        printf("bool aborted due to empty source or target");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        break;
      }

      if (doit) {
        auto fb1 = GeomToFBGeom(result);
        auto fb2 = GeomToFBGeom(secondGeom);

        result = FBGeomToGeom(fuzzybools::Subtract(fb1, fb2));
      }
    }
    results.push_back(result);
  }

  return flattenGeometry(results);
}

IfcGeometry ConwayGeometryProcessor::GetBooleanResult(
    ParamsGetBooleanResult parameters) {
  IfcGeometry resultGeometry =
      BoolSubtract(parameters.flatFirstMesh, parameters.flatSecondMesh);

  return resultGeometry;
}

IfcGeometry ConwayGeometryProcessor::GetHalfSpaceSolid(
    ParamsGetHalfspaceSolid parameters) {
  glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);

  if (parameters.flipWinding) {
    extrusionNormal *= -1;
  }

  double d =
      EXTRUSION_DISTANCE_HALFSPACE_M / parameters.optionalLinearScalingFactor;

  IfcProfile profile;
  profile.isConvex = false;
  profile.curve = GetRectangleCurve(d, d, glm::dmat3(1));

  auto geom = Extrude(profile, extrusionNormal, d);

  // @Refactor: duplicate of extrudedareasolid
  if (parameters.flipWinding) {
    for (uint32_t i = 0; i < geom.numFaces; i++) {
      uint32_t temp = geom.indexData[i * 3 + 0];
      temp = geom.indexData[i * 3 + 0];
      geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
      geom.indexData[i * 3 + 1] = temp;
    }
  }

  return geom;
}

IfcGeometry ConwayGeometryProcessor::GetPolygonalBoundedHalfspace(
    ParamsGetPolygonalBoundedHalfspace parameters) {
  if (!parameters.curve.IsCCW()) {
    parameters.curve.Invert();
  }

  glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);
  glm::dvec3 planeNormal = parameters.surface.transformation[2];
  glm::dvec3 planePosition = parameters.surface.transformation[3];

  glm::dmat4 invPosition = glm::inverse(parameters.position);
  glm::dvec3 localPlaneNormal = invPosition * glm::dvec4(planeNormal, 0);
  auto localPlanePos = invPosition * glm::dvec4(planePosition, 1);

  bool flipWinding = false;
  double extrudeDistance =
      EXTRUSION_DISTANCE_HALFSPACE_M / parameters.optionalLinearScalingFactor;
  bool halfSpaceInPlaneDirection = parameters.halfSpaceInPlaneDirection;
  bool extrudeInPlaneDirection =
      glm::dot(localPlaneNormal, extrusionNormal) > 0;
  bool ignoreDistanceInExtrude =
      (!halfSpaceInPlaneDirection && extrudeInPlaneDirection) ||
      (halfSpaceInPlaneDirection && !extrudeInPlaneDirection);

  if (ignoreDistanceInExtrude) {
    // spec says this should be * 0, but that causes issues for degenerate 0
    // volume pbhs hopefully we can get away by just inverting it
    extrudeDistance *= -1;
    flipWinding = true;
  }

  IfcProfile profile;
  profile.isConvex = false;
  profile.curve = parameters.curve;

  conway::geometry::IfcGeometry geom =
      Extrude(profile, extrusionNormal, extrudeDistance, localPlaneNormal,
              localPlanePos);
  // auto geom = Extrude(profile, surface.transformation, extrusionNormal,
  // EXTRUSION_DISTANCE_HALFSPACE);

  // @Refactor: duplicate of extrudedareasolid
  if (flipWinding) {
    for (uint32_t i = 0; i < geom.numFaces; i++) {
      uint32_t temp = geom.indexData[i * 3 + 0];
      temp = geom.indexData[i * 3 + 0];
      geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
      geom.indexData[i * 3 + 1] = temp;
    }
  }

#ifdef DUMP_CSG_MESHES
  DumpIfcGeometry(geom, L"pbhs.obj");
#endif

  return geom;
}

IfcGeometry ConwayGeometryProcessor::getBrep(ParamsGetBrep parameters) {
  IfcGeometry geometry;
  // set parameters
  ParamsAddFaceToGeometry paramsAddFaceToGeometry;
  paramsAddFaceToGeometry.boundsSize = parameters.boundsSize;
  paramsAddFaceToGeometry.indicesPerFace = parameters.indicesPerFace;
  paramsAddFaceToGeometry.boundsArray = parameters.boundsArray;
  paramsAddFaceToGeometry.advancedBrep = parameters.advancedBrep;
  paramsAddFaceToGeometry.surface = parameters.surface;

  for (size_t faceIndex = 0;
       faceIndex < parameters.numIndices / parameters.indicesPerFace;
       faceIndex++) {
    paramsAddFaceToGeometry.indices =
        &parameters.indices[faceIndex * parameters.indicesPerFace];

    AddFaceToGeometry(paramsAddFaceToGeometry, geometry);
  }

  return geometry;
}

void ConwayGeometryProcessor::AddFaceToGeometry(
    ParamsAddFaceToGeometry parameters, IfcGeometry &geometry) {
  if (!parameters.advancedBrep) {
    std::vector<IfcBound3D> bounds3D(parameters.boundsSize);

    for (int i = 0; i < parameters.boundsSize; i++) {
      bounds3D[i] = parameters.boundsArray[i];
    }

    TriangulateBounds(geometry, bounds3D);
  } else {
    std::vector<IfcBound3D> bounds3D(parameters.boundsSize);

    for (int i = 0; i < parameters.boundsSize; i++) {
      bounds3D[i] = parameters.boundsArray[i];
    }

    // auto surface = GetSurface(surfRef);

    if (parameters.surface == nullptr) {
      printf("surface was nullptr\n");
      return;
    }

    auto surface = parameters.surface[0];

    if (surface.BSplineSurface.Active) {
      TriangulateBspline(geometry, bounds3D, surface);
    } else if (surface.CylinderSurface.Active) {
      TriangulateCylindricalSurface(geometry, bounds3D, surface);
    } else if (surface.RevolutionSurface.Active) {
      TriangulateRevolution(geometry, bounds3D, surface);
    } else if (surface.ExtrusionSurface.Active) {
      TriangulateExtrusion(geometry, bounds3D, surface);
    } else {
      TriangulateBounds(geometry, bounds3D);
    }
  }
}

conway::geometry::ConwayGeometryProcessor::ResponseGeometryArray
ConwayGeometryProcessor::GetSurfaceModel(ParamsGetSurfaceModel parameters) {
  ResponseGeometryArray response;
  std::vector<IfcGeometry> geometryArray;

  for (uint32_t shellIndex = 0; shellIndex < parameters.numShellRefs;
       shellIndex++) {
    geometryArray.push_back(getBrep(parameters.shells[shellIndex]));
  }

  response.numRepresentations = parameters.numShellRefs;
  response.geometryArray = geometryArray.data();

  return response;
}

IfcSurface ConwayGeometryProcessor::GetSurface(ParamsGetSurface parameters) {
  if (parameters.isPlane) {
    IfcSurface surface;

    surface.transformation = parameters.transformation;

    return surface;
  } else if (parameters.isBsplineSurface) {
    IfcSurface surface;

    surface.BSplineSurface.Active = true;
    surface.BSplineSurface.UDegree = parameters.Udegree;
    surface.BSplineSurface.VDegree = parameters.Vdegree;
    surface.BSplineSurface.ControlPoints = parameters.ctrolPts;
    surface.BSplineSurface.ClosedU = parameters.closedU;
    surface.BSplineSurface.ClosedV = parameters.closedV;
    surface.BSplineSurface.CurveType = parameters.curveType;

    // TODO(nickcastel50): Old implementation wasn't returning a surface for this case.
    return surface;
  } else if (parameters.isBsplineSurfaceWithKnots ||
             parameters.isRationalBsplineSurfaceWithKnots) {
    IfcSurface surface;

    if (parameters.UKnots[parameters.UKnots.size() - 1] !=
        (int)parameters.UKnots[parameters.UKnots.size() - 1]) {
      for (uint32_t i = 0; i < parameters.UKnots.size(); i++) {
        parameters.UKnots[i] = parameters.UKnots[i] *
                               (parameters.UKnots.size() - 1) /
                               parameters.UKnots[parameters.UKnots.size() - 1];
      }
    }

    if (parameters.VKnots[parameters.VKnots.size() - 1] !=
        (int)parameters.VKnots[parameters.VKnots.size() - 1]) {
      for (uint32_t i = 0; i < parameters.VKnots.size(); i++) {
        parameters.VKnots[i] = parameters.VKnots[i] *
                               (parameters.VKnots.size() - 1) /
                               parameters.VKnots[parameters.VKnots.size() - 1];
      }
    }

    surface.BSplineSurface.Active = true;
    surface.BSplineSurface.UDegree = parameters.Udegree;
    surface.BSplineSurface.VDegree = parameters.Vdegree;
    surface.BSplineSurface.ControlPoints = parameters.ctrolPts;
    surface.BSplineSurface.UMultiplicity = parameters.UMultiplicity;
    surface.BSplineSurface.VMultiplicity = parameters.VMultiplicity;
    surface.BSplineSurface.UKnots = parameters.UKnots;
    surface.BSplineSurface.VKnots = parameters.VKnots;

    if (parameters.isRationalBsplineSurfaceWithKnots) {
      surface.BSplineSurface.WeightPoints = parameters.weightPts;
    }

    return surface;
  } else if (parameters.isCylindricalSurface) {
    IfcSurface surface;

    surface.transformation = parameters.transformation;
    surface.CylinderSurface.Active = true;
    surface.CylinderSurface.Radius = parameters.radius;

    return surface;
  } else if (parameters.isSurfaceOfRevolution) {
    IfcSurface surface;

    if (parameters.includeTransformation) {
      surface.transformation = parameters.transformation;
    }

    surface.RevolutionSurface.Active = true;
    surface.RevolutionSurface.Direction = parameters.revolutionDirection;
    surface.RevolutionSurface.Profile = parameters.revolutionProfile;

    return surface;
  } else if (parameters.isSurfaceOfLinearExtrusion) {
    IfcSurface surface;

    glm::dvec3 direction = parameters.extrusionDirection;

    double length = 0;
    if (parameters.customLength) {
      length = parameters.length;
    }

    surface.ExtrusionSurface.Active = true;
    surface.ExtrusionSurface.Length = length;
    surface.ExtrusionSurface.Profile = parameters.extrusionProfile;
    surface.ExtrusionSurface.Direction = direction;
    surface.transformation = parameters.transformation;

    return surface;
  }

  return IfcSurface();
}

glm::dmat3 ConwayGeometryProcessor::GetAxis2Placement2D(
    ParamsGetAxis2Placement2D parameters) {
  if (parameters.isAxis2Placement2D) {
    glm::dvec2 xAxis = glm::dvec2(1, 0);
    if (parameters.customAxis1Ref) {
      xAxis = glm::normalize(parameters.axis1Ref);
    }

    glm::dvec2 pos = parameters.position2D;

    glm::dvec2 yAxis = glm::normalize(glm::dvec2(xAxis.y, -xAxis.x));

    return glm::dmat3(glm::dvec3(xAxis, 0), glm::dvec3(yAxis, 0),
                      glm::dvec3(pos, 1));

  } else if (parameters.isCartesianTransformationOperator2D ||
             parameters.isCartesianTransformationOperator2DNonUniform) {
    double scale1 = 1.0;
    double scale2 = 1.0;

    glm::dvec2 Axis1(1, 0);
    glm::dvec2 Axis2(0, 1);

    if (parameters.customAxis1Ref) {
      Axis1 = glm::normalize(parameters.axis1Ref);
    }

    if (parameters.customAxis2Ref) {
      Axis2 = glm::normalize(parameters.axis2Ref);
    }

    glm::dvec2 pos = parameters.position2D;

    if (parameters.customScale) {
      scale1 = parameters.scale1;
    }

    if (parameters.isCartesianTransformationOperator2DNonUniform) {
      if (parameters.customScale2) {
        scale2 = parameters.scale2;
      }
    }

    if (parameters.isCartesianTransformationOperator2D) {
      scale2 = scale1;
    }

    return glm::dmat3(glm::dvec3(Axis1 * scale1, 0),
                      glm::dvec3(Axis2 * scale2, 0), glm::dvec3(pos, 1));
  }

  return glm::dmat3();
}

glm::dmat4 ConwayGeometryProcessor::GetAxis1Placement(
    ParamsAxis1Placement3D parameters) {
  glm::dvec3 zAxis(0, 0, 1);
  glm::dvec3 xAxis(1, 0, 0);

  if (parameters.normalizeZ) {
    zAxis = glm::normalize(parameters.zAxisRef);
  }

  glm::dvec3 pos = parameters.position;
  if (std::abs(glm::dot(xAxis, zAxis)) > 0.9) {
    xAxis = glm::dvec3(0, 1, 0);
  }

  glm::dvec3 yAxis = glm::normalize(glm::cross(zAxis, xAxis));
  xAxis = glm::normalize(glm::cross(zAxis, yAxis));

  glm::dmat4 result = glm::dmat4(glm::dvec4(xAxis, 0), glm::dvec4(yAxis, 0),
                                 glm::dvec4(zAxis, 0), glm::dvec4(pos, 1));

  return result;
}

glm::dmat4 ConwayGeometryProcessor::GetAxis2Placement3D(
    ParamsAxis2Placement3D parameters) {
  glm::dvec3 zAxis(0, 0, 1);
  glm::dvec3 xAxis(1, 0, 0);

  if (parameters.normalizeZ) {
    zAxis = glm::normalize(parameters.zAxisRef);
  }

  if (parameters.normalizeX) {
    xAxis = glm::normalize(parameters.xAxisRef);
  }

  glm::dvec3 pos = parameters.position;

  glm::dvec3 yAxis = glm::normalize(glm::cross(zAxis, xAxis));
  xAxis = glm::normalize(glm::cross(yAxis, zAxis));

  return glm::dmat4(glm::dvec4(xAxis, 0), glm::dvec4(yAxis, 0),
                    glm::dvec4(zAxis, 0), glm::dvec4(pos, 1));
}

glm::dmat4 ConwayGeometryProcessor::GetLocalPlacement(
    ParamsLocalPlacement parameters) {
  if (parameters.useRelPlacement) {
    glm::dmat4 result = parameters.relPlacement * parameters.axis2Placement;
    return result;
  } else {
    glm::dmat4 relPlacement(1);
    glm::dmat4 result = relPlacement * parameters.axis2Placement;
    return result;
  }
}

glm::dmat4 ConwayGeometryProcessor::GetCartesianTransformationOperator3D(
    ParamsCartesianTransformationOperator3D parameters) {
  double scale1 = 1.0;
  double scale2 = 1.0;
  double scale3 = 1.0;

  glm::dvec3 Axis1(1, 0, 0);
  glm::dvec3 Axis2(0, 1, 0);
  glm::dvec3 Axis3(0, 0, 1);

  if (parameters.normalizeAxis1) Axis1 = glm::normalize(parameters.axis1Ref);

  if (parameters.normalizeAxis2) Axis2 = glm::normalize(parameters.axis2Ref);

  glm::dvec3 pos = parameters.position;

  if (parameters.realScale) scale1 = parameters.scale1_;

  if (parameters.normalizeAxis3) Axis3 = glm::normalize(parameters.axis3Ref);

  if (parameters.nonUniform) {
    if (parameters.realScale) scale2 = parameters.scale2_;

    if (parameters.realScale) scale3 = parameters.scale3_;
  } else {
    scale2 = scale1;
    scale3 = scale1;
  }

  return glm::dmat4(glm::dvec4(Axis1 * scale1, 0),
                    glm::dvec4(Axis2 * scale2, 0),
                    glm::dvec4(Axis3 * scale3, 0), glm::dvec4(pos, 1));
}

IfcCurve ConwayGeometryProcessor::GetLoop(ParamsGetLoop parameters) {
  IfcCurve curve;
  if (!parameters.isEdgeLoop) {
    if (parameters.numPoints > 0) {
      curve.points.reserve(parameters.numPoints);

      glm::dvec3 prevPoint = parameters.points[0];

      curve.points.push_back(prevPoint);

      for (size_t index = 1; index < parameters.numPoints; index++) {
        glm::dvec3 currentPoint = parameters.points[index];
        // trim repeats
        if (currentPoint.x != prevPoint.x && currentPoint.y != prevPoint.y &&
            currentPoint.z != prevPoint.z) {
          curve.points.push_back(parameters.points[index]);
        }

        prevPoint = currentPoint;
      }
    }
  } else {
    // TODO(nickcastel50): Handle edge loop
    ;
    /*auto edges = _loader.GetSetArgument();
    int id = 0;

    for (auto &token : edges)
    {
            uint32_t edgeId = _loader.GetRefArgument(token);
            IfcCurve<3> edgeCurve = GetOrientedEdge(edgeId);

            // Important not to repeat the last point otherwise triangulation
    fails
            // if the list has zero points this is initial, no repetition is
    possible, otherwise we must check if (curve.points.size() == 0)
            {
                    for (auto &pt : edgeCurve.points)
                    {
                            curve.points.push_back(pt);
                            curve.indices.push_back(id);
                    }
            }
            else
            {
                    for (auto &pt : edgeCurve.points)
                    {
                            if (notPresent(pt, curve.points))
                            {
                                    curve.points.push_back(pt);
                                    curve.indices.push_back(id);
                            }
                    }
            }
            id++;
    }*/
  }

  return curve;
}

IfcBound3D ConwayGeometryProcessor::GetBound(ParamsGetBound parameters) {
  bool orient = parameters.orient;

  IfcBound3D bound;
  bound.curve = GetLoop(parameters.parametersGetLoop);
  bound.orientation = orient;
  bound.type = (parameters.isFaceOuterBound) ? IfcBoundType::OUTERBOUND
                                             : IfcBoundType::BOUND;

  if (!orient)
    std::reverse(bound.curve.points.begin(), bound.curve.points.end());

  return bound;
}

std::vector<IfcBound3D> ConwayGeometryProcessor::ReadIndexedPolygonalFace(
    ParamsReadIndexedPolygonalFace parameters) {
  std::vector<IfcBound3D> bounds;
  bounds.emplace_back();

  for (size_t index = 0; index < parameters.numIndices; index++) {
    uint32_t currentIndex = parameters.indices[index];

    glm::dvec3 point = parameters.points[currentIndex - 1];

    // I am not proud of this (I inherited this, will change - NC)
    bounds.back().curve.points.push_back(point);
  }

  if (!parameters.indexedPolygonalFaceWithVoids)
    return bounds;
  else
    // TODO(nickcastel50): handle case IFCINDEXEDPOLYGONALFACEWITHVOIDS
    ;

  return bounds;
}

conway::geometry::ConwayGeometryProcessor::ResultsGltf
ConwayGeometryProcessor::GeometryToGltf(conway::geometry::IfcGeometry geom,
                                        bool isGlb, bool outputDraco,
                                        std::string filePath, bool outputFile,
                                        glm::dmat4 transform) {
  ResultsGltf results;
  // Encode the geometry.
  draco::EncoderBuffer buffer;

  // create a mesh object
  std::unique_ptr<draco::Mesh> dracoMesh(new draco::Mesh());

  try {
    int32_t pos_att_id = -1;
    int32_t tex_att_id = -1;
    int32_t material_att_id = -1;

    // this internally populates the vertex float array, current storage type
    // is double
    geom.GetVertexData();

    if (outputDraco) {
      DracoOptions dracoOptions;

      // get number of faces
      uint32_t numFaces = 0;
      numFaces = geom.GetIndexDataSize() / 3;

      if (numFaces > 0) {
        // set number of faces
        dracoMesh->SetNumFaces(numFaces);

        // set number of indices
        dracoMesh->set_num_points(3 * numFaces);
      }

      uint32_t numPositions = 0;
      numPositions = geom.numPoints;

      // Add attributes if they are present in the input data.
      if (numPositions > 0) {
        draco::GeometryAttribute va;
        va.Init(draco::GeometryAttribute::POSITION, nullptr, 3,
                draco::DT_FLOAT32, false, sizeof(float) * 3, 0);
        pos_att_id = dracoMesh->AddAttribute(va, false, numPositions);
      }

      // TODO(nickcastel50): support multiple materials at some point when we add that
      int32_t numMaterials = 1;
      int32_t numTexCoords = 0;

      // populate position attribute
      float vertexVal[3];
      int32_t vertexCount = 0;

      for (uint32_t i = 0; i < geom.numPoints; i++) {
        glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
        vertexVal[0] = (float)t.x;
        vertexVal[1] = (float)t.y;
        vertexVal[2] = (float)t.z;
        dracoMesh->attribute(pos_att_id)
            ->SetAttributeValue(draco::AttributeValueIndex(i), vertexVal);
      }

      // no textures, just map vertices to face indices
      for (size_t triangleIndex = 0; triangleIndex < numFaces;
           ++triangleIndex) {
        const uint32_t *triangle = &(geom.indexData.data()[triangleIndex * 3]);

        uint32_t triangle0 = triangle[0];
        uint32_t triangle1 = triangle[1];
        uint32_t triangle2 = triangle[2];

        const draco::PointIndex vert_id_0(3 * triangleIndex);
        const draco::PointIndex vert_id_1(3 * triangleIndex + 1);
        const draco::PointIndex vert_id_2(3 * triangleIndex + 2);
        const int triangulated_index = 0;

        // map vertex to face index
        dracoMesh->attribute(pos_att_id)
            ->SetPointMapEntry(vert_id_0,
                               draco::AttributeValueIndex(triangle0));
        dracoMesh->attribute(pos_att_id)
            ->SetPointMapEntry(vert_id_1,
                               draco::AttributeValueIndex(triangle1));
        dracoMesh->attribute(pos_att_id)
            ->SetPointMapEntry(vert_id_2,
                               draco::AttributeValueIndex(triangle2));
      }

      // Add faces with identity mapping between vertex and corner indices.
      // Duplicate vertices will get removed below.
      draco::Mesh::Face face;
      for (draco::FaceIndex i(0); i < numFaces; ++i) {
        for (int c = 0; c < 3; ++c) {
          face[c] = 3 * i.value() + c;
        }
        dracoMesh->SetFace(i, face);
      }

#ifdef DRACO_ATTRIBUTE_VALUES_DEDUPLICATION_SUPPORTED
      if (dracoOptions.deduplicateInputValues) {
        dracoMesh->DeduplicateAttributeValues();
      }
#endif

#ifdef DRACO_ATTRIBUTE_INDICES_DEDUPLICATION_SUPPORTED
      dracoMesh->DeduplicatePointIds();
#endif

      // Convert compression level to speed (that 0 = slowest, 10 = fastest).
      const int32_t dracoCompressionSpeed = 10 - dracoOptions.compressionLevel;

      // set up Draco Encoder
      draco::Encoder encoder;

      // Setup encoder options.
      if (dracoOptions.posQuantizationBits > 0) {
        encoder.SetAttributeQuantization(
            draco::GeometryAttribute::POSITION,
            11);  // dracoOptions.posQuantizationBits );
      }
      if (dracoOptions.texCoordsQuantizationBits > 0) {
        encoder.SetAttributeQuantization(
            draco::GeometryAttribute::TEX_COORD,
            8);  // dracoOptions.texCoordsQuantizationBits );
      }
      if (dracoOptions.normalsQuantizationBits > 0) {
        encoder.SetAttributeQuantization(
            draco::GeometryAttribute::NORMAL,
            8);  // dracoOptions.normalsQuantizationBits );
      }
      if (dracoOptions.genericQuantizationBits > 0) {
        encoder.SetAttributeQuantization(
            draco::GeometryAttribute::GENERIC,
            8);  // dracoOptions.genericQuantizationBits );
      }

      encoder.SetSpeedOptions(dracoCompressionSpeed, dracoCompressionSpeed);

      // Convert to ExpertEncoder that allows us to set per-attribute options.
      std::unique_ptr<draco::ExpertEncoder> expert_encoder;
      expert_encoder.reset(new draco::ExpertEncoder(*dracoMesh));
      expert_encoder->Reset(encoder.CreateExpertEncoderOptions(*dracoMesh));

      // set up timer
      draco::CycleTimer timer;

      timer.Start();

      const draco::Status status = expert_encoder->EncodeToBuffer(&buffer);

      timer.Stop();

      if (!status.ok()) {
        printf("Failed to encode the mesh: %s", status.error_msg());
        results.success = false;
        return results;
      }

      if (outputFile) {
        printf("Encoded To Draco in %lld ms\n", timer.GetInMs());
      }
    }

    // The Document instance represents the glTF JSON manifest
    Microsoft::glTF::Document document;

    if (outputDraco) {
      document.extensionsRequired.insert("KHR_draco_mesh_compression");
      document.extensionsUsed.insert("KHR_draco_mesh_compression");
    }

    std::string base_filename =
        filePath.substr(filePath.find_last_of("/\\") + 1);

    // Create a Buffer - it will be the 'current' Buffer that all the
    // BufferViews created by this BufferBuilder will automatically reference

    // Pass the absolute path, without the filename, to the stream writer
    std::unique_ptr<StreamWriter> streamWriter = std::make_unique<StreamWriter>(
        std::filesystem::path(filePath).parent_path().c_str());
    std::unique_ptr<Microsoft::glTF::ResourceWriter> resourceWriter;

    if (!isGlb) {
      resourceWriter = std::make_unique<Microsoft::glTF::GLTFResourceWriter>(
          std::move(streamWriter));
    } else {
      resourceWriter = std::make_unique<Microsoft::glTF::GLBResourceWriter>(
          std::move(streamWriter));
    }

    std::string accessorIdIndices;
    std::string accessorIdPositions;
    std::string accessorIdUVs;

    // Use the BufferBuilder helper class to simplify the process of
    // constructing valid glTF Buffer, BufferView and Accessor entities
    Microsoft::glTF::BufferBuilder bufferBuilder(std::move(resourceWriter));

    // if Draco
    Microsoft::glTF::KHR::MeshPrimitives::DracoMeshCompression
        dracoMeshCompression;

    // no Draco
    //  Create all the resource data (e.g. triangle indices and
    //  vertex positions) that will be written to the binary buffer
    const char *bufferId = nullptr;

    // Specify the 'special' GLB buffer ID. This informs the GLBResourceWriter
    // that it should use the GLB container's binary chunk (usually the
    // desired buffer location when creating GLBs)
    if (isGlb) {
      bufferId = Microsoft::glTF::GLB_BUFFER_ID;
    } else {
      bufferId = base_filename.c_str();

      std::string bufferIdCopy = std::string(bufferId);
      bufferIdCopy += ".";
      bufferIdCopy += Microsoft::glTF::BUFFER_EXTENSION;
      results.bufferUris.push_back(bufferIdCopy);
    }

    // Create a Buffer - it will be the 'current' Buffer that all the
    // BufferViews created by this BufferBuilder will automatically reference
    bufferBuilder.AddBuffer(bufferId);

    // Add an Accessor for the indices and positions
    std::vector<float> positions;
    positions.resize(geom.numPoints * 3);

    std::vector<float> minValues(3U, std::numeric_limits<float>::max());
    std::vector<float> maxValues(3U, std::numeric_limits<float>::lowest());

    const size_t positionCount = positions.size();

    for (uint32_t i = 0; i < geom.numPoints; i++) {
      glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
      positions[3 * i + 0] = (float)t.x;
      positions[3 * i + 1] = (float)t.y;
      positions[3 * i + 2] = (float)t.z;
    }

    // Accessor min/max properties must be set for vertex position data so
    // calculate them here
    for (size_t i = 0U, j = 0U; i < positionCount; ++i, j = (i % 3U)) {
      minValues[j] = std::min(positions[i], minValues[j]);
      maxValues[j] = std::max(positions[i], maxValues[j]);
    }

    if (outputDraco) {
      Microsoft::glTF::Accessor positionsAccessor;
      positionsAccessor.min = minValues;
      positionsAccessor.max = maxValues;
      positionsAccessor.count = dracoMesh->num_points();
      positionsAccessor.componentType =
          Microsoft::glTF::ComponentType::COMPONENT_FLOAT;
      positionsAccessor.type = Microsoft::glTF::AccessorType::TYPE_VEC3;
      positionsAccessor.id = "0";
      accessorIdPositions = positionsAccessor.id;
      dracoMeshCompression.attributes.emplace(
          "POSITION", dracoMesh->attribute(pos_att_id)->unique_id());

      Microsoft::glTF::Accessor indicesAccessor;
      indicesAccessor.count = dracoMesh->num_faces() * 3;
      indicesAccessor.componentType = Microsoft::glTF::COMPONENT_UNSIGNED_INT;
      indicesAccessor.type = Microsoft::glTF::TYPE_SCALAR;
      indicesAccessor.id = "1";
      accessorIdIndices = indicesAccessor.id;

      // add position accessor first here
      document.accessors.Append(
          positionsAccessor, Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty);
      document.accessors.Append(
          indicesAccessor, Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty);
    } else {
      // Create a BufferView with a target of ELEMENT_ARRAY_BUFFER (as it will
      // reference index data) - it will be the 'current' BufferView that all
      // the Accessors created by this BufferBuilder will automatically
      // reference
      bufferBuilder.AddBufferView(
          Microsoft::glTF::BufferViewTarget::ELEMENT_ARRAY_BUFFER);

      // Copy the Accessor's id - subsequent calls to AddAccessor may
      // invalidate the returned reference
      accessorIdIndices =
          bufferBuilder
              .AddAccessor(geom.indexData,
                           {Microsoft::glTF::TYPE_SCALAR,
                            Microsoft::glTF::COMPONENT_UNSIGNED_INT})
              .id;

      // Create a BufferView with target ARRAY_BUFFER (as it will reference
      // vertex attribute data)
      bufferBuilder.AddBufferView(
          Microsoft::glTF::BufferViewTarget::ARRAY_BUFFER);

      // Add positions accessor
      accessorIdPositions =
          bufferBuilder
              .AddAccessor(
                  positions,
                  {Microsoft::glTF::TYPE_VEC3, Microsoft::glTF::COMPONENT_FLOAT,
                   false, std::move(minValues), std::move(maxValues)})
              .id;
    }

    // Construct a Material
    Microsoft::glTF::Material material;
    material.doubleSided = true;

    Microsoft::glTF::TextureInfo textureInfo;

    // Add it to the Document and store the generated ID
    auto materialId =
        document.materials
            .Append(std::move(material),
                    Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
            .id;

    // Construct a MeshPrimitive. Unlike most types in glTF, MeshPrimitives
    // are direct children of their parent Mesh entity rather than being
    // children of the Document. This is why they don't have an ID member.
    Microsoft::glTF::MeshPrimitive meshPrimitive;
    meshPrimitive.materialId = materialId;

    if (outputDraco) {
      // Finally put the encoded data in place.
      auto bufferView =
          bufferBuilder.AddBufferView(buffer.data(), buffer.size());
      dracoMeshCompression.bufferViewId = bufferView.id;

      // Add all of the Buffers, BufferViews and Accessors that were created
      // using BufferBuilder to the Document. Note that after this point, no
      // further calls should be made to BufferBuilder
      bufferBuilder.Output(document);

      const auto extensionSerializer =
          Microsoft::glTF::KHR::GetKHRExtensionSerializer();
      std::string strDracoMesh =
          Microsoft::glTF::KHR::MeshPrimitives::SerializeDracoMeshCompression(
              dracoMeshCompression, document, extensionSerializer);
      meshPrimitive.extensions["KHR_draco_mesh_compression"] = strDracoMesh;
      meshPrimitive.indicesAccessorId = accessorIdIndices;
      meshPrimitive.attributes[Microsoft::glTF::ACCESSOR_POSITION] =
          accessorIdPositions;

    } else {
      meshPrimitive.indicesAccessorId = accessorIdIndices;
      meshPrimitive.attributes[Microsoft::glTF::ACCESSOR_POSITION] =
          accessorIdPositions;

      // Add all of the Buffers, BufferViews and Accessors that were created
      // using BufferBuilder to the Document. Note that after this point, no
      // further calls should be made to BufferBuilder
      bufferBuilder.Output(document);
    }

    // Construct a Mesh and add the MeshPrimitive as a child
    Microsoft::glTF::Mesh mesh;
    mesh.primitives.push_back(std::move(meshPrimitive));
    // Add it to the Document and store the generated ID
    auto meshId = document.meshes
                      .Append(std::move(mesh),
                              Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
                      .id;

    // Construct a Node adding a reference to the Mesh
    Microsoft::glTF::Node node;
    node.meshId = meshId;

    // Add it to the Document and store the generated ID
    auto nodeId = document.nodes
                      .Append(std::move(node),
                              Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
                      .id;

    // Construct a Scene
    Microsoft::glTF::Scene scene;
    scene.nodes.push_back(nodeId);
    // Add it to the Document, using a utility method that also sets the Scene
    // as the Document's default
    document.SetDefaultScene(std::move(scene),
                             Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty);

    std::string manifest;

    try {
      // Serialize the glTF Document into a JSON manifest
      manifest = Serialize(document, Microsoft::glTF::SerializeFlags::Pretty);
    } catch (const Microsoft::glTF::GLTFException &ex) {
      std::stringstream ss;
      printf("Microsoft::glTF::Serialize failed: %s", ex.what());
    }

    auto &genericResourceWriter = bufferBuilder.GetResourceWriter();

    if (isGlb) {
      auto glbResourceWriter =
          static_cast<Microsoft::glTF::GLBResourceWriter *>(
              &genericResourceWriter);

      filePath += ".glb";

      glbResourceWriter->Flush(
          manifest,
          filePath);  // A GLB container isn't created until the
                      // GLBResourceWriter::Flush member function is called

      // get GLB string from ostringstream
      std::shared_ptr<std::ostream> stream_ =
          glbResourceWriter->GetOutputBuffer(filePath);
      std::ostream *streamPtr_ = stream_.get();
      std::ostringstream *ossPtr_ = (std::ostringstream *)(streamPtr_);

      std::string str = ossPtr_->str();
      std::string fileNameOnly =
          std::filesystem::path(filePath).filename().string();
      results.bufferUris.push_back(fileNameOnly);

      std::vector<uint8_t> glbBuffer(str.length());
      memcpy(glbBuffer.data(), str.c_str(), str.length());
      results.buffers.push_back(glbBuffer);

      // write file
      if (outputFile) {
        std::ofstream outfile(filePath);  // Open the file for writing

        if (outfile.is_open()) {
          // Write the data from the stream to the file
          outfile << str;

          // Check if the file was written to successfully
          if (outfile.fail()) {
            throw std::runtime_error("Failed to write to file " + filePath);
          } else {
            std::cout << "Data written to file " << filePath << " successfully."
                      << std::endl;
          }

          // Close the file
          outfile.close();
        } else {
          printf("Unable to open file %s", filePath.c_str());
        }
      }
    } else {
      filePath += ".gltf";

      // with the new StreamWriter, we must grab each output buffer
      // (ostringstream) and write it to the disk
      for (int bufferUrlIndex = 0; bufferUrlIndex < results.bufferUris.size();
           bufferUrlIndex++) {
        std::string bufferUri_ = results.bufferUris[bufferUrlIndex];

        std::shared_ptr<std::ostream> stream_ =
            genericResourceWriter.GetOutputBuffer(bufferUri_);
        std::ostream *streamPtr_ = stream_.get();
        std::ostringstream *ossPtr_ = (std::ostringstream *)(streamPtr_);

        std::string str = ossPtr_->str();
        std::vector<uint8_t> gltfBuffer(str.length());
        memcpy(gltfBuffer.data(), str.c_str(), str.length());

        results.buffers.push_back(gltfBuffer);

        if (outputFile) {
          std::ofstream outfile(bufferUri_);  // Open the file for writing

          if (outfile.is_open()) {
            printf("Writing: %s\n", bufferUri_.c_str());

            // Write the data from the stream to the file
            outfile << str;

            // Check if the file was written to successfully
            if (outfile.fail()) {
              throw std::runtime_error("Failed to write to file " + bufferUri_);
            } else {
              std::cout << "Data written to file " << bufferUri_
                        << " successfully." << std::endl;
            }

            // Close the file
            outfile.close();
          } else {
            throw std::runtime_error("Unable to open file " + filePath);
          }
        }
      }

      // push back the name of the manifest (.gltf)
      std::string fileNameOnly =
          std::filesystem::path(filePath).filename().string();
      results.bufferUris.push_back(fileNameOnly);

      std::vector<uint8_t> gltfManifest(manifest.length());
      memcpy(gltfManifest.data(), manifest.c_str(), manifest.length());
      results.buffers.push_back(gltfManifest);

      // write file
      if (outputFile) {
        printf("Writing: %s\n", filePath.c_str());

        std::ofstream outfile(filePath);  // Open the file for writing

        if (outfile.is_open()) {
          // Write the data from the stream to the file
          outfile << manifest.c_str();

          // Check if the file was written to successfully
          if (outfile.fail()) {
            throw std::runtime_error("Failed to write to file " + filePath);
          } else {
            std::cout << "Data written to file " << filePath << " successfully."
                      << std::endl;
          }

          // Close the file
          outfile.close();
        } else {
          printf("Unable to open file %s", filePath.c_str());
        }
      }
    }
    results.success = true;
    return results;
  } catch (const std::exception &ex) {
    printf("Couldn't write GLB file: %s", ex.what());
    results.success = false;
    return results;
  }
}

std::string ConwayGeometryProcessor::GeometryToObj(
    const conway::geometry::IfcGeometry &geom, size_t &offset,
    glm::dmat4 transform) {
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

IfcGeometry ConwayGeometryProcessor::getPolygonalFaceSetGeometry(
    ParamsPolygonalFaceSet parameters) {
  IfcGeometry geom;
  std::vector<IfcBound3D> bounds;

  ParamsReadIndexedPolygonalFace readIndexedPolygonalFaceParameters;
  readIndexedPolygonalFaceParameters.numIndices = parameters.indicesPerFace;
  readIndexedPolygonalFaceParameters.indexedPolygonalFaceWithVoids =
      parameters.indexedPolygonalFaceWithVoids;
  readIndexedPolygonalFaceParameters.points = parameters.points.data();

  for (size_t faceIndex = 0;
       faceIndex < parameters.numIndices / parameters.indicesPerFace;
       faceIndex++) {
    readIndexedPolygonalFaceParameters.indices =
        &parameters.indices[faceIndex * parameters.indicesPerFace];
    bounds = ReadIndexedPolygonalFace(readIndexedPolygonalFaceParameters);

    TriangulateBounds(geom, bounds);

    bounds.clear();
  }

  return geom;
}

}  // namespace conway::geometry