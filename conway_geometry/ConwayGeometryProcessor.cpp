#include "ConwayGeometryProcessor.h"

#include <glm/glm.hpp>

#include "fuzzy/fuzzy-bools.h"
// #include "legacy/math/bool-mesh-mesh.h"
// #include "legacy/math/intersect-mesh-mesh.h"
#include "operations/curve-utils.h"
#include "operations/geometryutils.h"
#include "operations/mesh_utils.h"
#include "representation/geometry.h"
#include "../logging/Logger.h"

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
    const fuzzybools::Face &f = geom.GetFace(i);

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

double normalDiff(glm::dvec3 extents) {
  double a = extents.x;

  if (a < extents.y) {
    a = extents.y;
  }
  if (a < extents.z) {
    a = extents.z;
  }

  double b = extents.x;

  if (b > extents.y) {
    b = extents.y;
  }
  if (b > extents.z) {
    b = extents.z;
  }

  if (a > 0) {
    return b / a;
  } else {
    return 0;
  }
}

/*IfcGeometry BoolJoinLegacy(const std::vector<IfcGeometry> &Geoms) {
  IfcGeometry result;

  if (Geoms.size() == 0) {
    return result;
  } else if (Geoms.size() == 1) {
    return Geoms[0];
  } else {
    bool first = true;
    for (auto &geom : Geoms) {
      if (first) {
        first = false;
        result = geom;
      } else {
        glm::dvec3 center;
        glm::dvec3 extents;
        result.GetCenterExtents(center, extents);

        glm::dvec3 s_center;
        glm::dvec3 s_extents;
        geom.GetCenterExtents(s_center, s_extents);

        if (normalDiff(extents) < EPS_BIG) {
          result = geom;
        } else if (normalDiff(s_extents) < EPS_BIG) {
          // result = result;
        } else if (result.numFaces > 0 && geom.numFaces > 0) {
          auto first = result.Normalize(center, extents);
          auto second = geom.Normalize(center, extents);

          IfcGeometry r1;
          IfcGeometry r2;

          intersectMeshMesh(first, second, r1, r2);

          result = (boolJoin(r1, r2)).DeNormalize(center, extents);
        }
      }
    }
    return result;
  }
}*/

/*IfcGeometry ConwayGeometryProcessor::BoolSubtractLegacy(
    const std::vector<IfcGeometry> &firstGeoms,
    std::vector<IfcGeometry> &secondGeoms) {
  IfcGeometry firstGeom = firstGeoms[0];//BoolJoinLegacy(firstGeoms);
  IfcGeometry secondGeom = secondGeoms[0];//BoolJoinLegacy(secondGeoms);

  IfcGeometry result;

  if (firstGeom.numFaces == 0 || secondGeoms.size() == 0) {
    printf("bool aborted due to empty source or target\n");
    printf("firstGeom.numFaces: %i\n", firstGeom.numFaces);
    printf("secondGeoms.size(): %i\n", secondGeoms.size());

    // bail out because we will get strange meshes
    // if this happens, probably there's an issue parsing the mesh that occurred
    // earlier
    return firstGeom;
  }

  result = boolMultiOp_CSGJSCPP(firstGeom, secondGeoms);

  return result;
}*/

IfcCurve ConwayGeometryProcessor::GetCShapeCurve(
    ParamsGetCShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetCShapedCurve(parameters.width, parameters.depth, parameters.girth,
                          parameters.thickness, parameters.hasFillet,
                          parameters.filletRadius, parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetIShapeCurve(
    ParamsGetIShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetIShapedCurve(parameters.width, parameters.depth,
                          parameters.webThickness, parameters.flangeThickness,
                          parameters.hasFillet, parameters.filletRadius,
                          parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetLShapeCurve(
    ParamsGetLShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetLShapedCurve(parameters.width, parameters.depth,
                          parameters.thickness, parameters.hasFillet,
                          parameters.filletRadius, parameters.edgeRadius,
                          parameters.legSlope, parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetTShapeCurve(
    ParamsGetTShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetTShapedCurve(parameters.width, parameters.depth,
                          parameters.webThickness, parameters.hasFillet,
                          parameters.filletRadius, parameters.flangeEdgeRadius,
                          parameters.flangeScope, parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetUShapeCurve(
    ParamsGetUShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetUShapedCurve(parameters.depth, parameters.flangeWidth,
                          parameters.webThickness, parameters.flangeThickness,
                          parameters.filletRadius, parameters.edgeRadius,
                          parameters.flangeScope, parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetZShapeCurve(
    ParamsGetZShapeCurve parameters) {
  IfcCurve curve;
  glm::dmat3 placement(1);
  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetZShapedCurve(parameters.depth, parameters.flangeWidth,
                          parameters.webThickness, parameters.flangeThickness,
                          parameters.filletRadius, parameters.edgeRadius,
                          parameters.placement);

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetRectangleProfileCurve(
    ParamsGetRectangleProfileCurve parameters) {
  IfcCurve curve;
  double xdim = parameters.xDim;
  double ydim = parameters.yDim;

  if (parameters.hasPlacement) {
    curve = GetRectangleCurve(xdim, ydim, parameters.matrix);
  } else {
    glm::dmat3 placement = glm::dmat3(glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0),
                                      glm::dvec3(0, 0, 1));
    curve = GetRectangleCurve(xdim, ydim, placement);
  }

  if (!curve.IsCCW()) {
    curve.Invert();
  }

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetRectangleHollowProfileHole(
    ParamsGetRectangleProfileCurve parameters) {
  IfcCurve curve;
  double xdim = parameters.xDim;
  double ydim = parameters.yDim;

  curve = GetRectangleCurve(xdim - parameters.thickness,
                            ydim - parameters.thickness, parameters.matrix);

  curve.Invert();

  return curve;
}

glm::dvec3 GetOrigin(const IfcGeometry &geometry) {
  return (geometry.min + geometry.max) * 0.5;
}

glm::dvec3 CalculateCentroid(const IfcGeometry &geometry) {
  glm::dvec3 centroid(0.0);
  uint32_t numVertices = 0;

  // I'm assuming the IfcGeometry class provides a method to iterate over all
  // vertices. If it doesn't, you might need to iterate over faces and access
  // vertices via face indices.
  for (uint32_t i = 0; i < geometry.numFaces; i++) {
    fuzzybools::Face f = geometry.GetFace(i);
    centroid += geometry.GetPoint(f.i0);
    centroid += geometry.GetPoint(f.i1);
    centroid += geometry.GetPoint(f.i2);
    numVertices += 3;
  }

  if (numVertices > 0) {
    centroid /= static_cast<double>(numVertices);
  }

  return centroid;
}

IfcGeometry ConwayGeometryProcessor::BoolSubtract(
    const std::vector<IfcGeometry> &firstGeoms,
    std::vector<IfcGeometry> &secondGeoms) {
  IfcGeometry finalResult;

  for (auto &firstGeom : firstGeoms) {
    fuzzybools::Geometry result = firstGeom;
    for (auto &secondGeom : secondGeoms) {
      bool doit = true;
      if (secondGeom.numFaces == 0) {
        Logger::logWarning("bool aborted due to empty source or target\n");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        doit = false;
      }

      if (result.numFaces == 0) {
        Logger::logWarning("bool aborted due to empty source or target\n");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        break;
      }

      if (doit) {
        if (secondGeom.halfSpace) {
          glm::dvec3 origin = secondGeom.halfSpaceOrigin;
          glm::dvec3 x = secondGeom.halfSpaceX - origin;
          glm::dvec3 y = secondGeom.halfSpaceY - origin;
          glm::dvec3 z = secondGeom.halfSpaceZ - origin;
          glm::dmat4 trans =
              glm::dmat4(glm::dvec4(x, 0), glm::dvec4(y, 0), glm::dvec4(z, 0),
                         glm::dvec4(0, 0, 0, 1));
          IfcGeometry newSecond;

          double scaleX = 1;
          double scaleY = 1;
          double scaleZ = 1;

          for (uint32_t i = 0; i < result.numPoints; i++) {
            glm::dvec3 p = result.GetPoint(i);
            glm::dvec3 vec = (p - origin);
            double dx = glm::dot(vec, x);
            double dy = glm::dot(vec, y);
            double dz = glm::dot(vec, z);
            if (glm::abs(dx) > scaleX) {
              scaleX = glm::abs(dx);
            }
            if (glm::abs(dy) > scaleY) {
              scaleY = glm::abs(dy);
            }
            if (glm::abs(dz) > scaleZ) {
              scaleZ = glm::abs(dz);
            }
          }
          newSecond.AddGeometry(secondGeom, trans, scaleX * 2, scaleY * 2,
                                scaleZ * 2, secondGeom.halfSpaceOrigin);
          result = fuzzybools::Subtract(result, newSecond);
        } else {
          result = fuzzybools::Subtract(result, secondGeom);
        }
      }
    }

    IfcGeometry newResult;
    newResult.AddGeometry(result);
    finalResult.AddPart(newResult);
    finalResult.AddGeometry(result);
  }

  return finalResult;
}

IfcGeometry ConwayGeometryProcessor::RelVoidSubtract(
    ParamsRelVoidSubtract parameters) {
  IfcGeometry resultGeometry;
  if (parameters.flatFirstMesh.size() <= 0) {
    Logger::logWarning("first mesh zero\n");
    return resultGeometry;
  }

  if (parameters.flatSecondMesh.size() <= 0) {
    Logger::logWarning("second mesh zero\n");
    return resultGeometry;
  }

  glm::dvec3 originFirstMesh = GetOrigin(parameters.flatFirstMesh[0]);
  // get origin
  if (parameters.flatFirstMesh[0].numFaces) {
    for (uint32_t i = 0; i < parameters.flatFirstMesh[0].numFaces; i++) {
      fuzzybools::Face f = parameters.flatFirstMesh[0].GetFace(i);
      originFirstMesh = parameters.flatFirstMesh[0].GetPoint(f.i0);
      break;
    }
  }

  // TODO: clean this up, remove origin translation
  auto normalizeMat = glm::translate(-originFirstMesh);
  glm::dmat4 newMatrix = normalizeMat;  // * parameters.transformationFirstMesh;
  bool transformationBreaksWinding = MatrixFlipsTriangles(newMatrix);
  IfcGeometry newGeomFirstMesh;

  for (uint32_t i = 0; i < parameters.flatFirstMesh[0].numFaces; i++) {
    fuzzybools::Face f = parameters.flatFirstMesh[0].GetFace(i);
    glm::dvec3 a =
        newMatrix * glm::dvec4(parameters.flatFirstMesh[0].GetPoint(f.i0), 1);
    glm::dvec3 b =
        newMatrix * glm::dvec4(parameters.flatFirstMesh[0].GetPoint(f.i1), 1);
    glm::dvec3 c =
        newMatrix * glm::dvec4(parameters.flatFirstMesh[0].GetPoint(f.i2), 1);

    if (transformationBreaksWinding) {
      newGeomFirstMesh.AddFace(b, a, c);
    } else {
      newGeomFirstMesh.AddFace(a, b, c);
    }
  }

  parameters.flatFirstMesh[0] = newGeomFirstMesh;
  transformationBreaksWinding = MatrixFlipsTriangles(newMatrix);
  IfcGeometry newGeomSecondMesh;

  for (uint32_t i = 0; i < parameters.flatSecondMesh[0].numFaces; i++) {
    fuzzybools::Face f = parameters.flatSecondMesh[0].GetFace(i);
    glm::dvec3 a =
        newMatrix * glm::dvec4(parameters.flatSecondMesh[0].GetPoint(f.i0), 1);
    glm::dvec3 b =
        newMatrix * glm::dvec4(parameters.flatSecondMesh[0].GetPoint(f.i1), 1);
    glm::dvec3 c =
        newMatrix * glm::dvec4(parameters.flatSecondMesh[0].GetPoint(f.i2), 1);

    if (transformationBreaksWinding) {
      newGeomSecondMesh.AddFace(b, a, c);
    } else {
      newGeomSecondMesh.AddFace(a, b, c);
    }
  }

  parameters.flatSecondMesh[0] = newGeomSecondMesh;

  resultGeometry =
      BoolSubtract(parameters.flatFirstMesh, parameters.flatSecondMesh);

  // BoolSubtractLegacy(parameters.flatFirstMesh, parameters.flatSecondMesh);

  glm::dmat4 combinedMatrix =
      glm::inverse(parameters.parentMatrix) * glm::translate(originFirstMesh);

  resultGeometry.ApplyTransform(combinedMatrix);

  /*resultGeometry.ApplyTransform(glm::translate(originFirstMesh));

  glm::dmat4 inverseParentMatrix = glm::inverse(parameters.parentMatrix);

  resultGeometry.ApplyTransform(inverseParentMatrix);*/

  return resultGeometry;
}

IfcGeometry ConwayGeometryProcessor::GetBooleanResult(
    ParamsGetBooleanResult *parameters) {
  IfcGeometry resultGeometry;
  if (parameters->flatFirstMesh.size() <= 0) {
    Logger::logWarning("first mesh zero\n");
    return resultGeometry;
  }

  if (parameters->flatSecondMesh.size() <= 0) {
    Logger::logWarning("second mesh zero\n");
    return resultGeometry;
  }

  resultGeometry =
      BoolSubtract(parameters->flatFirstMesh, parameters->flatSecondMesh);

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

/*IfcGeometry ConwayGeometryProcessor::getBrep(ParamsGetBrep parameters) {
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
}*/

void ConwayGeometryProcessor::AddFaceToGeometry(
    ParamsAddFaceToGeometry& parameters, IfcGeometry &geometry) {
  if (!parameters.advancedBrep) {
    if (parameters.boundsArray.size() > 0) {
      TriangulateBounds(geometry, parameters.boundsArray);
    }
  } else {
    if (parameters.boundsArray.size() > 0) {
      auto surface = parameters.surface;

      if (surface.BSplineSurface.Active) {
        TriangulateBspline(geometry, parameters.boundsArray, surface,
                           parameters.scaling);
      } else if (surface.CylinderSurface.Active) {
        TriangulateCylindricalSurface(geometry, parameters.boundsArray,
                                      surface);
      } else if (surface.RevolutionSurface.Active) {
        TriangulateRevolution(geometry, parameters.boundsArray, surface);
      } else if (surface.ExtrusionSurface.Active) {
        TriangulateExtrusion(geometry, parameters.boundsArray, surface);
      } else {
        TriangulateBounds(geometry, parameters.boundsArray);
      }
    }
  }
}

/*std::vector<IfcGeometry> ConwayGeometryProcessor::GetSurfaceModel(
    ParamsGetSurfaceModel parameters) {
  std::vector<IfcGeometry> geometryArray;

  for (uint32_t shellIndex = 0; shellIndex < parameters.numShellRefs;
       shellIndex++) {
    geometryArray.push_back(getBrep(parameters.shells[shellIndex]));
  }

  return geometryArray;
}*/

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

    // TODO(nickcastel50): Old implementation wasn't returning a surface for
    // this case.
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

IfcCurve ConwayGeometryProcessor::getPolyCurve(
    const ParamsGetPolyCurve &parameters) {
  IfcCurve curve;
  const float *points = reinterpret_cast<float *>(parameters.points_);

  if (parameters.dimensions == 2) {
    for (int i = 0; i < parameters.pointsLength * parameters.dimensions;
         i += parameters.dimensions) {
      glm::dvec3 point;
      point.x = points[i];
      point.y = points[i + 1];
      point.z = 0;

      curve.Add3d(point);
    }
  } else if (parameters.dimensions == 3) {
    for (int i = 0; i < parameters.pointsLength * parameters.dimensions;
         i += parameters.dimensions) {
      glm::dvec3 point;
      point.x = points[i];
      point.y = points[i + 1];
      point.z = points[i + 2];

      curve.Add3d(point);
    }
  }

  return curve;
}

IfcProfile ConwayGeometryProcessor::transformProfile(
    ParamsTransformProfile *parameters) {
  if (!parameters->profile.isComposite) {
    for (uint32_t i = 0; i < parameters->profile.curve.points.size(); i++) {
      parameters->profile.curve.points[i] =
          parameters->transformation *
          glm::dvec3(parameters->profile.curve.points[i].x,
                     parameters->profile.curve.points[i].y, 1);
    }
  } else {
    for (uint32_t j = 0; j < parameters->profile.profiles.size(); j++) {
      for (uint32_t i = 0;
           i < parameters->profile.profiles[j].curve.points.size(); i++) {
        parameters->profile.profiles[j].curve.points[i] =
            parameters->transformation *
            glm::dvec3(parameters->profile.profiles[j].curve.points[i].x,
                       parameters->profile.profiles[j].curve.points[i].y, 1);
      }
    }
  }

  return parameters->profile;
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
    const ParamsAxis2Placement3D &parameters) {
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
    const ParamsLocalPlacement &parameters) {
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
    const ParamsCartesianTransformationOperator3D &parameters) {
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
  if (parameters.edges.size() == 0) {
    if (parameters.points.size() > 0) {
      curve.points = parameters.points;
    }
  } else {
    int id = 0;
    for (int edgeIndex = 0; edgeIndex < parameters.edges.size(); ++edgeIndex) {
      // Important not to repeat the last point otherwise triangulation fails
      // if the list has zero points this is initial, no repetition is possible,
      // otherwise we must check
      if (curve.points.size() == 0) {
        for (auto &pt : parameters.edges[edgeIndex].points) {
          curve.points.push_back(pt);
          curve.indices.push_back(id);
        }
      } else {
        for (int i = 0; i < parameters.edges[edgeIndex].points.size(); ++i)
        // for (auto &pt : parameters.edges[edgeIndex].points)
        {
          glm::dvec3 pt = parameters.edges[edgeIndex].points[i];
          /*if (i == 0) {
            if (!notPresent(pt, curve.points)) {
              printf("point is PRESENT! - Point: x: %.3f, y: %.3f, z: %.3f\n",
                     pt.x, pt.y, pt.z);
            }
          } else {*/
          if (notPresent(pt, curve.points)) {
            curve.points.push_back(pt);
            curve.indices.push_back(id);
          }
          // }
        }
      }
      id++;
    }
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
    const ParamsReadIndexedPolygonalFace &parameters) {
  std::vector<IfcBound3D> bounds;

  bounds.emplace_back();

  // calculate loop upper bound
  size_t faceIndexUpperBound = 0;
  if (parameters.face.face_starts.size() > 1) {
    faceIndexUpperBound =
        parameters.face.face_starts[1] - parameters.face.face_starts[0];
  } else {
    faceIndexUpperBound = parameters.face.indices.size();
  }

  for (size_t index = 0; index < faceIndexUpperBound; index++) {
    uint32_t currentIndex = parameters.face.indices[index];

    glm::dvec3 point = parameters.points[currentIndex - 1];

    // I am not proud of this (I inherited this, will change - NC)
    bounds.back().curve.points.push_back(point);
  }

  if (parameters.face.face_starts.size() <= 1) {
    return bounds;
  } else {
    // TODO(nickcastel50): handle case IFCINDEXEDPOLYGONALFACEWITHVOIDS
    for (size_t i = 1; i < parameters.face.face_starts.size(); i++) {
      bounds.emplace_back();
      size_t startIdx = parameters.face.face_starts[i];
      size_t endIdx = 0;
      size_t faceVoidsSize = 0;
      if (i + 1 >= parameters.face.face_starts.size()) {
        endIdx = parameters.face.indices.size();
        faceVoidsSize = endIdx - startIdx;
      } else {
        endIdx = parameters.face.face_starts[i + 1];
        faceVoidsSize = endIdx - startIdx;
      }

      for (size_t index = startIdx; index < endIdx; index++) {
        uint32_t currentIndex = parameters.face.indices[index];

        glm::dvec3 point = parameters.points[currentIndex - 1];

        // I am not proud of this (I inherited this, will change - NC)
        bounds.back().curve.points.push_back(point);
      }
    }

    ;
  }

  return bounds;
}

conway::geometry::ConwayGeometryProcessor::ResultsGltf
ConwayGeometryProcessor::GeometryToGltf(
    std::span<conway::geometry::IfcGeometryCollection> geoms,
    std::span<conway::geometry::Material> materials, bool isGlb,
    bool outputDraco, std::string filePath, bool outputFile,
    glm::dmat4 transform) {
  ResultsGltf results;

  try {
    // Construct a Mesh and add the MeshPrimitive as a child
    Microsoft::glTF::Mesh mesh;
    // Construct a Scene
    Microsoft::glTF::Scene scene;
    // The Document instance represents the glTF JSON manifest
    Microsoft::glTF::Document document;

    std::vector<std::string> materialIds;

    int uniqueIdDraco = 0;

    for (conway::geometry::Material &conwayMaterial : materials) {
      Microsoft::glTF::Material material;

      material.doubleSided = conwayMaterial.doubleSided;
      material.metallicRoughness.metallicFactor = conwayMaterial.metallic;
      material.metallicRoughness.roughnessFactor = conwayMaterial.roughness;
      material.metallicRoughness.baseColorFactor =
          Microsoft::glTF::Color4(static_cast<float>(conwayMaterial.base.r),
                                  static_cast<float>(conwayMaterial.base.g),
                                  static_cast<float>(conwayMaterial.base.b),
                                  static_cast<float>(conwayMaterial.base.a));

      /* TODO - add support for specular colour and IOR here  */

      switch (conwayMaterial.alphaMode) {
        case BLEND_MODE::BLEND:

          material.alphaMode = Microsoft::glTF::AlphaMode::ALPHA_BLEND;
          break;

        case BLEND_MODE::MASK:

          material.alphaMode = Microsoft::glTF::AlphaMode::ALPHA_MASK;
          break;
      }

      // TODO - specular extension and IOR extension. - CS

      std::string materialId =
          document.materials
              .Append(std::move(material),
                      Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
              .id;

      materialIds.push_back(materialId);
    }

    Microsoft::glTF::Material defaultMaterial;

    defaultMaterial.doubleSided = true;

    std::optional<std::string> defaultMaterialId;

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

    // Use the BufferBuilder helper class to simplify the process of
    // constructing valid glTF Buffer, BufferView and Accessor entities
    Microsoft::glTF::BufferBuilder bufferBuilder(std::move(resourceWriter));

    std::string base_filename =
        filePath.substr(filePath.find_last_of("/\\") + 1);

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

    if (outputDraco) {
      document.extensionsRequired.insert("KHR_draco_mesh_compression");
      document.extensionsUsed.insert("KHR_draco_mesh_compression");
    }

    std::vector<Microsoft::glTF::KHR::MeshPrimitives::DracoMeshCompression>
        dracoPrimitives;

    for (conway::geometry::IfcGeometryCollection &geom : geoms) {
      // create a mesh object
      std::unique_ptr<draco::Mesh> dracoMesh;

      // Encode the geometry.
      draco::EncoderBuffer buffer;

      // Convert to ExpertEncoder that allows us to set per-attribute options.
      std::unique_ptr<draco::ExpertEncoder> expert_encoder;

      DracoOptions dracoOptions;

      // set up Draco Encoder
      draco::Encoder encoder;

      // if Draco
      Microsoft::glTF::KHR::MeshPrimitives::DracoMeshCompression
          dracoMeshCompression;

      int32_t pos_att_id = -1;
      int32_t tex_att_id = -1;
      int32_t material_att_id = -1;

      // this internally populates the vertex float array, current storage type
      // is double
      // geom.GetVertexData();

      if (outputDraco) {
        dracoMesh.reset(new draco::Mesh());

        // set the number of positions and faces
        uint32_t numPositions = 0;
        uint32_t numIndices = 0;

        for (conway::geometry::IfcGeometry *componentPtr : geom.components) {
          conway::geometry::IfcGeometry &component = *componentPtr;

          numPositions += component.numPoints;
          numIndices += component.GetIndexDataSize();
        }

        uint32_t numFaces = numIndices / 3;

        if (numFaces > 0) {
          // set number of faces
          dracoMesh->SetNumFaces(numFaces);

          // set number of indices
          dracoMesh->set_num_points(3 * numFaces);
        }

        // Add attributes if they are present in the input data.
        if (numPositions > 0) {
          draco::GeometryAttribute va;

          va.set_unique_id(uniqueIdDraco++);

          va.Init(draco::GeometryAttribute::POSITION, nullptr, 3,
                  draco::DT_FLOAT32, false, sizeof(float) * 3, 0);
          pos_att_id = dracoMesh->AddAttribute(va, false, numPositions);
        }

        // TODO(nickcastel50): support multiple materials at some point when we
        // add that
        int32_t numMaterials = 1;
        int32_t numTexCoords = 0;

        // populate position attribute

        uint32_t vertexCount = 0;
        uint32_t totalPoints = 0;
        uint32_t totalIndices = 0;

        // populate geometry components position attribute
        for (uint32_t geometryComponentIndex = 0;
             geometryComponentIndex < geom.components.size();
             ++geometryComponentIndex) {
          IfcGeometry &component = *geom.components[geometryComponentIndex];
          glm::mat4 geomTransform =
              transform * geom.transforms[geometryComponentIndex];

          for (const double *where = component.vertexData.data(),
                            *end = where + component.vertexData.size();
               where < end; where += VERTEX_FORMAT_SIZE_FLOATS) {
            glm::dvec4 t =
                geomTransform * glm::dvec4(where[0], where[1], where[2], 1);

            float vertexVal[3] = {static_cast<float>(t.x),
                                  static_cast<float>(t.y),
                                  static_cast<float>(t.z)};

            dracoMesh->attribute(pos_att_id)
                ->SetAttributeValue(draco::AttributeValueIndex(vertexCount++),
                                    vertexVal);
          }

          const uint32_t *triangles = component.indexData.data();

          // no textures, just map vertices to face indices
          for (size_t triangleIndex = 0, end = (component.indexData.size() / 3);
               triangleIndex < end; ++triangleIndex) {
            uint32_t indiceIndex = triangleIndex * 3;
            const uint32_t *triangle = triangles + indiceIndex;

            uint32_t triangle0 = triangle[0] + totalPoints;
            uint32_t triangle1 = triangle[1] + totalPoints;
            uint32_t triangle2 = triangle[2] + totalPoints;

            uint32_t compositeIndiceIndex = totalIndices + indiceIndex;

            const draco::PointIndex vert_id_0(compositeIndiceIndex);
            const draco::PointIndex vert_id_1(compositeIndiceIndex + 1);
            const draco::PointIndex vert_id_2(compositeIndiceIndex + 2);
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

          totalIndices += component.indexData.size();
          totalPoints += component.numPoints;
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
        const int32_t dracoCompressionSpeed =
            10 - dracoOptions.compressionLevel;

        // Setup encoder options.
        if (dracoOptions.posQuantizationBits > 0) {
          encoder.SetAttributeQuantization(
              draco::GeometryAttribute::POSITION,
              14);  // dracoOptions.posQuantizationBits );
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

        expert_encoder.reset(new draco::ExpertEncoder(*dracoMesh));

        expert_encoder->Reset(encoder.CreateExpertEncoderOptions(*dracoMesh));
        // set up timer
        draco::CycleTimer timer;

        timer.Start();

        const draco::Status status = expert_encoder->EncodeToBuffer(&buffer);

        timer.Stop();

        if (!status.ok()) {
          Logger::logError("Failed to encode the mesh: %s\n", status.error_msg());
          results.success = false;
          return results;
        }

        if (outputFile) {
          Logger::logInfo("Encoded To Draco in %lld ms\n", timer.GetInMs());
        }
      }

      std::string accessorIdIndices;
      std::string accessorIdPositions;
      std::string accessorIdUVs;

      uint32_t numPoints = 0;
      uint32_t numIndices = 0;

      for (conway::geometry::IfcGeometry *componentPtr : geom.components) {
        conway::geometry::IfcGeometry &component = *componentPtr;

        numPoints += component.numPoints;
        numIndices += component.GetIndexDataSize();
      }

      // printf("numPoints: %i\n", numPoints);
      // printf("numIndices: %i\n", numIndices);

      // Add an Accessor for the indices and positions
      // std::unique_ptr< std::vector< float > > positionsPtr    =
      //   std::make_unique< std::vector< float > >();
      // std::unique_ptr< std::vector< uint32_t > > indexDataPtr =
      //   std::make_unique< std::vector< uint32_t > >();

      // positionOwner.push_back( std::vector< float >( numPoints * 3 ) );
      // indiceOwner.push_back( std::vector< uint32_t >() );

      std::vector<float> positions;
      std::vector<uint32_t> indexData;

      positions.resize(numPoints * 3);
      indexData.reserve(numIndices);

      uint32_t numFaces = numIndices / 3;

      std::vector<float> minValues(3U, std::numeric_limits<float>::max());
      std::vector<float> maxValues(3U, std::numeric_limits<float>::lowest());

      const size_t positionCount = positions.size();
      float *positionOutputCursor = positions.data();

      uint32_t pointOffset = 0;
      size_t transformCursor = 0;

      for (conway::geometry::IfcGeometry *componentPtr : geom.components) {
        conway::geometry::IfcGeometry &component = *componentPtr;
        glm::dmat4 geomTransform =
            transform * geom.transforms[transformCursor++];

        for (const double *where = component.vertexData.data(),
                          *end = where + component.vertexData.size();
             where < end; where += VERTEX_FORMAT_SIZE_FLOATS) {
          glm::dvec4 t =
              geomTransform * glm::dvec4(where[0], where[1], where[2], 1);

          *(positionOutputCursor++) = (float)t.x;
          *(positionOutputCursor++) = (float)t.y;
          *(positionOutputCursor++) = (float)t.z;
        }

        size_t indexDataOffset = indexData.size();

        indexData.insert(indexData.end(), component.indexData.begin(),
                         component.indexData.end());

        for (uint32_t *where = indexData.data() + indexDataOffset,
                      *end = indexData.data() + indexData.size();
             where < end; where += 3) {
          where[0] += pointOffset;
          where[1] += pointOffset;
          where[2] += pointOffset;
        }

        pointOffset += component.numPoints;
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
        dracoMeshCompression.attributes.emplace(
            "POSITION", dracoMesh->attribute(pos_att_id)->unique_id());

        Microsoft::glTF::Accessor indicesAccessor;
        indicesAccessor.count = dracoMesh->num_faces() * 3;
        indicesAccessor.componentType = Microsoft::glTF::COMPONENT_UNSIGNED_INT;
        indicesAccessor.type = Microsoft::glTF::TYPE_SCALAR;

        // add position accessor first here
        accessorIdPositions =
            document.accessors
                .Append(positionsAccessor,
                        Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
                .id;

        accessorIdIndices =
            document.accessors
                .Append(indicesAccessor,
                        Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
                .id;

      } else {
        // Create a BufferView with a target of ELEMENT_ARRAY_BUFFER (as it will
        // reference index data) - it will be the 'current' BufferView that all
        // the Accessors created by this BufferBuilder will automatically
        // reference
        bufferBuilder.AddBufferView(
            Microsoft::glTF::BufferViewTarget::ELEMENT_ARRAY_BUFFER);

        // // Copy the Accessor's id - subsequent calls to AddAccessor may
        // // invalidate the returned reference
        accessorIdIndices =
            bufferBuilder
                .AddAccessor(indexData,
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
                .AddAccessor(positions,
                             {Microsoft::glTF::TYPE_VEC3,
                              Microsoft::glTF::COMPONENT_FLOAT, false,
                              std::move(minValues), std::move(maxValues)})
                .id;
      }

      Microsoft::glTF::TextureInfo textureInfo;

      // Construct a MeshPrimitive. Unlike most types in glTF, MeshPrimitives
      // are direct children of their parent Mesh entity rather than being
      // children of the Document. This is why they don't have an ID member.
      Microsoft::glTF::MeshPrimitive meshPrimitive;

      if (!geom.hasDefaultMaterial &&
          size_t(geom.materialIndex) < materialIds.size()) {
        meshPrimitive.materialId = materialIds[geom.materialIndex];
      } else {
        if (!defaultMaterialId.has_value()) {
          // Add it to the Document and store the generated ID
          defaultMaterialId =
              document.materials
                  .Append(std::move(defaultMaterial),
                          Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty)
                  .id;
        }

        meshPrimitive.materialId = *defaultMaterialId;
      }

      if (outputDraco) {
        // Finally put the encoded data in place.
        auto bufferView =
            bufferBuilder.AddBufferView(buffer.data(), buffer.size());
        meshPrimitive.indicesAccessorId = accessorIdIndices;
        meshPrimitive.attributes[Microsoft::glTF::ACCESSOR_POSITION] =
            accessorIdPositions;

        dracoMeshCompression.bufferViewId = bufferView.id;

        dracoPrimitives.push_back(dracoMeshCompression);

      } else {
        meshPrimitive.indicesAccessorId = accessorIdIndices;
        meshPrimitive.attributes[Microsoft::glTF::ACCESSOR_POSITION] =
            accessorIdPositions;
      }

      mesh.primitives.push_back(std::move(meshPrimitive));
    }

    // Add all of the Buffers, BufferViews and Accessors that were created
    // using BufferBuilder to the Document. Note that after this point, no
    // further calls should be made to BufferBuilder
    bufferBuilder.Output(document);

    if (outputDraco) {
      for (int i = 0; i < mesh.primitives.size(); i++) {
        const auto extensionSerializer =
            Microsoft::glTF::KHR::GetKHRExtensionSerializer();
        std::string strDracoMesh =
            Microsoft::glTF::KHR::MeshPrimitives::SerializeDracoMeshCompression(
                dracoPrimitives[i], document, extensionSerializer);
        mesh.primitives[i].extensions["KHR_draco_mesh_compression"] =
            strDracoMesh;
      }
    }

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
      Logger::logError("Microsoft::glTF::Serialize failed: %s\n", ex.what());
    } catch (...) {
      Logger::logError("Microsoft::glTF::Serialize failed: <unknown>\n");
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
          Logger::logError("Unable to open file %s\n", filePath.c_str());
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
        Logger::logInfo("Writing: %s\n", filePath.c_str());

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
    Logger::logError("Couldn't write GLB file: %s\n", ex.what());
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

IfcGeometry ConwayGeometryProcessor::getTriangulatedFaceSetGeometry(
    const ParamsGetTriangulatedFaceSetGeometry &parameters) {
  const uint32_t *indices =
      reinterpret_cast<uint32_t *>(parameters.indicesArray_);
  const float *points = reinterpret_cast<float *>(parameters.pointsArray_);

  IfcGeometry geom;

  for (size_t i = 0; i < parameters.indicesArrayLength; i += 3) {
    int i1 = indices[i + 0] - 1;
    int i2 = indices[i + 1] - 1;
    int i3 = indices[i + 2] - 1;

    geom.AddFace(points[i1], points[i2], points[i3]);
    // printf("adding face %i: x: %.3f, y: %.3f, z: %.3f\n", i, points[i1],
    // points[i2], points[i3]);
  }

  return geom;
}

IfcGeometry ConwayGeometryProcessor::getPolygonalFaceSetGeometry(
    const ParamsGetPolygonalFaceSetGeometry &parameters) {
  IfcGeometry geom;
  std::vector<IfcBound3D> bounds;

  for (size_t faceIndex = 0; faceIndex < parameters.faces.size(); faceIndex++) {
    ParamsReadIndexedPolygonalFace params(parameters.points,
                                          parameters.faces[faceIndex]);

    bounds = ReadIndexedPolygonalFace(params);

    TriangulateBounds(geom, bounds);

    bounds.clear();
  }

  return geom;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getIndexedPolyCurve(
    const ParamsGetIfcIndexedPolyCurve &parameters) {
  IfcCurve curve;

  if (parameters.dimensions == 2) {
    if (parameters.segments.size() > 0) {
      for (auto &sg : parameters.segments) {
        if (!sg.isArcType) {
          auto pts = parameters.points;
          for (auto &pt : sg.indices) {
            curve.Add2d(pts[pt - 1]);
          }
        }
        if (sg.isArcType) {
          auto pts = parameters.points;
          IfcCurve arc =
              BuildArc3Pt(pts[sg.indices[0] - 1], pts[sg.indices[1] - 1],
                          pts[sg.indices[2] - 1], CIRCLE_SEGMENTS_MEDIUM);
          for (auto &pt : arc.points) {
            curve.Add2d(pt);
          }
        }
      }
    } else {
      auto pts = parameters.points;
      for (auto &pt : pts) {
        curve.Add2d(pt);
      }
    }
  } else {
    Logger::logWarning("Parsing ifcindexedpolycurve in 3D is not possible\n");
  }

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getEllipseCurve(
    const ParamsGetEllipseCurve &parameters) {
  IfcCurve curve;

  double radiusX = parameters.radiusX;
  double radiusY = parameters.radiusY;

  if (parameters.hasPlacement) {
    curve = GetEllipseCurve(radiusX, radiusY, parameters.circleSegments,
                            parameters.placement);
  } else {
    glm::dmat3 placement = glm::dmat3(glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0),
                                      glm::dvec3(0, 0, 1));
    curve =
        GetEllipseCurve(radiusX, radiusY, parameters.circleSegments, placement);
  }

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getCircleCurve(
    const ParamsGetCircleCurve &parameters) {
  IfcCurve curve;

  double radius = parameters.radius;

  glm::dmat3 placement(1);

  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetCircleCurve(radius, ConwayGeometryProcessor::CIRCLE_SEGMENTS_HIGH,
                         placement);

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getCircleHoleCurve(
    const ParamsGetCircleCurve &parameters) {
  IfcCurve curve;

  double radius = parameters.radius - parameters.thickness;

  glm::dmat3 placement(1);

  if (parameters.hasPlacement) {
    placement = parameters.placement;
  }

  curve = GetCircleCurve(radius, ConwayGeometryProcessor::CIRCLE_SEGMENTS_HIGH,
                         placement);

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getBSplineCurve(
    const ParamsGetBSplineCurve &parameters) {
  IfcCurve curve;

  // bool condition = false;
  /*if (edge) {
    condition = !condition;
  }*/

  curve.points.clear();

  int dimensions = parameters.dimensions;

  int degree = parameters.degree;

  if (dimensions == 2) {
    // printf("parameters.points2 size: %i\n", parameters.points2.size());
    // printf("degree: %i\n", degree);
    // printf("parameters.knots:\n");

    // for (int i = 0; i < parameters.knots.size(); ++i) {
    //   printf("knot %i: %.3f\n", i, parameters.knots[i]);
    // }

    // printf("parameters.weights:\n");

    // for (int i = 0; i < parameters.weights.size(); ++i) {
    //   printf("weight %i: %.3f\n", i, parameters.weights[i]);
    // }

    std::vector<glm::dvec2> tempPoints = GetRationalBSplineCurveWithKnots(
        degree, parameters.points2, parameters.knots, parameters.weights);
    for (size_t i = 0; i < tempPoints.size(); i++) {
      // printf("Point %i: x: %.3f, y: %.3f\n", i, tempPoints[i].x,
      //        tempPoints[i].y);
      curve.Add2d(tempPoints[i]);
    }
  } else if (dimensions == 3) {
    // printf("parameters.points3 size: %i\n", parameters.points3.size());
    // printf("degree: %i\n", degree);
    // printf("parameters.knots:\n");

    // for (int i = 0; i < parameters.knots.size(); ++i) {
    //   printf("knot %i: %.3f\n", i, parameters.knots[i]);
    // }

    // printf("parameters.weights:\n");

    // for (int i = 0; i < parameters.weights.size(); ++i) {
    //   printf("weight %i: %.3f\n", i, parameters.weights[i]);
    // }
    
    std::vector<glm::dvec3> tempPoints = GetRationalBSplineCurveWithKnots(
        degree, parameters.points3, parameters.knots, parameters.weights);
    
    // printf( "Output points: %zu\n", tempPoints.size() );

     for (size_t i = 0; i < tempPoints.size(); i++) {
    //   printf("Point %i: x: %.3f, y: %.3f, z: %.3f\n", i, tempPoints[i].x,
    //          tempPoints[i].y, tempPoints[i].z);
      curve.Add3d(tempPoints[i]);
    }
  }

  /*if (condition) {
    std::reverse(curve.points.begin(), curve.points.end());
  }*/

  return curve;
}

enum IfcTrimmingPreference { CARTESIAN = 0, PARAMETER = 1, UNSPECIFIED = 2 };

conway::geometry::IfcCurve ConwayGeometryProcessor::getIfcCircle(
    const ParamsGetIfcCircle &parameters) {
  conway::geometry::IfcCurve curve;

  double radius = parameters.radius;

  double startDegrees = 0;
  double endDegrees = 360;

  if (parameters.paramsGetIfcTrimmedCurve.dimensions != 0) {
    if (parameters.paramsGetIfcTrimmedCurve.masterRepresentation ==
        IfcTrimmingPreference::PARAMETER) {
      startDegrees = parameters.paramsGetIfcTrimmedCurve.trim1Double;
      endDegrees = parameters.paramsGetIfcTrimmedCurve.trim2Double;
    } else {
      if (parameters.dimensions == 2) {
        glm::dmat3 placement = parameters.axis2Placement2D;
        double xx = placement[2].x -
                    parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.x;
        double yy = placement[2].y -
                    parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.y;
        startDegrees = VectorToAngle(xx, yy);
        xx = placement[2].x -
             parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D.x;
        yy = placement[2].y -
             parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D.y;
        endDegrees = VectorToAngle(xx, yy);
      } else if (parameters.dimensions == 3) {
        glm::dmat4 placement = parameters.axis2Placement3D;
        glm::dvec4 vecX = placement[0];
        glm::dvec4 vecY = placement[1];
        glm::dvec4 vecZ = placement[2];

        glm::dvec3 v1 =
            glm::dvec3(parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D.x -
                           placement[3].x,
                       parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D.y -
                           placement[3].y,
                       parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D.z -
                           placement[3].z);
        glm::dvec3 v2 =
            glm::dvec3(parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D.x -
                           placement[3].x,
                       parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D.y -
                           placement[3].y,
                       parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D.z -
                           placement[3].z);

        double dxS = vecX.x * v1.x + vecX.y * v1.y + vecX.z * v1.z;
        double dyS = vecY.x * v1.x + vecY.y * v1.y + vecY.z * v1.z;
        // double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

        double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
        double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
        // double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

        endDegrees = VectorToAngle(dxS, dyS) - 90;
        startDegrees = VectorToAngle(dxE, dyE) - 90;
      }
    }
  }

  double startRad = startDegrees / 180 * CONST_PI;
  double endRad = endDegrees / 180 * CONST_PI;
  double lengthDegrees = endDegrees - startDegrees;

  // unset or true
  // TODO(nickcastel50): incorporate trimSense / sameSense as parameters
  uint32_t trimSense = -1;
  uint32_t sameSense = -1;
  if (trimSense == 1 || trimSense == -1) {
    if (lengthDegrees < 0) {
      lengthDegrees += 360;
    }
  } else {
    if (lengthDegrees > 0) {
      lengthDegrees -= 360;
    }
  }

  double lengthRad = lengthDegrees / 180 * CONST_PI;

  size_t startIndex = curve.points.size();

  for (int i = 0; i < CIRCLE_SEGMENTS_MEDIUM; i++) {
    double ratio = static_cast<double>(i) / (CIRCLE_SEGMENTS_MEDIUM - 1);
    double angle = startRad + ratio * lengthRad;

    if (sameSense == 0) {
      angle = endRad - ratio * lengthRad;
    }

    if (parameters.dimensions == 2) {
      glm::dvec2 vec(0);
      vec[0] = radius * std::cos(angle);
      vec[1] =
          -radius *
          std::sin(angle);  // not sure why we need this, but we apparently do
      glm::dvec2 pos = parameters.axis2Placement2D * glm::dvec3(vec, 1);
      curve.Add2d(pos);
    } else {
      glm::dvec3 vec(0);
      vec[0] = radius * std::cos(angle);
      vec[1] = -radius * std::sin(angle);  // negative or not???
      glm::dvec3 pos =
          parameters.axis2Placement3D * glm::dvec4(glm::dvec3(vec), 1);
      curve.Add3d(pos);
    }
  }

  // without a trim, we close the circle
  if (parameters.paramsGetIfcTrimmedCurve.dimensions == 0) {
    if (parameters.dimensions == 2) {
      curve.Add2d(curve.points[startIndex]);
    } else if (parameters.dimensions == 3) {
      curve.Add3d(curve.points[startIndex]);
    }
  }

  return curve;
}

conway::geometry::IfcGeometry ConwayGeometryProcessor::getExtrudedAreaSolid(
    const ParamsGetExtrudedAreaSolid &parameters) {
  conway::geometry::IfcGeometry geom;
  double depth = parameters.depth;

  conway::geometry::IfcProfile profile = parameters.profile;

  if (!profile.isComposite) {
    if (profile.curve.points.empty()) {
      printf("empty curve, returning...\n");
      return geom;
    }
  } else {
    for (uint32_t i = 0; i < profile.profiles.size(); i++) {
      if (profile.profiles[i].curve.points.empty()) {
        printf("empty curve, returning...\n");
        return geom;
      }
    }
  }

  glm::dvec3 dir = parameters.dir;

  double dirDot = glm::dot(dir, glm::dvec3(0, 0, 1));
  // TODO(nickcastel50): I believe it can be now...
  // https://standards.buildingsmart.org/IFC/RELEASE/IFC2x3/TC1/HTML/ifcgeometricmodelresource/lexical/ifcextrudedareasolid.htm
  bool flipWinding = dirDot < 0;  // can't be perp according to spec

// TODO: correct dump in case of compositeProfile
#ifdef CSG_DEBUG_OUTPUT
  io::DumpSVGCurve(profile.curve.points, "IFCEXTRUDEDAREASOLID_curve.html");
#endif

  if (!profile.isComposite) {
    geom = Extrude(profile, dir, depth);
    if (flipWinding) {
      for (uint32_t i = 0; i < geom.numFaces; i++) {
        uint32_t temp = geom.indexData[i * 3 + 0];
        temp = geom.indexData[i * 3 + 0];
        geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
        geom.indexData[i * 3 + 1] = temp;
      }
    }
  } else {
    for (uint32_t i = 0; i < profile.profiles.size(); i++) {
      IfcGeometry geom_t = Extrude(profile.profiles[i], dir, depth);
      if (flipWinding) {
        for (uint32_t k = 0; k < geom_t.numFaces; k++) {
          uint32_t temp = geom_t.indexData[k * 3 + 0];
          temp = geom_t.indexData[k * 3 + 0];
          geom_t.indexData[k * 3 + 0] = geom_t.indexData[k * 3 + 1];
          geom_t.indexData[k * 3 + 1] = temp;
        }
      }
      geom.AddPart(geom_t);
      geom.AppendGeometry(geom_t);
    }
  }

// TODO: correct dump in case of compositeProfile
#ifdef CSG_DEBUG_OUTPUT
  io::DumpIfcGeometry(geom, "IFCEXTRUDEDAREASOLID_geom.obj");
#endif

  return geom;
}

}  // namespace conway::geometry
