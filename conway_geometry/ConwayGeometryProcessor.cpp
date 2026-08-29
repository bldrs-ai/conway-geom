#include "ConwayGeometryProcessor.h"

#include <utility>
#include <cmath>

#include <glm/glm.hpp>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wreturn-type"
//#include "fuzzy/fuzzy-bools.h"

#pragma clang diagnostic pop

#include "operations/curve_utils.h"
#include "operations/geometry_utils.h"
#include "operations/mesh_utils.h"
#include "representation/Geometry.h"
#include "Logger.h"

#include "csg/csg.h"

#include "structures/alloc_telemetry.h"
#include "structures/thread_pool.h"
#include "structures/vertex_welder.h"

namespace conway::geometry {

double normalDiff( const glm::dvec3& extents ) {
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


IfcCurve ConwayGeometryProcessor::GetCShapeCurve(
    const ParamsGetCShapeCurve& parameters) {
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
    const ParamsGetIShapeCurve& parameters) {
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
    const ParamsGetLShapeCurve& parameters) {
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
    const ParamsGetTShapeCurve& parameters) {
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
    const ParamsGetUShapeCurve& parameters) {
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
    const ParamsGetZShapeCurve& parameters) {
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
    const ParamsGetRectangleProfileCurve& parameters) {
  IfcCurve curve;
  double xdim = parameters.xDim;
  double ydim = parameters.yDim;

  if (parameters.hasPlacement) {

    if (parameters.hasRoundingRadius && parameters.roundingRadius > 0.0) {
      curve = GetRoundedRectangleCurve(xdim, ydim, parameters.roundingRadius, parameters.circleSegments, parameters.matrix);
    } else {
      curve = GetRectangleCurve(xdim, ydim, parameters.matrix);
    }
  } else {
    glm::dmat3 placement = glm::dmat3(glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0),
                                      glm::dvec3(0, 0, 1));

    if (parameters.hasRoundingRadius && parameters.roundingRadius > 0.0) {
      curve = GetRoundedRectangleCurve(xdim, ydim, parameters.roundingRadius, parameters.circleSegments, placement);
    } else {
      curve = GetRectangleCurve(xdim, ydim, placement);
    }
  }

  if (!curve.IsCCW()) {
    curve.Invert();
  }

  return curve;
}

IfcCurve ConwayGeometryProcessor::GetRectangleHollowProfileHole(
    const ParamsGetRectangleProfileCurve& parameters) {
  IfcCurve curve;
  double xdim = parameters.xDim;
  double ydim = parameters.yDim;

  curve = GetRectangleCurve(xdim - parameters.thickness,
                            ydim - parameters.thickness, parameters.matrix);

  curve.Invert();

  return curve;
}

glm::dvec3 GetOrigin(const Geometry &geometry) {
  
  box3 box = geometry.GetAABB();

  return box.centre();
}

glm::dvec3 CalculateCentroid( const Geometry &geometry ) {

  return centroid( geometry.vertices, geometry.triangles );
}

namespace {

  CSG csg;

}

Geometry ConwayGeometryProcessor::BoolSubtract(
    std::vector<Geometry> &firstGeoms,
    std::vector<Geometry> &secondGeoms,
    bool isSubtractOperand ) {

  // The single entry point for CSG/boolean composition (RelVoidSubtract and
  // GetBooleanResult both funnel here), and the other half of the call graph
  // conway#637 asks about. The unit is one boolean composition -- which may
  // run several csg.run() passes internally -- so peak-vs-retained here is per
  // composition, not per operand.
  AllocTelemetryScope telemetryScope( AllocSite::CsgBoolean );

  //printf( "%s\n", isSubtractOperand ? "sop" : "nop" );

  if ( firstGeoms.size() == 1 && secondGeoms.size() == 1 && !secondGeoms[ 0 ].halfSpace ) {
    
    Geometry resultMesh;

    {
      // Operand conditioning is separated from the kernel because the two have
      // opposite memory shapes: Cleanup()'s weld/dedupe structures are pure
      // transients an arena could back, while the kernel's output is retained.
      // Deciding whether to extend the arena here needs them apart.
      AllocTagScope prepTag( AllocSite::CsgOperandPrep );

      firstGeoms[ 0 ].Cleanup();
      secondGeoms[ 0 ].Cleanup( !isSubtractOperand );
    }

    {
      AllocTagScope kernelTag( AllocSite::CsgKernel );

      csg.run( CSG::Operation::SUBTRACTION, firstGeoms[ 0 ], secondGeoms[ 0 ], resultMesh, 0 );
    }

    return resultMesh;
  }

  Geometry finalResult;

  for ( auto &firstGeom : firstGeoms ) {

    {
      AllocTagScope prepTag( AllocSite::CsgOperandPrep );

      firstGeom.Cleanup();
    }

    Geometry result = firstGeom;

    for ( auto &secondGeom : secondGeoms ) {

      {
        AllocTagScope prepTag( AllocSite::CsgOperandPrep );

        result.Cleanup();
      }

      if ( result.IsEmpty() ) {

        Logger::logWarning("bool aborted due to empty first operand");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        break;
      }

      {
        AllocTagScope prepTag( AllocSite::CsgOperandPrep );

        secondGeom.Cleanup( !isSubtractOperand );
      }

      if ( secondGeom.IsEmpty() ) {
        
        Logger::logWarning("bool aborted due to empty second operand");

        // bail out because we will get strange meshes
        // if this happens, probably there's an issue parsing the mesh that
        // occurred earlier
        continue;
      }

      if ( secondGeom.halfSpace ) {

        glm::dvec3 origin = secondGeom.halfSpaceOrigin;
        glm::dvec3 x = secondGeom.halfSpaceX - origin;
        glm::dvec3 y = secondGeom.halfSpaceY - origin;
        glm::dvec3 z = secondGeom.halfSpaceZ - origin;
        glm::dmat4 trans =
            glm::dmat4(glm::dvec4(x, 0), glm::dvec4(y, 0), glm::dvec4(z, 0),
                        glm::dvec4(0, 0, 0, 1));
        Geometry newSecond;

        double scaleX = 1;
        double scaleY = 1;
        double scaleZ = 1;

        for (uint32_t i = 0; i < result.vertices.size(); i++) {

          const glm::dvec3& p = result.vertices[ i ];

          glm::dvec3 vec = ( p - origin );

          double dx = glm::dot( vec, x );
          double dy = glm::dot( vec, y );
          double dz = glm::dot( vec, z );

          if ( glm::abs( dx ) > scaleX ) {
            scaleX = glm::abs( dx );
          }

          if ( glm::abs( dy ) > scaleY ) {
            scaleY = glm::abs( dy );
          }

          if ( glm::abs( dz ) > scaleZ ) {
            scaleZ = glm::abs( dz );
          }
        }

        newSecond.AppendWithScalingTransform(
          secondGeom,
          trans,
          scaleX * 2,
          scaleY * 2,
          scaleZ * 2,
          secondGeom.halfSpaceOrigin );

        {
          AllocTagScope prepTag( AllocSite::CsgOperandPrep );

          newSecond.Cleanup( !isSubtractOperand );
        }

        Geometry newResult;

        {
          AllocTagScope kernelTag( AllocSite::CsgKernel );

          csg.run( CSG::Operation::SUBTRACTION, result, newSecond, newResult, 0 );
        }

        result = std::move( newResult );

      } else {

        Geometry newResult;

        {
          AllocTagScope kernelTag( AllocSite::CsgKernel );

          csg.run( CSG::Operation::SUBTRACTION, result, secondGeom, newResult, 0 );
        }

        result = std::move( newResult );
      }
    }

    if ( firstGeoms.size() == 1 ) {

      return result;
    }

    finalResult.AppendGeometry( result );
  }

  return finalResult;
}

Geometry ConwayGeometryProcessor::RelVoidSubtract(
      ParamsRelVoidSubtract& parameters) {
  Geometry resultGeometry;

  if (parameters.flatFirstMesh.size() <= 0) {
    Logger::logWarning("first mesh zero\n");
    return resultGeometry;
  }

  if (parameters.flatSecondMesh.size() <= 0) {
    Logger::logWarning("second mesh zero\n");
    return resultGeometry;
  }

  resultGeometry = BoolSubtract( parameters.flatFirstMesh, parameters.flatSecondMesh, false );

  glm::dmat4 combinedMatrix =
       glm::inverse(parameters.parentMatrix);// * glm::translate(originFirstMesh);

  resultGeometry.ApplyTransform( combinedMatrix );

  return resultGeometry;
}

Geometry ConwayGeometryProcessor::GetBlock(
    const ParamsGetBlock &parameters) {
  Geometry resultGeometry;

  if (parameters.xLength <= 0 || parameters.yLength <= 0 ||
      parameters.zLength <= 0) {
    Logger::logWarning("block size zero\n");
    return resultGeometry;
  }

  const double xLength = parameters.xLength;
  const double yLength = parameters.yLength;
  const double zLength = parameters.zLength;  
  glm::dmat4 placement = parameters.placement;

  //make vertices
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(xLength, yLength, zLength, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(0, yLength, zLength, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(0, 0, zLength, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(xLength, 0, zLength, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(xLength, yLength, 0, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(0, yLength, 0, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(0, 0, 0, 1)));
  resultGeometry.MakeVertex(glm::dvec3(placement * glm::dvec4(xLength, 0, 0, 1)));

  //make triangles
  resultGeometry.MakeTriangle(0, 1, 2);
  resultGeometry.MakeTriangle( 2, 3, 0 );

  resultGeometry.MakeTriangle( 7, 6, 5 );
  resultGeometry.MakeTriangle( 5, 4, 7 );

  resultGeometry.MakeTriangle( 0, 4, 5 );
  resultGeometry.MakeTriangle( 5, 1, 0 );

  resultGeometry.MakeTriangle( 1, 5, 6 );
  resultGeometry.MakeTriangle( 6, 2, 1 );

  resultGeometry.MakeTriangle( 2, 6, 7 );
  resultGeometry.MakeTriangle( 7, 3, 2 );

  resultGeometry.MakeTriangle( 3, 7, 4 );
  resultGeometry.MakeTriangle( 4, 0, 3 );

  if (MatrixFlipsTriangles(placement)) {
    resultGeometry.ReverseFaces();
  }

  return resultGeometry;
}

Geometry ConwayGeometryProcessor::GetBooleanResult(
    ParamsGetBooleanResult *parameters) {
  Geometry resultGeometry;
  if (parameters->flatFirstMesh.size() <= 0) {
    Logger::logWarning("first mesh zero\n");
    return resultGeometry;
  }

  if (parameters->flatSecondMesh.size() <= 0) {
    Logger::logWarning("second mesh zero\n");
    return resultGeometry;
  }

  resultGeometry =
    BoolSubtract(parameters->flatFirstMesh, parameters->flatSecondMesh, parameters->isSubtractOperand );

  return resultGeometry;
}

Geometry ConwayGeometryProcessor::GetHalfSpaceSolid(
    const ParamsGetHalfspaceSolid& parameters) {
  glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);

  if (parameters.flipWinding) {
    extrusionNormal *= -1;
  }

  double d = 1;

  IfcProfile profile;
  profile.isConvex = false;
  profile.curve = GetRectangleCurve(d, d, glm::dmat3(1));

  auto geom = Extrude(profile, extrusionNormal, d);
  geom.halfSpace = true;

  // @Refactor: duplicate of extrudedareasolid
  if ( parameters.flipWinding ) {
    
    geom.ReverseFaces();
  }

  return geom;
}

Geometry ConwayGeometryProcessor::GetPolygonalBoundedHalfspace(
    const ParamsGetPolygonalBoundedHalfspace& parameters) {

  glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);
  glm::dvec3 planeNormal = parameters.surface.transformation[2];
  glm::dvec3 planePosition = parameters.surface.transformation[3];


  glm::dmat4 invPosition = glm::inverse(parameters.position);
  glm::dvec3 localPlaneNormal = invPosition * glm::dvec4(planeNormal, 0);
  auto localPlanePos = invPosition * glm::dvec4(planePosition, 1);

  bool flipWinding = false;
  double extrudeDistance =
      EXTRUSION_DISTANCE_HALFSPACE_M / parameters.scaleFactor;
  bool halfSpaceInPlaneDirection = !parameters.agreement;
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

  conway::geometry::Geometry geom = 
      Extrude(
        profile,
        extrusionNormal,
        extrudeDistance,
        localPlaneNormal,
        localPlanePos );
  // auto geom = Extrude(
  //   profile,
  //   surface.transformation,
  //   extrusionNormal,
  //   EXTRUSION_DISTANCE_HALFSPACE);

  // @Refactor: duplicate of extrudedareasolid
  if (flipWinding) {
    
    geom.ReverseFaces();
  }

#ifdef DUMP_CSG_MESHES
  DumpIfcGeometry(geom, L"pbhs.obj");
#endif

  return geom;
}



namespace {

  /**
   * Name the surface a face was built on, for diagnostics.
   *
   * @param surface The face's surface.
   * @return A static, human-readable surface kind.
   */
  const char* surfaceKindName( const IfcSurface& surface ) {

    if ( surface.BSplineSurface.Active )   { return "b-spline"; }
    if ( surface.CylinderSurface.Active )  { return "cylindrical"; }
    if ( surface.SphericalSurface.Active ) { return "spherical"; }
    if ( surface.ToroidalSurface.Active )  { return "toroidal"; }
    if ( surface.ConicalSurface.Active )   { return "conical"; }
    if ( surface.RevolutionSurface.Active ){ return "revolution"; }
    if ( surface.ExtrusionSurface.Active ) { return "extrusion"; }

    return "planar/bounds";
  }

  /**
   * Warn that a face added nothing to its geometry.
   *
   * Reports the surface kind and bound count because the two causes want
   * different fixes: no bounds at all is an extraction-side gap, while
   * bounds that produced nothing is a triangulator bailing on geometry it
   * was handed.
   *
   * Reaches the host from pool threads as well as the main thread, but only
   * because FinalizeStagedFaces opens a Logger::DeferredScope around its
   * parallel_for. Without one, a bare EM_ASM evaluates in the calling
   * worker's own scope, where no sink is installed, and the warning vanishes
   * — non-deterministically, since parallel_for runs serially below its
   * threshold and the calling thread joins the work. Any future parallel
   * region that can log needs the same scope (conway#466).
   *
   * @param surfaceKind Human-readable surface kind.
   * @param boundCount How many bounds the face was given.
   */
  void warnFaceAddedNothing( const char* surfaceKind, size_t boundCount ) {

    Logger::logWarning(
      "Face contributed no geometry (%s surface, %zu bound(s))",
      surfaceKind,
      boundCount );
  }

  /**
   * Roll a geometry back to the state it was in before a face began
   * appending to it.
   *
   * A triangulator that throws part-way has usually already pushed some of
   * the face's vertices and triangles. Keeping them renders a torn fragment
   * of a face, which is worse than the missing face: it looks like geometry,
   * so nothing downstream treats it as absent. This is the same all-or-
   * nothing rule FinalizeStagedFaces applies with Geometry::Clear on its
   * per-face local geometry, expressed as a truncation because the immediate
   * path appends into a geometry that already holds earlier faces.
   *
   * Connectivity, the BVH and the reification are derived caches keyed on
   * the arrays being truncated, so they are dropped rather than fixed up.
   * None is built during face triangulation today; dropping them keeps that
   * from becoming a silent correctness bug if one ever is. The reification
   * matters most: GetVertexDataSize/GetIndexDataSize report its lengths
   * straight to the embind callers, which turn them into typed-array views
   * over the wasm heap, so a stale one is not a wrong picture but a
   * RangeError or an out-of-bounds read.
   *
   * Truncation is sound because MakeVertex is a plain push_back — there is
   * no vertex weld map to desync — and triangles only ever reference
   * vertices at or below their own index.
   *
   * @param geometry      Geometry to truncate.
   * @param vertexCount   vertices.size() before the face ran.
   * @param triangleCount triangles.size() before the face ran.
   */
  void rollbackFace(
    Geometry& geometry,
    size_t    vertexCount,
    size_t    triangleCount ) {

    if ( geometry.triangles.size() > triangleCount ) {
      geometry.triangles.resize( triangleCount );
    }

    if ( geometry.vertices.size() > vertexCount ) {
      geometry.vertices.resize( vertexCount );
    }

    geometry.ClearConnectivity();
    geometry.ClearReification();
    geometry.bvh.reset();
  }
}

void ConwayGeometryProcessor::AddFaceToGeometrySimple(
  ParamsAddFaceToGeometrySimple& parameters, Geometry &geometry) {
  AllocTelemetryScope telemetryScope( AllocSite::AdvancedFace );

  // Instrumented for the same reason as the advanced path: this one drops
  // just as silently, and it is the one IFC calls directly, so leaving it
  // out would have made coverage depend on which mode a model happened to
  // take — with today's default reporting less than the staged path.
  const size_t trianglesBefore = geometry.triangles.size();
  const size_t verticesBefore  = geometry.vertices.size();
  const size_t boundCount = static_cast< size_t >( parameters.boundsArray.size() );

  // Guarded for the same reason as the advanced path, which documents it:
  // this is also bound directly (`addFaceToGeometrySimple`) and also called
  // by the IFC front end, so an unguarded throw here costs the whole product.
  try {

    if (parameters.boundsArray.size() > 0) {

      TriangulateBounds(geometry, parameters.boundsArray);
    }

  } catch ( const std::exception& error ) {

    rollbackFace( geometry, verticesBefore, trianglesBefore );

    Logger::logError(
      "Face threw while triangulating (simple/bounds surface, %zu bound(s)): %s",
      boundCount,
      error.what() );

    return;

  } catch ( ... ) {

    rollbackFace( geometry, verticesBefore, trianglesBefore );

    Logger::logError(
      "Face threw an unknown exception while triangulating "
      "(simple/bounds surface, %zu bound(s))",
      boundCount );

    return;
  }

  if ( geometry.triangles.size() == trianglesBefore ) {

    warnFaceAddedNothing( "simple/bounds", boundCount );
  }
}


void ConwayGeometryProcessor::AddFaceToGeometry(
    ParamsAddFaceToGeometry& parameters, Geometry &geometry) {
  AllocTelemetryScope telemetryScope( AllocSite::AdvancedFace );

  // A face that contributes nothing is recoverable information, and until
  // now it was thrown away: every triangulator returns void, and one that
  // bails leaves a model quietly missing a face — or, when the face was
  // the whole body, missing entirely. A whole-sphere STEP loaded as zero
  // meshes with not one line of output, and a torus lost exactly half its
  // surface the same way; both took a 68-part corpus to notice rather
  // than a single load. See https://github.com/bldrs-ai/conway/issues/461.
  const size_t trianglesBefore = geometry.triangles.size();
  const size_t verticesBefore  = geometry.vertices.size();

  // A face that THROWS is the other half of that gap, and until now it cost
  // far more than the face: this is the path IFC uses (`addFaceToGeometry`
  // binds here; the staged path is not on the IFC front end), it had no
  // guard, and the exception propagated out through embind to the extractor's
  // per-product catch — so one bad face lost the entire product, reported as
  // a bare "CDT Exception at Triangulation.hpp:983" with nothing naming the
  // surface or even the face. ISSUE_159 shipped that way for months
  // (conway#473). Catch per face, name the surface kind and bound count the
  // way the added-nothing path does, and let the rest of the product build.
  //
  // Deliberately swallowed, not rethrown. Losing one face of a product beats
  // losing the product, and the loss is not silent: every occurrence logs an
  // error, which reaches errors.csv in the regression harness, and a model
  // that loses ALL its faces this way lands in zero_geometry.csv (conway#478).
  const char* const surfaceKind =
    parameters.advancedBrep ? surfaceKindName( parameters.surface ) : "non-advanced";
  const size_t boundCount = static_cast< size_t >( parameters.boundsArray.size() );

  try {

    if (!parameters.advancedBrep) {
      if (parameters.boundsArray.size() > 0) {
        conway::AllocTagScope tag( conway::AllocSite::TriBounds );
        TriangulateBounds(geometry, parameters.boundsArray);
      }
    } else {
      if (parameters.boundsArray.size() > 0) {
        auto surface = parameters.surface;

        if (surface.BSplineSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriBspline );
          TriangulateBspline(geometry, parameters.boundsArray, surface,
                             parameters.scaling, parameters.representationExtent);
        } else if (surface.CylinderSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriCylinder );
          TriangulateCylindricalSurface(geometry, parameters.boundsArray,
                                        surface, parameters.representationExtent);
        } else if (surface.SphericalSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriSphere );
          TriangulateSphericalSurface(geometry, parameters.boundsArray,
                                      surface, parameters.representationExtent);
        } else if (surface.ToroidalSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriToroidal );
          TriangulateToroidalSurface(geometry, parameters.boundsArray,
                                     surface, parameters.representationExtent);
        } else if (surface.ConicalSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriConical );
          TriangulateConicalSurface(geometry, parameters.boundsArray, surface,
                                    parameters.representationExtent);
        } else if (surface.RevolutionSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriRevolution );
          TriangulateRevolution(geometry, parameters.boundsArray, surface,
                                parameters.representationExtent);
        } else if (surface.ExtrusionSurface.Active) {
          conway::AllocTagScope tag( conway::AllocSite::TriExtrusion );
          TriangulateExtrusion(geometry, parameters.boundsArray, surface,
                               parameters.representationExtent);
        } else {
          conway::AllocTagScope tag( conway::AllocSite::TriBounds );
          TriangulateBounds(geometry, parameters.boundsArray);
        }
      }
    }

  } catch ( const std::exception& error ) {

    rollbackFace( geometry, verticesBefore, trianglesBefore );

    Logger::logError(
      "Face threw while triangulating (%s surface, %zu bound(s)): %s",
      surfaceKind,
      boundCount,
      error.what() );

    return;

  } catch ( ... ) {

    rollbackFace( geometry, verticesBefore, trianglesBefore );

    Logger::logError(
      "Face threw an unknown exception while triangulating "
      "(%s surface, %zu bound(s))",
      surfaceKind,
      boundCount );

    return;
  }

  if ( geometry.triangles.size() == trianglesBefore ) {

    warnFaceAddedNothing( surfaceKind, boundCount );
  }
}


namespace {

  /**
   * Auto-flush threshold for staged face jobs. Large enough to keep the
   * thread pool saturated per flush, small enough to bound the transient
   * memory held by staged bounds/surfaces on big models.
   */
  constexpr size_t STAGED_FACE_FLUSH_THRESHOLD = 4096;
}

void ConwayGeometryProcessor::StageFaceToGeometry(
    ParamsAddFaceToGeometry& parameters, Geometry& geometry) {

  stagedFaceJobs_.push_back(
      StagedFaceJob{ std::move( parameters ), &geometry } );

  if ( stagedFaceJobs_.size() >= STAGED_FACE_FLUSH_THRESHOLD ) {

    FinalizeStagedFaces();
  }
}

void ConwayGeometryProcessor::StageFaceToGeometrySimple(
    ParamsAddFaceToGeometrySimple& parameters, Geometry& geometry) {

  ParamsAddFaceToGeometry promoted;

  promoted.boundsArray  = std::move( parameters.boundsArray );
  promoted.advancedBrep = false;
  promoted.scaling      = parameters.scaling;
  promoted.representationExtent = parameters.representationExtent;

  // Delegate so the queue append and flush policy live in one place.
  StageFaceToGeometry( promoted, geometry );
}

void ConwayGeometryProcessor::FinalizeStagedFaces() {

  size_t jobCount = stagedFaceJobs_.size();

  if ( jobCount == 0 ) {
    return;
  }

  // Tessellate every staged face into its own local geometry in parallel,
  // then append the results in staging order so output ordering (and
  // therefore output data) matches the immediate AddFaceToGeometry path.
  std::vector< Geometry > results( jobCount );

  // Everything logged below happens on pool threads, where a bare EM_ASM
  // evaluates in the worker's own scope and is silently dropped - including
  // the catch handlers here, which are the only report a caller gets that a
  // face threw. Buffer for the duration and flush on this thread when the
  // scope closes. See conway#466.
  Logger::DeferredScope deferLogs;

  ThreadPool::instance().parallel_for( 0, jobCount, [&]( size_t where ) {

    StagedFaceJob& job = stagedFaceJobs_[ where ];

    try {

      AddFaceToGeometry( job.parameters, results[ where ] );

    } catch ( const std::exception& e ) {

      results[ where ].Clear();
      Logger::logError( "Exception tessellating staged face: %s", e.what() );

    } catch ( ... ) {

      results[ where ].Clear();
      Logger::logError( "Unknown exception tessellating staged face." );
    }
  } );

  for ( size_t where = 0; where < jobCount; ++where ) {

    const Geometry& result = results[ where ];

    if ( result.vertices.empty() && result.triangles.empty() ) {
      continue;
    }

    stagedFaceJobs_[ where ].target->AppendGeometry( result );
  }

  stagedFaceJobs_.clear();
}

IfcSurface ConwayGeometryProcessor::GetSurface(const ParamsGetSurface& parameters) {
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

    surface.BSplineSurface.UKnots = parameters.UKnots;
    surface.BSplineSurface.VKnots = parameters.VKnots;

    if (surface.BSplineSurface.UKnots[parameters.UKnots.size() - 1] !=
        (int)surface.BSplineSurface.UKnots[surface.BSplineSurface.UKnots.size() - 1]) {
      for (uint32_t i = 0; i < surface.BSplineSurface.UKnots.size(); i++) {
        surface.BSplineSurface.UKnots[i] =
          surface.BSplineSurface.UKnots[i] *
          (surface.BSplineSurface.UKnots.size() - 1) /
          surface.BSplineSurface.UKnots[surface.BSplineSurface.UKnots.size() - 1];
      }
    }

    if (surface.BSplineSurface.VKnots[surface.BSplineSurface.VKnots.size() - 1] !=
        (int)surface.BSplineSurface.VKnots[surface.BSplineSurface.VKnots.size() - 1]) {
      for (uint32_t i = 0; i < surface.BSplineSurface.VKnots.size(); i++) {
        surface.BSplineSurface.VKnots[i] =
          surface.BSplineSurface.VKnots[i] *
          (surface.BSplineSurface.VKnots.size() - 1) /
          surface.BSplineSurface.VKnots[surface.BSplineSurface.VKnots.size() - 1];
      }
    }

    surface.BSplineSurface.Active = true;
    surface.BSplineSurface.UDegree = parameters.Udegree;
    surface.BSplineSurface.VDegree = parameters.Vdegree;
    surface.BSplineSurface.ControlPoints = parameters.ctrolPts;
    surface.BSplineSurface.UMultiplicity = parameters.UMultiplicity;
    surface.BSplineSurface.VMultiplicity = parameters.VMultiplicity;

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
  } else if (parameters.isSphericalSurface) {
    IfcSurface surface;

    surface.transformation = parameters.transformation;
    surface.SphericalSurface.Active = true;
    surface.SphericalSurface.Radius = parameters.radius;

    return surface; 
  } else if (parameters.isToroidalSurface ) {
    IfcSurface surface;

    surface.transformation = parameters.transformation;
    surface.ToroidalSurface.Active = true;
    surface.ToroidalSurface.MajorRadius = parameters.radius;
    surface.ToroidalSurface.MinorRadius = parameters.radius2;

    return surface;
  } else if (parameters.isConicalSurface) {
    IfcSurface surface;

    surface.transformation = parameters.transformation;
    surface.ConicalSurface.Active = true;
    surface.ConicalSurface.Radius = parameters.radius;
    surface.ConicalSurface.SemiAngle = parameters.semiAngle;

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
  const double *points = reinterpret_cast<double *>(parameters.points_);

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
    const ParamsGetAxis2Placement2D& parameters) {
  if (parameters.isAxis2Placement2D) {
    glm::dvec2 xAxis = glm::dvec2(1, 0);
    if (parameters.customAxis1Ref) {
      xAxis = glm::normalize(parameters.axis1Ref);
    }

    glm::dvec2 pos = parameters.position2D;

    glm::dvec2 yAxis = glm::normalize(glm::dvec2(-xAxis.y, xAxis.x));

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
    const ParamsAxis1Placement3D& parameters) {
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

IfcCurve ConwayGeometryProcessor::GetLoop(const ParamsGetLoop& parameters) {
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
      //    if (notPresent(pt, curve.points)) {
          if ( curve.Add3d(pt) ) {
            curve.indices.push_back(id);
          }
      //    }
          // }
        }
      }
      id++;
    }
  }

  return curve;
}

IfcBound3D ConwayGeometryProcessor::GetBound(const ParamsGetBound& parameters) {
  bool orient = parameters.orient;

  IfcBound3D bound;
  bound.curve = GetLoop(parameters.parametersGetLoop);
  bound.orientation = orient;
  bound.type = (parameters.isFaceOuterBound) ? IfcBoundType::OUTERBOUND
                                             : IfcBoundType::BOUND;

  if (!orient) {
    std::reverse(bound.curve.points.begin(), bound.curve.points.end());
  }

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
      if (i + 1 >= parameters.face.face_starts.size()) {
        endIdx = parameters.face.indices.size();
      } else {
        endIdx = parameters.face.face_starts[i + 1];
      }

      for (size_t index = startIdx; index < endIdx; index++) {
        uint32_t currentIndex = parameters.face.indices[index];

        glm::dvec3 point = parameters.points[currentIndex - 1];

        // I am not proud of this (I inherited this, will change - NC)
        bounds.back().curve.points.push_back(point);
      }
    }
  }

  return bounds;
}

conway::geometry::ConwayGeometryProcessor::ResultsGltf
ConwayGeometryProcessor::GeometryToGltf(
    std::span<conway::geometry::GeometryCollection> geoms,
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

        default:

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
        std::filesystem::path(filePath).parent_path());
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
      results.bufferUris.push_back( filePath + "." + Microsoft::glTF::BUFFER_EXTENSION );
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

    bool hasFirstPoint = false;

    glm::dvec3 positionBias;

    for (conway::geometry::GeometryCollection &geom : geoms) {

      if ( geom.components.empty() ) {
        continue;
      }

      bool hasTriangles = false;

      for ( conway::geometry::Geometry *componentPtr : geom.components ) {
        conway::geometry::Geometry &component = *componentPtr;

        if ( !component.triangles.empty() ) {
          hasTriangles = true;
          break;
        }
      }

      if ( !hasTriangles ) {
        continue;
      }

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
      // this internally populates the vertex float array, current storage type
      // is double
      // geom.GetVertexData();

      if (outputDraco) {
        dracoMesh.reset(new draco::Mesh());

        // set the number of positions and faces
        uint32_t numPositions = 0;
        uint32_t numFaces = 0;

        for ( conway::geometry::Geometry *componentPtr : geom.components ) {

          conway::geometry::Geometry &component = *componentPtr;

          numPositions += component.vertices.size();
          numFaces     += component.triangles.size();
        }

        uint32_t numIndices = numFaces * 3;

        if (numFaces > 0) {
          // set number of faces
          dracoMesh->SetNumFaces(numFaces);

          // set number of indices
          dracoMesh->set_num_points(numIndices);
        }

        // Add attributes if they are present in the input data.
        if (numPositions > 0) {
          draco::GeometryAttribute va;

          va.set_unique_id(uniqueIdDraco++);

          va.Init(draco::GeometryAttribute::POSITION, nullptr, 3,
                  draco::DT_FLOAT32, false, sizeof(float) * 3, 0);
          pos_att_id = dracoMesh->AddAttribute(va, false, numPositions);
        }

        // populate position attribute

        uint32_t vertexCount = 0;
        uint32_t totalPoints = 0;
        uint32_t totalIndices = 0;

        // populate geometry components position attribute
        for (uint32_t geometryComponentIndex = 0;
             geometryComponentIndex < geom.components.size();
             ++geometryComponentIndex) {
          Geometry &component = *geom.components[geometryComponentIndex];
          glm::dmat4 geomTransform =
              transform * geom.transforms[geometryComponentIndex];

          for ( const glm::dvec3 *where = component.vertices.data(),
                            *end = where + component.vertices.size();
               where < end; ++where ) {
            glm::dvec3 t = glm::dvec3(
                geomTransform * glm::dvec4( *where, 1) );

            if ( !hasFirstPoint ) {

              hasFirstPoint = true;
              positionBias = t;
            }

            t = t - positionBias;

            float vertexVal[3] = {static_cast<float>(t.x),
                                  static_cast<float>(t.y),
                                  static_cast<float>(t.z)};

            dracoMesh->attribute( pos_att_id )
                ->SetAttributeValue( draco::AttributeValueIndex( vertexCount++ ),
                                    vertexVal );
          }

          const Triangle *triangles = component.triangles.data();

          // no textures, just map vertices to face indices
          for (
            size_t triangleIndex = 0, end = component.triangles.size();
            triangleIndex < end; 
            ++triangleIndex ) {

            const Triangle& triangle = triangles[ triangleIndex ];

            uint32_t triangle0 = triangle.vertices[ 0 ] + totalPoints;
            uint32_t triangle1 = triangle.vertices[ 1 ] + totalPoints;
            uint32_t triangle2 = triangle.vertices[ 2 ] + totalPoints;

            uint32_t compositeIndiceIndex = totalIndices + triangleIndex * 3;

            const draco::PointIndex vert_id_0(compositeIndiceIndex);
            const draco::PointIndex vert_id_1(compositeIndiceIndex + 1);
            const draco::PointIndex vert_id_2(compositeIndiceIndex + 2);

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

          totalIndices += static_cast< uint32_t >( component.triangles.size() * 3 );
          totalPoints  += static_cast< uint32_t >( component.vertices.size() );
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

      for (conway::geometry::Geometry *componentPtr : geom.components) {
        conway::geometry::Geometry &component = *componentPtr;

        numPoints  += static_cast< uint32_t >( component.vertices.size() );
        numIndices += component.triangles.size() * 3;
      }

      std::vector<float> positions;
      std::vector<uint32_t> indexData;

      positions.resize( numPoints * 3 );
      indexData.reserve( numIndices );

      std::vector<float> minValues( 3U, std::numeric_limits< float >::max() );
      std::vector<float> maxValues( 3U, std::numeric_limits< float >::lowest() );

      const size_t positionCount = positions.size();
      float *positionOutputCursor = positions.data();

      uint32_t pointOffset = 0;
      size_t transformCursor = 0;

      for (conway::geometry::Geometry *componentPtr : geom.components) {
        conway::geometry::Geometry &component = *componentPtr;
        glm::dmat4 geomTransform =
            transform * geom.transforms[transformCursor++];

        for (const glm::dvec3 *where = component.vertices.data(),
                              *end   = where + component.vertices.size();
             where < end;
             ++where ) {

          glm::dvec3 t =
              glm::dvec3( geomTransform * glm::dvec4( *where, 1 ) );

          if ( !hasFirstPoint ) {

            hasFirstPoint = true;
            positionBias = glm::dvec3( t );
          }

          t -= positionBias;

          *(positionOutputCursor++) = (float)t.x;
          *(positionOutputCursor++) = (float)t.y;
          *(positionOutputCursor++) = (float)t.z;
        }

        size_t indexDataOffset = indexData.size();

        indexData.insert(
          indexData.end(),
          &component.triangles[ 0 ].vertices[ 0 ],
          &component.triangles[ component.triangles.size() ].vertices[ 0 ] );

        for (uint32_t *where = indexData.data() + indexDataOffset,
                      *end = indexData.data() + indexData.size();
             where < end; where += 3) {
          where[0] += pointOffset;
          where[1] += pointOffset;
          where[2] += pointOffset;
        }

        pointOffset += component.vertices.size();
      }

      // Accessor min/max properties must be set for vertex position data so
      // calculate them here. Non-finite values (degenerate upstream geometry)
      // must be skipped: a NaN here poisons the accessor bounds and makes the
      // JSON serializer throw mid-manifest, truncating the whole GLB.
      for (size_t i = 0U, j = 0U; i < positionCount; ++i, j = (i % 3U)) {
        if (!std::isfinite(positions[i])) {
          continue;
        }
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

    // Construct a Node adding a reference to the Mesh. The vertex data was
    // rebased onto the first transformed vertex (positionBias) for float32
    // precision, so that offset must come back as the node translation —
    // otherwise every chunked GLB is displaced by its own (content-dependent)
    // first vertex and multi-chunk models no longer assemble.
    Microsoft::glTF::Node node;
    node.meshId = meshId;

    if (hasFirstPoint) {
      node.translation = Microsoft::glTF::Vector3(
          static_cast<float>(positionBias.x),
          static_cast<float>(positionBias.y),
          static_cast<float>(positionBias.z));
    }

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
      results.buffers.push_back( glbBuffer );

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
    const conway::geometry::Geometry &geom, size_t &offset,
    glm::dmat4 transform ) {
  std::string obj;
  obj.reserve( geom.vertices.size() * 32 + geom.triangles.size() * 32);  // preallocate memory

  const char *vFormat = "v %.6f %.6f %.6f\n";
  const char *fFormat = "f %zu// %zu// %zu//\n";

  for (uint32_t i = 0; i < geom.vertices.size(); ++i) {
    glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
    char vBuffer[64];
    snprintf(vBuffer, sizeof(vBuffer), vFormat, t.x, t.y, t.z);
    obj.append(vBuffer);
  }

  for ( const Triangle& triangle : geom.triangles ) {


    size_t f1 = triangle.vertices[ 0 ] + 1 + offset;
    size_t f2 = triangle.vertices[ 1 ] + 1 + offset;
    size_t f3 = triangle.vertices[ 2 ] + 1 + offset;
    char fBuffer[64];
    snprintf(fBuffer, sizeof(fBuffer), fFormat, f1, f2, f3);
    obj.append(fBuffer);
  }

  offset += geom.vertices.size();

  return obj;
}

Geometry ConwayGeometryProcessor::getTriangulatedFaceSetGeometry(
    const ParamsGetTriangulatedFaceSetGeometry &parameters) {
  const uint32_t *indices =
      reinterpret_cast<uint32_t *>(parameters.indicesArray_);

  Geometry geom;

  for ( const double* points = reinterpret_cast<double *>( parameters.pointsArray_ ),
                    * end    = reinterpret_cast<double *>( parameters.pointsArray_ ) + parameters.pointsArrayLength;
        points < end;
        points += 3 ) {

    geom.MakeVertex( points[ 0 ], points[ 1 ], points[ 2 ] );
  }

  for (size_t i = 0; i < parameters.indicesArrayLength; i += 3) {

    uint32_t i1 = ( indices[ i + 0 ] - 1 ) * 3;
    uint32_t i2 = ( indices[ i + 1 ] - 1 ) * 3;
    uint32_t i3 = ( indices[ i + 2 ] - 1 ) * 3;

    geom.MakeTriangle( indices[ i ] - 1, indices[ i + 1 ] - 1, indices[ i + 2 ] - 1 );
    // printf("adding face %i: x: %.3f, y: %.3f, z: %.3f\n", i, points[i1],
    // points[i2], points[i3]);
  }

  return geom;
}

Geometry ConwayGeometryProcessor::getPolygonalFaceSetGeometry(
    const ParamsGetPolygonalFaceSetGeometry &parameters) {
  Geometry geom;
  std::vector<IfcBound3D> bounds;

  for (size_t faceIndex = 0; faceIndex < parameters.faces.size(); faceIndex++) {
    ParamsReadIndexedPolygonalFace params(parameters.points,
                                          parameters.faces[faceIndex]);

    bounds = ReadIndexedPolygonalFace(params);

    // Unlike the staged AddFaceToGeometry path (whose FinalizeStagedFaces
    // wraps every face in a per-face try/catch), this unstaged loop calls
    // TriangulateBounds directly. Guard it so a single degenerate face — e.g.
    // the non-finite planar-projection case the tesselatePlane guard throws on
    // — is skipped rather than aborting the whole polygonal face set (and
    // propagating out across the embind boundary). This matches tesselatePlane's
    // own per-face CDT-error returns; byte-identical on any face that already
    // triangulates.
    try {
      TriangulateBounds(geom, bounds);
    } catch (const std::exception& e) {
      Logger::logError("Skipping polygonal face in face set: %s", e.what());
    }

    bounds.clear();
  }

  return geom;
}

Geometry ConwayGeometryProcessor::getPolygonalFaceSetGeometryPacked(
    const std::vector< glm::dvec3 >& points,
    const uint32_t* indices,
    const uint32_t* faceOffsets,
    uint32_t faceOffsetsLength,
    const uint32_t* startIndices,
    const uint32_t* startOffsets ) {

  Geometry geom;
  std::vector< IfcBound3D > bounds;
  bool loggedBadIndex = false;

  if ( faceOffsetsLength == 0 ) {
    return geom;
  }

  // Same contract as buildIndexedPolygonalFaceVector: the caller hands
  // the offset array's length (N+1). Deriving N here makes passing the
  // JS .length the only representable choice.
  const uint32_t faceCount = faceOffsetsLength - 1;
  const uint32_t pointCount = static_cast< uint32_t >( points.size() );

  for ( uint32_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {

    const uint32_t indicesBegin = faceOffsets[ faceIndex ];
    const uint32_t indicesEnd   = faceOffsets[ faceIndex + 1 ];
    const uint32_t startsBegin  = startOffsets[ faceIndex ];
    const uint32_t startsEnd    = startOffsets[ faceIndex + 1 ];
    const uint32_t startCount   = startsEnd - startsBegin;

    if ( startCount == 0 ) {
      continue;
    }

    // Keep the IfcBound3D objects (and their points-vector capacity)
    // across faces. bounds.clear() would ~IfcBound3D each face and
    // force a fresh allocation chain on the next — on PSB that's ~9.1M
    // point-vector allocations, the same order the packed path exists
    // to retire. resize() down after a rare voided face drops the extra
    // slots; the 1-bound case that dominates is a no-op + points.clear().
    // TriangulateBounds may swap slots, so reuse is capacity-only.
    bounds.resize( startCount );

    for ( uint32_t ring = 0; ring < startCount; ++ring ) {
      bounds[ ring ].curve.points.clear();
    }

    // Same layout ReadIndexedPolygonalFace walks: face_starts[0] is 0
    // for a plain face (the whole index slice), and extra starts mark
    // void rings. Indices are 1-based into `points`.
    const uint32_t outerEnd = startCount > 1 ?
        startIndices[ startsBegin + 1 ] : ( indicesEnd - indicesBegin );

    bool faceOk = true;

    auto appendIndex = [ & ]( uint32_t index, uint32_t ring ) -> bool {

      if ( index == 0 || index > pointCount ) {

        if ( !loggedBadIndex ) {

          Logger::logWarning(
              "Skipping polygonal face with CoordIndex %u outside 1..%u",
              index,
              pointCount );
          loggedBadIndex = true;
        }

        return false;
      }

      bounds[ ring ].curve.points.push_back( points[ index - 1 ] );
      return true;
    };

    for ( uint32_t where = 0; where < outerEnd; ++where ) {

      if ( !appendIndex( indices[ indicesBegin + where ], 0 ) ) {
        faceOk = false;
        break;
      }
    }

    for ( uint32_t ring = 1; faceOk && ring < startCount; ++ring ) {

      const uint32_t ringBegin = startIndices[ startsBegin + ring ];
      const uint32_t ringEnd   = ring + 1 < startCount ?
          startIndices[ startsBegin + ring + 1 ] : ( indicesEnd - indicesBegin );

      for ( uint32_t where = ringBegin; where < ringEnd; ++where ) {

        if ( !appendIndex( indices[ indicesBegin + where ], ring ) ) {
          faceOk = false;
          break;
        }
      }
    }

    if ( !faceOk ) {
      continue;
    }

    try {

      TriangulateBounds( geom, bounds );

    } catch ( const std::exception& e ) {

      Logger::logError( "Skipping polygonal face in face set: %s", e.what() );
    }
  }

  return geom;
}

conway::geometry::Geometry ConwayGeometryProcessor::getSweptDiskSolid(
  const ParamsGetSweptDiskSolid &parameters) {

    IfcCurve directrix = parameters.directrix;

    IfcProfile profile;

    profile.curve = GetCircleCurve(parameters.radius, parameters.circleSegments);

    Geometry geom = SweepCircular( parameters.scalingFactor, parameters.closed, profile, parameters.radius, directrix );

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
                          pts[sg.indices[2] - 1]);

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
  }

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getIndexedPolyCurve3D(
  const ParamsGetIfcIndexedPolyCurve3D &parameters) {
IfcCurve curve;

if (parameters.dimensions == 3) {
  if (parameters.segments.size() > 0) {
    for (auto &sg : parameters.segments) {
      if (!sg.isArcType) {
        auto pts = parameters.points;
        for (auto &pt : sg.indices) {
          curve.Add3d(pts[pt - 1]);
        }
      }
      if (sg.isArcType) {

        auto pts = parameters.points;

        IfcCurve arc =
            Build3DArc3Pt(pts[sg.indices[0] - 1], pts[sg.indices[1] - 1],
                        pts[sg.indices[2] - 1]);

        for (auto &pt : arc.points) {
          curve.Add3d(pt);
        }
      }
    }
  } else {
    auto pts = parameters.points;
    for (auto &pt : pts) {
      curve.Add3d(pt);
    }
  }
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

  int segments = static_cast< int >( ceil( radius * CIRCLE_SEGMENT_TO_RADIUS_RATIO ) );

  segments = std::max( segments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_LOW );
  segments = std::min( segments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_HIGH );

 // printf( "Segments %d\n", segments );

  curve = GetCircleCurve( radius, segments, placement );

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

  int segments = static_cast< int >( ceil( radius * CIRCLE_SEGMENT_TO_RADIUS_RATIO ) );

  segments = std::max( segments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_LOW );
  segments = std::min( segments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_HIGH );

 // printf( "Segments %d\n", segments );

  curve = GetCircleCurve( radius, segments, placement );

  return curve;
}


enum IfcTrimmingPreference { CARTESIAN = 0, PARAMETER = 1, UNSPECIFIED = 2 };

conway::geometry::IfcCurve ConwayGeometryProcessor::getBSplineCurve(
    const ParamsGetBSplineCurve &parameters) {
  IfcCurve curve;

  bool condition = !parameters.paramsGetIfcTrimmedCurve.senseAgreement;

  curve.points.clear();

  int dimensions = parameters.dimensions;

  int degree = parameters.degree;

  int32_t startOffset = 0;
  int32_t  endOffset   = 0;

  if (dimensions == 2) {
  
    double tStart = 0;
    double tEnd   = 1;

    if ( parameters.paramsGetIfcTrimmedCurve.trimExists ) {

      if ( parameters.paramsGetIfcTrimmedCurve.masterRepresentation ==
        IfcTrimmingPreference::PARAMETER ) {
        if (condition) {
          tStart = parameters.paramsGetIfcTrimmedCurve.trim1Double;
          tEnd   = parameters.paramsGetIfcTrimmedCurve.trim2Double;
        } else {
          tStart = parameters.paramsGetIfcTrimmedCurve.trim2Double;
          tEnd = parameters.paramsGetIfcTrimmedCurve.trim1Double;
        }

      } else {

        // InterpolateRationalBSplineCurveWithKnots takes a NORMALISED
        // parameter: its first act is tPrime = t * (high - low) + low, where
        // low/high are the ends of the curve's own knot domain. So the
        // interval to search for a cartesian trim is [0, 1], always — and
        // never the raw knot values, which is what this used to read.
        //
        // Reading the knot array here was correct only for an exporter that
        // happened to normalise its knots to [0, 1]; anything else searched a
        // wrong sub-interval, and a domain lying entirely outside [0, 1]
        // collapsed to a point. The clamp on the 3D branch turned that
        // collapse into a silent one rather than an out-of-range bisection —
        // it arrived with the trimming feature itself (d5ecd77) and was never
        // a fix for a specific model. On Orbiter_v1.1_Gear_7.5.step the thread
        // helices carry domains like [-1.418, -0.473] and [59.691, 65.974], so
        // tMin and tMax both clamped to 0 (or both to 1) and the edge came back
        // as its two endpoints — a chord where a helical arc belongs, and the
        // cause of both failure modes in bldrs-ai/conway#594.
        constexpr double tMin = 0.0;
        constexpr double tMax = 1.0;

        tStart = best_fit_param_bisection( [&]( double value ) {

          glm::dvec2 pointOnCurve = InterpolateRationalBSplineCurveWithKnots( 
            value,
            degree,
            parameters.points2,
            parameters.knots,
            parameters.weights ) - parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D;

          return glm::dot( pointOnCurve, pointOnCurve );

        }, tMin, tMax );

        tEnd = best_fit_param_bisection( [&]( double value ) {

          glm::dvec2 pointOnCurve = InterpolateRationalBSplineCurveWithKnots( 
            value,
            degree,
            parameters.points2,
            parameters.knots,
            parameters.weights ) - parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D;

          return glm::dot( pointOnCurve, pointOnCurve );

        }, tMin, tMax );

        startOffset = 1;
        endOffset   = -1;
        
        curve.Add2d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D );
      }
    }

    std::vector<glm::dvec2> tempPoints = GetRationalBSplineCurveWithKnots(
        degree, parameters.points2, parameters.knots, parameters.weights, tStart, tEnd);
    for (
      int32_t i = startOffset, end = static_cast< int32_t >( tempPoints.size() ) + endOffset;
      i < end;
      ++i ) {

      curve.Add2d(tempPoints[i]);
    }

    if ( endOffset != 0 ) {
      curve.Add2d( parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D );
    }

  } else if (dimensions == 3) {
        
    double tStart = 0;
    double tEnd   = 1;

    if ( parameters.paramsGetIfcTrimmedCurve.trimExists ) {

      if ( parameters.paramsGetIfcTrimmedCurve.masterRepresentation ==
        IfcTrimmingPreference::PARAMETER ) {
        if (condition) {
          tStart = parameters.paramsGetIfcTrimmedCurve.trim1Double;
          tEnd   = parameters.paramsGetIfcTrimmedCurve.trim2Double;
        } else {
          tStart = parameters.paramsGetIfcTrimmedCurve.trim2Double;
          tEnd = parameters.paramsGetIfcTrimmedCurve.trim1Double;
        }

      } else {

        // The normalised trim-search interval — see the 2D branch above for
        // why it is [0, 1] and not the curve's knot values. This branch also
        // used to clamp those values into [0, 1], which is what turned an
        // out-of-domain knot vector into a silently collapsed edge rather than
        // an out-of-range bisection (bldrs-ai/conway#594).
        constexpr double tMin = 0.0;
        constexpr double tMax = 1.0;

        tStart = best_fit_param_bisection( [&]( double value ) {

          glm::dvec3 pointOnCurve = InterpolateRationalBSplineCurveWithKnots( 
            value,
            degree,
            parameters.points3,
            parameters.knots,
            parameters.weights ) - parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D;

          return glm::dot( pointOnCurve, pointOnCurve );

        }, tMin, tMax );

        // dvec3, not dvec2. This is the 3D branch: the interpolation is over
        // points3 and the target is trim2Cartesian3D, so declaring the
        // difference as a dvec2 silently dropped z through GLM's implicit
        // vec2(vec3) conversion and minimised XY distance alone. The tStart
        // lambda above always had this right; only tEnd was truncated.
        //
        // It was dormant while this branch clamped tMin/tMax onto a single
        // point - the bisection had nothing to search. Searching [0, 1] makes
        // it decide the trim, and it decides it worst on exactly the geometry
        // this fix is for: the XY projection of a helix is a circle walked
        // once per turn, so every turn scores ~0 and the solve returns an
        // arbitrary one. See bldrs-ai/conway#594.
        tEnd = best_fit_param_bisection( [&]( double value ) {

          glm::dvec3 pointOnCurve = InterpolateRationalBSplineCurveWithKnots( 
            value,
            degree,
            parameters.points3,
            parameters.knots,
            parameters.weights ) - parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D;

          return glm::dot( pointOnCurve, pointOnCurve );

        }, tMin, tMax );
                
        startOffset = 1;
        endOffset   = -1;

        curve.Add3d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D );
      }
    }
    
    std::vector<glm::dvec3> tempPoints = GetRationalBSplineCurveWithKnots(
        degree, parameters.points3, parameters.knots, parameters.weights, tStart, tEnd );
    
    // printf( "Output points: %zu\n", tempPoints.size() );

    for (
      int32_t i = startOffset, end = static_cast< int32_t >( tempPoints.size() ) + endOffset;
      i < end;
      ++i ) {
   
      curve.Add3d( tempPoints[ i ] );
    }
        
    if ( endOffset != 0 ) {
      curve.Add3d( parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D );
    }
  }

  return curve;
}

conway::geometry::IfcCurve ConwayGeometryProcessor::getIfcLine(
  const ParamsGetIfcLine &parameters) {

  conway::geometry::IfcCurve curve;

  const ParamsGetIfcTrimmedCurve& trim = parameters.paramsGetIfcTrimmedCurve;

  const bool twoDimensional = parameters.dimensions == 2;

  const glm::dvec3 placement =
    twoDimensional ? glm::dvec3( parameters.cartesianPoint2D, 0 )
                   : parameters.cartesianPoint3D;

  const glm::dvec3 vector =
    parameters.vectorOrientation * parameters.vectorMagnitude;

  // Hoisted so the CARTESIAN arm and the UNSPECIFIED arm's Cartesian case
  // read one expression rather than two that have to be kept in step.
  const glm::dvec3 cartesianTrim1 =
    twoDimensional ? glm::dvec3( trim.trim1Cartesian2D, 0 )
                   : trim.trim1Cartesian3D;

  const glm::dvec3 cartesianTrim2 =
    twoDimensional ? glm::dvec3( trim.trim2Cartesian2D, 0 )
                   : trim.trim2Cartesian3D;

  glm::dvec3 p1;
  glm::dvec3 p2;

  if ( !trim.trimExists ) {

    // A bare, untrimmed LINE: the whole IFCVECTOR laid off from the point.
    // This arm used to sit *inside* the trimExists block (conway-geom#177),
    // where it could never run - with a trim the branches above took it,
    // without one control never entered the block - so an untrimmed line
    // silently produced a zero-point curve in both 2D and 3D.
    p1 = placement;
    p2 = placement + vector;

    if ( !trim.senseAgreement ) {
      std::swap( p1, p2 );
    }

  } else if ( trim.masterRepresentation == IfcTrimmingPreference::PARAMETER ) {

    p1 = placement + vector * trim.trim1Double;
    p2 = placement + vector * trim.trim2Double;

    if ( !trim.senseAgreement ) {
      std::swap( p1, p2 );
    }

  } else if ( trim.masterRepresentation == IfcTrimmingPreference::CARTESIAN ) {

    // Sense is deliberately NOT applied here. The Cartesian arm has never
    // swapped on senseAgreement while the parameter arm always has, and STEP
    // reaches this arm for nearly every trimmed line - AP214's extractor
    // emits only CARTESIAN or PARAMETER - so "fixing" the asymmetry would
    // move geometry across the corpus on no evidence. Left as-is on purpose;
    // conway-geom#182 tracks settling it with renders behind it.
    p1 = cartesianTrim1;
    p2 = cartesianTrim2;

  } else {

    // UNSPECIFIED (and any out-of-range value, which is a schema violation
    // and gets the same treatment). IFC does not mean "no trim given" by
    // this: it means EITHER representation may be used, so the reader has to
    // choose - which is what conway-geom#170 asks for, in its words, "use
    // the Cartesian trim points when both are present, else the parameter
    // values".
    //
    // Presence is the part C++ cannot see. ParamsGetIfcTrimmedCurve carries
    // no availability flags, and interface/parameters/
    // params_get_ifc_trimmed_curve.ts marks the four Cartesian members
    // optional while trim1Double/trim2Double are required - so an absent
    // Cartesian point arrives here as a DEFAULTED ZERO, not as an absence.
    // (It cannot arrive as an absence: embind's value_object rejects a
    // missing field outright - "Missing field: trim1Cartesian3D" - and
    // conway's own extractors zero-fill with `??= { x: 0, y: 0, z: 0 }`
    // before the call.)
    //
    // So use the only signal that survives the boundary: whether the two
    // Cartesian endpoints are distinct. A trimmed line whose endpoints
    // coincide carries no trim at all, and that is exactly the shape a
    // zero-defaulted pair takes. When real Cartesian points are supplied
    // this picks them and matches the CARTESIAN arm exactly; when they are
    // absent it falls through to the parameters instead of emitting a
    // zero-length segment at the origin.
    //
    // What this deliberately does NOT recover is a trim that MIXES the two
    // representations - Cartesian at one end, a parameter at the other. That
    // needs per-endpoint presence, and the sentinel does not support it: a
    // legitimate trim starting at the origin is indistinguishable from an
    // absent one, so a per-endpoint test would corrupt correct input to
    // rescue malformed input. The pair test only ever misfires on a trim
    // whose endpoints coincide, which carries no trim either way. Neither
    // conway extractor can express a mixed trim (both pick one
    // representation for both endpoints), so this is reachable only from a
    // direct embind caller; recovering it wants availability flags on
    // ParamsGetIfcTrimmedCurve, which is a cross-repo change.
    //
    // Unlike the arm above, sense IS applied. There is no prior behaviour to
    // preserve on this path - before conway-geom#170 it produced no points
    // at all in 2D and the untrimmed line in 3D - so the choice is between
    // matching the parameter arm, which reverses as IFC's SenseAgreement
    // requires, and reproducing the Cartesian arm's known asymmetry in a
    // brand-new path. It matches the parameter arm.
    const bool cartesianTrimIsDistinct = cartesianTrim1 != cartesianTrim2;

    if ( cartesianTrimIsDistinct ) {

      p1 = cartesianTrim1;
      p2 = cartesianTrim2;

    } else {

      p1 = placement + vector * trim.trim1Double;
      p2 = placement + vector * trim.trim2Double;
    }

    if ( !trim.senseAgreement ) {
      std::swap( p1, p2 );
    }
  }

  if ( twoDimensional ) {

    curve.Add2d( glm::dvec2( p1 ) );
    curve.Add2d( glm::dvec2( p2 ) );

  } else {

    curve.Add3d( p1 );
    curve.Add3d( p2 );
  }

  return curve;
}


conway::geometry::IfcCurve ConwayGeometryProcessor::getIfcCircle(
    const ParamsGetIfcCircle &parameters) {
  conway::geometry::IfcCurve curve;

  double radius1 = parameters.radius;
  double radius2 = parameters.radius2;

  double startDegrees = 0;
  double endDegrees = 360;

  bool byPos = false;

  int startOffset   = 0;
  int endOffset     = 0;
  bool isTrimmed360 = false;

  if (parameters.paramsGetIfcTrimmedCurve.trimExists)
  {
    if (parameters.paramsGetIfcTrimmedCurve.masterRepresentation ==
      IfcTrimmingPreference::PARAMETER) {
      startDegrees = parameters.paramsGetIfcTrimmedCurve.trim1Double;
      endDegrees = parameters.paramsGetIfcTrimmedCurve.trim2Double;
      
    }
    else
    {
      byPos = true;

      if (parameters.dimensions == 2)
      {
        glm::dmat3 placement = parameters.axis2Placement2D;
        double xx = parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.x - placement[2].x;
        double yy = parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.y - placement[2].y;
        startDegrees = VectorToAngle2D(xx, yy);
        xx = parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.x - placement[2].x;
        yy = parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D.y - placement[2].y;
        endDegrees = VectorToAngle2D(xx, yy);

        startOffset = 1;
        endOffset = -1;
        
        curve.Add2d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D );
      }
      else if (parameters.dimensions == 3)
      {
        glm::dmat4 placement = parameters.axis2Placement3D;
        glm::dvec3 vecX = placement[0];
        glm::dvec3 vecY = placement[1];

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

        double dxS = glm::dot(vecX, v1);
        double dyS = glm::dot(vecY, v1);
        // double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

        double dxE = glm::dot(vecX, v2);
        double dyE = glm::dot(vecY, v2);
        // double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

        endDegrees   = VectorToAngle(dxS, dyS) - 90;
        startDegrees = VectorToAngle(dxE, dyE) - 90;
        
        if ( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D == 
            parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D ) {

          isTrimmed360 = true;

          if ( parameters.paramsGetIfcTrimmedCurve.senseAgreement ) {
            endDegrees += 360;
          } else {
            startDegrees += 360;
          }
        }

        startOffset = 1;
        endOffset = -1;
            
        curve.Add3d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D );
      }
    }
  }

  while (startDegrees < 0)
  {
    startDegrees += 360;
  }

  while (endDegrees < 0)
  {
    endDegrees += 360;
  }

  while (startDegrees > 360)
  {
    startDegrees -= 360;
  }

  while (endDegrees > 360)
  {
    endDegrees -= 360;
  }

  double lengthDegrees = 0;      

  if (parameters.paramsGetIfcTrimmedCurve.senseAgreement)
  {
    if (startDegrees > endDegrees)
    {
      endDegrees += 360;
    }
    lengthDegrees = endDegrees - startDegrees; 

  }
  else
  {
    if (startDegrees < endDegrees)
    {
      startDegrees += 360;
    }
    lengthDegrees = endDegrees - startDegrees;
  }

  if (isTrimmed360) {
    lengthDegrees = 360.0;
  }
  
  if ( parameters.paramsGetIfcTrimmedCurve.trimExists &&
       parameters.paramsGetIfcTrimmedCurve.senseAgreement &&
       parameters.isEdge ) {
    // if this is an edge, we reverse the curve   
      
    lengthDegrees = -lengthDegrees;
    
    std::swap(startDegrees, endDegrees);
  }

  double startRad = degreesToRadians(startDegrees);
  double lengthRad = degreesToRadians(lengthDegrees);

  size_t startIndex = curve.points.size();

  for (int i = startOffset; i < ( CIRCLE_SEGMENTS_MEDIUM + endOffset ); i++)
  { 
    double ratio = static_cast<double>(i) / (CIRCLE_SEGMENTS_MEDIUM - 1);
    double angle = startRad + ratio * lengthRad;

    if (parameters.dimensions == 2)
    {
      glm::dvec2 vec(0);
      vec[0] = radius1 * std::cos(angle);
      vec[1] = radius2 * std::sin(angle); // not sure why we need this, but we apparently do
      glm::dmat3 dmat = parameters.axis2Placement2D;
      // If trimming by points no rotation is required
      if (byPos)
      {
        dmat[0] = glm::dvec3(1.0, 0.0, 0.0); // Assigning [1, 0, 0] to the X vector
        dmat[1] = glm::dvec3(0.0, 1.0, 0.0); // Assigning [0, 1, 0] to the Y vector
      }
      glm::dvec2 pos = dmat * glm::dvec3(vec, 1);
      curve.Add2d(pos);
    }
    else
    {
      glm::dvec3 vec(0);
      vec[0] = radius1 * std::cos(angle);
      vec[1] = -radius2 * std::sin(angle); // negative or not???
      glm::dvec3 pos = parameters.axis2Placement3D * glm::dvec4(glm::dvec3(vec), 1);
      curve.Add3d(pos);
    }
  }

  // without a trim, we close the circle
  if ( !parameters.paramsGetIfcTrimmedCurve.trimExists ) 
  {
    if (parameters.dimensions == 2) {
      curve.Add2d(curve.points[startIndex]);
    }
    else {
      curve.Add3d(curve.points[startIndex]);
    }
  } else if ( endOffset != 0 ) {

    if (parameters.dimensions == 2) {
      curve.Add2d(parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D);
    } else {
      curve.Add3d(parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D);
    }
  }

  if ( parameters.paramsGetIfcTrimmedCurve.trimExists && 
       parameters.paramsGetIfcTrimmedCurve.senseAgreement &&
       parameters.isEdge ) {
    // if this is an edge, we reverse the curve
    std::reverse(curve.points.begin(), curve.points.end());
  }

  return curve;
}


conway::geometry::IfcCurve ConwayGeometryProcessor::getAP214Circle(
    const ParamsGetIfcCircle &parameters) {
  conway::geometry::IfcCurve curve;

  double radius1 = parameters.radius;
  double radius2 = parameters.radius2;

  bool byPos = false;

  int startOffset   = 0;
  int endOffset     = 0;
  bool isTrimmed360 = false;

  double startDegrees  = 0.0;
  double lengthDegrees = 360.0;
  double endDegrees    = 360.0;

  if (parameters.paramsGetIfcTrimmedCurve.trimExists)
  {

    if (parameters.paramsGetIfcTrimmedCurve.masterRepresentation ==
      IfcTrimmingPreference::PARAMETER) {
      
      startDegrees = parameters.paramsGetIfcTrimmedCurve.trim1Double;      
      endDegrees   = parameters.paramsGetIfcTrimmedCurve.trim2Double;

    }
    else
    {
      byPos = true;

      if (parameters.dimensions == 2)
      {
        glm::dmat3 placement = parameters.axis2Placement2D;

        glm::dvec2 vecX = placement[0];
        glm::dvec2 vecY = placement[1];
        glm::dvec2 centre = placement[2];

        glm::dvec2 v1 = glm::normalize( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D - centre );
        glm::dvec2 v2 = glm::normalize( parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D - centre );

        double dxS = glm::dot(vecX, v1);
        double dyS = glm::dot(vecY, v1);
        // double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

        double dxE = glm::dot(vecX, v2);
        double dyE = glm::dot(vecY, v2);

        // Cartesian trims give geometric positions, but the sampling loop
        // below walks the *parametric* angle t of (x, y) = (r1 cos t, r2 sin t).
        // The polar angle of a point only matches t when r1 == r2, so for an
        // ellipse scale each local component by its semi-axis before atan2 —
        // otherwise an eccentric arc sweeps across the far side of the curve
        // (bldrs-ai/test-models#45). Circles skip the scaling to stay
        // bit-identical (the divisions would perturb atan2 by ULPs).
        if ( radius1 != radius2 && radius1 != 0.0 && radius2 != 0.0 ) {
          dxS /= radius1;
          dyS /= radius2;
          dxE /= radius1;
          dyE /= radius2;
        }

        startDegrees = 180.0 * atan2( dyS, dxS ) / M_PI;
        endDegrees   = 180.0 * atan2( dyE, dxE ) / M_PI;

        if ( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D ==
             parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D ) {
          
          if ( parameters.paramsGetIfcTrimmedCurve.senseAgreement ) {
            endDegrees = startDegrees + 360;
          } else {
            endDegrees = startDegrees - 360;
          }
        }


        startOffset = 1;
        endOffset = -1;
        
        curve.Add2d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian2D );
      }
      else if (parameters.dimensions == 3)
      {
        glm::dmat4 placement = parameters.axis2Placement3D;
        glm::dvec3 vecX = glm::normalize( placement[0] );
        glm::dvec3 vecY = glm::normalize( placement[1] );
        glm::dvec3 centre = placement[3];

        glm::dvec3 v1 = glm::normalize( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D - centre );
        glm::dvec3 v2 = glm::normalize( parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D - centre );

        double dxS = glm::dot(vecX, v1);
        double dyS = glm::dot(vecY, v1);
        // double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

        double dxE = glm::dot(vecX, v2);
        double dyE = glm::dot(vecY, v2);
        // double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

        // Same semi-axis scaling as the 2D branch above: convert the trim
        // points' polar angles to parametric angles so eccentric ellipse
        // edges keep to their trimmed span (bldrs-ai/test-models#45).
        if ( radius1 != radius2 && radius1 != 0.0 && radius2 != 0.0 ) {
          dxS /= radius1;
          dyS /= radius2;
          dxE /= radius1;
          dyE /= radius2;
        }

        startDegrees = ( 180.0 * atan2( dyS, dxS ) / M_PI );
        endDegrees   = ( 180.0 * atan2( dyE, dxE ) / M_PI );

        if ( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D == 
             parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D ) {

          if ( parameters.paramsGetIfcTrimmedCurve.senseAgreement ) {
            endDegrees = startDegrees + 360;
          } else {
            endDegrees = startDegrees - 360;
          }
        }

        startOffset = 1;
        endOffset = -1;
            
        curve.Add3d( parameters.paramsGetIfcTrimmedCurve.trim1Cartesian3D );
      }
    }
    
    lengthDegrees = endDegrees - startDegrees;

    if ( parameters.paramsGetIfcTrimmedCurve.senseAgreement ) {

      if ( lengthDegrees < 0.0 ) {
        lengthDegrees += 360.0;
      }

    } else  {

      // sense disagreement: sweep from trim1 to trim2 with DECREASING
      // parameter, i.e. a negative length (the sampling loop below walks
      // `startRad + ratio * lengthRad`, so a negative length lands exactly
      // on trim2). This used to negate back to a positive sweep, which
      // walked the complement arc CCW from trim1 and finished ~360deg away
      // from trim2 — the final Add3d(trim2) then stitched a chord across
      // the curve. The trim1 == trim2 full-circle case above stays correct
      // (start - 360 remains -360: one full CW revolution).
      if ( lengthDegrees > 0.0) {
        lengthDegrees -= 360.0;
      }
    }
  }

  lengthDegrees = std::clamp( lengthDegrees, -360.0, 360.0 );

  double startRad  = degreesToRadians(startDegrees);
  double lengthRad = degreesToRadians(lengthDegrees);

  size_t startIndex = curve.points.size();

  // Adapt the sample count to the swept arc length so micro-arcs (CAD
  // fillet slivers exported as fractions of a millimetre) don't emit a
  // fixed 23 interior points. Those samples land closer together than
  // ~100x the 2^-24 Add3d quantisation grid, and the quantisation wobble
  // on such short constraint segments is what produced the CDT
  // "Intersecting constraint edges" failures on otherwise-valid faces
  // (Onshape AmazingHand Right_Hand.step). The chord floor is expressed
  // in multiples of the quantum so it tracks the *actual* hazard —
  // quantisation noise — independent of the file's length unit. Arcs
  // longer than CIRCLE_SEGMENTS_MEDIUM chords keep the legacy fixed
  // density (bit-identical output for typical models).
  constexpr double LOOP_POINT_QUANTUM = 0x1p-24;
  constexpr double MIN_SAMPLE_CHORD   = 256.0 * LOOP_POINT_QUANTUM;

  double arcLength =
      std::abs( degreesToRadians( lengthDegrees ) ) *
      std::max( std::abs( radius1 ), std::abs( radius2 ) );

  // Trimmed arcs may collapse to their chord (endpoints are appended
  // exactly); an untrimmed full circle needs enough samples to stay a
  // ring, so keep the LOW floor there.
  int minimumSegments =
      parameters.paramsGetIfcTrimmedCurve.trimExists ? 2 : CIRCLE_SEGMENTS_LOW;

  int segments = std::clamp(
      static_cast< int >( std::ceil( arcLength / MIN_SAMPLE_CHORD ) ),
      minimumSegments,
      CIRCLE_SEGMENTS_MEDIUM );

  for (int i = startOffset; i < ( segments + endOffset ); i++)
  {
    double ratio = static_cast<double>(i) / (segments - 1);
    double angle = startRad + ratio * lengthRad;

    if (parameters.dimensions == 2)
    {
      glm::dvec2 vec(0);
      vec[0] = radius1 * std::cos(angle);
      vec[1] = radius2 * std::sin(angle); // not sure why we need this, but we apparently do
      glm::dmat3 dmat = parameters.axis2Placement2D;
      // If trimming by points no rotation is required
      if (byPos)
      {
        dmat[0] = glm::dvec3(1.0, 0.0, 0.0); // Assigning [1, 0, 0] to the X vector
        dmat[1] = glm::dvec3(0.0, 1.0, 0.0); // Assigning [0, 1, 0] to the Y vector
      }
      glm::dvec2 pos = dmat * glm::dvec3(vec, 1);
      curve.Add2d(pos);
    }
    else
    {
      glm::dvec3 vec(0);
      vec[0] = radius1 * std::cos(angle);
      vec[1] = radius2 * std::sin(angle); // negative or not???
      glm::dvec3 pos = parameters.axis2Placement3D * glm::dvec4(glm::dvec3(vec), 1);
      curve.Add3d(pos);
    }
  }

  // without a trim, we close the circle
  if ( !parameters.paramsGetIfcTrimmedCurve.trimExists ) 
  {
    if (parameters.dimensions == 2) {
      curve.Add2d(curve.points[startIndex]);
    }
    else {
      curve.Add3d(curve.points[startIndex]);
    }
  } else if ( endOffset != 0 ) {

    if (parameters.dimensions == 2) {
      curve.Add2d(parameters.paramsGetIfcTrimmedCurve.trim2Cartesian2D);
    } else {
      curve.Add3d(parameters.paramsGetIfcTrimmedCurve.trim2Cartesian3D);
    }
  }

  return curve;
}

conway::geometry::Geometry ConwayGeometryProcessor::getRevolvedAreaSolid(
    const ParamsGetRevolvedAreaSolid &parameters) {
  conway::geometry::Geometry geom;

  conway::geometry::IfcProfile profile = parameters.profile;
  glm::dmat4 placement = parameters.placement;
  glm::dvec3 axis = parameters.axis;

  double angle = parameters.angle;

  bool closed = false;

  glm::dvec3 pos = parameters.axisPosition;

  IfcCurve directrix = BuildArc(parameters.scalingFactor, pos, axis, angle, parameters.circleSegments);
  if(glm::distance(directrix.points[0], directrix.points[directrix.points.size() - 1]) < 1.0E-04) //EPS_BIG
  {
      closed = true;
  }

  if (!profile.isComposite)
  {
      geom = Sweep(parameters.scalingFactor, closed, profile, directrix, axis, false);
  }
  else
  {
      for (uint32_t i = 0; i < profile.profiles.size(); i++)
      {
          conway::geometry::Geometry geom_t = Sweep(parameters.scalingFactor, closed, profile.profiles[i], directrix, axis, false, false);
          // geom.AddPart(geom_t);
          geom.AppendGeometry(geom_t);
      }
  }

  return geom;
}

conway::geometry::Geometry ConwayGeometryProcessor::getExtrudedAreaSolid(
    const ParamsGetExtrudedAreaSolid &parameters) {
  conway::geometry::Geometry geom;
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
    
    if ( flipWinding ) {

      geom.ReverseFaces();
    }
  } else {
    for (uint32_t i = 0; i < profile.profiles.size(); i++) {
      Geometry geom_t = Extrude(profile.profiles[i], dir, depth);
      if (flipWinding) {
        
        geom.ReverseFaces();
      }

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
