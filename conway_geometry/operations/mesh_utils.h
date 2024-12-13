/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/operations/mesh_utils.h
 * */

#pragma once

#include <tinynurbs/tinynurbs.h>

#include <array>
#include <glm/glm.hpp>
#include <optional>
#include <unordered_map>
#include <vector>

#include "geometry_utils.h"
#include "tesselation_utils.h"

#define CONST_PI 3.141592653589793238462643383279502884L

namespace conway::geometry {

constexpr double MAX_DEFLECTION            = 0.001;
constexpr double MAX_TRIANGLE_AMPLIFACTION = 32;

// TODO: review and simplify
inline void TriangulateRevolution(Geometry &geometry,
                                  std::vector<IfcBound3D> &bounds,
                                  IfcSurface &surface) {
  // First we get the revolution data

  glm::dvec3 cent = surface.RevolutionSurface.Direction[3];
  glm::dvec3 vecX = glm::normalize(surface.RevolutionSurface.Direction[0]);
  glm::dvec3 vecY = glm::normalize(surface.RevolutionSurface.Direction[1]);
  glm::dvec3 vecZ = glm::normalize(surface.RevolutionSurface.Direction[2]);

  std::vector<std::vector<glm::dvec3>> newPoints;

  double numRots = 10;

  for (int r = 0; r < numRots; r++) {
    std::vector<glm::dvec3> newList;
    newPoints.push_back(newList);
  }

  std::vector<glm::dvec3> bounding;
  std::vector<double> angleVec;
  std::vector<double> angleDsp;

  // Now we construct the bounding box of the boundary ...
  // ... by adding the middle point of all curves

  for (size_t i = 0; i < bounds.size(); i++) {
    double xx = 0;
    double yy = 0;
    double zz = 0;
    double cc = 0;
    int lastTeam = bounds[i].curve.indices[0];
    for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {
      // If it is the first point of the group we close the previous group ...
      //  ... and create a new one. Else, the point is of the current group
      if (lastTeam != bounds[i].curve.indices[j] ||
          j == (bounds[i].curve.points.size() - 1)) {
        if (cc > 0) {
          xx /= cc;
          yy /= cc;
          zz /= cc;
          bounding.push_back(glm::dvec3(xx, yy, zz));
        }
        xx = bounds[i].curve.points[j].x;
        yy = bounds[i].curve.points[j].y;
        zz = bounds[i].curve.points[j].z;
        cc = 1;

        lastTeam = bounds[i].curve.indices[j];
      } else {
        xx += bounds[i].curve.points[j].x;
        yy += bounds[i].curve.points[j].y;
        zz += bounds[i].curve.points[j].z;
        cc++;
      }
    }
  }

  // There is a problem when points in the revolution are around 0 degrees
  // Numerical instabilities can make these points to jump from 0 to 360
  // It causes lots of trouble when drawing the boundaries in the revolution

  // The method presented here finds the angle of each point, measures the ...
  //  ... angular difference and then, if the difference is bigger than 180 ...
  //  ... corrects it to a lesser value. Finally it gets the first angle and ...
  //  ... adds the angular differences again, reconstructing a corrected
  //  boundary.

  // Now we find the angle of each point in the reference plane of the cylinder
  for (size_t j = 0; j < bounding.size(); j++) {
    double xx = bounding[j].x - cent.x;
    double yy = bounding[j].y - cent.y;
    double zz = bounding[j].z - cent.z;
    double dx = vecX.x * xx + vecX.y * yy + vecX.z * zz;
    double dy = vecY.x * xx + vecY.y * yy + vecY.z * zz;
    //				double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z *
    // zz;
    double temp = VectorToAngle(dx, dy);
    while (temp < 0) {
      temp += 360;
    }
    while (temp > 360) {
      temp -= 360;
    }
    angleVec.push_back(temp);
  }

  for (size_t i = 0; i < angleVec.size() - 1; i++) {
    if (angleVec[i] - angleVec[i + 1] > 180) {
      angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
    } else if (angleVec[i] - angleVec[i + 1] < -180) {
      angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
    } else {
      angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
    }
  }

  double startDegrees = angleVec[0];
  double endDegrees = angleVec[0];

  // Add angular differences starting from the first angle. We also correct the
  // start and end angles

  double temp = angleVec[0];
  for (size_t i = 0; i < angleDsp.size(); i++) {
    temp += angleDsp[i];
    if (endDegrees < temp) {
      endDegrees = temp;
    }
    if (startDegrees > temp) {
      startDegrees = temp;
    }
  }

  // Then we use the start and end angles as bounding boxes of the boundary ...
  //  ... we will represent this bounding box.

  double startRad = startDegrees / 180 * (double)CONST_PI;
  double endRad = endDegrees / 180 * (double)CONST_PI;
  double radSpan = endRad - startRad;
  double radStep = radSpan / (numRots - 1);

  for (size_t i = 0; i < surface.RevolutionSurface.Profile.curve.points.size();
       i++) {
    double xx = surface.RevolutionSurface.Profile.curve.points[i].x - cent.x;
    double yy = surface.RevolutionSurface.Profile.curve.points[i].y - cent.y;
    double zz = surface.RevolutionSurface.Profile.curve.points[i].z - cent.z;

    double dx = vecX.x * xx + vecX.y * yy + vecX.z * zz;
    double dy = vecY.x * xx + vecY.y * yy + vecY.z * zz;
    double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z * zz;
    double dd = sqrt(dx * dx + dy * dy);
    for (int r = 0; r < numRots; r++) {
      double angle = startRad + r * radStep;
      double dtempX = sin(angle) * dd;
      double dtempY = cos(angle) * dd;
      double newPx = dtempX * vecX.x + dtempY * vecY.x + dz * vecZ.x + cent.x;
      double newPy = dtempX * vecX.y + dtempY * vecY.y + dz * vecZ.y + cent.y;
      double newPz = dtempX * vecX.z + dtempY * vecY.z + dz * vecZ.z + cent.z;
      glm::dvec3 newPt = glm::dvec3(newPx, newPy, newPz);
      newPoints[r].push_back(newPt);
    }
  }
  for (int r = 0; r < numRots - 1; r++) {
    int r1 = r + 1;
    for (size_t s = 0; s < newPoints[r].size() - 1; s++) {

      uint32_t a = geometry.MakeVertex( newPoints[ r ][ s ] );
      uint32_t b = geometry.MakeVertex( newPoints[ r ][ s + 1 ] );
      uint32_t c = geometry.MakeVertex( newPoints[ r1 ][ s ] );
      uint32_t d = geometry.MakeVertex( newPoints[ r1 ][ s + 1 ] );

      geometry.MakeTriangle( a, b, c );
      geometry.MakeTriangle( c, b, d );
    }
  }
}

// TODO: review and simplify
inline void TriangulateCylindricalSurface(Geometry &geometry,
                                          const std::vector<IfcBound3D> &bounds,
                                          IfcSurface &surface) {
  // First we get the cylinder data

  if ( bounds.empty() ) {
    return;
  }

  double radius = surface.CylinderSurface.Radius;
  glm::dvec3 cent = surface.transformation[3];
  glm::dvec3 vecX = glm::normalize(surface.transformation[0]);
  glm::dvec3 vecY = glm::normalize(surface.transformation[1]);
  glm::dvec3 vecZ = glm::normalize(surface.transformation[2]);

  std::vector<std::vector<glm::dvec3>> newPoints;

  double minZ = DBL_MAX;
  double maxZ = -DBL_MAX;

  std::priority_queue< std::pair< double, size_t > > outsideMostBoundaries;

  // Find the relative coordinates of each curve point in the cylinder reference
  // plane Only retain the max and min relative Z
  for (size_t i = 0; i < bounds.size(); i++) {

    double localMaxZ = -DBL_MAX;
    double localMinZ = DBL_MAX;

    for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {
      glm::dvec3 vv = bounds[i].curve.points[j] - cent;
      //					double dx = glm::dot(vecX, vv);
      //					double dy = glm::dot(vecY, vv);
      double dz = glm::dot(vecZ, vv);
      
      localMaxZ = std::max( localMaxZ, dz );
      localMinZ = std::min( localMinZ, dz );
    }

    outsideMostBoundaries.push( std::make_pair( localMaxZ, i ) );

    maxZ = std::max( maxZ, localMaxZ );
    minZ = std::min( minZ, localMinZ );
  }

   using Point = std::array<double, 2>;
  std::vector<std::vector<Point>> uvBoundaryValues;
  std::vector<ParameterVertex> vertices;

  WingedEdgeMesh< ParameterVertex > mesh;

  double zBias = ( maxZ - minZ ) / radius;

  while ( !outsideMostBoundaries.empty() ) {

    std::vector<Point> points;

    size_t boundsIndex = outsideMostBoundaries.top().second;

    outsideMostBoundaries.pop();

    for (size_t j = 0; j < bounds[boundsIndex].curve.points.size(); j++) {

      glm::dvec3 pt = bounds[ boundsIndex ].curve.points[ j ];
      glm::dvec3 vv = pt - cent;

      double dx = glm::dot(vecX, vv);
      double dy = glm::dot(vecY, vv);
      double dz = glm::dot(vecZ, vv);

      glm::dvec2 pInv =
        glm::dvec2( dx, dy ) * ( zBias + ( dz - minZ ) );

      points.push_back({pInv.x, pInv.y});
      mesh.makeVertex( { pt, pInv } );
    }

    uvBoundaryValues.push_back( points );
  }

  // Triangulate projected boundary
  // Subdivide resulting triangles to increase definition
  // r indicates the level of subdivision, currently 3 you can increase it to
  // 5

  std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(uvBoundaryValues);

  for (size_t i = 0; i < indices.size(); i += 3) {

    mesh.makeTriangle( 
      indices[ i  + 0 ], 
      indices[ i  + 1 ], 
      indices[ i  + 2 ] );
  }

  tesselate(
    mesh,
    [&]( const glm::dvec3& point, [[maybe_unused]]const glm::dvec2& from ) { 
      
      glm::dvec3 vv       = point - cent;
      double     dx       = glm::dot(vecX, vv);
      double     dy       = glm::dot(vecY, vv);
      double     dz       = glm::dot(vecZ, vv);
      glm::dvec2 inPlane  = glm::dvec2( dx, dy );
      double     distance = glm::length( inPlane );

      glm::dvec2 newInPlane = inPlane * ( radius / distance );

      return cent + newInPlane.x * vecX + newInPlane.y * vecY + vecZ * dz;
    },
    mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
    MAX_DEFLECTION );

  appendMeshToGeometry( mesh, geometry );
}

// TODO: review and simplify
inline void TriangulateExtrusion(Geometry &geometry,
                                 std::vector<IfcBound3D> &bounds,
                                 IfcSurface &surface) {
  // NO EXAMPLE FILES ABOUT THIS CASE

  // THIS IS A SIMPLE EXTRUSION, NOT TRIMMED

  double len = surface.ExtrusionSurface.Length;
  glm::dvec3 dir = surface.ExtrusionSurface.Direction;

  if (!surface.ExtrusionSurface.Profile.isComposite) {
    for (size_t j = 0;
         j < surface.ExtrusionSurface.Profile.curve.points.size() - 1; j++) {
      int j2 = j + 1;

      double npx =
          surface.ExtrusionSurface.Profile.curve.points[j].x + dir.x * len;
      double npy =
          surface.ExtrusionSurface.Profile.curve.points[j].y + dir.y * len;
      double npz = dir.z * len;
      glm::dvec3 nptj1 = glm::dvec3(npx, npy, npz);
      npx = surface.ExtrusionSurface.Profile.curve.points[j2].x + dir.x * len;
      npy = surface.ExtrusionSurface.Profile.curve.points[j2].y + dir.y * len;
      npz = dir.z * len;
      glm::dvec3 nptj2 = glm::dvec3(npx, npy, npz);

      uint32_t nptj1i = geometry.MakeVertex( nptj1 );
      uint32_t nptj2i = geometry.MakeVertex( nptj2 );

      uint32_t a = geometry.MakeVertex(
        glm::dvec3(
          surface.ExtrusionSurface.Profile.curve.points[j].x,
          surface.ExtrusionSurface.Profile.curve.points[j].y,
          0 ) ); 

      uint32_t b = geometry.MakeVertex(
        glm::dvec3(
          surface.ExtrusionSurface.Profile.curve.points[j2].x,
          surface.ExtrusionSurface.Profile.curve.points[j2].y,
          0 ) ); 

      geometry.MakeTriangle(
        a,
        b,
        nptj1i );

      geometry.MakeTriangle(
        b,
        nptj2i,
        nptj1i );
    }
  } else {
    for (size_t i = 0; i < surface.ExtrusionSurface.Profile.profiles.size();
         i++) {
      for (size_t j = 0;
           j <
           surface.ExtrusionSurface.Profile.profiles[i].curve.points.size() - 1;
           j++) {
        int j2 = j + 1;

        double npx =
            surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].x +
            dir.x * len;
        double npy =
            surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].y +
            dir.y * len;
        double npz = dir.z * len;
        glm::dvec3 nptj1 = glm::dvec3(npx, npy, npz);
        npx = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].x +
              dir.x * len;
        npy = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].y +
              dir.y * len;
        npz = dir.z * len;
        glm::dvec3 nptj2 = glm::dvec3(npx, npy, npz);

        uint32_t nptj1i = geometry.MakeVertex( nptj1 );
        uint32_t nptj2i = geometry.MakeVertex( nptj2 );

      uint32_t a = geometry.MakeVertex(
        glm::dvec3(
                surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].x,
                surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].y,
                0) ); 
      
      uint32_t b = geometry.MakeVertex(
        glm::dvec3(
                surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].x,
                surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].y,
                0) ); 

        geometry.MakeTriangle(
          a,
          b,
          nptj1i);
        geometry.MakeTriangle(
          b,
          nptj2i,
          nptj1i);
      }
    }
  }
}

inline double InverseMethod(glm::dvec3 pt, const tinynurbs::RationalSurface3d& srf,
                            double pr, double rotations, double minError,
                            double maxError, double &fU, double &fV,
                            double &divisor, double maxDistance) {
  while (maxDistance > maxError && divisor < 10000) {
    for (double r = 1; r < 5; r++) {
      int round = 0;
      while (maxDistance > minError && round < 3) {
      //  printf("maxError: %.3f\n", maxError);
       // printf("minError: %.3f\n", minError);
      //  printf("round: %i\n", round);
        for (double i = 0; i < rotations; i++) {
          double rads = (i / rotations) * (double)CONST_PI * 2;
          double incU = glm::sin(rads) / (r * r * divisor);
          double incV = glm::cos(rads) / (r * r * divisor);
          if (pr > 1) {
            incV *= pr;
          } else {
            incU /= pr;
          }
          bool repeat = true;
          while (repeat) {
            double ffU = fU + incU;
            double ffV = fV + incV;
            glm::highp_dvec3 pt00 = tinynurbs::surfacePoint(srf, ffU, ffV);
            double di = glm::distance(pt00, pt);
            if (di < maxDistance) {
              maxDistance = di;
              fU = ffU;
              fV = ffV;
            } else {
              repeat = false;
            }
          }
        }
        round++;
      }
    }
    divisor *= 3;
    // printf("divisor: %.3f\n", divisor);
    // printf("maxError: %.3f\n", maxError);
    // printf("minError: %.3f\n", minError);
  }
  return maxDistance;
}

inline glm::dvec2 BSplineInverseEvaluation(glm::dvec3 pt,
                                           const tinynurbs::RationalSurface3d& srf,
                                           double scaling) {
  glm::highp_dvec3 ptc = tinynurbs::surfacePoint(srf, 0.0, 0.0);
  glm::highp_dvec3 pth = tinynurbs::surfacePoint(srf, 1.0, 0.0);
  glm::highp_dvec3 ptv = tinynurbs::surfacePoint(srf, 0.0, 1.0);

  double dh = glm::distance(ptc, pth);
  double dv = glm::distance(ptc, ptv);
  double pr = (dh + 1) / (dv + 1);

  double minError = 0.00001;
  double maxError = 0.001;
  double rotations = 6;

  double fU = 0.5;
  double fV = 0.5;
  double divisor = 100.0;
  double maxDistance = 1e+100;

  //printf("scaling: %.3f\n", scaling);
  maxDistance =
    InverseMethod(
        pt,
        srf,
        pr,
        rotations,
        minError / scaling,
        maxError / scaling,
        fU,
        fV,
        divisor,
        maxDistance );

  return glm::dvec2(fU, fV);
}

// TODO: review and simplify
inline void TriangulateBspline(Geometry &geometry,
                               const std::vector<IfcBound3D> &bounds,
                               IfcSurface &surface, double scaling) {

  tinynurbs::RationalSurface3d srf;
  srf.degree_u = surface.BSplineSurface.UDegree;
  srf.degree_v = surface.BSplineSurface.VDegree;
  size_t num_u = surface.BSplineSurface.ControlPoints.size();
  size_t num_v = surface.BSplineSurface.ControlPoints[0].size();

  std::vector<glm::dvec3> controlPoints;
  for (std::vector<glm::dvec3> row : surface.BSplineSurface.ControlPoints) {
    for (glm::dvec3 point : row) {
      controlPoints.push_back({point.x, point.y, point.z});
    }
  }
  srf.control_points = tinynurbs::array2(num_u, num_v, controlPoints);

  std::vector<double> weights;
  // for (std::vector<double> row : surface.BSplineSurface.WeightPoints) {
  for (std::vector<double> row : surface.BSplineSurface.Weights) {
    for (double weight : row) {
      weights.push_back(weight);
    }
  }

  if (weights.size() != num_u * num_v) {
    for (size_t i = 0; i < num_u * num_v; i++) {
      weights.push_back(1.0);
    }
  }
  
  srf.weights = tinynurbs::array2(num_u, num_v, weights);

  for (size_t i = 0; i < surface.BSplineSurface.UMultiplicity.size(); i++) {
    for (size_t r = 0; r < surface.BSplineSurface.UMultiplicity[i]; r++) {
      srf.knots_u.push_back(surface.BSplineSurface.UKnots[i]);
    }
  }

  for (size_t i = 0; i < surface.BSplineSurface.VMultiplicity.size(); i++) {
    for (size_t r = 0; r < surface.BSplineSurface.VMultiplicity[i]; r++) {
      srf.knots_v.push_back(surface.BSplineSurface.VKnots[i]);
    }
  }

  // If the NURBS surface is valid we continue

  if (tinynurbs::surfaceIsValid(srf)) {
    // Find projected boundary using NURBS inverse evaluation

    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> uvBoundaryValues;

    std::vector<ParameterVertex> vertices;

    WingedEdgeMesh< ParameterVertex > mesh;

    for ( size_t i = 0; i < bounds.size(); ++i ) {

      std::vector<Point> points;

      for (size_t j = 0; j < bounds[i].curve.points.size(); j++) {
        glm::dvec3 pt = bounds[i].curve.points[j];

        //hack 
        pt.x *= scaling;
        pt.y *= scaling;
        pt.z *= scaling;

        glm::dvec2 pInv = BSplineInverseEvaluation(pt, srf, 1.0f);

        // printf("[bounds[0]]: point %i, x: %.3f, y: %.3f, z: %.3f u: %.3f v: %.3f\n", j, pt.x,
        //         pt.y, pt.z, pInv.x, pInv.y);

        // pInv.x /= scaling;
        // pInv.y /= scaling;
        
        points.push_back({pInv.x, pInv.y});
        mesh.makeVertex( { pt, pInv } );
      }

      uvBoundaryValues.push_back(points);
    }

    // Triangulate projected boundary
    // Subdivide resulting triangles to increase definition
    // r indicates the level of subdivision, currently 3 you can increase it to
    // 5

    std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(uvBoundaryValues);

    for (size_t i = 0; i < indices.size(); i += 3) {

      mesh.makeTriangle( 
        indices[ i  + 0 ], 
        indices[ i  + 1 ], 
        indices[ i  + 2 ] );
    }

    tesselate(
      mesh,
      [&srf]( [[maybe_unused]]const glm::dvec3&, const glm::dvec2& from ) { 
        return tinynurbs::surfacePoint(srf, from.x, from.y); 
      },
      mesh.triangles.size() * MAX_TRIANGLE_AMPLIFACTION,
      MAX_DEFLECTION );

    appendMeshToGeometry( mesh, geometry );

  } else {
    printf("surface was not valid!\n");
  }
}
}  // namespace conway::geometry
