/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/operations/curve-utils.h
 * */

#pragma once

#include "../logging/Logger.h"
#include "representation/ConwayCurve.h"

namespace conway::geometry {

constexpr double MIN_T_DELTA = 0.01;
constexpr double MAX_CURVE_DEFLECTION = 0.001;
constexpr double MAX_CURVE_DEFLECTION2 =
    MAX_CURVE_DEFLECTION * MAX_CURVE_DEFLECTION;

inline bool isConvexOrColinear(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) >= 0;
}

inline IfcCurve Build3DArc3Pt(const glm::dvec3& p1, const glm::dvec3& p2,
                              const glm::dvec3& p3) {
  // Calculate the center of the circle
  glm::dvec3 v1 = p2 - p1;
  glm::dvec3 v2 = p3 - p1;
  glm::dvec3 center;

  if (glm::length(glm::cross(v1, v2)) < /*EPS_MINSIZE*/ 1.0E-12) {
    // Points are collinear, so there's no unique circle.
    // You can handle this case differently or return an error.
    // For simplicity, let's return an empty curve.
    return IfcCurve();
  } else {
    // Calculate the center of the circle
    double Cx = p2.x - p1.x;
    double Cy = p2.y - p1.y;
    double Cz = p2.z - p1.z;
    double Bx = p3.x - p1.x;
    double By = p3.y - p1.y;
    double Bz = p3.z - p1.z;
    double B2 = p1.x * p1.x - p3.x * p3.x + p1.y * p1.y - p3.y * p3.y +
                p1.z * p1.z - p3.z * p3.z;
    double C2 = p1.x * p1.x - p2.x * p2.x + p1.y * p1.y - p2.y * p2.y +
                p1.z * p1.z - p2.z * p2.z;

    double CByz = Cy * Bz - Cz * By;
    double CBxz = Cx * Bz - Cz * Bx;
    double CBxy = Cx * By - Cy * Bx;
    double ZZ1 = -(Bz - Cz * Bx / Cx) / (By - Cy * Bx / Cx);
    double Z01 = -(B2 - Bx / Cx * C2) / (2 * (By - Cy * Bx / Cx));
    double ZZ2 = -(ZZ1 * Cy + Cz) / Cx;
    double Z02 = -(2 * Z01 * Cy + C2) / (2 * Cx);

    double dz = -((Z02 - p1.x) * CByz - (Z01 - p1.y) * CBxz - p1.z * CBxy) /
                (ZZ2 * CByz - ZZ1 * CBxz + CBxy);
    double dx = ZZ2 * dz + Z02;
    double dy = ZZ1 * dz + Z01;

    center = glm::dvec3(dx, dy, dz);
  }

  // Calculate the radius
  double radius = glm::distance(center, p1);

  int circleSegments = static_cast<int>(
      ceil(radius * ConwayGeometryProcessor::CIRCLE_SEGMENT_TO_RADIUS_RATIO));

  circleSegments =
      std::max(circleSegments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_LOW);
  circleSegments =
      std::min(circleSegments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_MEDIUM);

  // Using geometrical subdivision to create points on the arc
  std::vector<glm::dvec3> pointList;
  pointList.push_back(p1);
  pointList.push_back(p2);
  pointList.push_back(p3);

  while (pointList.size() < (size_t)circleSegments) {
    std::vector<glm::dvec3> tempPointList;
    for (size_t j = 0; j < pointList.size() - 1; j++) {
      glm::dvec3 pt = (pointList[j] + pointList[j + 1]) / 2.0;
      glm::dvec3 vc = glm::normalize(pt - center);
      pt = center + vc * radius;
      tempPointList.push_back(pointList[j]);
      tempPointList.push_back(pt);
    }
    tempPointList.push_back(pointList.back());
    pointList = tempPointList;
  }

  IfcCurve curve;
  for (size_t j = 0; j < pointList.size(); j++) {
    curve.Add3d(pointList.at(j));
  }

  return curve;
}

inline IfcCurve BuildArc3Pt(const glm::dvec2& p1, const glm::dvec2& p2,
                            const glm::dvec2& p3) {
  double f1 = (p1.x * p1.x - p2.x * p2.x + p1.y * p1.y - p2.y * p2.y);
  double f2 = (p1.x * p1.x - p3.x * p3.x + p1.y * p1.y - p3.y * p3.y);
  double v =
      2 * (p1.x - p2.x) * (p1.y - p3.y) - 2 * (p1.x - p3.x) * (p1.y - p2.y);

  double cenX = ((p1.y - p3.y) * f1 - (p1.y - p2.y) * f2) / v;
  double cenYa = (f2 - 2 * cenX * (p1.x - p3.x)) / (2 * (p1.y - p3.y));
  double cenYb = (f1 - 2 * cenX * (p1.x - p2.x)) / (2 * (p1.y - p2.y));
  double cenY = cenYa;
  if (std::isnan(cenY)) {
    cenY = cenYb;
  }

  glm::dvec2 pCen;
  pCen.x = cenX;
  pCen.y = cenY;

  double radius = sqrt(pow(cenX - p1.x, 2) + pow(cenY - p1.y, 2));

  int circleSegments = static_cast<int>(
      ceil(radius * ConwayGeometryProcessor::CIRCLE_SEGMENT_TO_RADIUS_RATIO));

  circleSegments =
      std::max(circleSegments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_LOW);
  circleSegments =
      std::min(circleSegments, ConwayGeometryProcessor::CIRCLE_SEGMENTS_MEDIUM);

  // printf( "Segments %d\n", circleSegments );

  // Using geometrical subdivision to avoid complex calculus with angles

  std::vector<glm::dvec2> pointList;
  pointList.push_back(p1);
  pointList.push_back(p2);
  pointList.push_back(p3);

  while (pointList.size() < (size_t)circleSegments) {
    std::vector<glm::dvec2> tempPointList;
    for (uint32_t j = 0; j < pointList.size() - 1; j++) {
      glm::dvec2 pt = (pointList[j] + pointList[j + 1]);
      pt.x /= 2;
      pt.y /= 2;
      glm::dvec2 vc = glm::normalize(pt - pCen);
      pt = pCen + vc * radius;
      tempPointList.push_back(pointList[j]);
      tempPointList.push_back(pt);
    }
    tempPointList.push_back(pointList[pointList.size() - 1]);
    pointList = tempPointList;
  }
  IfcCurve curve;
  for (uint32_t j = 0; j < pointList.size(); j++) curve.Add2d(pointList.at(j));
  return curve;
}

inline glm::dvec3 InterpolateRationalBSplineCurveWithKnots(
    double t, int degree, const std::vector<glm::dvec3>& points,
    const std::vector<double>& knots, const std::vector<double>& weights) {
  glm::dvec3 point;

  int domainLow = degree;
  int domainHigh = knots.size() - 1 - degree;

  double low = knots[domainLow];
  double high = knots[domainHigh];

  double tPrime = t * (high - low) + low;
  if (tPrime < low || tPrime > high) {
    Logger::logWarning("BSpline tPrime out of bounds\n");
    return glm::dvec3(0, 0, 0);
  }

  // find s (the spline segment) for the [t] value provided
  int s = domainHigh - 1;

  for (int i = domainLow; i < domainHigh; i++) {
    if (knots[i] <= tPrime && tPrime < knots[i + 1]) {
      s = i;
      break;
    }
  }

  // TODO: this should be done before calling the function, instead of calling
  // it for each t convert points to homogeneous coordinates
  std::vector<glm::dvec4> homogeneousPoints;
  for (size_t i = 0; i < points.size(); i++) {
    glm::dvec3 p = points[i];
    glm::dvec4 h = glm::dvec4(p.x * weights[i], p.y * weights[i],
                              p.z * weights[i], weights[i]);
    homogeneousPoints.push_back(h);
  }

  // printf( "Degree %d HP %zu\n", degree, homogeneousPoints.size() );

  // l (level) goes from 1 to the curve degree + 1
  double alpha;
  for (int l = 1; l <= degree + 1; l++) {
    // build level l of the pyramid
    for (int i = s; i > s - degree - 1 + l; i--) {
      alpha = (tPrime - knots[i]) / (knots[i + degree + 1 - l] - knots[i]);

      // interpolate each component

      double x = (1 - alpha) * homogeneousPoints[i - 1].x +
                 alpha * homogeneousPoints[i].x;
      double y = (1 - alpha) * homogeneousPoints[i - 1].y +
                 alpha * homogeneousPoints[i].y;
      double z = (1 - alpha) * homogeneousPoints[i - 1].z +
                 alpha * homogeneousPoints[i].z;
      double w = (1 - alpha) * homogeneousPoints[i - 1].w +
                 alpha * homogeneousPoints[i].w;

      homogeneousPoints[i] = glm::dvec4(x, y, z, w);
    }
  }

  // convert back to cartesian and return
  point = glm::dvec3(homogeneousPoints[s].x / homogeneousPoints[s].w,
                     homogeneousPoints[s].y / homogeneousPoints[s].w,
                     homogeneousPoints[s].z / homogeneousPoints[s].w);

  return point;
}

inline glm::dvec2 InterpolateRationalBSplineCurveWithKnots(
    double t, int degree, const std::vector<glm::dvec2>& points,
    const std::vector<double>& knots, const std::vector<double>& weights) {
  glm::dvec2 point;

  int domainLow = degree;
  int domainHigh = knots.size() - 1 - degree;

  double low = knots[domainLow];
  double high = knots[domainHigh];

  double tPrime = t * (high - low) + low;
  if (tPrime < low || tPrime > high) {
    Logger::logWarning("BSpline tPrime out of bounds\n");
    return glm::dvec2(0, 0);
  }

  // find s (the spline segment) for the [t] value provided
  int s = domainHigh - 1;
  for (int i = domainLow; i < domainHigh; i++) {
    if (knots[i] <= tPrime && tPrime < knots[i + 1]) {
      s = i;
      break;
    }
  }

  // TODO: this should be done before calling the function, instead of calling
  // it for each t convert points to homogeneous coordinates
  std::vector<glm::dvec3> homogeneousPoints;
  for (size_t i = 0; i < points.size(); i++) {
    glm::dvec2 p = points[i];
    glm::dvec3 h = glm::dvec3(p.x * weights[i], p.y * weights[i], weights[i]);
    homogeneousPoints.push_back(h);
  }

  // l (level) goes from 1 to the curve degree + 1
  double alpha;
  for (int l = 1; l <= degree + 1; l++) {
    // build level l of the pyramid
    for (int i = s; i > s - degree - 1 + l; i--) {
      alpha = (tPrime - knots[i]) / (knots[i + degree + 1 - l] - knots[i]);

      // interpolate each component

      double x = (1 - alpha) * homogeneousPoints[i - 1].x +
                 alpha * homogeneousPoints[i].x;
      double y = (1 - alpha) * homogeneousPoints[i - 1].y +
                 alpha * homogeneousPoints[i].y;
      double w = (1 - alpha) * homogeneousPoints[i - 1].z +
                 alpha * homogeneousPoints[i].z;
      glm::dvec3 p = glm::dvec3(x, y, w);

      homogeneousPoints[i] = p;
    }
  }

  // convert back to cartesian and return
  point = glm::dvec2(homogeneousPoints[s].x / homogeneousPoints[s].z,
                     homogeneousPoints[s].y / homogeneousPoints[s].z);
  return point;
}

inline std::vector<glm::dvec3> GetRationalBSplineCurveWithKnots(
    int degree, const std::vector<glm::dvec3>& points,
    const std::vector<double>& knots, const std::vector<double>& weights) {
  std::vector<glm::dvec3> result;

  std::vector<std::pair<double, glm::dvec3> > spanStack;

  double acceptedT = 0;

  result.push_back(InterpolateRationalBSplineCurveWithKnots(0, degree, points,
                                                            knots, weights));

  spanStack.emplace_back(1, InterpolateRationalBSplineCurveWithKnots(
                                1, degree, points, knots, weights));

  while (!spanStack.empty()) {
    const glm::dvec3& acceptedPoint = result.back();

    double nextT = spanStack.back().first;
    const glm::dvec3& nextPoint = spanStack.back().second;

    double candidateT = acceptedT + ((nextT - acceptedT) * 0.5);
    double deltaT = (nextT - acceptedT) * 0.5;
    glm::dvec3 midPointLinear =
        (nextPoint - acceptedPoint) * 0.5 + acceptedPoint;
    glm::dvec3 midPointSpline = InterpolateRationalBSplineCurveWithKnots(
        acceptedT + deltaT, degree, points, knots, weights);

    glm::dvec3 deflectionVector = midPointSpline - midPointLinear;

    if (deltaT > MIN_T_DELTA &&
        glm::dot(deflectionVector, deflectionVector) > MAX_CURVE_DEFLECTION2) {
      spanStack.emplace_back(candidateT, midPointSpline);

    } else {
      result.push_back(nextPoint);
      spanStack.pop_back();
      acceptedT = nextT;
    }
  }

  return result;
}

inline std::vector<glm::dvec2> GetRationalBSplineCurveWithKnots(
    int degree, const std::vector<glm::dvec2>& points,
    const std::vector<double>& knots, const std::vector<double>& weights) {
  std::vector<glm::dvec2> result;

  std::vector<std::pair<double, glm::dvec2> > spanStack;

  double acceptedT = 0;

  result.push_back(InterpolateRationalBSplineCurveWithKnots(0, degree, points,
                                                            knots, weights));

  spanStack.emplace_back(1, InterpolateRationalBSplineCurveWithKnots(
                                1, degree, points, knots, weights));

  while (!spanStack.empty()) {
    const glm::dvec2& acceptedPoint = result.back();

    double nextT = spanStack.back().first;
    const glm::dvec2& nextPoint = spanStack.back().second;

    double candidateT = acceptedT + ((nextT - acceptedT) * 0.5);
    double deltaT = (nextT - acceptedT) * 0.5;
    glm::dvec2 midPointLinear =
        (nextPoint - acceptedPoint) * 0.5 + acceptedPoint;
    glm::dvec2 midPointSpline = InterpolateRationalBSplineCurveWithKnots(
        acceptedT + deltaT, degree, points, knots, weights);

    glm::dvec2 deflectionVector = midPointSpline - midPointLinear;

    if (deltaT > MIN_T_DELTA &&
        glm::dot(deflectionVector, deflectionVector) > MAX_CURVE_DEFLECTION2) {
      spanStack.emplace_back(candidateT, midPointSpline);

    } else {
      result.push_back(nextPoint);
      spanStack.pop_back();
      acceptedT = nextT;
    }
  }

  return result;
}

inline bool IsCurveConvex(IfcCurve& curve) {
  for (size_t i = 2; i < curve.points.size(); i++) {
    glm::dvec2 a = curve.points[i - 2];
    glm::dvec2 b = curve.points[i - 1];
    glm::dvec2 c = curve.points[i - 0];

    if (!isConvexOrColinear(a, b, c)) {
      return false;
    }
  }

  return true;
}

inline bool MatrixFlipsTriangles(const glm::dmat3& mat) {
  return glm::determinant(mat) < 0;
}

// TODO: support IfcAlignment types
inline IfcCurve GetEllipseCurve(double radiusX, double radiusY, int numSegments,
                                glm::dmat3 placement = glm::dmat3(1),
                                double startRad = 0, double endRad = M_PI * 2,
                                bool swap = true,
                                bool normalToCenterEnding = false) {
  IfcCurve c;
  if (normalToCenterEnding) {
    double sweep_angle = (endRad - startRad);

    double step = sweep_angle / (numSegments - 1);

    if (endRad > startRad) {
      startRad -= step / 2;
      endRad += step / 2;
    }
    if (endRad <= startRad) {
      startRad += step / 2;
      endRad -= step / 2;
    }

    for (int i = 0; i < numSegments; i++) {
      double ratio = static_cast<double>(i) / (numSegments - 1);
      double angle = startRad + ratio * (endRad - startRad);

      glm::dvec2 circleCoordinate;
      if (swap) {
        circleCoordinate =
            glm::dvec2(radiusX * std::cos(angle), radiusY * std::sin(angle));
      } else {
        circleCoordinate =
            glm::dvec2(radiusX * std::sin(angle), radiusY * std::cos(angle));
      }
      glm::dvec2 pos = placement * glm::dvec3(circleCoordinate, 1);
      c.points.push_back(glm::dvec3(pos, 0));
    }

    c.points[0] = (c.points[0] + c.points[1]) * 0.5;

    c.points[c.points.size() - 1] =
        (c.points[c.points.size() - 1] + c.points[c.points.size() - 2]) * 0.5;

    // check for a closed curve
    if (endRad == M_PI * 2 && startRad == 0) {
      c.points.push_back(c.points[0]);

      if (MatrixFlipsTriangles(placement)) {
        c.Invert();
      }
    }
  } else {
    for (int i = 0; i < numSegments; i++) {
      double ratio = static_cast<double>(i) / (numSegments - 1);
      double angle = startRad + ratio * (endRad - startRad);

      glm::dvec2 circleCoordinate;
      if (swap) {
        circleCoordinate =
            glm::dvec2(radiusX * std::cos(angle), radiusY * std::sin(angle));
      } else {
        circleCoordinate =
            glm::dvec2(radiusX * std::sin(angle), radiusY * std::cos(angle));
      }
      glm::dvec2 pos = placement * glm::dvec3(circleCoordinate, 1);
      c.points.push_back(glm::dvec3(pos, 0));
    }

    // check for a closed curve
    if (endRad == M_PI * 2 && startRad == 0) {
      c.points.push_back(c.points[0]);

      if (MatrixFlipsTriangles(placement)) {
        c.Invert();
      }
    }
  }
  return c;
}

inline IfcCurve GetCircleCurve(double radius, int numSegments,
                               glm::dmat3 placement = glm::dmat3(1)) {
  return GetEllipseCurve(radius, radius, numSegments, placement);
}

inline void AddArc2d(IfcCurve& c, const glm::dvec2& center, double radius,
                     double startAngle,     // in radians
                     double sweepAngle,     // in radians, positive CCW
                     uint16_t segments = 8  // number of line‐segments to use
) {
  for (uint16_t i = 0; i <= segments; ++i) {
    double t = double(i) / double(segments);
    double ang = startAngle + sweepAngle * t;
    glm::dvec2 p = center + glm::dvec2(std::cos(ang), std::sin(ang)) * radius;
    c.Add2d(p);
  }
}

inline IfcCurve GetRoundedRectangleCurve(double xdim, double ydim,
                                         double roundingRadius,
                                         uint32_t circleSegments,
                                         glm::dmat3 placement = glm::dmat3(1)) {
  double halfX = xdim * 0.5;
  double halfY = ydim * 0.5;
  const glm::dmat3& M = placement;

  IfcCurve c;
  // clamp radius so it never exceeds half‐size
  double r = std::min(roundingRadius, std::min(halfX, halfY));
  // build the four corner centers (in world space)
  glm::dvec2 ctrBR = glm::dvec2(M * glm::dvec3(halfX - r, -halfY + r, 1.0));
  glm::dvec2 ctrTR = glm::dvec2(M * glm::dvec3(halfX - r, halfY - r, 1.0));
  glm::dvec2 ctrTL = glm::dvec2(M * glm::dvec3(-halfX + r, halfY - r, 1.0));
  glm::dvec2 ctrBL = glm::dvec2(M * glm::dvec3(-halfX + r, -halfY + r, 1.0));

  // bottom-right: from -90° → 0°
  AddArc2d(c, ctrBR, r, -M_PI_2, M_PI_2, circleSegments);
  // top-right: from   0° →  90°
  AddArc2d(c, ctrTR, r, 0.0, M_PI_2, circleSegments);
  // top-left:  from  90° → 180°
  AddArc2d(c, ctrTL, r, M_PI_2, M_PI_2, circleSegments);
  // bottom-left: from 180° → 270°
  AddArc2d(c, ctrBL, r, M_PI, M_PI_2, circleSegments);

  // ensure a consistent CCW orientation in world‐space
  bool flip = MatrixFlipsTriangles(M);
  if (flip ^ !c.IsCCW()) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetRectangleCurve(double xdim, double ydim,
                                  glm::dmat3 placement = glm::dmat3(1)) {
  double halfX = xdim / 2;
  double halfY = ydim / 2;

  glm::dvec2 bl = placement * glm::dvec3(-halfX, -halfY, 1);
  glm::dvec2 br = placement * glm::dvec3(halfX, -halfY, 1);

  glm::dvec2 tl = placement * glm::dvec3(-halfX, halfY, 1);
  glm::dvec2 tr = placement * glm::dvec3(halfX, halfY, 1);

  IfcCurve c;
  c.Add2d(bl);
  c.Add2d(br);
  c.Add2d(tr);
  c.Add2d(tl);
  c.Add2d(bl);

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetIShapedCurve(double width, double depth, double webThickness,
                                double flangeThickness, bool hasFillet,
                                double filletRadius,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;

  double hw = width / 2;
  double hd = depth / 2;
  double hweb = webThickness / 2;

  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));  // TL
  c.points.push_back(placement * glm::dvec3(+hw, +hd, 1));  // TR
  c.points.push_back(placement *
                     glm::dvec3(+hw, +hd - flangeThickness, 1));  // TR knee

  if (hasFillet) {
    // TODO: interpolate
    c.points.push_back(placement * glm::dvec3(+hweb + filletRadius,
                                              +hd - flangeThickness,
                                              1));  // TR elbow start
    c.points.push_back(placement *
                       glm::dvec3(+hweb, +hd - flangeThickness - filletRadius,
                                  1));  // TR elbow end

    c.points.push_back(placement *
                       glm::dvec3(+hweb, -hd + flangeThickness + filletRadius,
                                  1));  // BR elbow start
    c.points.push_back(placement * glm::dvec3(+hweb + filletRadius,
                                              -hd + flangeThickness,
                                              1));  // BR elbow end
  } else {
    c.points.push_back(
        placement * glm::dvec3(+hweb, +hd - flangeThickness, 1));  // TR elbow
    c.points.push_back(
        placement * glm::dvec3(+hweb, -hd + flangeThickness, 1));  // BR elbow
  }

  c.points.push_back(placement *
                     glm::dvec3(+hw, -hd + flangeThickness, 1));  // BR knee
  c.points.push_back(placement * glm::dvec3(+hw, -hd, 1));        // BR

  c.points.push_back(placement * glm::dvec3(-hw, -hd, 1));  // BL
  c.points.push_back(placement *
                     glm::dvec3(-hw, -hd + flangeThickness, 1));  // BL knee

  if (hasFillet) {
    // TODO: interpolate
    c.points.push_back(placement * glm::dvec3(-hweb - filletRadius,
                                              -hd + flangeThickness,
                                              1));  // BL elbow start
    c.points.push_back(placement *
                       glm::dvec3(-hweb, -hd + flangeThickness + filletRadius,
                                  1));  // BL elbow end

    c.points.push_back(placement *
                       glm::dvec3(-hweb, +hd - flangeThickness - filletRadius,
                                  1));  // TL elbow start
    c.points.push_back(placement * glm::dvec3(-hweb - filletRadius,
                                              +hd - flangeThickness,
                                              1));  // TL elbow end
  } else {
    c.points.push_back(
        placement * glm::dvec3(-hweb, -hd + flangeThickness, 1));  // BL elbow
    c.points.push_back(
        placement * glm::dvec3(-hweb, +hd - flangeThickness, 1));  // TL elbow
  }

  c.points.push_back(placement *
                     glm::dvec3(-hw, +hd - flangeThickness, 1));  // TL knee
  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));        // TL

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetUShapedCurve(double depth, double flangeWidth,
                                double webThickness, double flangeThickness,
                                double filletRadius, double edgeRadius,
                                double flangeSlope,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;

  double hd = depth / 2;
  double hw = flangeWidth / 2;
  //		double hweb = webThickness / 2;
  double slopeOffsetRight = flangeSlope * hw;
  double slopeOffsetLeft = flangeSlope * (hw - webThickness);
  // double flangeReferencePointY = hd - flangeThickness;

  // TODO: implement the radius

  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));
  c.points.push_back(placement * glm::dvec3(+hw, +hd, 1));

  c.points.push_back(
      placement * glm::dvec3(+hw, +hd - flangeThickness + slopeOffsetRight, 1));

  c.points.push_back(placement * glm::dvec3(-hw + webThickness,
                                            +hd - flangeThickness,
                                            1 - slopeOffsetLeft));
  c.points.push_back(placement * glm::dvec3(-hw + webThickness,
                                            -hd + flangeThickness,
                                            1 + slopeOffsetLeft));

  c.points.push_back(
      placement * glm::dvec3(+hw, -hd + flangeThickness - slopeOffsetRight, 1));

  c.points.push_back(placement * glm::dvec3(+hw, -hd, 1));
  c.points.push_back(placement * glm::dvec3(-hw, -hd, 1));

  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetLShapedCurve(double width, double depth, double thickness,
                                bool hasFillet, double filletRadius,
                                double edgeRadius, double legSlope,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;

  double hw = width / 2;
  double hd = depth / 2;

  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));
  c.points.push_back(placement * glm::dvec3(-hw, -hd, 1));
  c.points.push_back(placement * glm::dvec3(+hw, -hd, 1));

  if (hasFillet) {
    // TODO: Create interpolation and sloped lines
    c.points.push_back(placement * glm::dvec3(+hw, -hd + thickness, 1));
    c.points.push_back(placement *
                       glm::dvec3(-hw + thickness, -hd + thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw + thickness, +hd, 1));
  } else {
    c.points.push_back(placement * glm::dvec3(+hw, -hd + thickness, 1));
    c.points.push_back(placement *
                       glm::dvec3(-hw + thickness, -hd + thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw + thickness, +hd, 1));
  }

  c.points.push_back(placement * glm::dvec3(-hw, +hd, 1));

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetTShapedCurve(double width, double depth, double thickness,
                                bool hasFillet, double filletRadius,
                                double edgeRadius, double legSlope,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;

  double hw = width / 2;
  double hd = depth / 2;
  double hweb = thickness / 2;

  c.points.push_back(placement * glm::dvec3(hw, hd, 1));

  if (hasFillet) {
    // TODO: Create interpolation and sloped lines
    c.points.push_back(placement * glm::dvec3(hw, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(hweb, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(hweb, -hd, 1));
    c.points.push_back(placement * glm::dvec3(-hweb, -hd, 1));
    c.points.push_back(placement * glm::dvec3(-hweb, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw, hd, 1));
  } else {
    c.points.push_back(placement * glm::dvec3(hw, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(hweb, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(hweb, -hd, 1));
    c.points.push_back(placement * glm::dvec3(-hweb, -hd, 1));
    c.points.push_back(placement * glm::dvec3(-hweb, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw, hd - thickness, 1));
    c.points.push_back(placement * glm::dvec3(-hw, hd, 1));
  }

  c.points.push_back(placement * glm::dvec3(hw, hd, 1));

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetCShapedCurve(double width, double depth, double girth,
                                double thickness, bool hasFillet,
                                double filletRadius,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;

  double hw = width / 2;
  double hd = depth / 2;
  // double hweb = thickness / 2;

  c.points.push_back(placement * glm::dvec3(-hw, hd, 1));
  c.points.push_back(placement * glm::dvec3(hw, hd, 1));

  if (hasFillet) {
    // TODO: Create interpolation and sloped lines
  } else {
    c.points.push_back(placement * glm::dvec3(hw, hd - girth, 1));
    c.points.push_back(placement * glm::dvec3(hw - thickness, hd - girth, 1));
    c.points.push_back(placement *
                       glm::dvec3(hw - thickness, hd - thickness, 1));
    c.points.push_back(placement *
                       glm::dvec3(-hw + thickness, hd - thickness, 1));
    c.points.push_back(placement *
                       glm::dvec3(-hw + thickness, -hd + thickness, 1));
    c.points.push_back(placement *
                       glm::dvec3(hw - thickness, -hd + thickness, 1));
    c.points.push_back(placement * glm::dvec3(hw - thickness, -hd + girth, 1));
    c.points.push_back(placement * glm::dvec3(hw, -hd + girth, 1));
    c.points.push_back(placement * glm::dvec3(hw, -hd, 1));
    c.points.push_back(placement * glm::dvec3(-hw, -hd, 1));
  }

  c.points.push_back(placement * glm::dvec3(-hw, hd, 1));

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetZShapedCurve(double depth, double flangeWidth,
                                double webThickness, double flangeThickness,
                                double filletRadius, double edgeRadius,
                                glm::dmat3 placement = glm::dmat3(1)) {
  IfcCurve c;
  double hw = flangeWidth / 2;
  double hd = depth / 2;
  double hweb = webThickness / 2;
  // double hfla = flangeThickness / 2;

  c.points.push_back(placement * glm::dvec3(-hw, hd, 1));
  c.points.push_back(placement * glm::dvec3(hweb, hd, 1));
  c.points.push_back(placement * glm::dvec3(hweb, -hd + flangeThickness, 1));
  c.points.push_back(placement * glm::dvec3(hw, -hd + flangeThickness, 1));
  c.points.push_back(placement * glm::dvec3(hw, -hd, 1));
  c.points.push_back(placement * glm::dvec3(-hweb, -hd, 1));
  c.points.push_back(placement * glm::dvec3(-hweb, hd - flangeThickness, 1));
  c.points.push_back(placement * glm::dvec3(-hw, hd - flangeThickness, 1));
  c.points.push_back(placement * glm::dvec3(-hw, hd, 1));

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve GetTrapeziumCurve(double bottomXDim, double topXDim,
                                  double yDim, double topXOffset,
                                  glm::dmat3 placement = glm::dmat3(1)) {
  double halfX = bottomXDim / 2;
  double halfY = yDim / 2;

  glm::dvec2 bl = placement * glm::dvec3(-halfX, -halfY, 1);
  glm::dvec2 br = placement * glm::dvec3(halfX, -halfY, 1);

  glm::dvec2 tl = placement * glm::dvec3(-halfX + topXOffset, halfY, 1);
  glm::dvec2 tr =
      placement * glm::dvec3(-halfX + topXOffset + topXDim, halfY, 1);

  IfcCurve c;
  c.Add2d(bl);
  c.Add2d(br);
  c.Add2d(tr);
  c.Add2d(tl);
  c.Add2d(bl);

  if (MatrixFlipsTriangles(placement)) {
    c.Invert();
  }

  return c;
}

inline IfcCurve BuildArc(double scale, const glm::dvec3& pos,
                         const glm::dvec3& axis, double angleRad,
                         uint16_t _circleSegments) {
  IfcCurve curve;

  // project pos onto axis

  double pdota = glm::dot(axis, pos);
  glm::dvec3 pproja = pdota * axis;

  glm::dvec3 right = -(pos - pproja);

  if (glm::length(right) == 0) {
    right = glm::dvec3(EPS_BIG2, 0, 0);
    glm::dvec3 up = glm::cross(axis, right);

    auto curve2D = GetEllipseCurve(1, 1, _circleSegments, glm::dmat3(1), 0,
                                   angleRad, true, true);

    for (auto& pt2D : curve2D.points) {
      glm::dvec3 pt3D = pos + pt2D.x * right + pt2D.y * up;
      curve.Add3d(pt3D);
    }
  } else {
    glm::dvec3 up = glm::cross(axis, right);

    auto curve2D = GetEllipseCurve(1, 1, _circleSegments, glm::dmat3(1), 0,
                                   angleRad, true);

    for (auto& pt2D : curve2D.points) {
      glm::dvec3 pt3D = pos + pt2D.x * right + pt2D.y * up;
      curve.Add3d(pt3D);
    }
  }

  return curve;
}

}  // namespace conway::geometry
