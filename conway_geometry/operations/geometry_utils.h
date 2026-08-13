/*
 * Decoupling:
 * https://github.com/nickcastel50/conway-geom/blob/59e9d56f6a19b5953186b78362de649437b46281/Decoupling.md
 * Ref:
 * https://github.com/IFCjs/web-ifc/blob/28681f5c4840b7ecf301e7888f98202f00adf306/src/wasm/geometry/operations/geometryutils.h
 * */

#pragma once

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#pragma clang diagnostic ignored "-Wunused-function"
#endif 

#include <mapbox/earcut.hpp>

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include "utility/LoaderError.h"
#include "representation/Geometry.h"
#include "representation/IfcGeometryReps.h"
#include "logging/Logger.h"
#include "manifold_utils.h"
#include "../structures/scratch_arena.h"

#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES 1
#endif

#include <math.h>
#include <algorithm>
#include <cmath>

namespace conway::geometry {

constexpr double EPS_BIG2 = 1e-3;

// areaOfTriangle2 lived here and had one caller, GetBasisFromCoplanarPoints,
// which is now on sinSquaredAt instead. Removed rather than left for reuse:
// the name says "area" but the value is |ab x ac|^2, i.e. (2 * area)^2, which
// scales as length^4 - and comparing that against a fixed constant, as its
// only caller did, is scale-dependent by a factor of 10^12 between a model in
// metres and the same model in millimetres. That was the bug. Anything
// needing a real area wants 0.5 * sqrt of it, and anything asking "are these
// three points usable as a basis" wants sinSquaredAt.

constexpr double DOUBLE_TO_RADIANS = M_PI / 180.0;

inline constexpr double degreesToRadians(double angle) {
  return angle * DOUBLE_TO_RADIANS;
}


inline glm::dvec3 projectOntoPlane(const glm::dvec3 &origin,
                                   const glm::dvec3 &normal,
                                   const glm::dvec3 &point,
                                   const glm::dvec3 &dir) {
  // project {et} onto the plane, following the extrusion normal
  double ldotn = glm::dot(dir, normal);
  if (ldotn == 0) {
    Logger::logWarning("0 direction in extrude");
    return glm::dvec3(0);
  } else {
    glm::dvec3 dpos = origin - glm::dvec3(point);
    double dist = glm::dot(dpos, normal) / ldotn;
    return point + dist * dir;
  }
}

inline glm::dvec3 computeNormal(const glm::dvec3 v1, const glm::dvec3 v2,
                                const glm::dvec3 v3) {
  glm::dvec3 v12(v2 - v1);
  glm::dvec3 v13(v3 - v1);

  glm::dvec3 norm = glm::cross(v12, v13);

  return glm::normalize(norm);
}

inline bool GetWindingOfTriangle(const glm::dvec3 &a, const glm::dvec3 &b,
                                 const glm::dvec3 &c) {
  auto norm = computeNormal(a, b, c);

  return /*norm.z > 0.0; */glm::dot(norm, glm::dvec3(0, 0, 1)) > 0.0;
}

//! This implementation generates much more vertices than needed, and does not
//! have smoothed normals
// TODO: Review rotate90 value, as it should be inferred from IFC but the source
// data had not been identified yet An arbitrary value has been added in
// IFCSURFACECURVESWEPTAREASOLID but this is a bad solution
inline Geometry Sweep(
    const double scaling, const bool closed, const IfcProfile &profile,
    const IfcCurve &directrix,
    const glm::dvec3 &initialDirectrixNormal = glm::dvec3(0),
    const bool rotate90 = false,
    const bool optimize = true) {
  Geometry geom;

  std::vector<glm::vec<3, glm::f64>> dpts;

  // Remove repeated points
  for (size_t i = 0; i < directrix.points.size(); i++)
  {
    if (i < directrix.points.size() - 1)
    {
      if (glm::distance(directrix.points[i], directrix.points[i + 1]) > EPS_BIG2 / scaling || !optimize)
      {
        dpts.push_back(directrix.points[i]);
      }
    }
    else
    {
      dpts.push_back(directrix.points[i]);
    }
  }

  if (closed)
  {
    glm::vec<3, glm::f64> dirStart = dpts[dpts.size() - 2] - dpts[dpts.size() - 1];
    glm::vec<3, glm::f64> dirEnd = dpts[1] - dpts[0];
    std::vector<glm::vec<3, glm::f64>> newDpts;
    newDpts.push_back(dpts[0] + dirStart);
    for (size_t i = 0; i < dpts.size(); i++)
    {
      newDpts.push_back(dpts[i]);
    }
    newDpts.push_back(dpts[dpts.size() - 1] + dirEnd);
    dpts = newDpts;
  }

  if (dpts.size() <= 1)
  {
      // nothing to sweep
    return geom;
  }

    // compute curve for each part of the directrix
  std::vector<IfcCurve> curves;
  std::vector<glm::dmat4> transforms;

  for (size_t i = 0; i < dpts.size(); i++)
  {
    IfcCurve segmentForCurve;

    glm::dvec3 planeNormal;
    glm::dvec3 directrixSegmentNormal;
    glm::dvec3 planeOrigin;

      if (i == 0) // start
      {
        planeNormal = glm::normalize(dpts[1] - dpts[0]);
        directrixSegmentNormal = planeNormal;
        planeOrigin = dpts[0];
      }
      else if (i == dpts.size() - 1) // end
      {
        planeNormal = glm::normalize(dpts[i] - dpts[i - 1]);
        directrixSegmentNormal = planeNormal;
        planeOrigin = dpts[i];
      }
      else // middle
      {
        // possibly the directrix is bad
        glm::dvec3 n1 = glm::normalize(dpts[i] - dpts[i - 1]);
        glm::dvec3 n2 = glm::normalize(dpts[i + 1] - dpts[i]);
        glm::dvec3 p = glm::normalize(glm::cross(n1, n2));

        // double prod = glm::dot(n1, n2);

        if (std::isnan(p.x))
        {
          // TODO: sometimes outliers cause the perp to become NaN!
          // this is bad news, as it nans the points added to the final mesh
          // also, it's hard to bail out now :/
          // see curve.add() for more info on how this is currently "solved"
          printf("NaN perp!\n");
        }

        glm::dvec3 u1 = glm::normalize(glm::cross(n1, p));
        glm::dvec3 u2 = glm::normalize(glm::cross(n2, p));

        // TODO: When n1 and n2 have similar direction but opposite side...
        // ... projection tend to infinity. -> glm::dot(n1, n2)
        // I implemented a bad solution to prevent projection to infinity
        if (glm::dot(n1, n2) < -0.9)
        {
          n2 = -n2;
          u2 = -u2;
        }

        glm::dvec3 au = glm::normalize(u1 + u2);
        planeNormal = glm::normalize(glm::cross(au, p));
        directrixSegmentNormal = n1; // n1 or n2 doesn't matter

        planeOrigin = dpts[i];
      }

      if (curves.empty())
      {
        // construct initial curve
        glm::dvec3 left;
        glm::dvec3 right;
        if (initialDirectrixNormal == glm::dvec3(0))
        {
          left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.y, directrixSegmentNormal.x, directrixSegmentNormal.z));
          if (left == glm::dvec3(0, 0, 0))
          {
            left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.x, directrixSegmentNormal.z, directrixSegmentNormal.y));
          }
          if (left == glm::dvec3(0, 0, 0))
          {
            left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.z, directrixSegmentNormal.y, directrixSegmentNormal.x));
          }
          right = glm::normalize(glm::cross(directrixSegmentNormal, left));
          left = glm::normalize(glm::cross(directrixSegmentNormal, right));
        }
        else
        {
          left = glm::cross(directrixSegmentNormal, initialDirectrixNormal);
          glm::dvec3 side = glm::normalize(initialDirectrixNormal);
          right = glm::normalize(glm::cross(directrixSegmentNormal, left));
          left = glm::normalize(glm::cross(directrixSegmentNormal, right));
          right *= side;
        }

        if (left == glm::dvec3(0, 0, 0))
        {
          printf("0 left vec in sweep!\n");
        }

        // project profile onto planeNormal, place on planeOrigin
        // TODO: look at holes
        auto &ppts = profile.curve.points;
        for (auto &pt2D : ppts)
        {				
          glm::dvec3 pt = -pt2D.x * left + -pt2D.y * right + planeOrigin;
          if(rotate90)
          {
            pt = -pt2D.x * right - pt2D.y * left + planeOrigin;
          }
          glm::dvec3 proj = projectOntoPlane(planeOrigin, planeNormal, pt, directrixSegmentNormal);
          
          segmentForCurve.Add3d(proj);
        }
      }
      else
      {
        // project previous curve onto the normal
        const IfcCurve &prevCurve = curves.back();

        auto &ppts = prevCurve.points;
        for (auto &pt : ppts)
        {
          glm::dvec3 proj = projectOntoPlane(planeOrigin, planeNormal, pt, directrixSegmentNormal);

          segmentForCurve.points.push_back( proj );
        }
      }

      if (!closed || (i != 0 && i != dpts.size() - 1))
      {
        curves.emplace_back( std::move( segmentForCurve ) );
      }
    }

    if (closed)
    {
      dpts.pop_back();
      dpts.erase(dpts.begin());
    }

    // connect the curves
    for (size_t i = 1; i < dpts.size(); i++) {

      const auto &c1 = curves[i - 1].points;
      const auto &c2 = curves[i].points;

      uint32_t capSize = static_cast< uint32_t >( c1.size() );
      for (size_t j = 1; j < capSize; j++) {
        uint32_t bli = geom.MakeVertex( c1[j - 1] );
        uint32_t bri = geom.MakeVertex( c1[j - 0] );

        uint32_t tli = geom.MakeVertex( c2[j - 1] );
        uint32_t tri = geom.MakeVertex( c2[j - 0] );

        geom.MakeTriangle( tli, bri, bli );
        geom.MakeTriangle( tli, tri, bri );
      }	
    }

    // io::DumpSVGCurve(directrix.points, "directrix.html");
    // DumpIfcGeometry(geom, "sweep.obj");

    return geom;
}

inline Geometry SweepCircular(
  const double scaling,
  const bool closed,
  const IfcProfile &profile,
  const double radius,
  const IfcCurve &directrix,
  const glm::dvec3 &initialDirectrixNormal = glm::dvec3(0),
  const bool rotate90 = false)
{
  Geometry geom;

  std::vector<glm::vec<3, glm::f64>> dpts;

  // Remove repeated points
  for (size_t i = 0; i < directrix.points.size(); i++)
  {
    if (i < directrix.points.size() - 1)
    {
      if (glm::distance(directrix.points[i], directrix.points[i + 1]) > EPS_BIG2 / scaling)
      {
        dpts.push_back(directrix.points[i]);
      }
    }
    else
    {
      dpts.push_back(directrix.points[i]);
    }
  }

  if (closed)
  {
    glm::vec<3, glm::f64> dirStart = dpts[dpts.size() - 2] - dpts[dpts.size() - 1];
    glm::vec<3, glm::f64> dirEnd = dpts[1] - dpts[0];
    std::vector<glm::vec<3, glm::f64>> newDpts;
    newDpts.push_back(dpts[0] + dirStart);
    for (size_t i = 0; i < dpts.size(); i++)
    {
      newDpts.push_back(dpts[i]);
    }
    newDpts.push_back(dpts[dpts.size() - 1] + dirEnd);
    dpts = newDpts;
  }

  if (dpts.size() <= 1)
  {
      // nothing to sweep
    return geom;
  }

    // compute curve for each part of the directrix
  std::vector<IfcCurve> curves;
  std::vector<glm::dmat4> transforms;

  for (size_t i = 0; i < dpts.size(); i++)
  {
    IfcCurve segmentForCurve;

    glm::dvec3 directrix2;
    glm::dvec3 planeNormal;
    glm::dvec3 directrixSegmentNormal;
    glm::dvec3 planeOrigin;

      if (i == 0) // start
      {
        planeNormal = glm::normalize(dpts[1] - dpts[0]);
        directrixSegmentNormal = planeNormal;
        planeOrigin = dpts[0];
        directrix2 = planeNormal;
      }
      else if (i == dpts.size() - 1) // end
      {
        planeNormal = glm::normalize(dpts[i] - dpts[i - 1]);
        directrixSegmentNormal = planeNormal;
        planeOrigin = dpts[i];
        directrix2 = planeNormal;
      }
      else // middle
      {
        // possibly the directrix is bad
        glm::dvec3 n1 = glm::normalize(dpts[i] - dpts[i - 1]);
        glm::dvec3 n2 = glm::normalize(dpts[i + 1] - dpts[i]);
        glm::dvec3 p = glm::normalize(glm::cross(n1, n2));
        directrix2 = -n1;

        // double prod = glm::dot(n1, n2);

        if (std::isnan(p.x))
        {
          // TODO: sometimes outliers cause the perp to become NaN!
          // this is bad news, as it nans the points added to the final mesh
          // also, it's hard to bail out now :/
          // see curve.add() for more info on how this is currently "solved"
          printf("NaN perp!\n");
        }

        glm::dvec3 u1 = glm::normalize(glm::cross(n1, p));
        glm::dvec3 u2 = glm::normalize(glm::cross(n2, p));

        // TODO: When n1 and n2 have similar direction but opposite side...
        // ... projection tend to infinity. -> glm::dot(n1, n2)
        // I implemented a bad solution to prevent projection to infinity
        if (glm::dot(n1, n2) < -0.9)
        {
          n2 = -n2;
          u2 = -u2;
        }

        glm::dvec3 au = glm::normalize(u1 + u2);
        planeNormal = glm::normalize(glm::cross(au, p));
        directrixSegmentNormal = n1; // n1 or n2 doesn't matter

        planeOrigin = dpts[i];
      }

              glm::dvec3 dz = glm::normalize(directrix2);
              glm::dvec3 dx = glm::dvec3(1, 0, 0);
              glm::dvec3 dy = glm::dvec3(0, 1, 0);

      double parallelZ = glm::abs(glm::dot(dz, glm::dvec3(0, 0, 1)));

      if(parallelZ > 1 - EPS_BIG2)
      {
        dx = glm::normalize(glm::cross(dz, glm::dvec3(0, 1, 0)));
      } else {
        dx = glm::normalize(glm::cross(dz, glm::dvec3(0, 0, 1)));
      }

      dy = glm::normalize(glm::cross(dz, dx));

              glm::dmat4 profileScale = glm::dmat4(
                  glm::dvec4(dx * radius, 0),
                  glm::dvec4(dy * radius, 0),
                  glm::dvec4(dz, 0),
                  glm::dvec4(planeOrigin, 1));

      transforms.push_back(profileScale);	

      if (curves.empty())
      {
        // construct initial curve
        glm::dvec3 left;
        glm::dvec3 right;
        if (initialDirectrixNormal == glm::dvec3(0))
        {
          left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.y, directrixSegmentNormal.x, directrixSegmentNormal.z));
          if (left == glm::dvec3(0, 0, 0))
          {
            left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.x, directrixSegmentNormal.z, directrixSegmentNormal.y));
          }
          if (left == glm::dvec3(0, 0, 0))
          {
            left = glm::cross(directrixSegmentNormal, glm::dvec3(directrixSegmentNormal.z, directrixSegmentNormal.y, directrixSegmentNormal.x));
          }
          right = glm::normalize(glm::cross(directrixSegmentNormal, left));
          left = glm::normalize(glm::cross(directrixSegmentNormal, right));
        }
        else
        {
          left = glm::cross(directrixSegmentNormal, initialDirectrixNormal);
          glm::dvec3 side = glm::normalize(initialDirectrixNormal);
          right = glm::normalize(glm::cross(directrixSegmentNormal, left));
          left = glm::normalize(glm::cross(directrixSegmentNormal, right));
          right *= side;
        }

        if (left == glm::dvec3(0, 0, 0))
        {
          printf("0 left vec in sweep!\n");
        }

        // project profile onto planeNormal, place on planeOrigin
        // TODO: look at holes
        auto &ppts = profile.curve.points;
        for (auto &pt2D : ppts)
        {				
          glm::dvec3 pt = -pt2D.x * left + -pt2D.y * right + planeOrigin;
          if(rotate90)
          {
            pt = -pt2D.x * right - pt2D.y * left + planeOrigin;
          }
          glm::dvec3 proj = projectOntoPlane(planeOrigin, planeNormal, pt, directrixSegmentNormal);
          
          segmentForCurve.Add3d(proj);
        }
      }
      else
      {
        // project previous curve onto the normal
        const IfcCurve &prevCurve = curves.back();

        auto &ppts = prevCurve.points;
        for (auto &pt : ppts)
        {
          glm::dvec3 proj = projectOntoPlane(planeOrigin, planeNormal, pt, directrixSegmentNormal);

          segmentForCurve.points.push_back( proj );
        }
      }

      if (!closed || (i != 0 && i != dpts.size() - 1))
      {
        curves.emplace_back( std::move( segmentForCurve ) );
      }
    }

    if (closed)
    {
      dpts.pop_back();
      dpts.erase(dpts.begin());
    }

    // connect the curves
    for (size_t i = 1; i < dpts.size(); i++)
    {
      glm::dvec3 p1 = dpts[i - 1];
      glm::dvec3 p2 = dpts[i];
      glm::dvec3 dir = p1 - p2;
      glm::dvec4 ddir = glm::dvec4(dir, 0);

      //Only segments smaller than 10 cm will be represented, those that are bigger will be standardized
      const auto &c1 = curves[ i - 1 ].points;
      const auto &c2 = curves[ i ].points;

      size_t capSize = c1.size();

      for (size_t j = 1; j < capSize; j++) {

        uint32_t bli = geom.MakeVertex( c1[ j - 1 ] );
        uint32_t bri = geom.MakeVertex( c1[ j - 0 ] );

        uint32_t tli = geom.MakeVertex( c2[ j - 1 ] );
        uint32_t tri = geom.MakeVertex( c2[ j - 0 ] );

        geom.MakeTriangle( tli, bri, bli );
        geom.MakeTriangle( tli, tri, bri );

      }
    }

    // DumpSVGCurve(directrix.points, glm::dvec3(), "directrix.html");
    // DumpIfcGeometry(geom, "sweep.obj");

    return geom;
  }

/**
 * Squared sine of the angle at `a` in triangle (a, b, c).
 *
 * areaOfTriangle2 returned |ab x ac|^2 - the SQUARED cross product magnitude,
 * so (2 * area)^2, scaling as length^4. Dividing by |ab|^2 |ac|^2 cancels the
 * scale entirely and leaves sin^2(theta), which is what "is this triangle
 * usable as a basis" actually asks. Being dimensionless is the point: one
 * threshold then means the same thing for a model in metres and the same
 * model in millimetres, where an absolute cutoff on a length^4 quantity is
 * off by 10^12 between the two.
 *
 * WHICH vertex matters. The value is symmetric in b and c but NOT in a:
 * sin at each corner shares the same 2 * area numerator and divides by that
 * corner's own two edge lengths, so the three differ by unbounded ratios.
 * Callers must measure at the corner whose edges they actually build from.
 *
 * @param a Vertex the angle is measured at.
 * @param b Second vertex.
 * @param c Third vertex.
 * @return {double} sin^2(theta), or 0 when either edge is degenerate.
 */
inline double sinSquaredAt( const glm::dvec3& a,
                            const glm::dvec3& b,
                            const glm::dvec3& c ) {

  const glm::dvec3 ab = b - a;
  const glm::dvec3 ac = c - a;

  const double abSqr = glm::dot( ab, ab );
  const double acSqr = glm::dot( ac, ac );

  if ( abSqr <= 0.0 || acSqr <= 0.0 ) {
    return 0.0;
  }

  const glm::dvec3 norm = glm::cross( ab, ac );

  return glm::dot( norm, norm ) / ( abSqr * acSqr );
}

/**
 * Pick three points of a bound that span a plane, for use as a basis.
 *
 * @param points The bound's points. Fewer than three cannot span a plane.
 * @param v1 First basis point.
 * @param v2 Second basis point.
 * @param v3 Third basis point.
 * @return {bool} True when a usable basis was found.
 */
inline bool GetBasisFromCoplanarPoints(std::vector<glm::dvec3> &points,
                                       glm::dvec3 &v1, glm::dvec3 &v2,
                                       glm::dvec3 &v3) {

  // One bar, and it is the noise floor: below this the caller's normalize() of
  // a cross product is rounding error rather than a direction, so the face is
  // refused and logged - which the caller already handles - rather than
  // propagating inf/NaN into every vertex it contributes, silently, since NaN
  // geometry is not obviously distinguishable downstream from ordinary output.
  //
  // The original EARLY_OUT had a second job, "this triple is good enough, stop
  // looking", and two attempts at a scale-free version of that bar were both
  // wrong, in the same way: there is no value clear of tessellated input.
  // Three consecutive vertices of a regular n-gon give sin^2(2*pi/n) at the
  // middle vertex, which runs 0.5 at n=8, 0.25 at n=12, 3.8e-2 at n=32 and
  // 9.6e-3 at n=64 - so any threshold in the useful range lands between two
  // ordinary tessellation densities and short-circuits one while searching the
  // next, deciding which basis a face gets on nothing more principled than its
  // segment count.
  //
  // So there is no early out. The search below is two O(n) passes over one
  // face's points and, with the first triple held and restored, can only
  // improve on it; measured across the public corpus the difference does not
  // clear the noise in a contended batch run. Paying it unconditionally buys
  // a basis that depends on the geometry rather than on where a constant fell.
  constexpr double MIN_SIN_SQUARED = 1e-16;

  // Three points are needed to span a plane. TriangulateBounds guarantees a
  // bound EXISTS, never that it has three points - a VERTEX_LOOP is one point
  // by design - so without this the reads below are out of bounds on a
  // std::vector, which is undefined behaviour rather than a caught error.
  if ( points.size() < 3 ) {
    return false;
  }

  v1 = points[0];
  v2 = points[1];
  v3 = points[2];

  // Measured AT v2, because that is the corner the caller builds from:
  //
  //   v12 = normalize( v3 - v2 ); v13 = normalize( v1 - v2 );
  //   n   = normalize( cross( v12, v13 ) );
  //
  // both edges leave v2, and crossing two unit vectors gives exactly
  // sin(angle at v2). Conditioning on the angle at v1 instead - which is what
  // areaOfTriangle2's argument order invites - can pass a triple whose angle
  // at v2 is arbitrarily small, since the two differ by the ratio of the edge
  // lengths and that ratio is unbounded.
  const double firstSinSquared = sinSquaredAt( v2, v1, v3 );

  // Held so the search cannot make things worse. It re-anchors v2 on the
  // farthest point, which is usually better but not always: a needle bound
  // whose points[0..2] are a perfectly serviceable corner can have a far
  // anchor from which every remaining point is nearly collinear. Without this
  // the face would be refused with a usable basis already in hand - a
  // regression against the old `bestArea > 0`, which never discarded one.
  const glm::dvec3 firstV1 = v1;
  const glm::dvec3 firstV2 = v2;
  const glm::dvec3 firstV3 = v3;

  // Farthest point from v1, to give the basis its longest first edge. This
  // one IS a squared distance, so it is compared against squared distances
  // and nothing else - the previous code shared a single constant between
  // this test and the length^4 one above, which cannot be right for both.
  double distanceSqr = 0;

  for (auto &p : points) {

    glm::dvec3 v1p = v1 - p;
    double candidate2 =  glm::dot( v1p, v1p );

    if ( candidate2 > distanceSqr ) {
      v2 = p;
      distanceSqr = candidate2;
    }
  }

  // Every point coincides with points[0], so there is no edge to build on.
  if ( distanceSqr <= 0.0 ) {
    return false;
  }

  // Ranked by distance from the LINE v1-v2, not by the angle at v2. The two
  // pick differently and the distance is the one that matters here: sin at v2
  // divides that distance by |p - v2|, so ranking on it is indifferent to how
  // long the second edge is and will happily choose a point a hair away from
  // v2 that happens to sit perpendicular. The basis is then built from a
  // near-zero difference vector, and any out-of-plane error in the bound tilts
  // the projection by roughly delta / |v3 - v2| - which on a sliver quad is
  // enough to fold the ring over itself and hand earcut a self-intersection.
  // The old areaOfTriangle2 ranking had this property, being base times
  // height; it is kept deliberately rather than lost to the change of measure.
  //
  // |cross| here is |v1 - v2| * perpendicular distance, and |v1 - v2| is fixed
  // across the loop, so it ranks identically to the distance without the
  // divide.
  const glm::dvec3 axis = v1 - v2;

  double bestPerpendicular = 0;

  for (auto &p : points) {

    const glm::dvec3 cross = glm::cross( axis, p - v2 );
    const double candidate = glm::dot( cross, cross );

    if ( candidate > bestPerpendicular ) {
      bestPerpendicular = candidate;
      v3 = p;
    }
  }

  // Ranking is one thing, acceptance another: whatever the search settled on,
  // what the caller needs is that ITS cross product is a direction and not
  // rounding error, and that is sin at v2.
  const double searchSinSquared = sinSquaredAt( v2, v1, v3 );

  if ( searchSinSquared > MIN_SIN_SQUARED ) {
    return true;
  }

  // Fall back to the first triple only when the search came back UNUSABLE, not
  // merely with a lower sin. The held triple exists to stop a face being
  // refused while a working basis is in hand, and nothing more: preferring it
  // on sin alone would reintroduce, at this exit, exactly the defect the v3
  // ranking above avoids. On a 1 m by 1 um sliver whose points[0..2] happen to
  // meet at a right angle across the short edge, sin there is ~1 and sin on
  // the searched long-edge basis is small - yet the long edge is the one worth
  // having, because the frame built from the 1 um edge amplifies the bound's
  // own out-of-plane deviation by the ratio of the two and folds the projected
  // ring.
  if ( firstSinSquared > MIN_SIN_SQUARED ) {

    v1 = firstV1;
    v2 = firstV2;
    v3 = firstV3;

    return true;
  }

  // Neither is usable. The only exit that can drop a face the old code kept,
  // which is why it is the noise floor and nothing higher: a poorly
  // conditioned basis still triangulates a sliver correctly enough to be worth
  // keeping.
  return false;
}

inline void TriangulateBounds(Geometry &geometry,
                              std::vector<IfcBound3D> &bounds) {
  if (bounds.size() == 1 && bounds[0].curve.points.size() == 3) {
    auto c = bounds[0].curve;

    geometry.MakeTriangle(

      geometry.MakeVertex( c.points[ 0 ] ),
      geometry.MakeVertex( c.points[ 1 ] ),
      geometry.MakeVertex( c.points[ 2 ] )

    );

    // size_t offset = geometry.numPoints;
  } else if (bounds.size() > 0 ) {
    // bound greater than 4 vertices or with holes, triangulate
    // TODO: modify to use glm::dvec2 with custom accessors
    using Point = std::array<double, 2>;

    // AFTP: the projected-point rings below are pure per-face scratch — built
    // here, fed to earcut, and only their *values* are copied into `geometry`.
    // Draw them from this thread's bump arena and rewind at function exit
    // (nestable checkpoint, so the several callers of TriangulateBounds are all
    // covered). No geometry changes: earcut reads the same doubles.
    ScratchArenaScope arenaScope;
    using ScratchPointVec =
        std::vector<Point, conway::ScratchAllocator<Point>>;

    std::vector<ScratchPointVec, conway::ScratchAllocator<ScratchPointVec>>
        polygon;

    uint32_t offset = geometry.vertices.size();

    // Whether any bound actually declared itself the outer one. Hoisted out of
    // the multi-bound block because the basis search below needs it: it is the
    // difference between "bounds[0] is the outer bound" and "bounds[0] is
    // whichever loop happened to be listed first".
    bool outerBoundKnown = bounds.size() == 1;

    // if more than one bound
    if (bounds.size() > 1) {
        // locate the outer bound index
      int outerIndex = -1;
      for (size_t i = 0; i < bounds.size(); i++) {
        if (bounds[i].type == IfcBoundType::OUTERBOUND) {
          outerIndex = i;
          break;
        }
      }

      if (outerIndex == -1) {
        Logger::logWarning( "Expected outer bound, using fallback tesselation." );
      } else {
        // swap the outer bound to the first position
        std::swap(bounds[0], bounds[outerIndex]);
        outerBoundKnown = true;
      }
    }


    glm::dvec3 v1, v2, v3;

    // Where no bound declared itself outer, whichever loop is listed first is
    // an ordering accident, so refusing the face because THAT one is
    // degenerate would discard every valid loop behind it. There, take the
    // basis from the first bound that can supply one. A too-short loop and a
    // collinear one are the same accident and get the same treatment.
    //
    // Where the outer bound IS known, bounds[0] is it and nothing may be
    // promoted over it. Everything else is a hole, and swapping a hole into
    // first position would hand the winding test the hole's winding and the
    // tesselatePlane fallback a ring to fill as solid - a face emitted inside
    // out with its hole plugged, silently, which is worse than the reported
    // drop this replaces.
    size_t usable = 0;
    const size_t searchEnd = outerBoundKnown ? 1 : bounds.size();

    while ( usable < searchEnd &&
            !GetBasisFromCoplanarPoints(
              bounds[usable].curve.points, v1, v2, v3 ) ) {
      ++usable;
    }

    if ( usable >= searchEnd ) {

      // Two categories, two messages, and no counts interpolated into either:
      // dedup keys on the message text, so a per-face number here would split
      // one family across a row per distinct count and stop `count` from
      // sizing it. The split that carries information is "the loop was never
      // long enough to be a face", whose fix is upstream in whatever built it,
      // versus "the loop is a full one that happens to be collinear" - which
      // the old single "No basis found for brep!" could not distinguish, and
      // that conflation is part of why the out-of-bounds read below went
      // unnoticed.
      //
      // GetBasisFromCoplanarPoints reads its argument's first three points and
      // nothing established there were three - the branch above only
      // guarantees a bound EXISTS - so that read ran off the end of a
      // std::vector, undefined behaviour rather than a caught error. Not
      // hypothetical: 63 bounds across three models of the public corpus have
      // one or two points.
      //
      // No trailing newline: the regression batch writes each distinct message
      // as one errors.csv cell, and one here makes every such row a quoted
      // multi-line field for no gain.
      // Scope matters as much as category. Only [0, searchEnd) was examined,
      // which is the outer bound alone whenever one was typed, so a message
      // about "no bound of this face" would be a false claim about the holes
      // on exactly the faces where they were never looked at.
      bool anyTooShort = false;
      bool anyCollinear = false;

      for ( size_t where = 0; where < searchEnd; ++where ) {

        if ( bounds[where].curve.points.size() < 3 ) {
          anyTooShort = true;
        } else {
          anyCollinear = true;
        }
      }

      // Both when the examined bounds are a mix, rather than one message
      // chosen by an OR. Picking one would be a false statement about the
      // other half, and since errors.csv dedups on message text it would fold
      // the two families back together - undoing the split this block exists
      // for. Only reachable with more than one bound examined, i.e. when no
      // outer bound was typed.
      if ( anyTooShort ) {

        Logger::logError( outerBoundKnown ?
          "Outer bound of this face has too few points to form it" :
          "This face has a bound with too few points to form it" );
      }

      if ( anyCollinear ) {

        Logger::logError( outerBoundKnown ?
          "Outer bound of this face is collinear, cannot form a basis" :
          "This face has a bound that is collinear, cannot form a basis" );
      }

      return;
    }

    if ( usable != 0 ) {
      std::swap( bounds[0], bounds[usable] );
    }

    glm::dvec3 v12(glm::normalize(v3 - v2));
    glm::dvec3 v13(glm::normalize(v1 - v2));
    glm::dvec3 n = glm::normalize(glm::cross(v12, v13));

    v12 = glm::cross(v13, n);
    
    // check winding of outer bound
    IfcCurve test;

    for (size_t i = 0; i < bounds[0].curve.points.size(); i++) {
      glm::dvec3 pt = bounds[0].curve.points[i];
      glm::dvec3 pt2 = pt - v1;

      glm::dvec2 proj(glm::dot(pt2, v12), glm::dot(pt2, v13));

      test.Add2d( proj );
    }

    // if the outer bound is clockwise under the current projection (v12,v13,n),
    // we invert the projection
    if (!test.IsCCW()) {
      n *= -1;
      std::swap(v12, v13);
    }

    // if the first bound is not an outer bound now, this is unexpected
    if ( bounds[0].type != IfcBoundType::OUTERBOUND && bounds.size() > 1 ) {
            
      if ( geometry.vertices.empty() && geometry.triangles.empty() ) {

        tesselatePlane(
          geometry, bounds, [&]( const glm::dvec3& vertex ) {

            glm::dvec3 pt2 = vertex - v1;

            return glm::dvec2(glm::dot(pt2, v12), glm::dot(pt2, v13));
          } );

      } else {

        Geometry intermediate;

        tesselatePlane(
          intermediate, bounds, [&]( const glm::dvec3& vertex ) {

            glm::dvec3 pt2 = vertex - v1;

            return glm::dvec2(glm::dot(pt2, v12), glm::dot(pt2, v13));
          } );

        geometry.AppendGeometry( intermediate );
      }

      return;
    }
    
    // if the first bound is not an outer bound now, this is unexpected
    if ( bounds[0].type != IfcBoundType::OUTERBOUND && bounds.size() > 1 ) {
      Logger::logWarning("Expected outer bound first!");
    }

    for (auto &bound : bounds) {

      ScratchPointVec points;

      for (size_t i = 0; i < bound.curve.points.size(); i++) {

        glm::dvec3 pt = bound.curve.points[i];

        geometry.MakeVertex( pt );

        // project pt onto plane of curve to obtain 2d coords
        glm::dvec3 pt2 = pt - v1;

        glm::dvec2 proj(glm::dot(pt2, v12), glm::dot(pt2, v13));

        points.push_back({proj.x, proj.y});
      }

      polygon.push_back(std::move(points));
    }

    std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);

    for (size_t i = 0; i < indices.size(); i += 3) {
      geometry.MakeTriangle(
        offset + indices[ i + 0 ],
        offset + indices[ i + 1 ],
        offset + indices[ i + 2 ] );
    }
  } else {

    // Reached only when bounds is EMPTY - the two branches above cover
    // size() == 1 with three points and size() > 0 - so the message this used
    // to print, which read bounds[0].curve.points.size(), indexed an empty
    // vector every time it fired. Same defect class as the one in
    // GetBasisFromCoplanarPoints above, found alongside it.
    Logger::logError( "Face has no bounds" );
  }
}

inline Geometry SectionedSurface(
    IfcCrossSections profiles) {
  Geometry geom;

  // Iterate over each profile, and create a surface by connecting the
  // corresponding points with faces.
  for (size_t i = 0; i < profiles.curves.size() - 1; i++) {
    IfcCurve &profile1 = profiles.curves[i];
    IfcCurve &profile2 = profiles.curves[i + 1];

    // Check that the profiles have the same number of points
    if (profile1.points.size() != profile2.points.size()) {
      Logger::logWarning(
          "profiles must have the same number of points in SectionedSurface");
    }

    std::vector<uint32_t> indices;

    // Create faces by connecting corresponding points from the two profiles
    for (size_t j = 0; j < profile1.points.size(); j++) {
      glm::dvec3 &p1 = profile1.points[j];
      int j2 = 0;
      if (profile1.points.size() > 1) {
        double pr = (double)j / (double)(profile1.points.size() - 1);
        j2 = static_cast< int >( pr * (profile2.points.size() - 1) );
      }
      glm::dvec3 &p2 = profile2.points[j2];

      // glm::dvec3 normal = glm::dvec3(0.0, 0.0, 1.0);

      // if (glm::distance(p1, p2) > 1E-5) {
      //   normal = glm::normalize(glm::cross(
      //       p2 - p1, glm::cross(p2 - p1, glm::dvec3(0.0, 0.0, 1.0))));
      // }

      indices.push_back( geom.MakeVertex( p1 ) );
      indices.push_back( geom.MakeVertex( p2 ) );
    }

    // Create the faces
    if (indices.size() > 0) {
      for (size_t j = 0; j < indices.size() - 2; j += 4) {
        geom.MakeTriangle( indices[ j ], indices[ j + 1 ], indices[ j + 2 ] );
        geom.MakeTriangle( indices[ j + 2 ], indices[j + 1], indices[j + 3] );
      }
    }
  }

  return geom;
}

inline Geometry Extrude(
  IfcProfile profile,
  const glm::dvec3& dir,
  double distance,
  const glm::dvec3& cuttingPlaneNormal = glm::dvec3(0),
  const glm::dvec3& cuttingPlanePos = glm::dvec3(0)) {
  
  Geometry geom;
  std::vector<bool> holesIndicesHash;

  // build the caps
  {
    using Point = std::array<double, 2>;
    int polygonCount = 1 + profile.holes.size();  // Main profile + holes
    std::vector<std::vector<Point>> polygon(polygonCount);

//    glm::dvec3 normal = dir;

    for (size_t i = 0; i < profile.curve.points.size(); i++) {
      glm::dvec2 pt = profile.curve.points[i];
      glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0) + dir * distance, 1);

      geom.MakeVertex( et );
      polygon[0].push_back({pt.x, pt.y});
    }

    for (size_t i = 0; i < profile.curve.points.size(); i++) {
      holesIndicesHash.push_back(false);
    }

    for (size_t i = 0; i < profile.holes.size(); i++) {
      IfcCurve hole = profile.holes[i];
      int pointCount = hole.points.size();

      bool firstPoint = true;

      for (int j = 0; j < pointCount; j++) {

        glm::dvec2 pt = hole.points[j];
        glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0) + dir * distance, 1);

        if ( profile.curve.Add2d(pt) ) {
          holesIndicesHash.push_back( firstPoint );
          firstPoint = false;

          geom.MakeVertex( et );

          polygon[i + 1].push_back(
              {pt.x, pt.y});  // Index 0 is main profile; see earcut reference
        }
      }
    }

    std::vector<uint32_t> indices = mapbox::earcut<uint32_t>( polygon );

    if (indices.size() < 3) {

      // Say WHICH degeneracy, and at the level it deserves. Measured over the
      // three worst models of the private corpus - 73 + 31 + 30 occurrences -
      // every single one is a three-point outer profile with no area at all,
      // i.e. collinear or coincident points in the source file. That extrudes
      // to a solid of zero volume, so nothing visible is lost and there is
      // nothing conway could recover: it is a note about the input, not a
      // failure of ours, and reporting it as an error put 134 rows of
      // unactionable noise in a baseline that is diffed for regressions.
      //
      // A polygon that has real area and still defeats earcut is a different
      // thing - self-intersecting, or a hole arrangement we mishandled - and
      // that IS potentially lost geometry, so it stays an error.
      //
      // The test is an UNSIGNED fan area about the ring's first point,
      // compared against the ring's own extent. Both properties are load
      // bearing, and a signed shoelace about the origin - the obvious version -
      // gets both wrong:
      //
      //   Unsigned, because signed areas cancel. A unit square traversed
      //   forward and then back, which is what a mis-stitched composite curve
      //   produces, sums to exactly zero signed area while bounding a real
      //   1x1 region - so the signed test calls the one case that must stay an
      //   error "zero area". Unsigned scores it 4.0 by the measure below.
      //
      //   Relative and recentred, because absolute zero is unreachable. Over
      //   random exactly-collinear rings the signed shoelace about the origin
      //   leaves residues up to 5e-11 of extent^2, growing with distance from
      //   the origin - so `== 0.0` would push most collinear profiles in a
      //   mm-unit model into the error branch this exists to keep them out of.
      //   Recentring first holds the residue under 4e-13 of extent^2.
      //
      // 1e-10 therefore sits ~250x above the worst collinear residue measured
      // and ~10 orders below a genuinely self-intersecting ring.
      constexpr double MIN_RELATIVE_AREA = 1e-10;

      const auto& ring = polygon[0];

      if ( ring.size() < 3 ) {

        Logger::logWarning(
          "Extruded profile has fewer than three points; nothing to extrude" );
        return geom;
      }

      const auto& origin = ring[0];

      double minX = 0, maxX = 0, minY = 0, maxY = 0;
      double unsignedDoubleArea = 0;

      for ( size_t where = 1; where < ring.size(); ++where ) {

        const double x = ring[where][0] - origin[0];
        const double y = ring[where][1] - origin[1];

        minX = std::min( minX, x );
        maxX = std::max( maxX, x );
        minY = std::min( minY, y );
        maxY = std::max( maxY, y );

        if ( where + 1 < ring.size() ) {

          const double nextX = ring[where + 1][0] - origin[0];
          const double nextY = ring[where + 1][1] - origin[1];

          unsignedDoubleArea += std::abs( x * nextY - nextX * y );
        }
      }

      const double extent = std::max( maxX - minX, maxY - minY );

      if ( unsignedDoubleArea <= MIN_RELATIVE_AREA * extent * extent ) {

        Logger::logWarning(
          "Extruded profile has no area; extrudes to zero volume" );
      } else {

        Logger::logError(
          "Extruded profile has area but could not be triangulated" );
      }

      return geom;
    }

    uint32_t offset = 0;

    bool winding = GetWindingOfTriangle( geom.vertices[ offset + indices[ 0 ] ],
                                         geom.vertices[ offset + indices[ 1 ] ],
                                         geom.vertices[ offset + indices[ 2 ] ] );
    bool flipWinding = !winding;

    for ( size_t i = 0, end = indices.size(); i < end; i += 3 ) {

      if (flipWinding) {

        geom.MakeTriangle(
          offset + indices[ i + 0 ],
          offset + indices[ i + 2 ],
          offset + indices[ i + 1 ] );
      } else {
        geom.MakeTriangle(
          offset + indices[ i + 0 ],
          offset + indices[ i + 1 ],
          offset + indices[ i + 2 ] );
      }
    }

    offset += geom.vertices.size();

   // normal = -dir;
   assert( offset <= profile.curve.points.size() );

    for (size_t i = 0; i < profile.curve.points.size(); i++) {
      glm::dvec2 pt = profile.curve.points[i];
      glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0), 1);

      if (cuttingPlaneNormal != glm::dvec3(0)) {
        et = glm::dvec4(glm::dvec3(pt, 0), 1);
        glm::dvec3 transDir = glm::dvec4(dir, 0);

        // project {et} onto the plane, following the extrusion normal
        double ldotn = glm::dot( transDir, cuttingPlaneNormal );

        if (ldotn == 0) {
          Logger::logWarning("0 direction in extrude");
        } else {
          glm::dvec3 dpos = cuttingPlanePos - glm::dvec3(et);
          double dist = glm::dot( dpos, cuttingPlaneNormal ) / ldotn;
          // we want to apply dist, even when negative
          et = et + glm::dvec4( dist * transDir, 1 );
        }
      }

      geom.MakeVertex( et );
    }

    for (size_t i = 0; i < indices.size(); i += 3) {

      if ( flipWinding ) {
        geom.MakeTriangle(
          offset + indices[ i + 0 ],
          offset + indices[ i + 1 ],
          offset + indices[ i + 2 ]);
      } else {
        geom.MakeTriangle(
          offset + indices[ i + 0 ],
          offset + indices[ i + 2 ],
          offset + indices[ i + 1 ] );
      }
    }
  }

  uint32_t capSize = profile.curve.points.size();

  for (size_t i = 1; i < capSize; i++) {
    // https://github.com/tomvandig/web-ifc/issues/5
    if (holesIndicesHash[i]) {
      continue;
    }

    uint32_t bl = i - 1;
    uint32_t br = i - 0;

    uint32_t tl = capSize + i - 1;
    uint32_t tr = capSize + i - 0;

    // this winding should be correct
    geom.MakeTriangle( tl, br, bl );
    geom.MakeTriangle( tl, tr, br );
  }

  return geom;
}

inline double VectorToAngle2D(double x, double y)
{
  double dd = sqrt(x * x + y * y);
  double xx = x / dd;
  double yy = y / dd;

  double angle = acos(xx);
  double cosv = cos(angle);
  double sinv = sin(angle);
  if (glm::abs(xx - cosv) > 1e-5 || glm::abs(yy - sinv) > 1e-5)
  {
    angle = asin(yy);
    sinv = sin(angle);
    cosv = cos(angle);
    if (glm::abs(xx - cosv) > 1e-5 || glm::abs(yy - sinv) > 1e-5)
    {
      angle = angle + (M_PI - angle) * 2;
      sinv = sin(angle);
      cosv = cos(angle);
      if (glm::abs(xx - cosv) > 1e-5 || glm::abs(yy - sinv) > 1e-5)
      {
        angle = angle + M_PI;
      }
    }
  }

  return (angle / (2 * M_PI)) * 360;
}

inline double VectorToAngle(double x, double y)
{
  double dd = sqrt(x * x + y * y);
  double xx = x / dd;
  double yy = y / dd;

  double angle = asin(xx);
  double cosv = cos(angle);

  if (glm::abs(yy - cosv) > 1e-5)
  {
    angle = acos(yy);
    double sinv = sin(angle);
    cosv = cos(angle);
    if (glm::abs(yy - cosv) > 1e-5 || glm::abs(xx - sinv) > 1e-5)
    {
      angle = angle + (M_PI - angle) * 2;
      sinv = sin(angle);
      cosv = cos(angle);
      if (glm::abs(yy - cosv) > 1e-5 || glm::abs(xx - sinv) > 1e-5)
      {
        angle = angle + M_PI;
      }
    }
  }

  return (angle / (2 * M_PI)) * 360;
}

inline bool MatrixFlipsTriangles(const glm::dmat4 &mat) {
  return glm::determinant(mat) < 0;
}

inline bool equals(glm::dvec3 A, glm::dvec3 B, double eps = 0) {
  return std::fabs(A.x - B.x) <= eps && std::fabs(A.y - B.y) <= eps &&
         std::fabs(A.z - B.z) <= eps;
}

inline double areaOfTriangle(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c) {
  glm::dvec3 ab = b - a;
  glm::dvec3 ac = c - a;

  glm::dvec3 norm = glm::cross(ab, ac);
  return glm::length(norm) / 2;
}

inline double cross2d(const glm::dvec2 &point1, const glm::dvec2 &point2) {
  return point1.x * point2.y - point1.y * point2.x;
}

inline double areaOfTriangle(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c) {
  glm::dvec2 ab = b - a;
  glm::dvec2 ac = c - a;

  double norm = cross2d(ab, ac) / 2;
  return std::fabs(norm);
}

inline double RandomDouble(double lo, double hi) {
  return lo + static_cast<double>(rand()) /
                  (static_cast<double>(RAND_MAX / (hi - lo)));
}

inline std::optional<glm::dvec3> GetOriginRec(
    IfcComposedMesh &mesh,
    std::unordered_map<uint32_t, Geometry> &geometryMap, glm::dmat4 mat) {
  glm::dmat4 newMat = mat * mesh.transformation;

  auto geomIt = geometryMap.find(mesh.expressID);

  if (geomIt != geometryMap.end()) {
    auto meshGeom = geomIt->second;

    if ( meshGeom.triangles.size() > 0 ) {
      for ( uint32_t i = 0, end = meshGeom.triangles.size(); i < end; ++i ) {
        
        const Triangle& f = meshGeom.triangles[ i ];

        glm::dvec3 a = newMat * glm::dvec4( meshGeom.vertices[ f.vertices[ 0 ] ], 1 );

        return a;
      }
    }
  }

  for (auto &c : mesh.children) {
    auto v = GetOriginRec(c, geometryMap, newMat);
    if (v.has_value()) {
      return v;
    }
  }

  return std::nullopt;
}

inline glm::dvec3 GetOrigin(
    IfcComposedMesh &mesh,
    std::unordered_map<uint32_t, Geometry> &geometryMap) {
  auto v = GetOriginRec(mesh, geometryMap, glm::dmat4(1));

  if (v.has_value()) {
    return *v;
  } else {
    return glm::dvec3(0);
  }
}

//TODO: on typescript geometry parser, must apply matrix transforms to halfspace in a different way 

inline std::array<double, 16> FlattenTransformation(
    const glm::dmat4 &transformation) {
  std::array<double, 16> flatTransformation;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      flatTransformation[i * 4 + j] = transformation[i][j];
    }
  }

  return flatTransformation;
}

inline bool notPresent( const glm::dvec3& pt, const std::vector<glm::dvec3>& points) {
  for (const auto &pt2 : points) {
    if (pt.x == pt2.x && pt.y == pt2.y && pt.z == pt2.z) {
      return false;
    }
  }
  return true;
}
}  // namespace conway::geometry
