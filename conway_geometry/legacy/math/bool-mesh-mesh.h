/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <glm/glm.hpp>

#include "../../operations/curve-utils.h"
#include "../../operations/geometryutils.h"
#include "../../operations/mesh_utils.h"
#include "./is-inside-mesh.h"
#define CSGJSCPP_REAL double
#define CSGJSCPP_IMPLEMENTATION
#include "mycsgjs.h"

namespace conway
{
    bool IsInsideCenterExtents(const glm::dvec3 &pt, const glm::dvec3 &center, const glm::dvec3 &extents)
	{
		glm::dvec3 delta = pt - center;
		delta = glm::abs(delta);
		glm::dvec3 offset = delta - extents;

		return offset.x < EPS_SMALL && offset.y < EPS_SMALL && offset.z < EPS_SMALL;
	}
    static void clipMesh(geometry::IfcGeometry& source, geometry::IfcGeometry& target, geometry::IfcGeometry& result, bool invert, bool flip, bool keepBoundary)
    {
        glm::dvec3 targetCenter;
        glm::dvec3 targetExtents;
        target.GetCenterExtents(targetCenter, targetExtents);

        for (uint32_t i = 0; i < source.numFaces; i++)
        {
            geometry::Face tri = source.GetFace(i);
            glm::dvec3 a = source.GetPoint(tri.i0);
            glm::dvec3 b = source.GetPoint(tri.i1);
            glm::dvec3 c = source.GetPoint(tri.i2);

            glm::dvec3 n = geometry::computeNormal(a, b, c);

            glm::dvec3 triCenter = (a + b + c) * 1.0 / 3.0;

            auto isInsideTarget = MeshLocation::INSIDE;

            if (IsInsideCenterExtents(triCenter, targetCenter, targetExtents))
            {
                isInsideTarget = isInsideMesh(triCenter, n, target);
            }
            else
            {
                isInsideTarget = MeshLocation::OUTSIDE;
            }

            if ((isInsideTarget == MeshLocation::INSIDE && !invert) || (isInsideTarget == MeshLocation::OUTSIDE && invert) || (isInsideTarget == MeshLocation::BOUNDARY && keepBoundary))
            {
                // emit triangle
                if (flip)
                {
                    result.AddFace(a, c, b);
                }
                else
                {
                    result.AddFace(a, b, c);
                }
            }
        }
    }

    static geometry::IfcGeometry boolIntersect(geometry::IfcGeometry& mesh1, geometry::IfcGeometry& mesh2)
    {
        geometry::IfcGeometry resultingMesh;

        clipMesh(mesh1, mesh2, resultingMesh, false, false, true);
        clipMesh(mesh2, mesh1, resultingMesh, false, false, false);

        return resultingMesh;
    }

    static geometry::IfcGeometry boolJoin(geometry::IfcGeometry& mesh1, geometry::IfcGeometry& mesh2)
    {
        geometry::IfcGeometry resultingMesh;

        clipMesh(mesh1, mesh2, resultingMesh, true, false, true);
        clipMesh(mesh2, mesh1, resultingMesh, true, false, false);

        return resultingMesh;
    }

    geometry::IfcGeometry boolSubtract(geometry::IfcGeometry& mesh1, geometry::IfcGeometry& mesh2)
    {

        geometry::IfcGeometry resultingMesh;

        clipMesh(mesh1, mesh2, resultingMesh, true, false, false);
        clipMesh(mesh2, mesh1, resultingMesh, false, true, false);

        return resultingMesh;
    }

    // TODO: I don't think XOR works right now...
    static geometry::IfcGeometry boolXOR(geometry::IfcGeometry& mesh1, geometry::IfcGeometry& mesh2)
    {
        geometry::IfcGeometry resultingMesh;

        clipMesh(mesh1, mesh2, resultingMesh, true, false, false);
        clipMesh(mesh2, mesh1, resultingMesh, true, false, false);

        return resultingMesh;
    }

    csgjscpp::Model IfcGeometryToCSGModel(const geometry::IfcGeometry& mesh1)
    {
        std::vector<csgjscpp::Polygon> polygons1;

        for (uint32_t i = 0; i < mesh1.numFaces; i++)
        {
            geometry::Face f = mesh1.GetFace(i);
            std::vector<csgjscpp::Vertex> verts;

            glm::dvec3 a = mesh1.GetPoint(f.i0);
            glm::dvec3 b = mesh1.GetPoint(f.i1);
            glm::dvec3 c = mesh1.GetPoint(f.i2);

            glm::dvec3 n = geometry::computeNormal(a, b, c);


            verts.push_back(csgjscpp::Vertex{ csgjscpp::Vector(a.x, a.y, a.z), csgjscpp::Vector(n.x, n.y, n.z), 0 });
            verts.push_back(csgjscpp::Vertex{ csgjscpp::Vector(b.x, b.y, b.z), csgjscpp::Vector(n.x, n.y, n.z), 0 });
            verts.push_back(csgjscpp::Vertex{ csgjscpp::Vector(c.x, c.y, c.z), csgjscpp::Vector(n.x, n.y, n.z), 0 });

            polygons1.push_back(csgjscpp::Polygon(verts));
        }

        return csgjscpp::modelfrompolygons(polygons1);
    }

    geometry::IfcGeometry boolMultiOp_CSGJSCPP(const geometry::IfcGeometry& firstGeom, const std::vector<geometry::IfcGeometry>& secondGeoms)
    {
        csgjscpp::Model model;

        for(auto& geom : secondGeoms)
        {
          auto m = IfcGeometryToCSGModel(geom);
          model = csgjscpp::csgunion(model, m);
        }

        auto firstModel = IfcGeometryToCSGModel(firstGeom);
        auto ModelResult = csgjscpp::csgsubtract(firstModel, model);;

        geometry::IfcGeometry result;

        for (uint32_t i = 0; i < ModelResult.indices.size(); i += 3)
        {
            uint32_t i0 = ModelResult.indices[i + 0];
            uint32_t i1 = ModelResult.indices[i + 1];
            uint32_t i2 = ModelResult.indices[i + 2];

            csgjscpp::Vertex va = ModelResult.vertices[i0];
            csgjscpp::Vertex vb = ModelResult.vertices[i1];
            csgjscpp::Vertex vc = ModelResult.vertices[i2];

            glm::dvec3 a(va.pos.x, va.pos.y, va.pos.z);
            glm::dvec3 b(vb.pos.x, vb.pos.y, vb.pos.z);
            glm::dvec3 c(vc.pos.x, vc.pos.y, vc.pos.z);

            result.AddFace(a, b, c);
        }

        return result;
    }
}