/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <string>
#include <vector>
#include <array>
#include <deque>
#include <unordered_map>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <mapbox/earcut.hpp>

#include "math/intersect-mesh-mesh.h"
#include "math/bool-mesh-mesh.h"

#include <tinynurbs/tinynurbs.h>

#include "ifc-schema.h"
#include "web-ifc.h"
#include "util.h"

const double EXTRUSION_DISTANCE_HALFSPACE_M = 50;

const bool DEBUG_DUMP_SVG = false;


//not sure if we'll need this anymore 
struct GeometryStatistics
{
	uint32_t meshCacheHits = 0;
	uint32_t meshCacheMisses = 0;
	uint32_t verticesCached = 0;
	uint32_t totalVertices = 0;

	double GetCacheRatio()
	{
		return static_cast<double>(meshCacheHits) / (meshCacheHits + meshCacheMisses);
	}
};

namespace webifc
{

    struct Geometry
	{
        uint32_t numPoints = 0;
		uint32_t numFaces = 0;
		std::vector<float> fvertexData;
		std::vector<double> vertexData;
		std::vector<uint32_t> indexData;
		glm::dvec3 min = glm::dvec3(DBL_MAX, DBL_MAX, DBL_MAX);
		glm::dvec3 max = glm::dvec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
		bool normalized = false;

		glm::dvec3 GetExtent() const
		{
			return max - min;
		}

		// set all vertices relative to min
		void Normalize()
		{
			for (size_t i = 0; i < vertexData.size(); i += 6)
			{
				vertexData[i + 0] = vertexData[i + 0] - min.x;
				vertexData[i + 1] = vertexData[i + 1] - min.y;
				vertexData[i + 2] = vertexData[i + 2] - min.z;
			}

			normalized = true;
		}

		inline void AddPoint(glm::dvec4 &pt, glm::dvec3 &n)
		{
			glm::dvec3 p = pt;
			AddPoint(p, n);
		}

		inline void AddPoint(glm::dvec3 &pt, glm::dvec3 &n)
		{
			// vertexData.reserve((numPoints + 1) * VERTEX_FORMAT_SIZE_FLOATS);
			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 0] = pt.x;
			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 1] = pt.y;
			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 2] = pt.z;
			vertexData.push_back(pt.x);
			vertexData.push_back(pt.y);
			vertexData.push_back(pt.z);

			min = glm::min(min, pt);
			max = glm::max(max, pt);

			vertexData.push_back(n.x);
			vertexData.push_back(n.y);
			vertexData.push_back(n.z);

			if (std::isnan(pt.x) || std::isnan(pt.y) || std::isnan(pt.z))
			{
				printf("NaN in geom!\n");
			}

			if (std::isnan(n.x) || std::isnan(n.y) || std::isnan(n.z))
			{
				printf("NaN in geom!\n");
			}

			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 3] = n.x;
			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 4] = n.y;
			// vertexData[numPoints * VERTEX_FORMAT_SIZE_FLOATS + 5] = n.z;

			numPoints += 1;
		}

		inline void AddFace(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c)
		{
			glm::dvec3 normal;
			if (!computeSafeNormal(a, b, c, normal))
			{
				// bail out, zero area triangle
				printf("zero tri");
				return;
			}

			AddFace(numPoints + 0, numPoints + 1, numPoints + 2);

			AddPoint(a, normal);
			AddPoint(b, normal);
			AddPoint(c, normal);
		}

		inline void AddFace(uint32_t a, uint32_t b, uint32_t c)
		{
			// indexData.reserve((numFaces + 1) * 3);
			// indexData[numFaces * 3 + 0] = a;
			// indexData[numFaces * 3 + 1] = b;
			// indexData[numFaces * 3 + 2] = c;
			indexData.push_back(a);
			indexData.push_back(b);
			indexData.push_back(c);

			numFaces++;
		}

		inline Face GetFace(uint32_t index) const
		{
			Face f;
			f.i0 = indexData[index * 3 + 0];
			f.i1 = indexData[index * 3 + 1];
			f.i2 = indexData[index * 3 + 2];
			return f;
		}

		inline glm::dvec3 GetPoint(uint32_t index) const
		{
			return glm::dvec3(
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]);
		}

		void GetCenterExtents(glm::dvec3 &center, glm::dvec3 &extents) const
		{
			glm::dvec3 min(DBL_MAX, DBL_MAX, DBL_MAX);
			glm::dvec3 max(DBL_MIN, DBL_MIN, DBL_MIN);

			for (size_t i = 0; i < numPoints; i++)
			{
				auto pt = GetPoint(i);
				min = glm::min(min, pt);
				max = glm::max(max, pt);
			}

			extents = (max - min);
			center = min + extents / 2.0;
		}

		IfcGeometry Normalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			IfcGeometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z));

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto a = (GetPoint(face.i0) - center) / scale;
				auto b = (GetPoint(face.i1) - center) / scale;
				auto c = (GetPoint(face.i2) - center) / scale;

				newGeom.AddFace(a, b, c);
			}

			return newGeom;
		}

		IfcGeometry DeNormalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			IfcGeometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z));

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto a = GetPoint(face.i0) * scale + center;
				auto b = GetPoint(face.i1) * scale + center;
				auto c = GetPoint(face.i2) * scale + center;

				newGeom.AddFace(a, b, c);
			}

			return newGeom;
		}

		uint32_t GetVertexData()
		{
			// unfortunately webgl can't do doubles
			if (fvertexData.size() != vertexData.size())
			{
				fvertexData.resize(vertexData.size());
				for (size_t i = 0; i < vertexData.size(); i += 6)
				{
					fvertexData[i + 0] = vertexData[i + 0];
					fvertexData[i + 1] = vertexData[i + 1];
					fvertexData[i + 2] = vertexData[i + 2];

					fvertexData[i + 3] = vertexData[i + 3];
					fvertexData[i + 4] = vertexData[i + 4];
					fvertexData[i + 5] = vertexData[i + 5];
				}

				// cleanup
				// vertexData = {};
			}

			if (fvertexData.empty())
			{
				return 0;
			}

			return (uint32_t)(size_t)&fvertexData[0];
		}

		void AddGeometry(IfcGeometry geom)
		{
			uint32_t maxIndex = numPoints;
			numPoints += geom.numPoints;
			min = glm::min(min, geom.min);
			max = glm::max(max, geom.max);
			vertexData.insert(vertexData.end(), geom.vertexData.begin(), geom.vertexData.end());
			for (uint32_t k = 0; k < geom.numFaces; k++)
			{
				AddFace(
					maxIndex + geom.indexData[k * 3 + 0],
					maxIndex + geom.indexData[k * 3 + 1],
					maxIndex + geom.indexData[k * 3 + 2]);
			}
		}

		uint32_t GetVertexDataSize()
		{
			return (uint32_t)fvertexData.size();
		}

		uint32_t GetIndexData()
		{
			return (uint32_t)(size_t)&indexData[0];
		}

		uint32_t GetIndexDataSize()
		{
			return (uint32_t)indexData.size();
		}

		bool IsEmpty()
		{
			return vertexData.empty();
		}
	};

    struct ComposedMesh
	{
		glm::dvec4 color;
		glm::dmat4 transformation;
		uint32_t expressID;
		bool hasGeometry = false;
		bool hasColor = false;
        Geometry meshGeometry;
		std::vector<ComposedMesh> children;

		std::optional<glm::dvec4> GetColor()
		{
			if (hasColor)
			{
				return color;
			}
			else
			{
				for (auto &c : children)
				{
					auto col = c.GetColor();
					if (col.has_value())
					{
						return col;
					}
				}
			}

			return std::nullopt;
		}
	};

    /*
     * The purpose of this class is to process all of the geometry extraction / transform request queries as required by the IFC file format. 
     * This will later be extended to handle other 3D file formats. 
    */
	class ConwayGeometryProcessor
	{ 
        public:
		ConwayGeometryProcessor(IfcLoader &l)
		{
		}

        //case ifc::IFCMAPPEDITEM:
        ComposedMesh getMappedItem(uint32_t ifcPresentation, uint32_t localPlacement)
        {

        }

        private:

        //case ifc::IFCAXIS1PLACEMENT:
        struct ParamsAxis1Placement3D
        {
            glm::dvec3 position;
            glm::dvec3 xAxisRef;
            glm::dvec3 zAxisRef;
            bool normalizeZ;
        };
        glm::dmat4 GetAxis1Placement(ParamsAxis1Placement3D parameters)
        {
            glm::dvec3 zAxis(0, 0, 1);
            glm::dvec3 xAxis(1, 0, 0);

            if (parameters.normalizeZ)
            {
                zAxis = glm::normalize(parameters.zAxisRef);
            }

            glm::dvec3 pos = parameters.position;
            if (std::abs(glm::dot(xAxis, zAxis)) > 0.9)
            {
                xAxis = glm::dvec3(0, 1, 0);
            }

            glm::dvec3 yAxis = glm::normalize(glm::cross(zAxis, xAxis));
            xAxis = glm::normalize(glm::cross(zAxis, yAxis));

            glm::dmat4 result = glm::dmat4(
                glm::dvec4(xAxis, 0),
                glm::dvec4(yAxis, 0),
                glm::dvec4(zAxis, 0),
                glm::dvec4(pos, 1));

            return result;
        }

        //case ifc::IFCAXIS2PLACEMENT3D:
        struct ParamsAxis2Placement3D
        {
            glm::dvec3 position;
            glm::dvec3 xAxisRef;
            glm::dvec3 zAxisRef;
            bool normalizeZ;
            bool normalizeX;
        };
        glm::dmat4 GetAxis2Placement3D(ParamsAxis2Placement3D parameters)
        {
            glm::dvec3 zAxis(0, 0, 1);
            glm::dvec3 xAxis(1, 0, 0);

            if (parameters.normalizeZ)
            {
                zAxis = glm::normalize(parameters.zAxisRef);
            }

            if (parameters.normalizeX)
            {
                xAxis = glm::normalize(parameters.xAxisRef);
            }

            glm::dvec3 pos = parameters.position;

            glm::dvec3 yAxis = glm::normalize(glm::cross(zAxis, xAxis));
            xAxis = glm::normalize(glm::cross(yAxis, zAxis));

            return glm::dmat4(
                glm::dvec4(xAxis, 0),
                glm::dvec4(yAxis, 0),
                glm::dvec4(zAxis, 0),
                glm::dvec4(pos, 1));
        }

        //case ifc::IFCLOCALPLACEMENT:
        //This case just recursively calls GetLocalPlacement, not sure if needed. See GetLocalPlacement 
		struct ParamsLocalPlacement {
			bool useRelPlacement;
			glm::dmat4 axis2Placement;
			glm::dmat4 relPlacement;
		};

		glm::dmat4 GetLocalPlacement(ParamsLocalPlacement parameters)
		{
			if ( parameters.useRelPlacement )
			{
				glm::dmat4 result = parameters.relPlacement * parameters.axis2Placement;
				return result;
			} 
			else 
			{
				glm::dmat4 relPlacement(1);
				glm::dmat4 result = relPlacement * parameters.axis2Placement;
				return result;
			}
		}

        //case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3D:
		//case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM:
        struct ParamsCartesianTransformationOperator3D
        {
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

        glm::dmat4 GetCartesianTransformationOperator3D(ParamsCartesianTransformationOperator3D parameters)
        {
            double scale1 = 1.0;
            double scale2 = 1.0;
            double scale3 = 1.0;

            glm::dvec3 Axis1( 1, 0, 0 );
            glm::dvec3 Axis2( 0, 1, 0 );
            glm::dvec3 Axis3( 0, 0, 1 );

            if (parameters.normalizeAxis1)
            {
                Axis1 = glm::normalize(parameters.axis1Ref);
            }

            if (parameters.normalizeAxis2)
            {
                Axis2 = glm::normalize(parameters.axis2Ref);
            }

            glm::dvec3 pos = parameters.position;

            if (parameters.realScale)
            {
                scale1 = parameters.scale1_;
            }

            if ( parameters.normalizeAxis3 )
            {
                Axis3 = glm::normalize(parameters.axis3Ref);
            }

            if ( parameters.nonUniform )
            {
                if ( parameters.realScale )
                {
                    scale2 = parameters.scale2_;
                }

                if ( parameters.realScale )
                {
                    scale3 = parameters.scale3_;
                }
            }
            else
            {
                scale2 = scale1;
                scale3 = scale1;
            }

            return glm::dmat4(
                glm::dvec4(Axis1 * scale1, 0),
                glm::dvec4(Axis2 * scale2, 0),
                glm::dvec4(Axis3 * scale3, 0),
                glm::dvec4(pos, 1));
        }
    };
}