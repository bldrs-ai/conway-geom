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

		//case ifc::IFCAXIS2PLACEMENT2D:
		//case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2D:
		//case ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM:
		struct ParamsGetAxis2Placement2D
		{
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

		glm::dmat3 GetAxis2Placement2D(ParamsGetAxis2Placement2D parameters)
		{
			if ( parameters.isAxis2Placement2D)
			{
				glm::dvec2 xAxis = glm::dvec2(1, 0);
				if (parameters.customAxis1Ref)
				{
					xAxis = glm::normalize(parameters.axis1Ref);
				}

				glm::dvec2 pos = parameters.position2D;

				glm::dvec2 yAxis = glm::normalize(glm::dvec2(xAxis.y, -xAxis.x));

				return glm::dmat3(
					glm::dvec3(xAxis, 0),
					glm::dvec3(yAxis, 0),
					glm::dvec3(pos, 1));

			} else if (parameters.isCartesianTransformationOperator2D || parameters.isCartesianTransformationOperator2DNonUniform)
			{
				double scale1 = 1.0;
				double scale2 = 1.0;

				glm::dvec2 Axis1(1, 0);
				glm::dvec2 Axis2(0, 1);

				if ( parameters.customAxis1Ref )
				{
					Axis1 = glm::normalize(parameters.axis1Ref);
				}

				if ( parameters.customAxis2Ref )
				{
					Axis2 = glm::normalize(parameters.axis2Ref);
				}

				glm::dvec2 pos = parameters.position2D;

				if ( parameters.customScale )
				{
					scale1 = parameters.scale1;
				}

				if (parameters.isCartesianTransformationOperator2DNonUniform)
				{
					if (parameters.customScale2)
					{
						scale2 = parameters.scale2;
					}
				}

				if (parameters.isCartesianTransformationOperator2D)
				{
					scale2 = scale1;
				}

				return glm::dmat3(
					glm::dvec3(Axis1 * scale1, 0),
					glm::dvec3(Axis2 * scale2, 0),
					glm::dvec3(pos, 1));
			}

			return glm::dmat3();
		}

		template <uint32_t DIM>
		void ComputePolylineCurve(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			for (auto &point : points)
			{
				curve.Add(point);
			}

			if (edge)
			{
				if (sameSense == 1 || sameSense == -1)
				{
					std::reverse(curve.points.begin(), curve.points.end());
				}
			}
		}

		template <uint32_t DIM>
		void ComputeCompositeCurve(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			//todo: recursion 
		}

		template <uint32_t DIM>
		void ComputeCompositeCurveSegment(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			//todo: recursion 
		}

		struct TrimmingSelect
		{
			bool hasParam = false;
			bool hasPos = false;
			double param;
			glm::dvec2 pos;
			glm::dvec3 pos3D;
		};

		struct TrimmingArguments
		{
			bool exist = false;
			glm::dvec2 position;
			glm::dvec3 direction;
			TrimmingSelect start;
			TrimmingSelect end;
		};

		template <uint32_t DIM>
		void ComputeCurveLine(IfcCurve<DIM> &curve, TrimmingArguments trim, int sameSense = -1)
		{
			bool condition = sameSense == 1 || sameSense == -1;
			if (edge)
			{
				condition = !condition;
			}
			if constexpr (DIM == 2)
			{
				if (trim.start.hasPos && trim.end.hasPos)
				{

					if (condition)
					{
						curve.Add(trim.start.pos);
						curve.Add(trim.end.pos);
					}
					else
					{
						curve.Add(trim.end.pos);
						curve.Add(trim.start.pos);
					}
				}
				else if (trim.start.hasParam && trim.end.hasParam)
				{
					glm::dvec3 placement = glm::dvec3(trim.position, 0);
					glm::dvec3 vector = trim.direction;

					if (condition)
					{
						glm::dvec3 p1 = placement + vector * trim.start.param;
						glm::dvec3 p2 = placement + vector * trim.end.param;
						curve.Add(p1);
						curve.Add(p2);
					}
					else
					{
						glm::dvec3 p2 = placement + vector * trim.start.param;
						glm::dvec3 p1 = placement + vector * trim.end.param;
						curve.Add(p1);
						curve.Add(p2);
					}
				}
				else
				{
					printf("Unsupported trimmingselect IFCLINE\n");
				}
			}
			else if constexpr (DIM == 3)
			{
				if (trim.start.hasPos && trim.end.hasPos)
				{
					if (condition)
					{
						curve.Add(trim.start.pos3D);
						curve.Add(trim.end.pos3D);
					}
					else
					{
						curve.Add(trim.end.pos3D);
						curve.Add(trim.start.pos3D);
					}
				}
				else
				{
					printf("Unsupported trimmingselect IFCLINE\n");
				}
			}
		}

		template <uint32_t DIM>
		void ComputeTrimmedCurve(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			//todo: recursion 
		}

		enum CurveType 
		{
			lineIndex = 0,
			arcIndex = 1,
		};

		struct CurveSegmentIndex
		{
			uint32_t curveType;
			size_t indicesCount;
			uint32_t* indices;
		};

		template <uint32_t DIM>
		void ComputeIndexedPolycurve(IfcCurve<DIM> &curve, std::vector<glm::dvec2> points, uint32_t numSegments, CurveSegmentIndex* segments, bool curveLineIndex, bool selfIntersects)
		{
			if (selfIntersects)
			{
				sprintf("Self intersecting ifcindexedpolycurve");
				return;
			}

			if constexpr (DIM == 2)
			{
				for (size_t segmentIndex = 0; segmentIndex < numSegments; segmentIndex++)
					{
						CurveSegmentIndex segment = segments[segmentIndex];
						if ( segment.curveType == CurveType::lineIndex )
						{
							for (size_t currentIndex = 0; currentIndex < segment.indicesCount; currentIndex++ )
							{
								curve.Add( points[ segment.indices[currentIndex] - 1]  );
							}
						} else if ( segment.curveType == CurveType::arcIndex )
						{
							IfcCurve<2> arc = BuildArc3Pt(points[ segment.indices[ 0 ] - 1 ], points[ segment.indices[ 1 ] - 1 ], points[ segment.indices[ 2 ] - 1 ] );

							for (auto &pt : arc.points)
							{
								curve.Add(pt);
							}
						}
					}
			}
			else
			{
				sprintf("Parsing ifcindexedpolycurve in 3D is not possible");
			}
		}

		//TODO: Rework this parameter list, I really don't like it. Probably should split all 2D and 3D functions and get rid of the DIM template.
		template <uint32_t DIM>
		void ComputeCircleCurve(IfcCurve<DIM> &curve, glm::mat4 placementMat4, glm::mat3 placementMat3, double radius, bool selfIntersects, TrimmingArguments trim, int sameSense = -1)
		{
			glm::mat<DIM + 1, DIM + 1, glm::f64, glm::defaultp> placement;
			if constexpr (DIM == 2)
			{
				placement = placementMat3;
			}
			else
			{
				placement = placementMat4;
			}

			double startDegrees = 0;
			double endDegrees = 360;

			if (trim.exist)
			{
				if (trim.start.hasParam && trim.end.hasParam)
				{
					startDegrees = trim.start.param;
					endDegrees = trim.end.param;
				}
				else if (trim.start.hasPos && trim.end.hasPos)
				{
					if constexpr (DIM == 2)
					{
						double xx = placement[2].x - trim.start.pos.x;
						double yy = placement[2].y - trim.start.pos.y;
						startDegrees = VectorToAngle(xx, yy);
						xx = placement[2].x - trim.end.pos.x;
						yy = placement[2].y - trim.end.pos.y;
						endDegrees = VectorToAngle(xx, yy);
					}
					else if constexpr (DIM == 3)
					{
						glm::dvec4 vecX = placement[0];
						glm::dvec4 vecY = placement[1];
						glm::dvec4 vecZ = placement[2];

						glm::dvec3 v1 = glm::dvec3(
							trim.start.pos3D.x - placement[3].x,
							trim.start.pos3D.y - placement[3].y,
							trim.start.pos3D.z - placement[3].z);
						glm::dvec3 v2 = glm::dvec3(
							trim.end.pos3D.x - placement[3].x,
							trim.end.pos3D.y - placement[3].y,
							trim.end.pos3D.z - placement[3].z);

						double dxS = vecX.x * v1.x + vecX.y * v1.y + vecX.z * v1.z;
						double dyS = vecY.x * v1.x + vecY.y * v1.y + vecY.z * v1.z;
						double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

						double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
						double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
						double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

						endDegrees = VectorToAngle(dxS, dyS) - 90;
						startDegrees = VectorToAngle(dxE, dyE) - 90;
					}
				}
			}

			double startRad = startDegrees / 180 * CONST_PI;
			double endRad = endDegrees / 180 * CONST_PI;
			double lengthDegrees = endDegrees - startDegrees;

			// unset or true
			if (trimSense == 1 || trimSense == -1)
			{
				if (lengthDegrees < 0)
				{
					lengthDegrees += 360;
				}
			}
			else
			{
				if (lengthDegrees > 0)
				{
					lengthDegrees -= 360;
				}
			}

			double lengthRad = lengthDegrees / 180 * CONST_PI;

			size_t startIndex = curve.points.size();

			//const int numSegments = _loader.GetSettings().CIRCLE_SEGMENTS_HIGH;
			const int numSegments = LoaderSettings::CIRCLE_SEGMENTS_HIGH;

			for (int i = 0; i < numSegments; i++)
			{
				double ratio = static_cast<double>(i) / (numSegments - 1);
				double angle = startRad + ratio * lengthRad;

				if (sameSense == 0)
				{
					angle = endRad - ratio * lengthRad;
				}

				glm::vec<DIM, glm::f64> pos;
				if constexpr (DIM == 2)
				{
					glm::dvec2 vec(0);
					vec[0] = radius * std::cos(angle);
					vec[1] = -radius * std::sin(angle); // not sure why we need this, but we apparently do
					pos = placement * glm::dvec3(vec, 1);
				}
				else
				{
					glm::dvec3 vec(0);
					vec[0] = radius * std::cos(angle);
					vec[1] = -radius * std::sin(angle); // negative or not???
					pos = placement * glm::dvec4(glm::dvec3(vec), 1);
				}
				curve.Add(pos);
			}

			// without a trim, we close the circle
			if (!trim.exist)
			{
				curve.Add(curve.points[startIndex]);
			}
		}

		//TODO: Rework this parameter list, I really don't like it. Probably should split all 2D and 3D functions and get rid of the DIM template.
		template <uint32_t DIM>
		void ComputeEllipseCurve(IfcCurve<DIM> &curve, glm::mat4 placementMat4, glm::mat3 placementMat3, double radius1, double radius2, TrimmingArguments trim, int trimSense = -1)
		{
			glm::mat<DIM + 1, DIM + 1, glm::f64, glm::defaultp> placement;
			if constexpr (DIM == 2)
			{
				placement = placementMat3;
			}
			else
			{
				placement = placementMat4;
			}

			double startDegrees = 0;
			double endDegrees = 360;

			if (trim.exist)
			{
				if (trim.start.hasParam && trim.end.hasParam)
				{
					startDegrees = trim.start.param;
					endDegrees = trim.end.param;
				}
				else if (trim.start.hasPos && trim.end.hasPos)
				{
					if constexpr (DIM == 2)
					{
						double xx = placement[2].x - trim.start.pos.x;
						double yy = placement[2].y - trim.start.pos.y;
						startDegrees = VectorToAngle(xx, yy);
						xx = placement[2].x - trim.end.pos.x;
						yy = placement[2].y - trim.end.pos.y;
						endDegrees = VectorToAngle(xx, yy);
					}
					else if constexpr (DIM == 3)
					{
						glm::dvec4 vecX = placement[0];
						glm::dvec4 vecY = placement[1];
						glm::dvec4 vecZ = placement[2];

						glm::dvec3 v1 = glm::dvec3(
							trim.start.pos3D.x - placement[3].x,
							trim.start.pos3D.y - placement[3].y,
							trim.start.pos3D.z - placement[3].z);
						glm::dvec3 v2 = glm::dvec3(
							trim.end.pos3D.x - placement[3].x,
							trim.end.pos3D.y - placement[3].y,
							trim.end.pos3D.z - placement[3].z);

						double dxS = vecX.x * v1.x + vecX.y * v1.y + vecX.z * v1.z;
						double dyS = vecY.x * v1.x + vecY.y * v1.y + vecY.z * v1.z;
						double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

						double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
						double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
						double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

						endDegrees = VectorToAngle(dxS, dyS) - 90;
						startDegrees = VectorToAngle(dxE, dyE) - 90;
					}
				}
			}

			double startRad = startDegrees / 180 * CONST_PI;
			double endRad = endDegrees / 180 * CONST_PI;

			// TODO: Because this is an ellipse you need to correct the angles

			// startRad = atan((radius1 / radius2) * tan(startDegrees));
			// endRad = atan((radius1 / radius2) * tan(endDegrees));

			double lengthDegrees = endDegrees - startDegrees;

			// unset or true
			if (trimSense == 1 || trimSense == -1)
			{
				if (lengthDegrees < 0)
				{
					lengthDegrees += 360;
				}
			}
			else
			{
				if (lengthDegrees > 0)
				{
					lengthDegrees -= 360;
				}
			}

			double lengthRad = lengthDegrees / 180 * CONST_PI;

			size_t startIndex = curve.points.size();

			//const int numSegments = _loader.GetSettings().CIRCLE_SEGMENTS_HIGH;
			const int numSegments = LoaderSettings::CIRCLE_SEGMENTS_HIGH;

			for (int i = 0; i < numSegments; i++)
			{
				double ratio = static_cast<double>(i) / (numSegments - 1);
				double angle = startRad + ratio * lengthRad;
				if (sameSense == 0)
				{
					angle = endRad - ratio * lengthRad;
				}
				glm::vec<DIM, glm::f64> pos;
				if constexpr (DIM == 2)
				{
					glm::dvec2 vec(0);
					vec[0] = radius1 * std::cos(angle);
					vec[1] = -radius2 * std::sin(angle);
					pos = placement * glm::dvec3(vec, 1);
				}
				else
				{
					glm::dvec3 vec(0);
					vec[0] = radius1 * std::cos(angle);
					vec[1] = -radius2 * std::sin(angle); // negative or not???
					pos = placement * glm::dvec4(glm::dvec3(vec), 1);
				}
				curve.Add(pos);
			}

			// without a trim, we close the circle
			if (!trim.exist)
			{
				curve.Add(curve.points[startIndex]);
			}
		}

		/*struct ParamsComputeCurve 
		{
			//cases 
			bool isPolyLine;
			bool isCompositeCurve;
			bool isCompositeCurveSegment;
			bool isLine;
			bool isTrimmedCurve;
			bool isIndexedPolyCurve;
			bool isCircle;
			bool isEllipse;
			bool isBsplineCurve;
			bool isBsplineCurveWithKnots;
			bool isRationalBsplineCurveWithKnots;


			bool edge;
			int sameSense = -1;
			int trimSense = -1;
			
			//IfcCurve<DIM> &curve;
		};
		*/

		template <uint32_t DIM>
		void ComputeBsplineCurve(IfcCurve<DIM> &curve, std::vector<glm::vec<DIM, glm::f64>> ctrolPts, std::vector<glm::f64> knots, 
								 std::vector<glm::f64> weights, double degree, bool edge, int sameSense = -1)
		{
			bool condition = sameSense == 0;
			if (edge)
			{
				condition = !condition;
			}

			if (knots.size() != ctrolPts.size() + degree + 1)
			{
				std::cout << "Error: Knots and control points do not match" << std::endl;
			}

			std::vector<glm::vec<DIM, glm::f64>> tempPoints = GetRationalBSplineCurveWithKnots(degree, ctrolPts, knots, weights);
			for (int i = 0; i < tempPoints.size(); i++)
			{
				curve.Add(tempPoints[i]);
			}

			if (condition)
			{
				std::reverse(curve.points.begin(), curve.points.end());
			}
		}


		/*template <uint32_t DIM>
		void ComputeCurve(uint32_t expressID, IfcCurve<DIM> &curve, bool edge, int sameSense = -1, int trimSense = -1, IfcTrimmingArguments trim = {})
		{
			uint32_t lineID = _loader.ExpressIDToLineID(expressID);
			auto &line = _loader.GetLine(lineID);
			printf("%s\n", GetReadableNameFromTypeCode(line.ifcType));
			switch (line.ifcType)
			{
			case ifc::IFCPOLYLINE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto points = _loader.GetSetArgument();

				for (auto &token : points)
				{
					uint32_t pointId = _loader.GetRefArgument(token);
					curve.Add(GetCartesianPoint<DIM>(pointId));
				}

				if (edge)
				{
					if (sameSense == 1 || sameSense == -1)
					{
						std::reverse(curve.points.begin(), curve.points.end());
					}
				}

				break;
			}
			case ifc::IFCCOMPOSITECURVE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto segments = _loader.GetSetArgument();
				auto selfIntersects = _loader.GetStringArgument();

				if (selfIntersects == "T")
				{
					// TODO: this is probably bad news
					_loader.ReportError({LoaderErrorType::UNSPECIFIED, "Self intersecting composite curve", line.expressID});
				}

				for (auto &token : segments)
				{
					if (DEBUG_DUMP_SVG)
					{
						if constexpr (DIM == 2)
						{
							DumpSVGCurve(curve.points, L"partial_curve.html");
						}
					}

					uint32_t segmentId = _loader.GetRefArgument(token);

					ComputeCurve<DIM>(segmentId, curve, edge, sameSense, trimSense);
				}

				break;
			}
			case ifc::IFCCOMPOSITECURVESEGMENT:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto transition = _loader.GetStringArgument();
				auto sameSenseS = _loader.GetStringArgument();
				auto parentID = _loader.GetRefArgument();

				bool sameSense = sameSenseS == "T";

				ComputeCurve<DIM>(parentID, curve, edge, sameSense, trimSense);

				break;
			}

			// TODO: review and simplify
			case ifc::IFCLINE:
			{
				bool condition = sameSense == 1 || sameSense == -1;
				if (edge)
				{
					condition = !condition;
				}
				if constexpr (DIM == 2)
				{
					if (trim.start.hasPos && trim.end.hasPos)
					{

						if (condition)
						{
							curve.Add(trim.start.pos);
							curve.Add(trim.end.pos);
						}
						else
						{
							curve.Add(trim.end.pos);
							curve.Add(trim.start.pos);
						}
					}
					else if (trim.start.hasParam && trim.end.hasParam)
					{
						_loader.MoveToArgumentOffset(line, 0);
						auto positionID = _loader.GetRefArgument();
						auto vectorID = _loader.GetRefArgument();
						glm::dvec3 placement = glm::dvec3(GetCartesianPoint2D(positionID), 0);
						glm::dvec3 vector;
						vector = GetVector(vectorID);

						if (condition)
						{
							glm::dvec3 p1 = placement + vector * trim.start.param;
							glm::dvec3 p2 = placement + vector * trim.end.param;
							curve.Add(p1);
							curve.Add(p2);
						}
						else
						{
							glm::dvec3 p2 = placement + vector * trim.start.param;
							glm::dvec3 p1 = placement + vector * trim.end.param;
							curve.Add(p1);
							curve.Add(p2);
						}
					}
					else
					{
						_loader.ReportError({LoaderErrorType::UNSUPPORTED_TYPE, "Unsupported trimmingselect IFCLINE", line.expressID, line.ifcType});
					}
				}
				else if constexpr (DIM == 3)
				{
					if (trim.start.hasPos && trim.end.hasPos)
					{
						if (condition)
						{
							curve.Add(trim.start.pos3D);
							curve.Add(trim.end.pos3D);
						}
						else
						{
							curve.Add(trim.end.pos3D);
							curve.Add(trim.start.pos3D);
						}
					}
					else
					{
						_loader.ReportError({LoaderErrorType::UNSUPPORTED_TYPE, "Unsupported trimmingselect IFCLINE", line.expressID, line.ifcType});
					}
				}

				break;
			}
			case ifc::IFCTRIMMEDCURVE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto basisCurveID = _loader.GetRefArgument();
				auto trim1Set = _loader.GetSetArgument();
				auto trim2Set = _loader.GetSetArgument();

				auto senseAgreementS = _loader.GetStringArgument();
				auto trimmingPreference = _loader.GetStringArgument();

				auto trim1 = ParseTrimSelect(DIM, trim1Set);
				auto trim2 = ParseTrimSelect(DIM, trim2Set);

				IfcTrimmingArguments trim;

				trim.exist = true;
				trim.start = trim1;
				trim.end = trim2;

				bool senseAgreement = senseAgreementS == "T";

				if (trimSense == 0)
				{
					senseAgreement = !senseAgreement;
				}

				ComputeCurve<DIM>(basisCurveID, curve, edge, sameSense, senseAgreement, trim);

				break;
			}
			case ifc::IFCINDEXEDPOLYCURVE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto pts2DRef = _loader.GetRefArgument();

				_loader.MoveToArgumentOffset(line, 2);

				if (_loader.GetTokenType() != webifc::IfcTokenType::EMPTY)
				{
					_loader.Reverse();
					auto selfIntersects = _loader.GetStringArgument();

					if (selfIntersects == "T")
					{
						// TODO: this is probably bad news
						_loader.ReportError({LoaderErrorType::UNSPECIFIED, "Self intersecting ifcindexedpolycurve", line.expressID});
					}
				}

				if constexpr (DIM == 2)
				{
					_loader.MoveToArgumentOffset(line, 1);
					if (_loader.GetTokenType() != webifc::IfcTokenType::EMPTY)
					{
						auto pnSegment = ReadCurveIndices();
						for (auto &sg : pnSegment)
						{
							if (sg.type == "IFCLINEINDEX")
							{
								auto pts = ReadIfcCartesianPointList2D(pts2DRef);
								for (auto &pt : sg.indexs)
								{
									curve.Add(pts[pt - 1]);
								}
							}
							if (sg.type == "IFCARCINDEX")
							{
								auto pts = ReadIfcCartesianPointList2D(pts2DRef);
								IfcCurve<2> arc = BuildArc3Pt(pts[sg.indexs[0] - 1], pts[sg.indexs[1] - 1], pts[sg.indexs[2] - 1]);
								for (auto &pt : arc.points)
								{
									curve.Add(pt);
								}
							}
						}
					}
					else
					{
						auto pts = ReadIfcCartesianPointList2D(pts2DRef);
						for (auto &pt : pts)
						{
							curve.Add(pt);
						}
					}
				}
				else
				{
					_loader.ReportError({LoaderErrorType::UNSPECIFIED, "Parsing ifcindexedpolycurve in 3D is not possible", line.expressID});
				}

				break;
			}

			// TODO: review and simplify
			case ifc::IFCCIRCLE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto positionID = _loader.GetRefArgument();
				double radius = _loader.GetDoubleArgument();

				glm::mat<DIM + 1, DIM + 1, glm::f64, glm::defaultp> placement;
				if constexpr (DIM == 2)
				{
					placement = GetAxis2Placement2D(positionID);
				}
				else
				{
					placement = GetLocalPlacement(positionID);
				}

				double startDegrees = 0;
				double endDegrees = 360;

				if (trim.exist)
				{
					if (trim.start.hasParam && trim.end.hasParam)
					{
						startDegrees = trim.start.param;
						endDegrees = trim.end.param;
					}
					else if (trim.start.hasPos && trim.end.hasPos)
					{
						if constexpr (DIM == 2)
						{
							double xx = placement[2].x - trim.start.pos.x;
							double yy = placement[2].y - trim.start.pos.y;
							startDegrees = VectorToAngle(xx, yy);
							xx = placement[2].x - trim.end.pos.x;
							yy = placement[2].y - trim.end.pos.y;
							endDegrees = VectorToAngle(xx, yy);
						}
						else if constexpr (DIM == 3)
						{
							glm::dvec4 vecX = placement[0];
							glm::dvec4 vecY = placement[1];
							glm::dvec4 vecZ = placement[2];

							glm::dvec3 v1 = glm::dvec3(
								trim.start.pos3D.x - placement[3].x,
								trim.start.pos3D.y - placement[3].y,
								trim.start.pos3D.z - placement[3].z);
							glm::dvec3 v2 = glm::dvec3(
								trim.end.pos3D.x - placement[3].x,
								trim.end.pos3D.y - placement[3].y,
								trim.end.pos3D.z - placement[3].z);

							double dxS = vecX.x * v1.x + vecX.y * v1.y + vecX.z * v1.z;
							double dyS = vecY.x * v1.x + vecY.y * v1.y + vecY.z * v1.z;
							double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

							double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
							double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
							double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

							endDegrees = VectorToAngle(dxS, dyS) - 90;
							startDegrees = VectorToAngle(dxE, dyE) - 90;
						}
					}
				}

				double startRad = startDegrees / 180 * CONST_PI;
				double endRad = endDegrees / 180 * CONST_PI;
				double lengthDegrees = endDegrees - startDegrees;

				// unset or true
				if (trimSense == 1 || trimSense == -1)
				{
					if (lengthDegrees < 0)
					{
						lengthDegrees += 360;
					}
				}
				else
				{
					if (lengthDegrees > 0)
					{
						lengthDegrees -= 360;
					}
				}

				double lengthRad = lengthDegrees / 180 * CONST_PI;

				size_t startIndex = curve.points.size();

				const int numSegments = _loader.GetSettings().CIRCLE_SEGMENTS_HIGH;

				for (int i = 0; i < numSegments; i++)
				{
					double ratio = static_cast<double>(i) / (numSegments - 1);
					double angle = startRad + ratio * lengthRad;

					if (sameSense == 0)
					{
						angle = endRad - ratio * lengthRad;
					}

					glm::vec<DIM, glm::f64> pos;
					if constexpr (DIM == 2)
					{
						glm::dvec2 vec(0);
						vec[0] = radius * std::cos(angle);
						vec[1] = -radius * std::sin(angle); // not sure why we need this, but we apparently do
						pos = placement * glm::dvec3(vec, 1);
					}
					else
					{
						glm::dvec3 vec(0);
						vec[0] = radius * std::cos(angle);
						vec[1] = -radius * std::sin(angle); // negative or not???
						pos = placement * glm::dvec4(glm::dvec3(vec), 1);
					}
					curve.Add(pos);
				}

				// without a trim, we close the circle
				if (!trim.exist)
				{
					curve.Add(curve.points[startIndex]);
				}

				break;
			}

			// TODO: review and simplify
			case ifc::IFCELLIPSE:
			{
				_loader.MoveToArgumentOffset(line, 0);
				auto positionID = _loader.GetRefArgument();
				double radius1 = _loader.GetDoubleArgument();
				double radius2 = _loader.GetDoubleArgument();

				glm::mat<DIM + 1, DIM + 1, glm::f64, glm::defaultp> placement;
				if constexpr (DIM == 2)
				{
					placement = GetAxis2Placement2D(positionID);
				}
				else
				{
					placement = GetLocalPlacement(positionID);
				}

				double startDegrees = 0;
				double endDegrees = 360;

				if (trim.exist)
				{
					if (trim.start.hasParam && trim.end.hasParam)
					{
						startDegrees = trim.start.param;
						endDegrees = trim.end.param;
					}
					else if (trim.start.hasPos && trim.end.hasPos)
					{
						if constexpr (DIM == 2)
						{
							double xx = placement[2].x - trim.start.pos.x;
							double yy = placement[2].y - trim.start.pos.y;
							startDegrees = VectorToAngle(xx, yy);
							xx = placement[2].x - trim.end.pos.x;
							yy = placement[2].y - trim.end.pos.y;
							endDegrees = VectorToAngle(xx, yy);
						}
						else if constexpr (DIM == 3)
						{
							glm::dvec4 vecX = placement[0];
							glm::dvec4 vecY = placement[1];
							glm::dvec4 vecZ = placement[2];

							glm::dvec3 v1 = glm::dvec3(
								trim.start.pos3D.x - placement[3].x,
								trim.start.pos3D.y - placement[3].y,
								trim.start.pos3D.z - placement[3].z);
							glm::dvec3 v2 = glm::dvec3(
								trim.end.pos3D.x - placement[3].x,
								trim.end.pos3D.y - placement[3].y,
								trim.end.pos3D.z - placement[3].z);

							double dxS = vecX.x * v1.x + vecX.y * v1.y + vecX.z * v1.z;
							double dyS = vecY.x * v1.x + vecY.y * v1.y + vecY.z * v1.z;
							double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

							double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
							double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
							double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

							endDegrees = VectorToAngle(dxS, dyS) - 90;
							startDegrees = VectorToAngle(dxE, dyE) - 90;
						}
					}
				}

				double startRad = startDegrees / 180 * CONST_PI;
				double endRad = endDegrees / 180 * CONST_PI;

				// TODO: Because this is an ellipse you need to correct the angles

				// startRad = atan((radius1 / radius2) * tan(startDegrees));
				// endRad = atan((radius1 / radius2) * tan(endDegrees));

				double lengthDegrees = endDegrees - startDegrees;

				// unset or true
				if (trimSense == 1 || trimSense == -1)
				{
					if (lengthDegrees < 0)
					{
						lengthDegrees += 360;
					}
				}
				else
				{
					if (lengthDegrees > 0)
					{
						lengthDegrees -= 360;
					}
				}

				double lengthRad = lengthDegrees / 180 * CONST_PI;

				size_t startIndex = curve.points.size();

				const int numSegments = _loader.GetSettings().CIRCLE_SEGMENTS_HIGH;

				for (int i = 0; i < numSegments; i++)
				{
					double ratio = static_cast<double>(i) / (numSegments - 1);
					double angle = startRad + ratio * lengthRad;
					if (sameSense == 0)
					{
						angle = endRad - ratio * lengthRad;
					}
					glm::vec<DIM, glm::f64> pos;
					if constexpr (DIM == 2)
					{
						glm::dvec2 vec(0);
						vec[0] = radius1 * std::cos(angle);
						vec[1] = -radius2 * std::sin(angle);
						pos = placement * glm::dvec3(vec, 1);
					}
					else
					{
						glm::dvec3 vec(0);
						vec[0] = radius1 * std::cos(angle);
						vec[1] = -radius2 * std::sin(angle); // negative or not???
						pos = placement * glm::dvec4(glm::dvec3(vec), 1);
					}
					curve.Add(pos);
				}

				// without a trim, we close the circle
				if (!trim.exist)
				{
					curve.Add(curve.points[startIndex]);
				}

				break;
			}
			case ifc::IFCBSPLINECURVE:
			{
				bool condition = sameSense == 0;
				if (edge)
				{
					condition = !condition;
				}

				std::vector<glm::vec<DIM, glm::f64>> ctrolPts;
				std::vector<glm::f64> knotMultiplicities;
				std::vector<glm::f64> distinctKnots;
				std::vector<glm::f64> knots;
				std::vector<glm::f64> indexes;
				std::vector<glm::f64> weights;

				_loader.MoveToArgumentOffset(line, 0);
				double degree = _loader.GetDoubleArgument();
				auto points = _loader.GetSetArgument();
				auto curveType = _loader.GetStringArgument();
				auto closed = _loader.GetStringArgument();
				auto selfIntersect = _loader.GetStringArgument();

				for (auto &token : points)
				{
					uint32_t pointId = _loader.GetRefArgument(token);
					ctrolPts.push_back(GetCartesianPoint<DIM>(pointId));
				}

				// build default knots
				for (int k = 0; k < points.size() + degree + 1; k++)
				{
					knots.push_back(k);
				}

				if (knots.size() != ctrolPts.size() + degree + 1)
				{
					std::cout << "Error: Knots and control points do not match" << std::endl;
				}

				// build default weights vector
				for (int i = 0; i < ctrolPts.size(); i++)
				{
					weights.push_back(1);
				}

				std::vector<glm::vec<DIM, glm::f64>> tempPoints = GetRationalBSplineCurveWithKnots(degree, ctrolPts, knots, weights);
				for (int i = 0; i < tempPoints.size(); i++)
				{
					curve.Add(tempPoints[i]);
				}

				if (condition)
				{
					std::reverse(curve.points.begin(), curve.points.end());
				}

				break;
			}
			case ifc::IFCBSPLINECURVEWITHKNOTS:
			{
				bool condition = sameSense == 0;
				if (edge)
				{
					condition = !condition;
				}

				// The IfcBSplineCurveWithKnots is a spline curve parameterized by spline functions for which the knot values are explicitly given.
				// This is the type of b-spline curve for which the knot values are explicitly given. This subtype shall be used to represent non-uniform B-spline curves and may be used for other knot types.
				// Let L denote the number of distinct values amongst the d+k+2 knots in the knot list; L will be referred to as the ‘upper index on knots’. Let mj denote the multiplicity (i.e., number of repetitions) of the _j_th distinct knot.

				// All knot multiplicities except the first and the last shall be in the range 1,...,d; the first and last may have a maximum value of d + 1. In evaluating the basis functions, a knot u of, e.g., multiplicity 3 is interpreted as a sequence u, u, u,; in the knot array.
				std::vector<glm::vec<DIM, glm::f64>> ctrolPts;
				std::vector<glm::f64> knotMultiplicities;
				std::vector<glm::f64> distinctKnots;
				std::vector<glm::f64> knots;
				std::vector<glm::f64> indexes;
				std::vector<glm::f64> weights;

				_loader.MoveToArgumentOffset(line, 0);
				double degree = _loader.GetDoubleArgument();
				auto points = _loader.GetSetArgument();
				auto curveType = _loader.GetStringArgument();
				auto closed = _loader.GetStringArgument();
				auto selfIntersect = _loader.GetStringArgument();
				auto knotMultiplicitiesSet = _loader.GetSetArgument(); // The multiplicities of the knots. This list defines the number of times each knot in the knots list is to be repeated in constructing the knot array.
				auto knotSet = _loader.GetSetArgument();			   // The list of distinct knots used to define the B-spline basis functions.

				for (auto &token : points)
				{
					uint32_t pointId = _loader.GetRefArgument(token);
					ctrolPts.push_back(GetCartesianPoint<DIM>(pointId));
				}

				for (auto &token : knotMultiplicitiesSet)
				{
					knotMultiplicities.push_back(_loader.GetDoubleArgument(token));
				}

				for (auto &token : knotSet)
				{
					distinctKnots.push_back(_loader.GetDoubleArgument(token));
				}

				for (int k = 0; k < distinctKnots.size(); k++)
				{
					double knot = distinctKnots[k];
					for (int i = 0; i < knotMultiplicities[k]; i++)
					{
						knots.push_back(knot);
					}
				}

				if (knots.size() != ctrolPts.size() + degree + 1)
				{
					std::cout << "Error: Knots and control points do not match" << std::endl;
				}

				// build default weights vector
				for (int i = 0; i < ctrolPts.size(); i++)
				{
					weights.push_back(1);
				}

				std::vector<glm::vec<DIM, glm::f64>> tempPoints = GetRationalBSplineCurveWithKnots(degree, ctrolPts, knots, weights);
				for (int i = 0; i < tempPoints.size(); i++)
				{
					curve.Add(tempPoints[i]);
				}

				if (condition)
				{
					std::reverse(curve.points.begin(), curve.points.end());
				}

				break;
			}
			case ifc::IFCRATIONALBSPLINECURVEWITHKNOTS:
			{

				bool condition = sameSense == 0;
				if (edge)
				{
					condition = !condition;
				}

				std::vector<glm::vec<DIM, glm::f64>> ctrolPts;
				std::vector<glm::f64> distinctKnots;
				std::vector<glm::u32> knotMultiplicities;
				std::vector<glm::f64> knots;
				std::vector<glm::f64> weights;
				_loader.MoveToArgumentOffset(line, 0);
				double degree = _loader.GetDoubleArgument();
				auto points = _loader.GetSetArgument();
				auto curveType = _loader.GetStringArgument();
				auto closed = _loader.GetStringArgument();
				auto selfIntersect = _loader.GetStringArgument();
				auto knotMultiplicitiesSet = _loader.GetSetArgument(); // The multiplicities of the knots. This list defines the number of times each knot in the knots list is to be repeated in constructing the knot array.
				auto knotSet = _loader.GetSetArgument();
				auto knotSpec = _loader.GetStringArgument(); // The description of the knot type. This is for information only.
				auto weightsSet = _loader.GetSetArgument();

				for (auto &token : points)
				{
					uint32_t pointId = _loader.GetRefArgument(token);
					ctrolPts.push_back(GetCartesianPoint<DIM>(pointId));
				}

				for (auto &token : knotMultiplicitiesSet)
				{
					knotMultiplicities.push_back(_loader.GetDoubleArgument(token));
				}

				for (auto &token : knotSet)
				{
					distinctKnots.push_back(_loader.GetDoubleArgument(token));
				}

				for (auto &token : weightsSet)
				{
					weights.push_back(_loader.GetDoubleArgument(token));
				}

				for (int k = 0; k < distinctKnots.size(); k++)
				{
					double knot = distinctKnots[k];
					for (int i = 0; i < knotMultiplicities[k]; i++)
					{
						knots.push_back(knot);
					}
				}

				if (knots.size() != ctrolPts.size() + degree + 1)
				{
					std::cout << "Error: Knots and control points do not match" << std::endl;
				}

				std::vector<glm::vec<DIM, glm::f64>> tempPoints = GetRationalBSplineCurveWithKnots(degree, ctrolPts, knots, weights);
				for (int i = 0; i < tempPoints.size(); i++)
				{
					curve.Add(tempPoints[i]);
				}

				if (condition)
				{
					std::reverse(curve.points.begin(), curve.points.end());
				}

				break;
			}
			default:
				_loader.ReportError({LoaderErrorType::UNSUPPORTED_TYPE, "Unsupported curve type", line.expressID, line.ifcType});
				break;
			}

			if (DEBUG_DUMP_SVG)
			{
				if constexpr (DIM == 2)
				{
					DumpSVGCurve(curve.points, L"partial_curve.html");
				}
			}
		}*/

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

		//case ifc::IFCCONNECTEDFACESET:
		//case ifc::IFCCLOSEDSHELL:
		//case ifc::IFCOPENSHELL:
		//These cases are handled by getBrep()
		struct ParamsGetBrep
        {
            size_t numFaces;
            uint32_t* indices;
        };

		Geometry getBrep(ParamsGetBrep parameters)
		{

		}

		bool notPresent(glm::dvec3 pt, std::vector<glm::dvec3> points)
		{
			for (auto &pt2 : points)
			{
				if (pt.x == pt2.x && pt.y == pt2.y && pt.z == pt2.z)
				{
					return false;
				}
			}
			return true;
		}

		//case ifc::IFCPOLYLOOP:
		//case ifc::IFCEDGELOOP:
		struct ParamsGetLoop {
			bool isEdgeLoop;
			size_t numPoints;
			glm::dvec3* points;
		};

		IfcCurve3D GetLoop(ParamsGetLoop parameters)
		{
			IfcCurve3D curve;

			if (!parameters.isEdgeLoop)
			{
				if (parameters.numPoints > 0)
				{
					curve.points.reserve(parameters.numPoints);

					glm::dvec3 prevPoint = parameters.points[0];

					curve.points.push_back(prevPoint);

					for (size_t index = 1; index < parameters.numPoints; index++)
					{
						glm::dvec3 currentPoint = parameters.points[index];
						//trim repeats 
						if (currentPoint.x != prevPoint.x && currentPoint.y != prevPoint.y && currentPoint.z != prevPoint.z)
						{
							curve.points.push_back(parameters.points[index]);
						}

						prevPoint = currentPoint;
					}
				}
			} else 
			{
				//TODO: Handle edge loop 
				/*auto edges = _loader.GetSetArgument();
				int id = 0;

				for (auto &token : edges)
				{
					uint32_t edgeId = _loader.GetRefArgument(token);
					IfcCurve<3> edgeCurve = GetOrientedEdge(edgeId);

					// Important not to repeat the last point otherwise triangulation fails
					// if the list has zero points this is initial, no repetition is possible, otherwise we must check
					if (curve.points.size() == 0)
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


		struct ParamsGetBound
		{
			bool isFaceOuterBound;
			bool orient;
			ParamsGetLoop parametersGetLoop;
		};

		IfcBound3D GetBound(ParamsGetBound parameters)
		{
			
			bool orient = parameters.orient;

			IfcBound3D bound;
			bound.curve = GetLoop(parameters.parametersGetLoop);
			bound.orientation = orient;
			bound.type = (parameters.isFaceOuterBound) ?  IfcBoundType::OUTERBOUND : IfcBoundType::BOUND;

			if (!orient)
			{
				std::reverse(bound.curve.points.begin(), bound.curve.points.end());
			}

			return bound;

		}

		//case ifc::IFCINDEXEDPOLYGONALFACEWITHVOIDS:
		//case ifc::IFCINDEXEDPOLYGONALFACE:
		struct ParamsReadIndexedPolygonalFace {
			size_t numPoints;
			size_t numIndices;
			bool indexedPolygonalFaceWithVoids;
			glm::dvec3* points;
			uint32_t* indices;
		};


		std::vector<IfcBound3D> ReadIndexedPolygonalFace(ParamsReadIndexedPolygonalFace parameters)
		{
			std::vector<IfcBound3D> bounds;
			bounds.emplace_back();

			for (size_t index = 0; index < parameters.numIndices; index++)
			{
				uint32_t currentIndex = parameters.indices[index];

				glm::dvec3 point = parameters.points[currentIndex - 1];

				// I am not proud of this (I inherited this, will change - NC)
				bounds.back().curve.points.push_back(point);

			}

			if (!parameters.indexedPolygonalFaceWithVoids)
			{
				return bounds;
			}


			//TODO: // handle case IFCINDEXEDPOLYGONALFACEWITHVOIDS


		}

		void TriangulateBounds(Geometry &geometry, std::vector<IfcBound3D> &bounds)
		{
			if (bounds.size() == 1 && bounds[0].curve.points.size() == 3)
			{
				auto c = bounds[0].curve;

				size_t offset = geometry.numPoints;

				geometry.AddFace(c.points[0], c.points[1], c.points[2]);
			}
			else if (bounds.size() > 0 && bounds[0].curve.points.size() >= 3)
			{
				// bound greater than 4 vertices or with holes, triangulate
				// TODO: modify to use glm::dvec2 with custom accessors
				using Point = std::array<double, 2>;
				std::vector<std::vector<Point>> polygon;

				uint32_t offset = geometry.numPoints;

				// if more than one bound
				if (bounds.size() > 1)
				{
					// locate the outer bound index
					int outerIndex = -1;
					for (int i = 0; i < bounds.size(); i++)
					{
						if (bounds[i].type == IfcBoundType::OUTERBOUND)
						{
							outerIndex = i;
							break;
						}
					}

					if (outerIndex == -1)
					{
						printf("Expected outer bound!\n");
					}
					else
					{
						// swap the outer bound to the first position
						std::swap(bounds[0], bounds[outerIndex]);
					}
				}

				// if the first bound is not an outer bound now, this is unexpected
				if (bounds[0].type != IfcBoundType::OUTERBOUND)
				{
					printf("Expected outer bound first!\n");
				}

				glm::dvec3 v1, v2, v3;
				if (!GetBasisFromCoplanarPoints(bounds[0].curve.points, v1, v2, v3))
				{
					// these points are on a line
					printf("No basis found for brep!\n");
					return;
				}

				glm::dvec3 v12(glm::normalize(v3 - v2));
				glm::dvec3 v13(glm::normalize(v1 - v2));
				glm::dvec3 n = glm::normalize(glm::cross(v12, v13));
				v12 = glm::cross(v13, n);

				// check winding of outer bound
				IfcCurve<2> test;
				for (int i = 0; i < bounds[0].curve.points.size(); i++)
				{
					glm::dvec3 pt = bounds[0].curve.points[i];
					glm::dvec3 pt2 = pt - v1;

					glm::dvec2 proj(
						glm::dot(pt2, v12),
						glm::dot(pt2, v13));

					test.Add(proj);
				}

				// if the outer bound is clockwise under the current projection (v12,v13,n), we invert the projection
				if (!test.IsCCW())
				{
					n *= -1;
					std::swap(v12, v13);
				}

				for (auto &bound : bounds)
				{
					std::vector<Point> points;
					for (int i = 0; i < bound.curve.points.size(); i++)
					{
						glm::dvec3 pt = bound.curve.points[i];
						geometry.AddPoint(pt, n);

						// project pt onto plane of curve to obtain 2d coords
						glm::dvec3 pt2 = pt - v1;

						glm::dvec2 proj(
							glm::dot(pt2, v12),
							glm::dot(pt2, v13));

						points.push_back({proj.x, proj.y});
					}

					polygon.push_back(points);
				}

				std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);

				for (int i = 0; i < indices.size(); i += 3)
				{
					geometry.AddFace(offset + indices[i + 0], offset + indices[i + 1], offset + indices[i + 2]);
				}
			}
			else
			{
				printf("bad bound\n");
			}
		}

		//case ifc::IFCPOLYGONALFACESET:
		struct ParamsPolygonalFaceSet {
			size_t numPoints;
			size_t numIndices;
			uint32_t indicesPerFace;
			bool indexedPolygonalFaceWithVoids;
			glm::dvec3* points;
			uint32_t* indices;
		};
		Geometry getPolygonalFaceSetGeometry(ParamsPolygonalFaceSet parameters)
		{
			Geometry geom;
			std::vector<IfcBound3D> bounds;

			ParamsReadIndexedPolygonalFace readIndexedPolygonalFaceParameters;
			readIndexedPolygonalFaceParameters.numIndices = parameters.indicesPerFace;
			readIndexedPolygonalFaceParameters.indexedPolygonalFaceWithVoids = parameters.indexedPolygonalFaceWithVoids;
			readIndexedPolygonalFaceParameters.points = parameters.points;

			for (size_t faceIndex = 0; faceIndex < parameters.numIndices / parameters.indicesPerFace; faceIndex++)
			{
				readIndexedPolygonalFaceParameters.indices = &parameters.indices[faceIndex * parameters.indicesPerFace];
				bounds = ReadIndexedPolygonalFace(readIndexedPolygonalFaceParameters);

				TriangulateBounds(geom, bounds);

				bounds.clear();
			}

			return geom;
		}

    };
}