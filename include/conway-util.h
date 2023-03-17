/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <optional>
#include <cstring>
#include <functional>

#include <glm/glm.hpp>

#include <tinynurbs/tinynurbs.h>

#define CONST_PI 3.141592653589793238462643383279502884L

namespace conway
{

    const double EPS_MINISCULE = 1e-12; // what?
	const double EPS_TINY = 1e-9;
	const double EPS_SMALL = 1e-6;
	const double EPS_BIG = 1e-4;

    struct Face
	{
		int i0;
		int i1;
		int i2;
	};

	struct Loop
	{
		bool hasOne;
		glm::dvec2 v1;
		glm::dvec2 v2;
	};

	constexpr int VERTEX_FORMAT_SIZE_FLOATS = 6;

	glm::dvec3 computeNormal(const glm::dvec3 v1, const glm::dvec3 v2, const glm::dvec3 v3)
	{
		glm::dvec3 v12(v2 - v1);
		glm::dvec3 v13(v3 - v1);

		glm::dvec3 norm = glm::cross(v12, v13);

		return glm::normalize(norm);
	}

	// just follow the ifc spec, damn
	bool computeSafeNormal(const glm::dvec3 v1, const glm::dvec3 v2, const glm::dvec3 v3, glm::dvec3 &normal, double eps = 0)
	{
		glm::dvec3 v12(v2 - v1);
		glm::dvec3 v13(v3 - v1);

		glm::dvec3 norm = glm::cross(v12, v13);

		double len = glm::length(norm);

		if (len <= eps)
		{
			return false;
		}

		normal = norm / len;

		return true;
	}

	bool IsInsideCenterExtents(const glm::dvec3 &pt, const glm::dvec3 &center, const glm::dvec3 &extents)
	{
		glm::dvec3 delta = pt - center;
		delta = glm::abs(delta);
		glm::dvec3 offset = delta - extents;

		return offset.x < EPS_SMALL && offset.y < EPS_SMALL && offset.z < EPS_SMALL;
	}
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

		Geometry Normalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

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

		Geometry DeNormalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

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

		void AddGeometry(Geometry geom)
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

    bool equals2d(glm::dvec2 A, glm::dvec2 B, double eps = 0)
	{
		return std::fabs(A.x - B.x) <= eps && std::fabs(A.y - B.y) <= eps;
	}

	bool equals(glm::dvec3 A, glm::dvec3 B, double eps = 0)
	{
		return std::fabs(A.x - B.x) <= eps && std::fabs(A.y - B.y) <= eps && std::fabs(A.z - B.z) <= eps;
	}

    template <uint32_t DIM>
	struct IfcCurve
	{
		std::vector<glm::vec<DIM, glm::f64>> points;
		inline void Add(const glm::vec<DIM, glm::f64> &pt)
		{
			if (points.empty())
			{
				points.push_back(pt);
			}
			else
			{
				bool add = false;
				if constexpr (DIM == 2)
				{
					add = !equals2d(pt, points.back(), EPS_TINY);
				}
				else
				{
					add = !equals(pt, points.back(), EPS_TINY);
				}

				if (add)
				{
					points.push_back(pt);
				}
				else
				{
					// TODO: we are now discarding these points, but they probably should not even be generated
				}
			}
		}

		void Invert()
		{
			std::reverse(points.begin(), points.end());
		}

		bool IsCCW()
		{
			double sum = 0;

			for (int i = 0; i < points.size(); i++)
			{
				glm::dvec2 pt1 = points[(i - 1) % points.size()];
				glm::dvec2 pt2 = points[i];

				sum += (pt2.x - pt1.x) * (pt2.y + pt1.y);
			}

			return sum < 0;
		}
	};

    struct IfcProfile
	{
		std::string type;
		IfcCurve<2> curve;
		std::vector<IfcCurve<2>> holes;
		bool isConvex;
		bool isComposite = false;
		std::vector<IfcProfile> profiles;
	};

	struct IfcProfile3D
	{
		std::string type;
		IfcCurve<3> curve;
		bool isConvex;
	};

	struct Cylinder
	{
		bool Active = false;
		double Radius;
	};

	struct BSpline
	{
		bool Active = false;
		double UDegree;
		double VDegree;
		std::string ClosedU;
		std::string ClosedV;
		std::string CurveType;
		std::vector<std::vector<double>> Weights;
		std::vector<std::vector<glm::dvec3>> ControlPoints;
		std::vector<glm::f64> UMultiplicity;
		std::vector<glm::f64> VMultiplicity;
		std::vector<glm::f64> UKnots;
		std::vector<glm::f64> VKnots;
		std::vector<std::vector<glm::f64>> WeightPoints;
	};

	struct Revolution
	{
		bool Active = false;
		glm::dmat4 Direction;
		IfcProfile3D Profile;
	};

	struct Extrusion
	{
		bool Active = false;
		glm::dvec3 Direction;
		IfcProfile Profile;
		double Length;
	};

    struct IfcSurface
	{
		glm::dmat4 transformation;
		BSpline BSplineSurface;
		Cylinder CylinderSurface;
		Revolution RevolutionSurface;
		Extrusion ExtrusionSurface;

		glm::dvec3 normal()
		{
			if (!CylinderSurface.Active && !BSplineSurface.Active && !RevolutionSurface.Active)
			{
				return transformation[2];
			}
			else
			{
				if (BSplineSurface.Active)
				{
					printf("Normal to bspline still not implemented\n");
				}
				if (CylinderSurface.Active)
				{
					printf("Normal to cylinder still not implemented\n");
				}
				if (RevolutionSurface.Active)
				{
					printf("Normal to revolution still not implemented\n");
				}
				return glm::dvec3(0);
			}
		}
	};

    double VectorToAngle(double x, double y)
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
				angle = angle + (CONST_PI - angle) * 2;
				sinv = sin(angle);
				cosv = cos(angle);
				if (glm::abs(yy - cosv) > 1e-5 || glm::abs(xx - sinv) > 1e-5)
				{
					angle = angle + CONST_PI;
				}
			}
		}

		return (angle / (2 * CONST_PI)) * 360;
	}

    // TODO: review and simplify
	glm::dvec2 BSplineInverseEvaluation(glm::dvec3 pt, tinynurbs::RationalSurface3d srf)
	{
		// Initial data

		glm::highp_dvec3 ptc = tinynurbs::surfacePoint(srf, 0.0, 0.0);
		glm::highp_dvec3 pth = tinynurbs::surfacePoint(srf, 1.0, 0.0);
		glm::highp_dvec3 ptv = tinynurbs::surfacePoint(srf, 0.0, 1.0);

		double dh = glm::distance(ptc, pth);
		double dv = glm::distance(ptc, ptv);
		double pr = (dh + 1) / (dv + 1);

		double step1 = 0.01;
		double minError = 0.0001;
		double maxError = 0.01;
		double rotacions = 6;
		double stepOld = step1;

		// First approximation

		double fU = 0.5;
		double fV = 0.5;
		double divisor = 100;
		double maxdi = 1e+100;
		double extension = 0;

		while (maxdi > maxError && divisor < 10000)
		{
			for (double r = 1; r < 5; r++)
			{
				int round = 0;
				while (maxdi > minError && round < 3)
				{
					for (double i = 0; i < rotacions; i++)
					{
						double rads = (i / rotacions) * CONST_PI * 2;
						double incU = glm::sin(rads) / (r * r * divisor);
						double incV = glm::cos(rads) / (r * r * divisor);
						if (pr > 1)
						{
							incV *= pr;
						}
						else
						{
							incU /= pr;
						}
						bool repeat = true;
						while (repeat)
						{
							double ffU = fU + incU;
							double ffV = fV + incV;
							glm::highp_dvec3 pt00 = tinynurbs::surfacePoint(srf, ffU, ffV);
							double di = glm::distance(pt00, pt);
							if (di < maxdi)
							{
								maxdi = di;
								fU = ffU;
								fV = ffV;
							}
							else
							{
								repeat = false;
							}
						}
					}
					round++;
				}
			}
			divisor *= 3;
		}

		// If first method fails to provide a precise solution we use second slow but reliable method

		double repetition = 0;
		double maxdis = maxdi;
		double fUs = fU;
		double fVs = fV;
		while (maxdi > maxError && repetition < 8)
		{
			double extension = 1;
			double repetitionTemp = repetition;
			while (repetitionTemp > 4)
			{
				repetitionTemp -= 3;
				extension++;
			}
			if (repetitionTemp == 0)
			{
				fU = extension;
				fV = 0;
			}
			if (repetitionTemp == 1)
			{
				fU = 0;
				fV = extension;
			}
			if (repetitionTemp == 2)
			{
				fU = -extension;
				fV = 0;
			}
			if (repetitionTemp == 3)
			{
				fU = 0;
				fV = -extension;
			}

			maxdi = 1e+100;
			divisor = 100;
			rotacions = 6;
			while (maxdi > maxError && divisor < 10000)
			{
				for (double r = 1; r < 5; r++)
				{
					int round = 0;
					while (maxdi > minError && round < 3)
					{
						for (double i = 0; i < rotacions; i++)
						{
							double rads = (i / rotacions) * CONST_PI * 2;
							double incU = glm::sin(rads) / (r * r * divisor);
							double incV = glm::cos(rads) / (r * r * divisor);
							if (pr > 1)
							{
								incV *= pr;
							}
							else
							{
								incU /= pr;
							}
							bool repeat = true;
							while (repeat)
							{
								double ffU = fU + incU;
								double ffV = fV + incV;
								glm::highp_dvec3 pt00 = tinynurbs::surfacePoint(srf, ffU, ffV);
								double di = glm::distance(pt00, pt);
								if (di < maxdi)
								{
									maxdi = di;
									fU = ffU;
									fV = ffV;
									if (di < maxdis)
									{
										maxdis = di;
										fUs = ffU;
										fVs = ffV;
									}
								}
								else
								{
									repeat = false;
								}
							}
						}
						round++;
					}
				}
				divisor *= 3;
			}
			repetition++;
		}

		// If the second method fails then we go to the third method
		while (maxdi > maxError * 3 && repetition < 32)
		{
			double extension = 1;
			double repetitionTemp = repetition;
			while (repetitionTemp > 7)
			{
				repetitionTemp -= 8;
				extension++;
			}
			if (repetitionTemp == 0)
			{
				fU = extension;
				fV = 0;
			}
			if (repetitionTemp == 1)
			{
				fU = 0;
				fV = extension;
			}
			if (repetitionTemp == 2)
			{
				fU = -extension;
				fV = 0;
			}
			if (repetitionTemp == 3)
			{
				fU = 0;
				fV = -extension;
			}

			if (repetitionTemp == 4)
			{
				fU = extension * 0.707;
				fV = extension * 0.707;
			}
			if (repetitionTemp == 5)
			{
				fU = -extension * 0.707;
				fV = extension * 0.707;
			}
			if (repetitionTemp == 6)
			{
				fU = extension * 0.707;
				fV = -extension * 0.707;
			}
			if (repetitionTemp == 7)
			{
				fU = -extension * 0.707;
				fV = -extension * 0.707;
			}

			maxdi = 1e+100;
			divisor = 100;
			rotacions = 6;
			while (maxdi > maxError && divisor < 10000)
			{
				for (double r = 1; r < 5; r++)
				{
					int round = 0;
					while (maxdi > minError && round < 3)
					{
						for (double i = 0; i < rotacions; i++)
						{
							double rads = (i / rotacions) * CONST_PI * 2;
							double incU = glm::sin(rads) / (r * r * divisor);
							double incV = glm::cos(rads) / (r * r * divisor);
							if (pr > 1)
							{
								incV *= pr;
							}
							else
							{
								incU /= pr;
							}
							bool repeat = true;
							while (repeat)
							{
								double ffU = fU + incU;
								double ffV = fV + incV;
								glm::highp_dvec3 pt00 = tinynurbs::surfacePoint(srf, ffU, ffV);
								double di = glm::distance(pt00, pt);
								if (di < maxdi)
								{
									maxdi = di;
									fU = ffU;
									fV = ffV;
									if (di < maxdis)
									{
										maxdis = di;
										fUs = ffU;
										fVs = ffV;
									}
								}
								else
								{
									repeat = false;
								}
							}
						}
						round++;
					}
				}
				divisor *= 3;
			}
			repetition++;
		}

		return glm::dvec2(fUs, fVs);
	}

    struct IfcTrimmingSelect
	{
		bool hasParam = false;
		bool hasPos = false;
		double param;
		glm::dvec2 pos;
		glm::dvec3 pos3D;
	};

	struct IfcTrimmingArguments
	{
		bool exist = false;
		IfcTrimmingSelect start;
		IfcTrimmingSelect end;
	};

	struct IfcCurve3D
	{
		std::vector<glm::dvec3> points;
		std::vector<int> indices;
	};

    bool GetBasisFromCoplanarPoints(std::vector<glm::dvec3> &points, glm::dvec3 &v1, glm::dvec3 &v2, glm::dvec3 &v3)
	{
		v1 = points[0];

		for (auto &p : points)
		{
			if (v1 != p)
			{
				v2 = p;
				break;
			}
		}

		glm::dvec3 normal;
		// multiple tries to find the best match
		for (double i = 0; i < 4; i++)
		{
			double EPS = EPS_SMALL;
			if (i == 0)
			{
				EPS = 100;
			}
			if (i == 1)
			{
				EPS = 1;
			}
			if (i == 2)
			{
				EPS = 0.01;
			}
			for (auto &p : points)
			{
				if (computeSafeNormal(v1, v2, p, normal, EPS))
				{
					v3 = p;
					return true;
				}
			}
		}

		return false;
	}

    enum class IfcBoundType
	{
		OUTERBOUND,
		BOUND
	};

	// TODO: IfcBound3D can probably be merged with IfcProfile
	struct IfcBound3D
	{
		IfcBoundType type;
		bool orientation;
		IfcCurve3D curve;
	};

    struct IfcComposedMesh
	{
		glm::dvec4 color;
		glm::dmat4 transformation;
		uint32_t expressID;
		bool hasGeometry = false;
		bool hasColor = false;
		std::vector<IfcComposedMesh> children;

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

	glm::dvec3 GetOrigin(Geometry &geometry, glm::dmat4 transformation)
	{
		glm::dmat4 mat = glm::dmat4(1);
		glm::dmat4 newMat = mat * transformation;

		if (geometry.numFaces)
		{
			Face f = geometry.GetFace(0);
			glm::dvec3 a = newMat * glm::dvec4(geometry.GetPoint(f.i0), 1);
			return a;
		}

		return glm::dvec3(0);;
	}
}