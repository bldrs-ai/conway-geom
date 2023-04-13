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

#include "math/intersect-mesh-mesh-conway.h"
#include "math/bool-mesh-mesh-conway.h"

#include <tinynurbs/tinynurbs.h>

#include "ifc-schema.h"
#include "web-ifc.h"
#include "conway-util.h"

//draco
#include <draco/compression/config/compression_shared.h>
#include <draco/compression/encode.h>
#include <draco/compression/expert_encode.h>
#include <draco/core/cycle_timer.h>
#include <draco/io/file_utils.h>
#include <draco/io/mesh_io.h>
#include <draco/io/point_cloud_io.h>

//GLTFSDK
#include <GLTFSDK/GLTF.h>
#include <GLTFSDK/BufferBuilder.h>
#include <GLTFSDK/GLTFResourceWriter.h>
#include <GLTFSDK/GLBResourceWriter.h>
#include <GLTFSDK/IStreamWriter.h>
#include <GLTFSDK/Serialize.h>
#include <GLTFSDK/ExtensionsKHR.h>

//Note - A new IStreamWriter will need to be written for emscripten builds. 
namespace
{
    // The glTF SDK is decoupled from all file I/O by the IStreamWriter (and IStreamReader)
    // interface(s) and the C++ stream-based I/O library. This allows the glTF SDK to be used in
    // sandboxed environments, such as WebAssembly modules and UWP apps, where any file I/O code
    // must be platform or use-case specific.
	class StreamWriter : public Microsoft::glTF::IStreamWriter
    {
      public:
        StreamWriter( std::filesystem::path pathBase )
            : m_pathBase( std::move( pathBase ) )
        {
            assert( m_pathBase.has_root_path() );
        }

		const std::vector<std::string> getUris() const
		{
			return uris;
		}

        // Resolves the relative URIs of any external resources declared in the glTF manifest
        std::shared_ptr< std::ostream > GetOutputStream( const std::string& filename ) const override
        {
            auto stream = std::make_shared<std::ostringstream>();
			uris.push_back(filename);
        	return stream;
        }

      private:
        std::filesystem::path m_pathBase;
		mutable std::vector<std::string> uris;
    };
} // namespace

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

    enum GeometryProcessorSettings : int
	{

		CIRCLE_SEGMENTS_LOW = 5,
		CIRCLE_SEGMENTS_MEDIUM = 8,
		CIRCLE_SEGMENTS_HIGH = 12,
		BOOL_ABORT_THRESHOLD = 10000 // 10k verts
	};

namespace conway
{

	 struct DracoOptions
        {
            bool        isPointCloud               = false;
            int         posQuantizationBits        = 11;
            int         texCoordsQuantizationBits  = 8;
            bool        texCoordsDeleted           = false;
            int         normalsQuantizationBits    = 8;
            bool        normalsDeleted             = false;
            int         genericQuantizationBits    = 8;
            bool        genericDeleted             = false;
            int         compressionLevel           = 7;
            bool        preservePolygons           = false;
            bool        useMetadata                = false;
            bool        deduplicateInputValues     = true;
        };

    struct ComposedMesh
	{
		glm::dvec4 color;
		glm::dmat4 transformation;
		uint32_t expressID;
		bool hasGeometry = false;
		bool hasColor = false;
        IfcGeometry meshGeometry;
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

	bool COORDINATE_TO_ORIGIN = false;
	bool USE_FAST_BOOLS = true; //TODO: This needs to be fixed in the future to rely on elalish/manifold
	bool DUMP_CSG_MESHES = false;
	int CIRCLE_SEGMENTS_LOW = 5;
	int CIRCLE_SEGMENTS_MEDIUM = 8;
	int CIRCLE_SEGMENTS_HIGH = 12;
	bool MESH_CACHE = false;
	int BOOL_ABORT_THRESHOLD = 10000; // 10k verts

    /*
     * The purpose of this class is to process all of the geometry extraction / transform request queries as required by the IFC file format. 
     * This will later be extended to handle other 3D file formats. 
    */
	class ConwayGeometryProcessor
	{ 
        public:
		ConwayGeometryProcessor()
		{
		}

        //case ifc::IFCMAPPEDITEM:
		struct ParamsGetMappedItem
		{
			glm::dmat4 transformation;
			IfcComposedMesh ifcPresentationMesh;
		};

        IfcComposedMesh getMappedItem(ParamsGetMappedItem parameters)
        {
			IfcComposedMesh mesh;
			mesh.transformation = parameters.transformation;
			mesh.children.push_back(parameters.ifcPresentationMesh);

			return mesh;
        }

		IfcGeometry BoolSubtract(const IfcGeometry &firstGeom, const std::vector<IfcGeometry> &secondGeoms)
		{
			IfcGeometry result;
			IfcGeometry secondGeom;

			if (USE_FAST_BOOLS/*LoaderSettings::USE_FAST_BOOLS*/)
			{
				for (auto geom : secondGeoms)
				{
					if (geom.numFaces != 0)
					{
						if (secondGeom.numFaces == 0)
						{
							secondGeom = geom;
						}
						else
						{
							secondGeom = boolJoin(secondGeom, geom);
						}

						if (DUMP_CSG_MESHES/*_loader.GetSettings().DUMP_CSG_MESHES*/)
						{
							DumpIfcGeometry(geom, L"geom.obj");
						}
					}
				}
				if (firstGeom.numFaces == 0 || secondGeom.numFaces == 0)
				{
					printf("bool aborted due to empty source or target");

					// bail out because we will get strange meshes
					// if this happens, probably there's an issue parsing the mesh that occurred earlier
					return firstGeom;
				}

				IfcGeometry r1;
				IfcGeometry r2;

				intersectMeshMesh(firstGeom, secondGeom, r1, r2);

				if (DUMP_CSG_MESHES/*_loader.GetSettings().DUMP_CSG_MESHES*/)
				{
					DumpIfcGeometry(r1, L"r1.obj");
					DumpIfcGeometry(r2, L"r2.obj");
				}
				result = boolSubtract(r1, r2);
			}
			else
			{
				const int threshold = BOOL_ABORT_THRESHOLD/*LoaderSettings().BOOL_ABORT_THRESHOLD*/;
				std::vector<IfcGeometry> seconds;

				for (auto &geom : secondGeoms)
				{
					if (geom.numPoints < threshold)
					{
						seconds.push_back(geom);
					}
					else
					{
						printf("complex bool aborted due to BOOL_ABORT_THRESHOLD");
					}

					if (DUMP_CSG_MESHES/*_loader.GetSettings().DUMP_CSG_MESHES*/)
					{
						DumpIfcGeometry(geom, L"geom.obj");
					}
				}

				if (firstGeom.numPoints > threshold)
				{
					printf("complex bool aborted due to BOOL_ABORT_THRESHOLD");

					// bail out because we expect this operation to take too long
					return firstGeom;
				}

				if (firstGeom.numFaces == 0 || seconds.size() == 0)
				{
					printf("bool aborted due to empty source or target");

					// bail out because we will get strange meshes
					// if this happens, probably there's an issue parsing the mesh that occurred earlier
					return firstGeom;
				}

				result = boolMultiOp_Manifold(firstGeom, seconds);
			}

			if (DUMP_CSG_MESHES/*_loader.GetSettings().DUMP_CSG_MESHES*/)
			{
				DumpIfcGeometry(firstGeom, L"first.obj");
				DumpIfcGeometry(secondGeom, L"second.obj");
				DumpIfcGeometry(result, L"result.obj");
			}

			return result;
		}

		//case ifc::IFCBOOLEANCLIPPINGRESULT:
		//case ifc::IFCBOOLEANRESULT:
		struct ParamsGetBooleanResult 
		{
			IfcGeometry& firstMesh;
			IfcGeometry& secondMesh;
			bool clippingResult;
		};

		IfcGeometry GetBooleanResult(ParamsGetBooleanResult parameters)
		{
			//glm::dvec3 origin = GetOrigin(parameters.firstMesh, parameters.firstMeshTranslation);
			//auto normalizeMat = glm::translate(-origin);
			std::vector<IfcGeometry> flatSecondGeoms;
			flatSecondGeoms.push_back(parameters.secondMesh);

			IfcGeometry resultGeometry = BoolSubtract(parameters.firstMesh, flatSecondGeoms);

			return resultGeometry;

		} 

		IfcGeometry Extrude(IfcProfile profile, glm::dvec3 dir, double distance, glm::dvec3 cuttingPlaneNormal = glm::dvec3(0), glm::dvec3 cuttingPlanePos = glm::dvec3(0))
		{
			IfcGeometry geom;
			std::vector<bool> holesIndicesHash;

			// build the caps
			{
				using Point = std::array<double, 2>;
				int polygonCount = 1 + profile.holes.size(); // Main profile + holes
				std::vector<std::vector<Point>> polygon(polygonCount);

				glm::dvec3 normal = dir;

				for (int i = 0; i < profile.curve.points.size(); i++)
				{
					glm::dvec2 pt = profile.curve.points[i];
					glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0) + dir * distance, 1);

					geom.AddPoint(et, normal);
					polygon[0].push_back({pt.x, pt.y});
				}

				for (int i = 0; i < profile.curve.points.size(); i++)
				{
					holesIndicesHash.push_back(false);
				}

				for (int i = 0; i < profile.holes.size(); i++)
				{
					IfcCurve<2> hole = profile.holes[i];
					int pointCount = hole.points.size();

					for (int j = 0; j < pointCount; j++)
					{
						holesIndicesHash.push_back(j == 0);

						glm::dvec2 pt = hole.points[j];
						glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0) + dir * distance, 1);

						profile.curve.Add(pt);
						geom.AddPoint(et, normal);
						polygon[i + 1].push_back({pt.x, pt.y}); // Index 0 is main profile; see earcut reference
					}
				}

				std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);

				if (indices.size() < 3)
				{
					// probably a degenerate polygon
					printf("degenerate polygon in extrude");
					return geom;
				}

				uint32_t offset = 0;
				bool winding = GetWindingOfTriangle(geom.GetPoint(offset + indices[0]), geom.GetPoint(offset + indices[1]), geom.GetPoint(offset + indices[2]));
				bool flipWinding = !winding;

				for (int i = 0; i < indices.size(); i += 3)
				{
					if (flipWinding)
					{
						geom.AddFace(offset + indices[i + 0], offset + indices[i + 2], offset + indices[i + 1]);
					}
					else
					{
						geom.AddFace(offset + indices[i + 0], offset + indices[i + 1], offset + indices[i + 2]);
					}
				}

				offset += geom.numPoints;

				normal = -dir;

				for (int i = 0; i < profile.curve.points.size(); i++)
				{
					glm::dvec2 pt = profile.curve.points[i];
					glm::dvec4 et = glm::dvec4(glm::dvec3(pt, 0), 1);

					if (cuttingPlaneNormal != glm::dvec3(0))
					{
						et = glm::dvec4(glm::dvec3(pt, 0), 1);
						glm::dvec3 transDir = glm::dvec4(dir, 0);

						// project {et} onto the plane, following the extrusion normal
						double ldotn = glm::dot(transDir, cuttingPlaneNormal);
						if (ldotn == 0)
						{
							printf("0 direction in extrude\n");
						}
						else
						{
							glm::dvec3 dpos = cuttingPlanePos - glm::dvec3(et);
							double dist = glm::dot(dpos, cuttingPlaneNormal) / ldotn;
							// we want to apply dist, even when negative
							et = et + glm::dvec4(dist * transDir, 1);
						}
					}

					geom.AddPoint(et, normal);
				}

				for (int i = 0; i < indices.size(); i += 3)
				{
					if (flipWinding)
					{
						geom.AddFace(offset + indices[i + 0], offset + indices[i + 1], offset + indices[i + 2]);
					}
					else
					{
						geom.AddFace(offset + indices[i + 0], offset + indices[i + 2], offset + indices[i + 1]);
					}
				}
			}

			uint32_t capSize = profile.curve.points.size();
			for (int i = 1; i < capSize; i++)
			{
				// https://github.com/tomvandig/web-ifc/issues/5
				if (holesIndicesHash[i])
				{
					continue;
				}

				uint32_t bl = i - 1;
				uint32_t br = i - 0;

				uint32_t tl = capSize + i - 1;
				uint32_t tr = capSize + i - 0;

				// this winding should be correct
				geom.AddFace(geom.GetPoint(tl),
							 geom.GetPoint(br),
							 geom.GetPoint(bl));

				geom.AddFace(geom.GetPoint(tl),
							 geom.GetPoint(tr),
							 geom.GetPoint(br));
			}

			return geom;
		}

		//case ifc::IFCHALFSPACESOLID:
		struct ParamsGetHalfspaceSolid
		{
			IfcSurface& surface;
			bool flipWinding;
			double optionalLinearScalingFactor = 1.0;
		};

		IfcGeometry GetHalfSpaceSolid(ParamsGetHalfspaceSolid parameters)
		{
			glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);

			if (parameters.flipWinding)
			{
				extrusionNormal *= -1;
			}

			double d = EXTRUSION_DISTANCE_HALFSPACE_M / parameters.optionalLinearScalingFactor;

			IfcProfile profile;
			profile.isConvex = false;
			profile.curve = GetRectangleCurve(d, d, glm::dmat3(1));

			auto geom = Extrude(profile, extrusionNormal, d);

			// @Refactor: duplicate of extrudedareasolid
			if (parameters.flipWinding)
			{
				for (uint32_t i = 0; i < geom.numFaces; i++)
				{
					uint32_t temp = geom.indexData[i * 3 + 0];
					temp = geom.indexData[i * 3 + 0];
					geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
					geom.indexData[i * 3 + 1] = temp;
				}
			}

			return geom;
		}
		
		//case ifc::IFCPOLYGONALBOUNDEDHALFSPACE
		struct ParamsGetPolygonalBoundedHalfspace
		{
			IfcSurface& surface; 
			glm::dmat4 position;
			conway::IfcCurve<2> curve;
			bool halfSpaceInPlaneDirection; //agreement != "T"
			double optionalLinearScalingFactor = 1.0;
		};

		IfcGeometry GetPolygonalBoundedHalfspace(ParamsGetPolygonalBoundedHalfspace parameters)
		{
			if (!parameters.curve.IsCCW())
			{
				parameters.curve.Invert();
			}

			glm::dvec3 extrusionNormal = glm::dvec3(0, 0, 1);
			glm::dvec3 planeNormal = parameters.surface.transformation[2];
			glm::dvec3 planePosition = parameters.surface.transformation[3];

			glm::dmat4 invPosition = glm::inverse(parameters.position);
			glm::dvec3 localPlaneNormal = invPosition * glm::dvec4(planeNormal, 0);
			auto localPlanePos = invPosition * glm::dvec4(planePosition, 1);

			bool flipWinding = false;
			double extrudeDistance = EXTRUSION_DISTANCE_HALFSPACE_M / parameters.optionalLinearScalingFactor;
			bool halfSpaceInPlaneDirection = parameters.halfSpaceInPlaneDirection;
			bool extrudeInPlaneDirection = glm::dot(localPlaneNormal, extrusionNormal) > 0;
			bool ignoreDistanceInExtrude = (!halfSpaceInPlaneDirection && extrudeInPlaneDirection) || (halfSpaceInPlaneDirection && !extrudeInPlaneDirection);

			if (ignoreDistanceInExtrude)
			{
				// spec says this should be * 0, but that causes issues for degenerate 0 volume pbhs
				// hopefully we can get away by just inverting it
				extrudeDistance *= -1;
				flipWinding = true;
			}

			IfcProfile profile;
			profile.isConvex = false;
			profile.curve = parameters.curve;

			conway::IfcGeometry geom = Extrude(profile, extrusionNormal, extrudeDistance, localPlaneNormal, localPlanePos);
			// auto geom = Extrude(profile, surface.transformation, extrusionNormal, EXTRUSION_DISTANCE_HALFSPACE);

			// @Refactor: duplicate of extrudedareasolid
			if (flipWinding)
			{
				for (uint32_t i = 0; i < geom.numFaces; i++)
				{
					uint32_t temp = geom.indexData[i * 3 + 0];
					temp = geom.indexData[i * 3 + 0];
					geom.indexData[i * 3 + 0] = geom.indexData[i * 3 + 1];
					geom.indexData[i * 3 + 1] = temp;
				}
			}

			if (DUMP_CSG_MESHES/*_loader.GetSettings().DUMP_CSG_MESHES*/)
			{
				DumpIfcGeometry(geom, L"pbhs.obj");
			}

			return geom;
		}

		//case ifc::IFCREPRESENTATIONMAP
		//TODO: see if this is needed


		//case ifc::IFCCONNECTEDFACESET:
		//case ifc::IFCCLOSEDSHELL:
		//case ifc::IFCOPENSHELL:
		//These cases are handled by getBrep()
		struct ParamsAddFaceToGeometry
		{
			uint32_t boundsSize;
			uint32_t* indices;
			uint32_t indicesPerFace;
			IfcBound3D* boundsArray;
			bool advancedBrep;
			IfcSurface* surface;
		};

		struct ParamsGetBrep
        {
			uint32_t boundsSize;
			uint32_t indicesPerFace;
            size_t numIndices;
            uint32_t* indices;
			IfcBound3D* boundsArray;
			bool advancedBrep;
			IfcSurface* surface;
        };

		IfcGeometry getBrep(ParamsGetBrep parameters)
		{
			IfcGeometry geometry;
			//set parameters 
			ParamsAddFaceToGeometry paramsAddFaceToGeometry;
			paramsAddFaceToGeometry.boundsSize = parameters.boundsSize;
			paramsAddFaceToGeometry.indicesPerFace = parameters.indicesPerFace;
			paramsAddFaceToGeometry.boundsArray = parameters.boundsArray;
			paramsAddFaceToGeometry.advancedBrep = parameters.advancedBrep;
			paramsAddFaceToGeometry.surface = parameters.surface;


			for (size_t faceIndex = 0; faceIndex < parameters.numIndices / parameters.indicesPerFace; faceIndex++)
			{
				paramsAddFaceToGeometry.indices = &parameters.indices[faceIndex * parameters.indicesPerFace];

				
				AddFaceToGeometry(paramsAddFaceToGeometry, geometry);
			}

			return geometry;
		}

		//case ifc::IFCFACE:
		//case ifc::IFCADVANCEDFACE:

		void AddFaceToGeometry(ParamsAddFaceToGeometry parameters, IfcGeometry &geometry)
		{
			if (!parameters.advancedBrep)
			{
				std::vector<IfcBound3D> bounds3D(parameters.boundsSize);

				for (int i = 0; i < parameters.boundsSize; i++)
				{
					bounds3D[i] = parameters.boundsArray[i];
				}

				TriangulateBounds(geometry, bounds3D);
			}
			else
			{
				std::vector<IfcBound3D> bounds3D(parameters.boundsSize);

				for (int i = 0; i < parameters.boundsSize; i++)
				{
					bounds3D[i] = parameters.boundsArray[i];
				}

				//auto surface = GetSurface(surfRef);

				if (parameters.surface == nullptr)
				{
					printf("surface was nullptr\n");
					return;
				}

				auto surface = parameters.surface[0];

				// TODO: place the face in the surface and tringulate
				if (surface.BSplineSurface.Active)
				{
					TriangulateBspline(geometry, bounds3D, surface);
				}
				else if (surface.CylinderSurface.Active)
				{
					TriangulateCylindricalSurface(geometry, bounds3D, surface);
				}
				else if (surface.RevolutionSurface.Active)
				{
					TriangulateRevolution(geometry, bounds3D, surface);
				}
				else if (surface.ExtrusionSurface.Active)
				{
					TriangulateExtrusion(geometry, bounds3D, surface);
				}
				else
				{
					TriangulateBounds(geometry, bounds3D);
				}
			}
		}


		//case ifc::IFCFACEBASEDSURFACEMODEL:
		//case ifc::IFCSHELLBASEDSURFACEMODEL:
		struct ParamsGetSurfaceModel
		{
			uint32_t numShellRefs;
			ParamsGetBrep* shells;
		};

		struct ResponseGeometryArray
		{
			uint32_t numRepresentations;
			IfcGeometry* geometryArray;
		};

		ResponseGeometryArray GetSurfaceModel(ParamsGetSurfaceModel parameters)
		{
			ResponseGeometryArray response;
			std::vector<IfcGeometry> geometryArray;

			for (uint32_t shellIndex = 0; shellIndex < parameters.numShellRefs; shellIndex++)
			{
				geometryArray.push_back(getBrep(parameters.shells[shellIndex]));
			}

			response.numRepresentations = parameters.numShellRefs;
			response.geometryArray = geometryArray.data();

			return response;

		}

		//case ifc::IFCPLANE:
		//case ifc::IFCBSPLINESURFACE:
		//case ifc::IFCBSPLINESURFACEWITHKNOTS:
		//case ifc::IFCRATIONALBSPLINESURFACEWITHKNOTS:
		struct ParamsGetSurface 
		{
			bool isPlane;
			glm::dmat4 transformation;
			bool isBsplineSurface;
			double Udegree;
			double Vdegree;
			//TODO: How do we pass these across? 
			std::vector<std::vector<glm::vec<3, glm::f64>>> ctrolPts;
			std::string curveType;
			std::string closedU;
			std::string closedV;
			std::string selfIntersect;
			bool isBsplineSurfaceWithKnots;
			std::vector<glm::f64> UMultiplicity;
			std::vector<glm::f64> VMultiplicity;
			std::vector<glm::f64> UKnots;
			std::vector<glm::f64> VKnots;
			bool isRationalBsplineSurfaceWithKnots;
			std::vector<std::vector<glm::f64>> weightPts;
			bool isCylindricalSurface;
			double radius;
			bool isSurfaceOfRevolution;
			glm::dmat4 revolutionDirection;
			IfcProfile3D revolutionProfile;
			bool includeTransformation;
			bool isSurfaceOfLinearExtrusion;
			glm::dvec3 extrusionDirection;
			IfcProfile extrusionProfile;
			bool customLength;
			double length;
		};
		IfcSurface GetSurface(ParamsGetSurface parameters)
		{
			if (parameters.isPlane)
			{
				IfcSurface surface;

				surface.transformation = parameters.transformation;

				return surface;
			}
			else if ( parameters.isBsplineSurface )
			{
				IfcSurface surface;

				surface.BSplineSurface.Active = true;
				surface.BSplineSurface.UDegree = parameters.Udegree;
				surface.BSplineSurface.VDegree = parameters.Vdegree;
				surface.BSplineSurface.ControlPoints = parameters.ctrolPts;
				surface.BSplineSurface.ClosedU = parameters.closedU;
				surface.BSplineSurface.ClosedV = parameters.closedV;
				surface.BSplineSurface.CurveType = parameters.curveType;

				//TODO: Old implementation wasn't returning a surface for this case.
				return surface;
			}
			else if (parameters.isBsplineSurfaceWithKnots || parameters.isRationalBsplineSurfaceWithKnots)
			{
				IfcSurface surface;

				if (parameters.UKnots[parameters.UKnots.size() - 1] != (int)parameters.UKnots[parameters.UKnots.size() - 1])
				{
					for (uint32_t i = 0; i < parameters.UKnots.size(); i++)
					{
						parameters.UKnots[i] = parameters.UKnots[i] * (parameters.UKnots.size() - 1) / parameters.UKnots[parameters.UKnots.size() - 1];
					}
				}

				if (parameters.VKnots[parameters.VKnots.size() - 1] != (int)parameters.VKnots[parameters.VKnots.size() - 1])
				{
					for (uint32_t i = 0; i < parameters.VKnots.size(); i++)
					{
						parameters.VKnots[i] = parameters.VKnots[i] * (parameters.VKnots.size() - 1) / parameters.VKnots[parameters.VKnots.size() - 1];
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

				if ( parameters.isRationalBsplineSurfaceWithKnots )
				{
					surface.BSplineSurface.WeightPoints = parameters.weightPts;
				}

				return surface;
			}
			else if (parameters.isCylindricalSurface)
			{
				IfcSurface surface;

				surface.transformation = parameters.transformation;
				surface.CylinderSurface.Active = true;
				surface.CylinderSurface.Radius = parameters.radius;

				return surface;
			}
			else if ( parameters.isSurfaceOfRevolution )
			{
				IfcSurface surface;

				if (parameters.includeTransformation)
				{
					surface.transformation = parameters.transformation;
				}

				surface.RevolutionSurface.Active = true;
				surface.RevolutionSurface.Direction = parameters.revolutionDirection;
				surface.RevolutionSurface.Profile = parameters.revolutionProfile;

				return surface;
			}
			else if (parameters.isSurfaceOfLinearExtrusion)
			{
				IfcSurface surface;

				glm::dvec3 direction = parameters.extrusionDirection;

				double length = 0;
				if (parameters.customLength)
				{
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

		//case ifc::IFCPOLYLINE
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

		//case ifc::IFCCOMPOSITECURVE
		template <uint32_t DIM>
		void ComputeCompositeCurve(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			//todo: recursion 
		}

		//case ifc::IFCCOMPOSITECURVESEGMENT
		template <uint32_t DIM>
		void ComputeCompositeCurveSegment(IfcCurve<DIM> &curve, glm::vec<DIM, glm::f64> points, bool edge, int sameSense = -1)
		{
			//todo: recursion 
		}

		IfcCurve<2> BuildArc3Pt(const glm::dvec2 &p1, const glm::dvec2 &p2, const glm::dvec2 &p3)
		{
			double f1 = (p1.x * p1.x - p2.x * p2.x + p1.y * p1.y - p2.y * p2.y);
			double f2 = (p1.x * p1.x - p3.x * p3.x + p1.y * p1.y - p3.y * p3.y);
			double v = 2 * (p1.x - p2.x) * (p1.y - p3.y) - 2 * (p1.x - p3.x) * (p1.y - p2.y);

			double cenX = ((p1.y - p3.y) * f1 - (p1.y - p2.y) * f2) / v;
			double cenYa = (f2 - 2 * cenX * (p1.x - p3.x)) / (2 * (p1.y - p3.y));
			double cenYb = (f1 - 2 * cenX * (p1.x - p2.x)) / (2 * (p1.y - p2.y));
			double cenY = cenYa;
			if (isnan(cenY))
			{
				cenY = cenYb;
			}

			glm::dvec2 pCen;
			pCen.x = cenX;
			pCen.y = cenY;

			double radius = sqrt(pow(cenX - p1.x, 2) + pow(cenY - p1.y, 2));

			//Using geometrical subdivision to avoid complex calculus with angles

			std::vector<glm::dvec2> pointList;
			pointList.push_back(p1);
			pointList.push_back(p2);
			pointList.push_back(p3);

			while(pointList.size() < GeometryProcessorSettings::CIRCLE_SEGMENTS_MEDIUM)
			{
				std::vector<glm::dvec2> tempPointList;
				for (uint32_t j = 0; j < pointList.size() - 1; j++)
				{
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
			IfcCurve<2> curve;
			curve.points = pointList;

			return curve;
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

		//case ifc::IFCLINE
		template <uint32_t DIM>
		void ComputeCurveLine(IfcCurve<DIM> &curve, 
		bool edge,
		TrimmingArguments trim, 
		int sameSense = -1)
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

		//case ifc::IFCTRIMMEDCURVE
		template <uint32_t DIM>
		void ComputeTrimmedCurve(IfcCurve<DIM> &curve, 
		glm::vec<DIM, glm::f64> points, 
		bool edge, 
		int sameSense = -1)
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


		//case ifc::IFCINDEXEDPOLYCURVE
		template <uint32_t DIM>
		void ComputeIndexedPolycurve(IfcCurve<DIM> &curve, 
		std::vector<glm::dvec2> points, 
		uint32_t numSegments, 
		CurveSegmentIndex* segments, 
		bool curveLineIndex, 
		bool selfIntersects)
		{
			if (selfIntersects)
			{
				printf("Self intersecting ifcindexedpolycurve");
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
				printf("Parsing ifcindexedpolycurve in 3D is not possible");
			}
		}

		//case ifc::IFCCIRCLE
		//TODO: Rework this parameter list, I really don't like it. Probably should split all 2D and 3D functions and get rid of the DIM template.
		template <uint32_t DIM>
		void ComputeCircleCurve(IfcCurve<DIM> &curve, 
		glm::mat4 placementMat4, 
		glm::mat3 placementMat3, 
		double radius, 
		bool selfIntersects, 
		TrimmingArguments trim, 
		int trimSense = -1,
		int sameSense = -1)
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
						//double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

						double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
						double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
						//double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

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
			const int numSegments = GeometryProcessorSettings::CIRCLE_SEGMENTS_HIGH;

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

		//case ifc::IFCELLIPSE
		//TODO: Rework this parameter list, I really don't like it. Probably should split all 2D and 3D functions and get rid of the DIM template.
		template <uint32_t DIM>
		void ComputeEllipseCurve(IfcCurve<DIM> &curve, 
		glm::mat4 placementMat4, 
		glm::mat3 placementMat3, 
		double radius1, 
		double radius2, 
		TrimmingArguments trim, 
		int trimSense = -1, 
		int sameSense = -1 )
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
						//double dzS = vecZ.x * v1.x + vecZ.y * v1.y + vecZ.z * v1.z;

						double dxE = vecX.x * v2.x + vecX.y * v2.y + vecX.z * v2.z;
						double dyE = vecY.x * v2.x + vecY.y * v2.y + vecY.z * v2.z;
						//double dzE = vecZ.x * v2.x + vecZ.y * v2.y + vecZ.z * v2.z;

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
			const int numSegments = GeometryProcessorSettings::CIRCLE_SEGMENTS_HIGH;

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


		//case ifc::IFCBSPLINECURVE
		template <uint32_t DIM>
		void ComputeBsplineCurve(IfcCurve<DIM> &curve, 
		std::vector<glm::vec<DIM, glm::f64>> ctrolPts, 
		std::vector<glm::f64> knots, 
		std::vector<glm::f64> weights, 
		double degree,
		bool edge, 
		int sameSense = -1)
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

		/*
		* The IfcBSplineCurveWithKnots is a spline curve parameterized by spline functions for which the knot values are explicitly given.
		* This is the type of b-spline curve for which the knot values are explicitly given. This subtype shall be used to represent non-uniform B-spline curves and may be used for other knot types.
		* Let L denote the number of distinct values amongst the d+k+2 knots in the knot list; L will be referred to as the ‘upper index on knots’. Let mj denote the multiplicity (i.e., number of repetitions) of the _j_th distinct knot.

		* All knot multiplicities except the first and the last shall be in the range 1,...,d; the first and last may have a maximum value of d + 1. In evaluating the basis functions, a knot u of, e.g., multiplicity 3 is interpreted as a sequence u, u, u,; in the knot array.
		* @param ctrolPts == list GetCartesianPoint<DIM>(pointId)
		* @param knotMultiplicities == The multiplicities of the knots. This list defines the number of times each knot in the knots list is to be repeated in constructing the knot array.
		* @param distinctKnots == The list of distinct knots used to define the B-spline basis functions.
		*/
		//case ifc::IFCBSPLINECURVEWITHKNOTS
		template <uint32_t DIM>
		void ComputeBsplineCurveWithKnots(IfcCurve<DIM> &curve, 
		std::vector<glm::vec<DIM, glm::f64>> ctrolPts, 
		std::vector<glm::f64> knotMultiplicities,
		std::vector<glm::f64> distinctKnots,
		double degree,
		bool edge,
		int sameSense = -1)
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
			std::vector<glm::f64> knots;
			std::vector<glm::f64> weights;

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
		}

		/*
		* The IfcBSplineCurveWithKnots is a spline curve parameterized by spline functions for which the knot values are explicitly given.
		* This is the type of b-spline curve for which the knot values are explicitly given. This subtype shall be used to represent non-uniform B-spline curves and may be used for other knot types.
		* Let L denote the number of distinct values amongst the d+k+2 knots in the knot list; L will be referred to as the ‘upper index on knots’. Let mj denote the multiplicity (i.e., number of repetitions) of the _j_th distinct knot.

		* All knot multiplicities except the first and the last shall be in the range 1,...,d; the first and last may have a maximum value of d + 1. In evaluating the basis functions, a knot u of, e.g., multiplicity 3 is interpreted as a sequence u, u, u,; in the knot array.
		* @param ctrolPts == list GetCartesianPoint<DIM>(pointId)
		* @param knotMultiplicities == The multiplicities of the knots. This list defines the number of times each knot in the knots list is to be repeated in constructing the knot array.
		* @param distinctKnots == The list of distinct knots used to define the B-spline basis functions.
		* @param weights == list double 
		*/
		//case ifc::IFCRATIONALBSPLINECURVEWITHKNOTS
		template <uint32_t DIM>
		void ComputeRationalBsplineCurveWithKnots(IfcCurve<DIM> &curve, 
		std::vector<glm::vec<DIM, glm::f64>> ctrolPts, 
		std::vector<glm::f64> knotMultiplicities,
		std::vector<glm::f64> distinctKnots,
		std::vector<glm::f64> weights,
		double degree,
		bool edge,
		int sameSense = -1)
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
			std::vector<glm::f64> knots;

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
		}

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
			glm::vec3* points;
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
			} else 
			{
				//TODO: handle case IFCINDEXEDPOLYGONALFACEWITHVOIDS
			}


			
			return bounds;

		}

		void TriangulateRevolution(IfcGeometry &geometry, std::vector<IfcBound3D> &bounds, conway::IfcSurface &surface)
		{
			// First we get the revolution data

			glm::dvec3 cent = surface.RevolutionSurface.Direction[3];
			glm::dvec3 vecX = glm::normalize(surface.RevolutionSurface.Direction[0]);
			glm::dvec3 vecY = glm::normalize(surface.RevolutionSurface.Direction[1]);
			glm::dvec3 vecZ = glm::normalize(surface.RevolutionSurface.Direction[2]);

			std::vector<std::vector<glm::dvec3>> newPoints;

			double numRots = 10;

			for (int r = 0; r < numRots; r++)
			{
				std::vector<glm::dvec3> newList;
				newPoints.push_back(newList);
			}

			std::vector<glm::dvec3> bounding;
			std::vector<double> angleVec;
			std::vector<double> angleDsp;

			// Now we construct the bounding box of the boundary ...
			// ... by adding the middle point of all curves

			for (int i = 0; i < bounds.size(); i++)
			{
				double xx = 0;
				double yy = 0;
				double zz = 0;
				double cc = 0;
				int lastTeam = bounds[i].curve.indices[0];
				for (int j = 0; j < bounds[i].curve.points.size(); j++)
				{
					// If it is the first point of the group we close the previous group ...
					//  ... and create a new one. Else, the point is of the current group
					if (lastTeam != bounds[i].curve.indices[j] || j == (bounds[i].curve.points.size() - 1))
					{
						if (cc > 0)
						{
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
					}
					else
					{
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
			//  ... adds the angular differences again, reconstructing a corrected boundary.

			// Now we find the angle of each point in the reference plane of the cylinder
			for (int j = 0; j < bounding.size(); j++)
			{
				double xx = bounding[j].x - cent.x;
				double yy = bounding[j].y - cent.y;
				double zz = bounding[j].z - cent.z;
				double dx = vecX.x * xx + vecX.y * yy + vecX.z * zz;
				double dy = vecY.x * xx + vecY.y * yy + vecY.z * zz;
				//double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z * zz;
				double temp = VectorToAngle(dx, dy);
				while (temp < 0)
				{
					temp += 360;
				}
				while (temp > 360)
				{
					temp -= 360;
				}
				angleVec.push_back(temp);
			}

			for (int i = 0; i < angleVec.size() - 1; i++)
			{
				if (angleVec[i] - angleVec[i + 1] > 180)
				{
					angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
				}
				else if (angleVec[i] - angleVec[i + 1] < -180)
				{
					angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
				}
				else
				{
					angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
				}
			}

			double startDegrees = angleVec[0];
			double endDegrees = angleVec[0];

			// Add angular differences starting from the first angle. We also correct the start and end angles

			double temp = angleVec[0];
			for (int i = 0; i < angleDsp.size(); i++)
			{
				temp += angleDsp[i];
				if (endDegrees < temp)
				{
					endDegrees = temp;
				}
				if (startDegrees > temp)
				{
					startDegrees = temp;
				}
			}

			// Then we use the start and end angles as bounding boxes of the boundary ...
			//  ... we will represent this bounding box.

			double startRad = startDegrees / 180 * CONST_PI;
			double endRad = endDegrees / 180 * CONST_PI;
			double radSpan = endRad - startRad;
			double radStep = radSpan / (numRots - 1);

			for (int i = 0; i < surface.RevolutionSurface.Profile.curve.points.size(); i++)
			{
				double xx = surface.RevolutionSurface.Profile.curve.points[i].x - cent.x;
				double yy = surface.RevolutionSurface.Profile.curve.points[i].y - cent.y;
				double zz = surface.RevolutionSurface.Profile.curve.points[i].z - cent.z;

				double dx = vecX.x * xx + vecX.y * yy + vecX.z * zz;
				double dy = vecY.x * xx + vecY.y * yy + vecY.z * zz;
				double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z * zz;
				double dd = sqrt(dx * dx + dy * dy);
				for (int r = 0; r < numRots; r++)
				{
					double angle = startRad + r * radStep;
					double dtempX = sin(angle) * dd;
					double dtempY = cos(angle) * dd;
					double newPx = dtempX * vecX.x + dtempY * vecY.x + dz * vecZ.x + cent.x;
					double newPy = dtempX * vecX.y + dtempY * vecY.y + dz * vecZ.y + cent.y;
					double newPz = dtempX * vecX.z + dtempY * vecY.z + dz * vecZ.z + cent.z;
					glm::dvec3 newPt = glm::dvec3(
						newPx,
						newPy,
						newPz);
					newPoints[r].push_back(newPt);
				}
			}
			for (int r = 0; r < numRots - 1; r++)
			{
				int r1 = r + 1;
				for (int s = 0; s < newPoints[r].size() - 1; s++)
				{
					geometry.AddFace(newPoints[r][s], newPoints[r][s + 1], newPoints[r1][s]);
					geometry.AddFace(newPoints[r1][s], newPoints[r][s + 1], newPoints[r1][s + 1]);
				}
			}
		}

		void TriangulateCylindricalSurface(IfcGeometry &geometry, std::vector<IfcBound3D> &bounds, conway::IfcSurface &surface)
		{
			// First we get the cylinder data

			double radius = surface.CylinderSurface.Radius;
			glm::dvec3 cent = surface.transformation[3];
			glm::dvec3 vecX = glm::normalize(surface.transformation[0]);
			glm::dvec3 vecY = glm::normalize(surface.transformation[1]);
			glm::dvec3 vecZ = glm::normalize(surface.transformation[2]);

			std::vector<std::vector<glm::dvec3>> newPoints;

			double numRots = 10;
			double minZ = 1e+10;
			double maxZ = -1e+10;

			// Find the relative coordinates of each curve point in the cylinder reference plane
			// Only retain the max and min relative Z
			for (int i = 0; i < bounds.size(); i++)
			{
				for (int j = 0; j < bounds[i].curve.points.size(); j++)
				{
					glm::dvec3 vv = bounds[i].curve.points[j] - cent;
					//double dx = glm::dot(vecX, vv);
					//double dy = glm::dot(vecY, vv);
					double dz = glm::dot(vecZ, vv);
					if (maxZ < dz)
					{
						maxZ = dz;
					}
					if (minZ > dz)
					{
						minZ = dz;
					}
				}
			}

			for (int r = 0; r < numRots; r++)
			{
				std::vector<glm::dvec3> newList;
				newPoints.push_back(newList);
			}

			std::vector<glm::dvec3> bounding;
			std::vector<double> angleVec;
			std::vector<double> angleDsp;

			// Find the max. curve index in the boundary

			int maxTeam = 0;
			for (int i = 0; i < bounds.size(); i++)
			{
				for (int j = 0; j < bounds[i].curve.indices.size(); j++)
				{
					if (bounds[i].curve.indices[j] > maxTeam)
					{
						maxTeam = bounds[i].curve.indices[j];
					}
				}
			}

			std::vector<std::vector<glm::dvec3>> boundingGroups;

			// We group each point with their boundary

			for (int r = 0; r < maxTeam; r++)
			{
				std::vector<glm::dvec3> boundingTemp = std::vector<glm::dvec3>();
				for (int i = 0; i < bounds.size(); i++)
				{
					for (int j = 0; j < bounds[i].curve.points.size(); j++)
					{
						if (bounds[i].curve.indices[j] == r)
						{
							boundingTemp.push_back(bounds[i].curve.points[j]);
						}
					}
				}
				boundingGroups.push_back(boundingTemp);
			}

			int repeats = 0;
			bool start = false;
			bool end = false;
			int id = 0;

			// In the case of boundary lines having only 2 endings...
			//... we omit these lines and add solely curves having > 2 points...
			//... starting from a 2 point line, by doing it this way we don't have repeated points
			while (!end && repeats < maxTeam * 3)
			{
				if (id >= boundingGroups.size())
				{
					id = 0;
				}
				if (boundingGroups[id].size() < 3)
				{
					if (!start)
					{
						start = true;
					}
					else
					{
						break;
					}
				}
				if (boundingGroups[id].size() > 2 && start)
				{
					for (int i = 0; i < boundingGroups[id].size(); i++)
					{
						bounding.push_back(boundingGroups[id][i]);
					}
				}
				id++;
				repeats++;
			}

			// If the previous method finds nothing, then we don't have straight lines ...
			//... then we add all boundary points directly
			if (bounding.size() == 0)
			{
				for (int j = 0; j < boundingGroups.size(); j++)
				{
					for (int i = 0; i < boundingGroups[j].size(); i++)
					{
						bounding.push_back(boundingGroups[j][i]);
					}
				}
			}

			double startDegrees = 0;
			double endDegrees = 360;

			// Now we project the points in the cylinder surface
			// There is a problem when points in the cylinder are around 0 degrees
			// Numerical instabilities can make these points to jump from 0 to 360
			// It causes lots of trouble when drawing the boundaries in the cylinder

			// The method presented here finds the angle of each point, measures the ...
			//  ... angular difference and then, if the difference is bigger than 180 ...
			//  ... corrects it to a lesser value. Finally it gets the first angle and ...
			//  ... adds the angular differences again, reconstructing a corrected boundary.

			// Now we find the angle of each point in the reference plane of the cylinder
			for (int j = 0; j < bounding.size(); j++)
			{
				glm::dvec3 vv = bounding[j] - cent;
				double dx = glm::dot(vecX, vv);
				double dy = glm::dot(vecY, vv);
				//double dz = glm::dot(vecZ, vv);
				double temp = VectorToAngle(dx, dy);
				while (temp < 0)
				{
					temp += 360;
				}
				while (temp > 360)
				{
					temp -= 360;
				}
				angleVec.push_back(temp);
			}

			// Then we find the angular difference
			for (int i = 0; i < angleVec.size() - 1; i++)
			{
				if (angleVec[i] - angleVec[i + 1] > 180)
				{
					angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
				}
				else if (angleVec[i] - angleVec[i + 1] < -180)
				{
					angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
				}
				else
				{
					angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
				}
			}

			startDegrees = angleVec[0];
			endDegrees = angleVec[0];

			// Add angular differences starting from the first angle. We also correct the start and end angles

			double temp = angleVec[0];
			for (int i = 0; i < angleDsp.size(); i++)
			{
				temp += angleDsp[i];
				if (endDegrees < temp)
				{
					endDegrees = temp;
				}
				if (startDegrees > temp)
				{
					startDegrees = temp;
				}
			}

			// Then we use the start and end angles as bounding boxes of the boundary ...
			//  ... we will represent this bounding box.

			while (startDegrees < -360)
			{
				startDegrees += 360;
			}
			double startRad = startDegrees / 180 * CONST_PI;
			double endRad = endDegrees / 180 * CONST_PI;
			double radSpan = endRad - startRad;
			double radStep = radSpan / (numRots - 1);

			for (int r = 0; r < numRots; r++)
			{
				double angle = startRad + r * radStep;
				double dtempX = sin(angle) * radius;
				double dtempY = cos(angle) * radius;
				double newPx = dtempX * vecX.x + dtempY * vecY.x + minZ * vecZ.x + cent.x;
				double newPy = dtempX * vecX.y + dtempY * vecY.y + minZ * vecZ.y + cent.y;
				double newPz = dtempX * vecX.z + dtempY * vecY.z + minZ * vecZ.z + cent.z;
				glm::dvec3 newPt = glm::dvec3(
					newPx,
					newPy,
					newPz);
				newPoints[r].push_back(newPt);
			}
			for (int r = 0; r < numRots; r++)
			{
				double angle = startRad + r * radStep;
				double dtempX = sin(angle) * radius;
				double dtempY = cos(angle) * radius;
				double newPx = dtempX * vecX.x + dtempY * vecY.x + maxZ * vecZ.x + cent.x;
				double newPy = dtempX * vecX.y + dtempY * vecY.y + maxZ * vecZ.y + cent.y;
				double newPz = dtempX * vecX.z + dtempY * vecY.z + maxZ * vecZ.z + cent.z;
				glm::dvec3 newPt = glm::dvec3(
					newPx,
					newPy,
					newPz);
				newPoints[r].push_back(newPt);
			}

			for (int r = 0; r < numRots - 1; r++)
			{
				int r1 = r + 1;
				for (int s = 0; s < newPoints[r].size() - 1; s++)
				{
					geometry.AddFace(newPoints[r][s], newPoints[r][s + 1], newPoints[r1][s]);
					geometry.AddFace(newPoints[r1][s], newPoints[r][s + 1], newPoints[r1][s + 1]);
				}
			}
		}

		// TODO: review and simplify
		void TriangulateExtrusion(IfcGeometry &geometry, std::vector<IfcBound3D> &bounds, conway::IfcSurface &surface)
		{
			// NO EXAMPLE FILES ABOUT THIS CASE

			// THIS IS A SIMPLE EXTRUSION, NOT TRIMMED

			double len = surface.ExtrusionSurface.Length;
			glm::dvec3 dir = surface.ExtrusionSurface.Direction;

			if (!surface.ExtrusionSurface.Profile.isComposite)
			{
				for (int j = 0; j < surface.ExtrusionSurface.Profile.curve.points.size() - 1; j++)
				{
					int j2 = j + 1;

					double npx = surface.ExtrusionSurface.Profile.curve.points[j].x + dir.x * len;
					double npy = surface.ExtrusionSurface.Profile.curve.points[j].y + dir.y * len;
					double npz = dir.z * len;
					glm::dvec3 nptj1 = glm::dvec3(
						npx,
						npy,
						npz);
					npx = surface.ExtrusionSurface.Profile.curve.points[j2].x + dir.x * len;
					npy = surface.ExtrusionSurface.Profile.curve.points[j2].y + dir.y * len;
					npz = dir.z * len;
					glm::dvec3 nptj2 = glm::dvec3(
						npx,
						npy,
						npz);
					geometry.AddFace(
						glm::dvec3(surface.ExtrusionSurface.Profile.curve.points[j], 0),
						glm::dvec3(surface.ExtrusionSurface.Profile.curve.points[j2], 0),
						nptj1);
					geometry.AddFace(
						glm::dvec3(surface.ExtrusionSurface.Profile.curve.points[j2], 0),
						nptj2,
						nptj1);
				}
			}
			else
			{
				for (uint32_t i = 0; i < surface.ExtrusionSurface.Profile.profiles.size(); i++)
				{
					for (int j = 0; j < surface.ExtrusionSurface.Profile.profiles[i].curve.points.size() - 1; j++)
					{
						int j2 = j + 1;

						double npx = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].x + dir.x * len;
						double npy = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j].y + dir.y * len;
						double npz = dir.z * len;
						glm::dvec3 nptj1 = glm::dvec3(
							npx,
							npy,
							npz);
						npx = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].x + dir.x * len;
						npy = surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2].y + dir.y * len;
						npz = dir.z * len;
						glm::dvec3 nptj2 = glm::dvec3(
							npx,
							npy,
							npz);
						geometry.AddFace(
							glm::dvec3(surface.ExtrusionSurface.Profile.profiles[i].curve.points[j], 0),
							glm::dvec3(surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2], 0),
							nptj1);
						geometry.AddFace(
							glm::dvec3(surface.ExtrusionSurface.Profile.profiles[i].curve.points[j2], 0),
							nptj2,
							nptj1);
					}
				}
			}
		}

		// TODO: review and simplify
		void TriangulateBspline(IfcGeometry &geometry, std::vector<IfcBound3D> &bounds, conway::IfcSurface &surface)
		{
			//double limit = 1e-4;

			// First: We define the Nurbs surface

			tinynurbs::RationalSurface3d srf;
			srf.degree_u = surface.BSplineSurface.UDegree;
			srf.degree_v = surface.BSplineSurface.VDegree;
			size_t num_u = surface.BSplineSurface.ControlPoints.size();
			size_t num_v = surface.BSplineSurface.ControlPoints[0].size();

			std::vector<glm::dvec3> controlPoints;
			for (std::vector<glm::dvec3> row : surface.BSplineSurface.ControlPoints)
			{
				for (glm::dvec3 point : row)
				{
					controlPoints.push_back({point.x, point.y, point.z});
				}
			}
			srf.control_points = tinynurbs::array2(num_u, num_v, controlPoints);

			std::vector<double> weights;
			for (std::vector<double> row : surface.BSplineSurface.Weights)
			{
				for (double weight : row)
				{
					weights.push_back(weight);
				}
			}
			if (weights.size() != num_u * num_v)
			{
				for (int i = 0; i < num_u * num_v; i++)
				{
					weights.push_back(1.0);
				}
			}
			srf.weights = tinynurbs::array2(num_u, num_v, weights);

			for (int i = 0; i < surface.BSplineSurface.UMultiplicity.size(); i++)
			{
				for (int r = 0; r < surface.BSplineSurface.UMultiplicity[i]; r++)
				{
					srf.knots_u.push_back(surface.BSplineSurface.UKnots[i]);
				}
			}

			for (int i = 0; i < surface.BSplineSurface.VMultiplicity.size(); i++)
			{
				for (int r = 0; r < surface.BSplineSurface.VMultiplicity[i]; r++)
				{
					srf.knots_v.push_back(surface.BSplineSurface.VKnots[i]);
				}
			}

			// If the NURBS surface is valid we continue

			if (tinynurbs::surfaceIsValid(srf))
			{

				// Find projected boundary using NURBS inverse evaluation

				using Point = std::array<double, 2>;
				std::vector<std::vector<Point>> uvBoundaryValues;

				std::vector<Point> points;
				for (int j = 0; j < bounds[0].curve.points.size(); j++)
				{
					glm::dvec3 pt = bounds[0].curve.points[j];
					glm::dvec2 pInv = BSplineInverseEvaluation(pt, srf);
					points.push_back({pInv.x, pInv.y});
				}
				uvBoundaryValues.push_back(points);

				// Triangulate projected boundary
				// Subdivide resulting triangles to increase definition
				// r indicates the level of subdivision, currently 3 you can increase it to 5

				std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(uvBoundaryValues);

				for (int r = 0; r < 3; r++)
				{
					std::vector<uint32_t> newIndices;
					std::vector<Point> newUVPoints;

					for (int i = 0; i < indices.size(); i += 3)
					{
						Point p0 = uvBoundaryValues[0][indices[i + 0]];
						Point p1 = uvBoundaryValues[0][indices[i + 1]];
						Point p2 = uvBoundaryValues[0][indices[i + 2]];

						newUVPoints.push_back(p0);
						newUVPoints.push_back(p1);
						newUVPoints.push_back(p2);

						Point p3 = {(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2};
						Point p4 = {(p0[0] + p2[0]) / 2, (p0[1] + p2[1]) / 2};
						Point p5 = {(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2};

						newUVPoints.push_back(p3);
						newUVPoints.push_back(p4);
						newUVPoints.push_back(p5);

						int offset = newUVPoints.size() - 6;

						newIndices.push_back(offset + 0);
						newIndices.push_back(offset + 3);
						newIndices.push_back(offset + 4);

						newIndices.push_back(offset + 3);
						newIndices.push_back(offset + 5);
						newIndices.push_back(offset + 4);

						newIndices.push_back(offset + 3);
						newIndices.push_back(offset + 1);
						newIndices.push_back(offset + 5);

						newIndices.push_back(offset + 4);
						newIndices.push_back(offset + 5);
						newIndices.push_back(offset + 2);
					}

					uvBoundaryValues[0] = newUVPoints;
					indices = newIndices;
				}

				for (int i = 0; i < indices.size(); i += 3)
				{
					Point p0 = uvBoundaryValues[0][indices[i + 0]];
					Point p1 = uvBoundaryValues[0][indices[i + 1]];
					Point p2 = uvBoundaryValues[0][indices[i + 2]];
					glm::dvec3 pt00 = tinynurbs::surfacePoint(srf, p0[0], p0[1]);
					glm::dvec3 pt01 = tinynurbs::surfacePoint(srf, p1[0], p1[1]);
					glm::dvec3 pt10 = tinynurbs::surfacePoint(srf, p2[0], p2[1]);
					geometry.AddFace(pt00, pt01, pt10);
				}
			}
		}

		void TriangulateBounds(IfcGeometry &geometry, std::vector<IfcBound3D> &bounds)
		{
			if (bounds.size() == 1 && bounds[0].curve.points.size() == 3)
			{
				auto c = bounds[0].curve;

				//size_t offset = geometry.numPoints;

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

		struct ResultsGltf {
			bool success = false;
			std::vector<std::string> bufferUris;
			std::vector<std::vector<uint8_t>> buffers;
		};
		ResultsGltf GeometryToGltf(conway::IfcGeometry geom, bool isGlb, bool outputDraco, std::string filePath, bool outputFile, glm::dmat4 transform = glm::dmat4(1))
		{

			ResultsGltf results;
			// Encode the geometry.
			draco::EncoderBuffer buffer;

			//create a mesh object
			std::unique_ptr< draco::Mesh > dracoMesh( new draco::Mesh() );

			try 
			{
				int32_t              pos_att_id      = -1;
				int32_t              tex_att_id      = -1;
				int32_t              material_att_id = -1;

				//this internally populates the vertex float array, current storage type is double 
				geom.GetVertexData();

				if ( outputDraco )
				{
					DracoOptions dracoOptions;

					//get number of faces
					uint32_t numFaces = 0;
					numFaces = geom.GetIndexDataSize() / 3;

					if ( numFaces > 0 )
					{
						//set number of faces
						dracoMesh->SetNumFaces( numFaces );

						//set number of indices
						dracoMesh->set_num_points( 3 * numFaces );
					}

					uint32_t numPositions = 0;
					numPositions          = geom.numPoints;

					// Add attributes if they are present in the input data.
					if ( numPositions > 0 )
					{
						draco::GeometryAttribute va;
						va.Init( draco::GeometryAttribute::POSITION, nullptr, 3, draco::DT_FLOAT32, false, sizeof( float ) * 3, 0 );
						pos_att_id = dracoMesh->AddAttribute( va, false, numPositions );
					}

					//TODO: support multiple materials at some point when we add that
					int32_t numMaterials = 1;
					int32_t numTexCoords = 0;

					//populate position attribute
					float   vertexVal[ 3 ];
					int32_t vertexCount = 0;

					for (uint32_t i = 0; i < geom.numPoints; i++)
					{
						glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
						vertexVal[0] = (float)t.x;
						vertexVal[1] = (float)t.y;
						vertexVal[2] = (float)t.z;
						dracoMesh->attribute( pos_att_id )->SetAttributeValue( draco::AttributeValueIndex( i ), vertexVal );
					}

					//no textures, just map vertices to face indices
					for ( size_t triangleIndex = 0; triangleIndex < numFaces; ++triangleIndex )
					{
						const uint32_t* triangle = &( geom.indexData.data()[ triangleIndex * 3 ] );

						uint32_t triangle0 = triangle[ 0 ];
						uint32_t triangle1 = triangle[ 1 ];
						uint32_t triangle2 = triangle[ 2 ];

						const draco::PointIndex vert_id_0( 3 * triangleIndex );
						const draco::PointIndex vert_id_1( 3 * triangleIndex + 1 );
						const draco::PointIndex vert_id_2( 3 * triangleIndex + 2 );
						const int               triangulated_index = 0;

						//map vertex to face index
						dracoMesh->attribute( pos_att_id )
							->SetPointMapEntry( vert_id_0, draco::AttributeValueIndex( triangle0 ) );
						dracoMesh->attribute( pos_att_id )
							->SetPointMapEntry( vert_id_1, draco::AttributeValueIndex( triangle1 ) );
						dracoMesh->attribute( pos_att_id )
							->SetPointMapEntry( vert_id_2, draco::AttributeValueIndex( triangle2 ) );
					}

					// Add faces with identity mapping between vertex and corner indices.
					// Duplicate vertices will get removed below.
					draco::Mesh::Face face;
					for ( draco::FaceIndex i( 0 ); i < numFaces; ++i )
					{
						for ( int c = 0; c < 3; ++c )
						{
							face[ c ] = 3 * i.value() + c;
						}
						dracoMesh->SetFace( i, face );
					}

		#ifdef DRACO_ATTRIBUTE_VALUES_DEDUPLICATION_SUPPORTED
					if ( dracoOptions.deduplicateInputValues )
					{
						dracoMesh->DeduplicateAttributeValues();
					}
		#endif

		#ifdef DRACO_ATTRIBUTE_INDICES_DEDUPLICATION_SUPPORTED
					dracoMesh->DeduplicatePointIds();
		#endif

					// Convert compression level to speed (that 0 = slowest, 10 = fastest).
					const int32_t dracoCompressionSpeed = 10 - dracoOptions.compressionLevel;

					//set up Draco Encoder
					draco::Encoder encoder;

					// Setup encoder options.
					if ( dracoOptions.posQuantizationBits > 0 )
					{
						encoder.SetAttributeQuantization( draco::GeometryAttribute::POSITION, 11);//dracoOptions.posQuantizationBits );
					}
					if ( dracoOptions.texCoordsQuantizationBits > 0 )
					{
						encoder.SetAttributeQuantization( draco::GeometryAttribute::TEX_COORD, 8);//dracoOptions.texCoordsQuantizationBits );
					}
					if ( dracoOptions.normalsQuantizationBits > 0 )
					{
						encoder.SetAttributeQuantization( draco::GeometryAttribute::NORMAL, 8);//dracoOptions.normalsQuantizationBits );
					}
					if ( dracoOptions.genericQuantizationBits > 0 )
					{
						encoder.SetAttributeQuantization( draco::GeometryAttribute::GENERIC, 8);//dracoOptions.genericQuantizationBits );
					}

					encoder.SetSpeedOptions( dracoCompressionSpeed, dracoCompressionSpeed );

					// Convert to ExpertEncoder that allows us to set per-attribute options.
					std::unique_ptr< draco::ExpertEncoder > expert_encoder;
					expert_encoder.reset( new draco::ExpertEncoder( *dracoMesh ) );
					expert_encoder->Reset( encoder.CreateExpertEncoderOptions( *dracoMesh ) );

					//set up timer
					draco::CycleTimer timer;

					timer.Start();

					const draco::Status status = expert_encoder->EncodeToBuffer( &buffer );

					timer.Stop();

					if ( !status.ok() )
					{
						printf( "Failed to encode the mesh: %s", status.error_msg() );
						results.success = false;
						return results;
					}

					printf("Encoded To Draco in %lld ms\n", timer.GetInMs());
				}


				//The Document instance represents the glTF JSON manifest
				Microsoft::glTF::Document document;

				if ( outputDraco )
				{
					document.extensionsRequired.insert( "KHR_draco_mesh_compression" );
					document.extensionsUsed.insert( "KHR_draco_mesh_compression" );
				}

				std::string base_filename = filePath.substr(filePath.find_last_of("/\\") + 1);

				//Create a Buffer - it will be the 'current' Buffer that all the BufferViews
				//created by this BufferBuilder will automatically reference

				// Pass the absolute path, without the filename, to the stream writer
				std::unique_ptr< StreamWriter >                    streamWriter = std::make_unique< StreamWriter >( std::filesystem::path(filePath).parent_path().c_str() );
				std::unique_ptr< Microsoft::glTF::ResourceWriter > resourceWriter;

				if ( !isGlb )
				{
					resourceWriter = std::make_unique< Microsoft::glTF::GLTFResourceWriter >( std::move( streamWriter ) );
				}
				else
				{
					resourceWriter = std::make_unique< Microsoft::glTF::GLBResourceWriter >( std::move( streamWriter ) );
				}

				std::string accessorIdIndices;
				std::string accessorIdPositions;
				std::string accessorIdUVs;

				// Use the BufferBuilder helper class to simplify the process of
				// constructing valid glTF Buffer, BufferView and Accessor entities
				Microsoft::glTF::BufferBuilder bufferBuilder( std::move( resourceWriter ) );

				//if Draco
				Microsoft::glTF::KHR::MeshPrimitives::DracoMeshCompression dracoMeshCompression;

				//no Draco
				// Create all the resource data (e.g. triangle indices and
				// vertex positions) that will be written to the binary buffer
				const char* bufferId = nullptr;

				// Specify the 'special' GLB buffer ID. This informs the GLBResourceWriter that it should use
				// the GLB container's binary chunk (usually the desired buffer location when creating GLBs)
				if ( isGlb )
				{
					bufferId = Microsoft::glTF::GLB_BUFFER_ID;
				}
				else
				{
					bufferId = base_filename.c_str();

					std::string bufferIdCopy = std::string(bufferId);
					bufferIdCopy +=  ".";
					bufferIdCopy += Microsoft::glTF::BUFFER_EXTENSION;
					results.bufferUris.push_back(bufferIdCopy);
				}

				// Create a Buffer - it will be the 'current' Buffer that all the BufferViews
				// created by this BufferBuilder will automatically reference
				bufferBuilder.AddBuffer( bufferId );

				// Add an Accessor for the indices and positions
				std::vector< float >    positions;
				positions.resize(geom.numPoints * 3);

				std::vector< float > minValues( 3U, std::numeric_limits< float >::max() );
				std::vector< float > maxValues( 3U, std::numeric_limits< float >::lowest() );

				const size_t positionCount = positions.size();

				for (uint32_t i = 0; i < geom.numPoints; i++)
				{
					glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
					positions[3 * i + 0] = (float)t.x;
					positions[3 * i + 1] = (float)t.y;
					positions[3 * i + 2] = (float)t.z;
				}


				// Accessor min/max properties must be set for vertex position data so calculate them here
				for ( size_t i = 0U, j = 0U; i < positionCount; ++i, j = ( i % 3U ) )
				{
					minValues[ j ] = std::min( positions[ i ], minValues[ j ] );
					maxValues[ j ] = std::max( positions[ i ], maxValues[ j ] );
				}

				if (outputDraco)
				{
					Microsoft::glTF::Accessor positionsAccessor;
					positionsAccessor.min           = minValues;
					positionsAccessor.max           = maxValues;
					positionsAccessor.count         = dracoMesh->num_points();
					positionsAccessor.componentType = Microsoft::glTF::ComponentType::COMPONENT_FLOAT;
					positionsAccessor.type          = Microsoft::glTF::AccessorType::TYPE_VEC3;
					positionsAccessor.id            = "0";
					accessorIdPositions             = positionsAccessor.id;
					dracoMeshCompression.attributes.emplace( "POSITION", dracoMesh->attribute( pos_att_id )->unique_id() );

					Microsoft::glTF::Accessor indicesAccessor;
					indicesAccessor.count         = dracoMesh->num_faces() * 3;
					indicesAccessor.componentType = Microsoft::glTF::COMPONENT_UNSIGNED_INT;
					indicesAccessor.type          = Microsoft::glTF::TYPE_SCALAR;
					indicesAccessor.id            = "1";
					accessorIdIndices             = indicesAccessor.id;

					//add position accessor first here
					document.accessors.Append( positionsAccessor, Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty );
					document.accessors.Append( indicesAccessor, Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty );
				}
				else
				{
					// Create a BufferView with a target of ELEMENT_ARRAY_BUFFER (as it will reference index
					// data) - it will be the 'current' BufferView that all the Accessors created by this
					// BufferBuilder will automatically reference
					bufferBuilder.AddBufferView( Microsoft::glTF::BufferViewTarget::ELEMENT_ARRAY_BUFFER );
				

					// Copy the Accessor's id - subsequent calls to AddAccessor may invalidate the returned reference
					accessorIdIndices = bufferBuilder.AddAccessor( geom.indexData, { Microsoft::glTF::TYPE_SCALAR, Microsoft::glTF::COMPONENT_UNSIGNED_INT } ).id;

					// Create a BufferView with target ARRAY_BUFFER (as it will reference vertex attribute data)
					bufferBuilder.AddBufferView( Microsoft::glTF::BufferViewTarget::ARRAY_BUFFER );

					//Add positions accessor
					accessorIdPositions = bufferBuilder.AddAccessor( positions, { Microsoft::glTF::TYPE_VEC3, Microsoft::glTF::COMPONENT_FLOAT, false, std::move( minValues ), std::move( maxValues ) } ).id;
				}

				// Construct a Material
				Microsoft::glTF::Material    material;
				material.doubleSided = true;

				Microsoft::glTF::TextureInfo textureInfo;

				// Add it to the Document and store the generated ID
				auto materialId = document.materials.Append( std::move( material ), Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty ).id;

				// Construct a MeshPrimitive. Unlike most types in glTF, MeshPrimitives are direct children
				// of their parent Mesh entity rather than being children of the Document. This is why they
				// don't have an ID member.
				Microsoft::glTF::MeshPrimitive meshPrimitive;
				meshPrimitive.materialId = materialId;

				if (outputDraco)
				{
					// Finally put the encoded data in place.
					auto bufferView                   = bufferBuilder.AddBufferView( buffer.data(), buffer.size() );
					dracoMeshCompression.bufferViewId = bufferView.id;

					// Add all of the Buffers, BufferViews and Accessors that were created using BufferBuilder to
					// the Document. Note that after this point, no further calls should be made to BufferBuilder
					bufferBuilder.Output( document );

					//dracoMeshCompression.bufferViewId = bufferView.id;
					const auto  extensionSerializer                                = Microsoft::glTF::KHR::GetKHRExtensionSerializer();
					std::string strDracoMesh                                       = Microsoft::glTF::KHR::MeshPrimitives::SerializeDracoMeshCompression( dracoMeshCompression, document, extensionSerializer );
					meshPrimitive.extensions[ "KHR_draco_mesh_compression" ]       = strDracoMesh;
					meshPrimitive.indicesAccessorId                                = accessorIdIndices;
					meshPrimitive.attributes[ Microsoft::glTF::ACCESSOR_POSITION ] = accessorIdPositions;

				}
				else 
				{
					meshPrimitive.indicesAccessorId                                = accessorIdIndices;
					meshPrimitive.attributes[ Microsoft::glTF::ACCESSOR_POSITION ] = accessorIdPositions;

					// Add all of the Buffers, BufferViews and Accessors that were created using BufferBuilder to
					// the Document. Note that after this point, no further calls should be made to BufferBuilder
					bufferBuilder.Output( document );
				}

				// Construct a Mesh and add the MeshPrimitive as a child
				Microsoft::glTF::Mesh mesh;
				mesh.primitives.push_back( std::move( meshPrimitive ) );
				// Add it to the Document and store the generated ID
				auto meshId = document.meshes.Append( std::move( mesh ), Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty ).id;

				// Construct a Node adding a reference to the Mesh
				Microsoft::glTF::Node node;
				node.meshId = meshId;

				// Add it to the Document and store the generated ID
				auto nodeId = document.nodes.Append( std::move( node ), Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty ).id;

				// Construct a Scene
				Microsoft::glTF::Scene scene;
				scene.nodes.push_back( nodeId );
				// Add it to the Document, using a utility method that also sets the Scene as the Document's default
				document.SetDefaultScene( std::move( scene ), Microsoft::glTF::AppendIdPolicy::GenerateOnEmpty );

				std::string manifest;

				try
				{
					// Serialize the glTF Document into a JSON manifest
					manifest = Serialize( document, Microsoft::glTF::SerializeFlags::Pretty );
				}
				catch ( const Microsoft::glTF::GLTFException& ex )
				{
					std::stringstream ss;

					printf( "Microsoft::glTF::Serialize failed: %s", ex.what() );
				}

				auto& genericResourceWriter = bufferBuilder.GetResourceWriter();

				if ( isGlb )
				{
					
					auto glbResourceWriter = static_cast< Microsoft::glTF::GLBResourceWriter* >( &genericResourceWriter );

					filePath += ".glb";

					glbResourceWriter->Flush( manifest, filePath ); // A GLB container isn't created until the GLBResourceWriter::Flush member function is called

					//get GLB string from ostringstream
					std::shared_ptr<std::ostream> stream_ = glbResourceWriter->GetOutputBuffer(filePath);
					std::ostream* streamPtr_ = stream_.get();
					std::ostringstream* ossPtr_ = (std::ostringstream*)(streamPtr_);

					std::string str = ossPtr_->str();
					std::string fileNameOnly = std::filesystem::path(filePath).filename().string();
					results.bufferUris.push_back(fileNameOnly);

					std::vector<uint8_t> glbBuffer(str.length());
					memcpy(glbBuffer.data(), str.c_str(), str.length());
					results.buffers.push_back(glbBuffer);

					//write file 
					if (outputFile)
					{
						std::ofstream outfile(filePath); // Open the file for writing

						if (outfile.is_open()) 
						{
							// Write the data from the stream to the file
							outfile << str;

							// Check if the file was written to successfully
							if (outfile.fail()) 
							{
								throw std::runtime_error("Failed to write to file " + filePath);
							} else
							{
								std::cout << "Data written to file " << filePath << " successfully." << std::endl;
							}

							// Close the file
							outfile.close();
						} 
						else 
						{
							printf("Unable to open file %s", filePath.c_str());
						}
					}
				}
				else
				{
					filePath += ".gltf";
					

					//with the new StreamWriter, we must grab each output buffer (ostringstream) and write
					//it to the disk
					for (int bufferUrlIndex = 0; bufferUrlIndex < results.bufferUris.size(); bufferUrlIndex++)
					{
						std::string bufferUri_ = results.bufferUris[bufferUrlIndex];

						std::shared_ptr<std::ostream> stream_ = genericResourceWriter.GetOutputBuffer(bufferUri_);
						std::ostream* streamPtr_ = stream_.get();
						std::ostringstream* ossPtr_ = (std::ostringstream*)(streamPtr_);
						
						std::string str = ossPtr_->str();
						std::vector<uint8_t> gltfBuffer(str.length());
						memcpy(gltfBuffer.data(), str.c_str(), str.length());

						results.buffers.push_back(gltfBuffer);

						if (outputFile)
						{
							std::ofstream outfile(bufferUri_); // Open the file for writing

							if (outfile.is_open()) 
							{

								printf("Writing: %s\n", bufferUri_.c_str());

								// Write the data from the stream to the file
								outfile << str;

								// Check if the file was written to successfully
								if (outfile.fail()) 
								{
									throw std::runtime_error("Failed to write to file " + bufferUri_);
								} else
								{
									std::cout << "Data written to file " << bufferUri_ << " successfully." << std::endl;
								}

								// Close the file
								outfile.close();
							} 
							else 
							{
								throw std::runtime_error("Unable to open file " + filePath);
							}
						}
					}

					//push back the name of the manifest (.gltf)
					std::string fileNameOnly = std::filesystem::path(filePath).filename().string();
					results.bufferUris.push_back(fileNameOnly);

					std::vector<uint8_t> gltfManifest(manifest.length());
					memcpy(gltfManifest.data(), manifest.c_str(), manifest.length());
					results.buffers.push_back(gltfManifest);

					//write file 
					if (outputFile)
					{
						printf("Writing: %s\n", filePath.c_str());

						std::ofstream outfile(filePath); // Open the file for writing

						if (outfile.is_open()) 
						{
							// Write the data from the stream to the file
							outfile << manifest.c_str();

							// Check if the file was written to successfully
							if (outfile.fail()) 
							{
								throw std::runtime_error("Failed to write to file " + filePath);
							} else
							{
								std::cout << "Data written to file " << filePath << " successfully." << std::endl;
							}

							// Close the file
							outfile.close();
						} 
						else 
						{
							printf("Unable to open file %s", filePath.c_str());
						}
					}
				}
				results.success = true;
				return results;
			}
			catch ( const std::exception& ex )
			{
				printf( "Couldn't write GLB file: %s", ex.what() );
				results.success = false;
				return results;
			}
		}

		std::string GeometryToObj(const IfcGeometry &geom, size_t &offset, glm::dmat4 transform = glm::dmat4(1))
		{	
			std::stringstream obj;

			double scale = 1.0;

			for (uint32_t i = 0; i < geom.numPoints; i++)
			{
				glm::dvec4 t = transform * glm::dvec4(geom.GetPoint(i), 1);
				obj << "v " << t.x * scale << " " << t.y * scale << " " << t.z * scale << "\n";
			}

			for (uint32_t i = 0; i < geom.numFaces; i++)
			{
				Face f = geom.GetFace(i);
				obj << "f " << (f.i0 + 1 + offset) << "// " << (f.i1 + 1 + offset) << "// " << (f.i2 + 1 + offset) << "//\n";
			}

			offset += geom.numPoints;

			return obj.str();
		}

		//case ifc::IFCPOLYGONALFACESET:
		struct ParamsPolygonalFaceSet {
			size_t numPoints;
			size_t numIndices;
			uint32_t indicesPerFace;
			bool indexedPolygonalFaceWithVoids;
			std::vector<glm::vec3> points;
			std::vector<uint32_t> indices;
		};
		IfcGeometry getPolygonalFaceSetGeometry(ParamsPolygonalFaceSet parameters)
		{
			IfcGeometry geom;
			std::vector<IfcBound3D> bounds;

			ParamsReadIndexedPolygonalFace readIndexedPolygonalFaceParameters;
			readIndexedPolygonalFaceParameters.numIndices = parameters.indicesPerFace;
			readIndexedPolygonalFaceParameters.indexedPolygonalFaceWithVoids = parameters.indexedPolygonalFaceWithVoids;
			readIndexedPolygonalFaceParameters.points = parameters.points.data();

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