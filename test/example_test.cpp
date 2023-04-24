#include <TinyCppTest.hpp>
#include <glm/glm.hpp>
#include "../conway_geometry/ConwayGeometryProcessor.h"

TEST (BasicTest)
{
	conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

	// example struct
	conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
		parametersAxis2Placement3D;
	parametersAxis2Placement3D.zAxisRef.x = 0.000;
	parametersAxis2Placement3D.zAxisRef.y = 0.000;
	parametersAxis2Placement3D.zAxisRef.z = 1.000;
	parametersAxis2Placement3D.normalizeZ = true;
	parametersAxis2Placement3D.xAxisRef.x = 1.000;
	parametersAxis2Placement3D.xAxisRef.y = 0.000;
	parametersAxis2Placement3D.xAxisRef.z = 0.000;
	parametersAxis2Placement3D.normalizeX = true;
	parametersAxis2Placement3D.position.x = 0.000;
	parametersAxis2Placement3D.position.y = 0.000;
	parametersAxis2Placement3D.position.z = 0.000;

	glm::dmat4 localPlacementMatrix =
		conwayGeometryProcessor.GetAxis2Placement3D(parametersAxis2Placement3D);
	ASSERT_EQ_EPS (parametersAxis2Placement3D.zAxisRef.x, 0.0, conway::geometry::EPS_TINY);
}

