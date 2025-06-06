#include "../conway_geometry/ConwayGeometryProcessor.h"

#include <TinyCppTest.hpp>
#include <glm/glm.hpp>

void testAddFaceToGeometry() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsAddFaceToGeometry
      paramsAddFaceToGeometry;
  conway::geometry::Geometry geometry;
  conwayGeometryProcessor.AddFaceToGeometry(paramsAddFaceToGeometry, geometry);
}

void testBoolSubtract() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  std::vector<conway::geometry::Geometry> geometryArr1;
  std::vector<conway::geometry::Geometry> geometryArr2;

  conwayGeometryProcessor.BoolSubtractLegacy(geometryArr1, geometryArr2);
}

void testGeometryToGltf() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::Geometry ifcGeometry;
  conway::geometry::Material ifcMaterial;

  std::vector<conway::geometry::Material> materials;
  materials.push_back(ifcMaterial);

  glm::dmat4 identityMatrix(1.0);
  std::vector<conway::geometry::IfcGeometryCollection> geometryCollection;
  conway::geometry::IfcGeometryCollection geometryCollectionSingle;
  geometryCollectionSingle.AddComponentWithTransform(&ifcGeometry,
                                                     identityMatrix);
  geometryCollection.push_back(geometryCollectionSingle);

  std::string testFilePath = "./test";
  conwayGeometryProcessor.GeometryToGltf(geometryCollection, materials, false, false,
                                         testFilePath, false);
}

void testGeometryToGlb() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::Geometry ifcGeometry;
  conway::geometry::Material ifcMaterial;

  std::vector<conway::geometry::Material> materials;
  materials.push_back(ifcMaterial);

  glm::dmat4 identityMatrix(1.0);
  std::vector<conway::geometry::IfcGeometryCollection> geometryCollection;
  conway::geometry::IfcGeometryCollection geometryCollectionSingle;
  geometryCollectionSingle.AddComponentWithTransform(&ifcGeometry,
                                                     identityMatrix);
  geometryCollection.push_back(geometryCollectionSingle);

  std::string testFilePath = "./test";
  conwayGeometryProcessor.GeometryToGltf(geometryCollection, materials, true, false,
                                         testFilePath, false);
}

void testGeometryToGltfDraco() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::Geometry ifcGeometry;
  conway::geometry::Material ifcMaterial;

  std::vector<conway::geometry::Material> materials;
  materials.push_back(ifcMaterial);

  std::string testFilePath = "./test";
  glm::dmat4 identityMatrix(1.0);
  std::vector<conway::geometry::IfcGeometryCollection> geometryCollection;
  conway::geometry::IfcGeometryCollection geometryCollectionSingle;
  geometryCollectionSingle.AddComponentWithTransform(&ifcGeometry,
                                                     identityMatrix);
  geometryCollection.push_back(geometryCollectionSingle);
  conwayGeometryProcessor.GeometryToGltf(geometryCollection, materials, false,
                                         true, testFilePath, false);
}

void testGeometryToGlbDraco() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::Geometry ifcGeometry;
  conway::geometry::Material ifcMaterial;

  std::vector<conway::geometry::Material> materials;
  materials.push_back(ifcMaterial);

  glm::dmat4 identityMatrix(1.0);
  std::vector<conway::geometry::IfcGeometryCollection> geometryCollection;
  conway::geometry::IfcGeometryCollection geometryCollectionSingle;
  geometryCollectionSingle.AddComponentWithTransform(&ifcGeometry,
                                                     identityMatrix);
  geometryCollection.push_back(geometryCollectionSingle);

  std::string testFilePath = "./test";
  conwayGeometryProcessor.GeometryToGltf(geometryCollection, materials, true, true,
                                         testFilePath, false);
}

void testGetAxis1Placement() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsAxis1Placement3D
      paramsAxis1Placement3D;
  conwayGeometryProcessor.GetAxis1Placement(paramsAxis1Placement3D);
}

void testGetAxis2Placement2D() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetAxis2Placement2D
      paramsAxis2Placement2D;
  conwayGeometryProcessor.GetAxis2Placement2D(paramsAxis2Placement2D);
}

void testGetAxis2Placement3D() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsAxis2Placement3D
      paramsAxis2Placement3D;
  conwayGeometryProcessor.GetAxis2Placement3D(paramsAxis2Placement3D);
}

void testGetBooleanResult() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetBooleanResult
      paramsGetBooleanResult;
  conwayGeometryProcessor.GetBooleanResult(&paramsGetBooleanResult);
}

void testGetBound() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetBound paramsGetBound;
  conwayGeometryProcessor.GetBound(paramsGetBound);
}

void testGetCartesianTransformationOperator3D() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::
      ParamsCartesianTransformationOperator3D
          paramsCartesianTransformationOperator3D;
  conwayGeometryProcessor.GetCartesianTransformationOperator3D(
      paramsCartesianTransformationOperator3D);
}

void testGetHalfSpaceSolid() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetHalfspaceSolid
      paramsGetHalfSpaceSolid;
  conwayGeometryProcessor.GetHalfSpaceSolid(paramsGetHalfSpaceSolid);
}

void testGetLocalPlacement() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsLocalPlacement
      paramsLocalPlacement;
  conwayGeometryProcessor.GetLocalPlacement(paramsLocalPlacement);
}

void testGetLoop() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetLoop paramsGetLoop;
  conwayGeometryProcessor.GetLoop(paramsGetLoop);
}

void testGetMappedItem() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetMappedItem
      paramsGetMappedItem;
  conwayGeometryProcessor.getMappedItem(paramsGetMappedItem);
}

void testGetPolygonalBoundedHalfspace() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalBoundedHalfspace
      paramsGetPolygonalBoundedHalfspace;
  conwayGeometryProcessor.GetPolygonalBoundedHalfspace(
      paramsGetPolygonalBoundedHalfspace);
}

void testGetPolygonalFaceSetGeometry() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
      paramsPolygonalFaceSet;
  conwayGeometryProcessor.getPolygonalFaceSetGeometry(paramsPolygonalFaceSet);
}

void testGetIfcSurface() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();

  conway::geometry::ConwayGeometryProcessor::ParamsGetSurface paramsSurface;
  conwayGeometryProcessor.GetSurface(paramsSurface);
}

void testReadIndexedPolygonalFace() {
  conway::geometry::ConwayGeometryProcessor conwayGeometryProcessor =
      conway::geometry::ConwayGeometryProcessor();
  std::vector<glm::vec3> points;
  conway::geometry::ConwayGeometryProcessor::IndexedPolygonalFace face;
  conway::geometry::ConwayGeometryProcessor::ParamsReadIndexedPolygonalFace
      paramsReadIndexedPolygonalFace(points, face);

  conwayGeometryProcessor.ReadIndexedPolygonalFace(
      paramsReadIndexedPolygonalFace);
}

TEST(AddFaceToGeometryFunctionalTest) { testAddFaceToGeometry(); }

TEST(BoolSubtractFunctionalTest) { testBoolSubtract(); }

TEST(GeometryToGltfFunctionalTest) { testGeometryToGltf(); }

TEST(GeometryToGlbFunctionalTest) { testGeometryToGlb(); }

TEST(GeometryToGltfDracoFunctionalTest) { testGeometryToGltfDraco(); }

TEST(GeometryToGlbDracoFunctionalTest) { testGeometryToGlbDraco(); }

TEST(GetAxis1PlacementFunctionalTest) { testGetAxis1Placement(); }

TEST(GetAxis2Placement2DFunctionalTest) { testGetAxis2Placement2D(); }

TEST(GetAxis2Placement3DFunctionalTest) { testGetAxis2Placement3D(); }

TEST(GetBooleanResultFunctionalTest) { testGetBooleanResult(); }

TEST(GetBoundFunctionalTest) { testGetBound(); }

TEST(GetCartesianTransformationOperator3DFunctionalTest) {
  testGetCartesianTransformationOperator3D();
}

TEST(GetHalfSpaceSolidFunctionalTest) { testGetHalfSpaceSolid(); }

TEST(GetLocalPlacementFunctionalTest) { testGetLocalPlacement(); }

TEST(GetLoopFunctionalTest) { testGetLoop(); }

TEST(GetMappedItemFunctionalTest) { testGetMappedItem(); }

TEST(GetPolygonalBoundedHalfspaceFunctionalTest) {
  testGetPolygonalBoundedHalfspace();
}

TEST(GetPolygonalFaceSetGeometryFunctionalTest) {
  testGetPolygonalFaceSetGeometry();
}

TEST(GetSurfaceFunctionalTest) { testGetIfcSurface(); }

TEST(ReadIndexedPolygonalFaceFunctionalTest) { testReadIndexedPolygonalFace(); }