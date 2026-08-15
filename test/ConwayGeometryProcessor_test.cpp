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

void assertGeometriesEqual(
    const conway::geometry::Geometry& unpacked,
    const conway::geometry::Geometry& packed ) {

  ASSERT_EQ( unpacked.GetVertexCount(), packed.GetVertexCount() );
  ASSERT_EQ( unpacked.GetTriangleCount(), packed.GetTriangleCount() );

  for ( uint32_t where = 0; where < unpacked.GetVertexCount(); ++where ) {

    const glm::dvec3& left  = unpacked.GetPoint( where );
    const glm::dvec3& right = packed.GetPoint( where );

    ASSERT_EQ( left.x, right.x );
    ASSERT_EQ( left.y, right.y );
    ASSERT_EQ( left.z, right.z );
  }

  for ( uint32_t where = 0; where < unpacked.GetTriangleCount(); ++where ) {

    ASSERT_EQ(
        unpacked.triangles[ where ].vertices[ 0 ],
        packed.triangles[ where ].vertices[ 0 ] );
    ASSERT_EQ(
        unpacked.triangles[ where ].vertices[ 1 ],
        packed.triangles[ where ].vertices[ 1 ] );
    ASSERT_EQ(
        unpacked.triangles[ where ].vertices[ 2 ],
        packed.triangles[ where ].vertices[ 2 ] );
  }
}

void testPackedFacesetMatchesUnpacked() {
  conway::geometry::ConwayGeometryProcessor processor;

  std::vector< glm::dvec3 > points = {
      { 0, 0, 0 }, { 4, 0, 0 }, { 4, 4, 0 }, { 0, 4, 0 },
      { 1, 1, 0 }, { 3, 1, 0 }, { 3, 3, 0 }, { 1, 3, 0 },
  };

  conway::geometry::ConwayGeometryProcessor::ParamsGetPolygonalFaceSetGeometry
      unpacked;
  unpacked.points = points;
  unpacked.faces.resize( 2 );

  // Plain quad: one ring, face_starts = [0].
  unpacked.faces[ 0 ].indices     = { 1, 2, 3, 4 };
  unpacked.faces[ 0 ].face_starts = { 0 };

  // Voided quad: outer then hole, face_starts = [0, 4].
  unpacked.faces[ 1 ].indices     = { 1, 2, 3, 4, 5, 6, 7, 8 };
  unpacked.faces[ 1 ].face_starts = { 0, 4 };

  const uint32_t indices[] = {
      1, 2, 3, 4,
      1, 2, 3, 4, 5, 6, 7, 8,
  };
  const uint32_t faceOffsets[]  = { 0, 4, 12 };
  const uint32_t startIndices[] = { 0, 0, 4 };
  const uint32_t startOffsets[] = { 0, 1, 3 };

  const conway::geometry::Geometry fromUnpacked =
      processor.getPolygonalFaceSetGeometry( unpacked );
  const conway::geometry::Geometry fromPacked =
      processor.getPolygonalFaceSetGeometryPacked(
          points,
          indices,
          faceOffsets,
          3,
          startIndices,
          startOffsets );

  ASSERT_EQ( fromUnpacked.IsEmpty(), false );
  assertGeometriesEqual( fromUnpacked, fromPacked );
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

TEST(PackedFacesetMatchesUnpackedTest) { testPackedFacesetMatchesUnpacked(); }