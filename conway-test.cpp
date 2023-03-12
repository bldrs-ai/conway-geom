/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <iostream>
#include <fstream>
#include <filesystem>

#include "include/web-ifc.h"
#include "include/conway-geometry.h"
//#include "include/web-ifc-geometry.h"
#include "include/math/triangulate-with-boundaries.h"
#include "include/ifc-schema.h"

std::string ReadFile(std::string filename)
{
    std::ifstream t(filename);
    std::stringstream buffer;
    buffer << t.rdbuf();
    return buffer.str();
}

void DumpRefs(std::unordered_map<uint32_t, std::vector<uint32_t>> &refs)
{
    std::ofstream of("refs.txt");

    int32_t prev = 0;
    for (auto &it : refs)
    {
        if (!it.second.empty())
        {
            for (auto &i : it.second)
            {
                of << (((int32_t)i) - (prev));
                prev = i;
            }
        }
    }
}

struct BenchMarkResult
{
    std::string file;
    long long timeMS;
    long long sizeBytes;
};

void Benchmark()
{
    std::vector<BenchMarkResult> results;
    std::string path = "../../../benchmark/ifcfiles";
    for (const auto &entry : std::filesystem::directory_iterator(path))
    {
        if (entry.path().extension().string() != ".ifc")
        {
            continue;
        }

        std::string filePath = entry.path().string();
        std::string filename = entry.path().filename().string();

        std::string content = ReadFile(filePath);

        webifc::IfcLoader loader;
        auto start = webifc::ms();
        {
            // loader.LoadFile(content);
        }
        auto time = webifc::ms() - start;

        BenchMarkResult result;
        result.file = filename;
        result.timeMS = time;
        result.sizeBytes = entry.file_size();
        results.push_back(result);

        std::cout << "Reading " << result.file << " took " << time << "ms" << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl;

    double avgMBsec = 0;
    for (auto &result : results)
    {
        double MBsec = result.sizeBytes / 1000.0 / result.timeMS;
        avgMBsec += MBsec;
        std::cout << result.file << ": " << MBsec << " MB/sec" << std::endl;
    }

    avgMBsec /= results.size();

    std::cout << std::endl;
    std::cout << "Average: " << avgMBsec << " MB/sec" << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
}

void TestTriangleDecompose()
{
    const int NUM_TESTS = 100;
    const int PTS_PER_TEST = 100;
    const int EDGE_PTS_PER_TEST = 10;

    const double scaleX = 650;
    const double scaleY = 1;

    glm::dvec2 a(0, 0);
    glm::dvec2 b(scaleX, 0);
    glm::dvec2 c(0, scaleY);

    for (int i = 0; i < NUM_TESTS; i++)
    {
        srand(i);

        std::vector<glm::dvec2> points;

        // random points
        for (int j = 0; j < PTS_PER_TEST; j++)
        {
            points.push_back({webifc::RandomDouble(0, scaleX),
                              webifc::RandomDouble(0, scaleY)});
        }

        // points along the edges
        for (int j = 0; j < EDGE_PTS_PER_TEST; j++)
        {
            glm::dvec2 e1 = b - a;
            glm::dvec2 e2 = c - a;
            glm::dvec2 e3 = b - c;

            points.push_back(a + e1 * webifc::RandomDouble(0, 1));
            points.push_back(a + e2 * webifc::RandomDouble(0, 1));
            points.push_back(c + e3 * webifc::RandomDouble(0, 1));
        }

        std::vector<webifc::Loop> loops;

        for (auto &pt : points)
        {
            // if (pt.x > scaleX / 2)
            {
                webifc::Loop l;
                l.hasOne = true;
                l.v1 = pt;
                loops.push_back(l);
            }
        }

        std::cout << "Start test " << i << std::endl;

        bool swapped = false;
        auto triangles = webifc::triangulate(a, b, c, loops, swapped);

        // webifc::IsValidTriangulation(triangles, points);

        std::vector<webifc::Point> pts;

        for (auto &pt : points)
        {
            webifc::Point p;
            p.x = pt.x;
            p.y = pt.y;
            pts.push_back(p);
        }

        webifc::DumpSVGTriangles(triangles, webifc::Point(), webifc::Point(), L"triangles.svg", pts);
    }
}

int main()
{
    std::cout << "Hello Conway-Geom test!\n";


    std::vector<webifc::Geometry> geometryVec;

    //taken from web ifc obj dump code 
    glm::dmat4 NormalizeMat(
		glm::dvec4(1, 0, 0, 0),
		glm::dvec4(0, 0, -1, 0),
		glm::dvec4(0, 1, 0, 0),
		glm::dvec4(0, 0, 0, 1));
    
    //from outside debugger 
    //std::string content = ReadFile("../../../index.ifc");

    auto start = webifc::ms();

    webifc::ConwayGeometryProcessor conwayGeometryProcessor = webifc::ConwayGeometryProcessor();

    //generated from webifc_geom output 
    webifc::ConwayGeometryProcessor::ParamsAxis2Placement3D parametersAxis2Placement3D;
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

    glm::dmat4 localPlacementMatrix = conwayGeometryProcessor.GetAxis2Placement3D(parametersAxis2Placement3D);

    printf("localPlacementMatrix[0][1] = %.3f, localPlacementMatrix[0][2] = %.3f, localPlacementMatrix[0][3] = %.3f, localPlacementMatrix[0][4] = %.3f,\n", localPlacementMatrix[0][0], localPlacementMatrix[0][1], localPlacementMatrix[0][2], localPlacementMatrix[0][3]);
    printf("localPlacementMatrix[1][1] = %.3f, localPlacementMatrix[1][2] = %.3f, localPlacementMatrix[1][3] = %.3f, localPlacementMatrix[1][4] = %.3f,\n", localPlacementMatrix[1][0], localPlacementMatrix[1][1], localPlacementMatrix[1][2], localPlacementMatrix[1][3]);
    printf("localPlacementMatrix[2][1] = %.3f, localPlacementMatrix[2][2] = %.3f, localPlacementMatrix[2][3] = %.3f, localPlacementMatrix[2][4] = %.3f,\n", localPlacementMatrix[2][0], localPlacementMatrix[2][1], localPlacementMatrix[2][2], localPlacementMatrix[2][3]);
    printf("localPlacementMatrix[3][1] = %.3f, localPlacementMatrix[3][2] = %.3f, localPlacementMatrix[3][3] = %.3f, localPlacementMatrix[3][4] = %.3f,\n", localPlacementMatrix[3][0], localPlacementMatrix[3][1], localPlacementMatrix[3][2], localPlacementMatrix[3][3]);
    

    webifc::ConwayGeometryProcessor::ParamsPolygonalFaceSet parametersPolygonalFaceset;
    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 76.0000;
    parametersPolygonalFaceset.points[0].y = -11.4504;
    parametersPolygonalFaceset.points[0].z = 0.0000;

    parametersPolygonalFaceset.points[1].x = 76.0000;
    parametersPolygonalFaceset.points[1].y = -11.4504;
    parametersPolygonalFaceset.points[1].z = 15.0000;

    parametersPolygonalFaceset.points[2].x = 76.0000;
    parametersPolygonalFaceset.points[2].y = 0.0000;
    parametersPolygonalFaceset.points[2].z = 15.0000;

    parametersPolygonalFaceset.points[3].x = 76.0000;
    parametersPolygonalFaceset.points[3].y = 0.0000;
    parametersPolygonalFaceset.points[3].z = 0.0000;

    parametersPolygonalFaceset.points[4].x = 86.0000;
    parametersPolygonalFaceset.points[4].y = -11.4504;
    parametersPolygonalFaceset.points[4].z = 15.0000;

    parametersPolygonalFaceset.points[5].x = 86.0000;
    parametersPolygonalFaceset.points[5].y = -11.4504;
    parametersPolygonalFaceset.points[5].z = 0.0000;

    parametersPolygonalFaceset.points[6].x = 86.0000;
    parametersPolygonalFaceset.points[6].y = 0.0000;
    parametersPolygonalFaceset.points[6].z = 15.0000;

    parametersPolygonalFaceset.points[7].x = 86.0000;
    parametersPolygonalFaceset.points[7].y = 0.0000;
    parametersPolygonalFaceset.points[7].z = 0.0000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 5;
    parametersPolygonalFaceset.indices[5] = 2;
    parametersPolygonalFaceset.indices[6] = 1;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 2;
    parametersPolygonalFaceset.indices[9] = 5;
    parametersPolygonalFaceset.indices[10] = 7;
    parametersPolygonalFaceset.indices[11] = 3;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 6;
    parametersPolygonalFaceset.indices[17] = 1;
    parametersPolygonalFaceset.indices[18] = 4;
    parametersPolygonalFaceset.indices[19] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 5;
    parametersPolygonalFaceset.indices[21] = 6;
    parametersPolygonalFaceset.indices[22] = 8;
    parametersPolygonalFaceset.indices[23] = 7;

    webifc::Geometry geometry = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;

    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 48.0000;
    parametersPolygonalFaceset.points[0].y = -11.4504;
    parametersPolygonalFaceset.points[0].z = 0.0000;

    parametersPolygonalFaceset.points[1].x = 48.0000;
    parametersPolygonalFaceset.points[1].y = -11.4504;
    parametersPolygonalFaceset.points[1].z = 30.0000;

    parametersPolygonalFaceset.points[2].x = 48.0000;
    parametersPolygonalFaceset.points[2].y = 0.0000;
    parametersPolygonalFaceset.points[2].z = 30.0000;

    parametersPolygonalFaceset.points[3].x = 48.0000;
    parametersPolygonalFaceset.points[3].y = 0.0000;
    parametersPolygonalFaceset.points[3].z = 0.0000;

    parametersPolygonalFaceset.points[4].x = 58.0000;
    parametersPolygonalFaceset.points[4].y = -11.4504;
    parametersPolygonalFaceset.points[4].z = 0.0000;

    parametersPolygonalFaceset.points[5].x = 58.0000;
    parametersPolygonalFaceset.points[5].y = -11.4504;
    parametersPolygonalFaceset.points[5].z = 30.0000;

    parametersPolygonalFaceset.points[6].x = 58.0000;
    parametersPolygonalFaceset.points[6].y = 0.0000;
    parametersPolygonalFaceset.points[6].z = 30.0000;

    parametersPolygonalFaceset.points[7].x = 58.0000;
    parametersPolygonalFaceset.points[7].y = 0.0000;
    parametersPolygonalFaceset.points[7].z = 0.0000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 2;
    parametersPolygonalFaceset.indices[5] = 1;
    parametersPolygonalFaceset.indices[6] = 5;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 7;
    parametersPolygonalFaceset.indices[9] = 3;
    parametersPolygonalFaceset.indices[10] = 2;
    parametersPolygonalFaceset.indices[11] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 4;
    parametersPolygonalFaceset.indices[17] = 8;
    parametersPolygonalFaceset.indices[18] = 5;
    parametersPolygonalFaceset.indices[19] = 1;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 5;
    parametersPolygonalFaceset.indices[21] = 8;
    parametersPolygonalFaceset.indices[22] = 7;
    parametersPolygonalFaceset.indices[23] = 6;

    webifc::Geometry geometry2 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry2);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;


    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 0.0000;
    parametersPolygonalFaceset.points[0].y = -11.4504;
    parametersPolygonalFaceset.points[0].z = 0.0000;

    parametersPolygonalFaceset.points[1].x = 0.0000;
    parametersPolygonalFaceset.points[1].y = -11.4504;
    parametersPolygonalFaceset.points[1].z = 30.0000;

    parametersPolygonalFaceset.points[2].x = 0.0000;
    parametersPolygonalFaceset.points[2].y = 0.0000;
    parametersPolygonalFaceset.points[2].z = 30.0000;

    parametersPolygonalFaceset.points[3].x = 0.0000;
    parametersPolygonalFaceset.points[3].y = 0.0000;
    parametersPolygonalFaceset.points[3].z = 0.0000;

    parametersPolygonalFaceset.points[4].x = 10.0000;
    parametersPolygonalFaceset.points[4].y = -11.4504;
    parametersPolygonalFaceset.points[4].z = 0.0000;

    parametersPolygonalFaceset.points[5].x = 10.0000;
    parametersPolygonalFaceset.points[5].y = -11.4504;
    parametersPolygonalFaceset.points[5].z = 30.0000;

    parametersPolygonalFaceset.points[6].x = 10.0000;
    parametersPolygonalFaceset.points[6].y = 0.0000;
    parametersPolygonalFaceset.points[6].z = 30.0000;

    parametersPolygonalFaceset.points[7].x = 10.0000;
    parametersPolygonalFaceset.points[7].y = 0.0000;
    parametersPolygonalFaceset.points[7].z = 0.0000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 2;
    parametersPolygonalFaceset.indices[5] = 1;
    parametersPolygonalFaceset.indices[6] = 5;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 7;
    parametersPolygonalFaceset.indices[9] = 3;
    parametersPolygonalFaceset.indices[10] = 2;
    parametersPolygonalFaceset.indices[11] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 4;
    parametersPolygonalFaceset.indices[17] = 8;
    parametersPolygonalFaceset.indices[18] = 5;
    parametersPolygonalFaceset.indices[19] = 1;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 8;
    parametersPolygonalFaceset.indices[21] = 7;
    parametersPolygonalFaceset.indices[22] = 6;
    parametersPolygonalFaceset.indices[23] = 5;

    webifc::Geometry geometry3 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry3);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;

    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 0.00232305;
    parametersPolygonalFaceset.points[0].y = -12.647637;
    parametersPolygonalFaceset.points[0].z = 0.000000;

    parametersPolygonalFaceset.points[1].x = 0.00232305;
    parametersPolygonalFaceset.points[1].y = -24.098042;
    parametersPolygonalFaceset.points[1].z = 0.000000;

    parametersPolygonalFaceset.points[2].x = 0.00232305;
    parametersPolygonalFaceset.points[2].y = -24.098042;
    parametersPolygonalFaceset.points[2].z = 15.000000;

    parametersPolygonalFaceset.points[3].x = 0.00232305;
    parametersPolygonalFaceset.points[3].y = -12.647637;
    parametersPolygonalFaceset.points[3].z = 15.000000;

    parametersPolygonalFaceset.points[4].x = 10.002323;
    parametersPolygonalFaceset.points[4].y = -12.647637;
    parametersPolygonalFaceset.points[4].z = 0.000000;

    parametersPolygonalFaceset.points[5].x = 10.002323;
    parametersPolygonalFaceset.points[5].y = -24.098042;
    parametersPolygonalFaceset.points[5].z = 0.000000;

    parametersPolygonalFaceset.points[6].x = 10.002323;
    parametersPolygonalFaceset.points[6].y = -24.098042;
    parametersPolygonalFaceset.points[6].z = 15.000000;

    parametersPolygonalFaceset.points[7].x = 10.002323;
    parametersPolygonalFaceset.points[7].y = -12.647637;
    parametersPolygonalFaceset.points[7].z = 15.000000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 2;
    parametersPolygonalFaceset.indices[5] = 1;
    parametersPolygonalFaceset.indices[6] = 5;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 2;
    parametersPolygonalFaceset.indices[9] = 6;
    parametersPolygonalFaceset.indices[10] = 7;
    parametersPolygonalFaceset.indices[11] = 3;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 1;
    parametersPolygonalFaceset.indices[17] = 4;
    parametersPolygonalFaceset.indices[18] = 8;
    parametersPolygonalFaceset.indices[19] = 5;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 7;
    parametersPolygonalFaceset.indices[21] = 6;
    parametersPolygonalFaceset.indices[22] = 5;
    parametersPolygonalFaceset.indices[23] = 8;

    webifc::Geometry geometry4 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry4);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;


    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 24.0000;
    parametersPolygonalFaceset.points[0].y = -11.4504;
    parametersPolygonalFaceset.points[0].z = 0.0000;

    parametersPolygonalFaceset.points[1].x = 24.0000;
    parametersPolygonalFaceset.points[1].y = -11.4504;
    parametersPolygonalFaceset.points[1].z = 30.0000;

    parametersPolygonalFaceset.points[2].x = 24.0000;
    parametersPolygonalFaceset.points[2].y = 0.0000;
    parametersPolygonalFaceset.points[2].z = 30.0000;

    parametersPolygonalFaceset.points[3].x = 24.0000;
    parametersPolygonalFaceset.points[3].y = 0.0000;
    parametersPolygonalFaceset.points[3].z = 0.0000;

    parametersPolygonalFaceset.points[4].x = 34.0000;
    parametersPolygonalFaceset.points[4].y = -11.4504;
    parametersPolygonalFaceset.points[4].z = 30.0000;

    parametersPolygonalFaceset.points[5].x = 34.0000;
    parametersPolygonalFaceset.points[5].y = -11.4504;
    parametersPolygonalFaceset.points[5].z = 0.0000;

    parametersPolygonalFaceset.points[6].x = 34.0000;
    parametersPolygonalFaceset.points[6].y = 0.0000;
    parametersPolygonalFaceset.points[6].z = 30.0000;

    parametersPolygonalFaceset.points[7].x = 34.0000;
    parametersPolygonalFaceset.points[7].y = 0.0000;
    parametersPolygonalFaceset.points[7].z = 0.0000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 5;
    parametersPolygonalFaceset.indices[5] = 2;
    parametersPolygonalFaceset.indices[6] = 1;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 2;
    parametersPolygonalFaceset.indices[9] = 5;
    parametersPolygonalFaceset.indices[10] = 7;
    parametersPolygonalFaceset.indices[11] = 3;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 6;
    parametersPolygonalFaceset.indices[17] = 1;
    parametersPolygonalFaceset.indices[18] = 4;
    parametersPolygonalFaceset.indices[19] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 5;
    parametersPolygonalFaceset.indices[21] = 6;
    parametersPolygonalFaceset.indices[22] = 8;
    parametersPolygonalFaceset.indices[23] = 7;

    webifc::Geometry geometry5 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);


    geometryVec.push_back(geometry5);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;


    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 47.859639;
    parametersPolygonalFaceset.points[0].y = 0.973380;
    parametersPolygonalFaceset.points[0].z = 0.000000;

    parametersPolygonalFaceset.points[1].x = 47.859639;
    parametersPolygonalFaceset.points[1].y = 0.973380;
    parametersPolygonalFaceset.points[1].z = 15.000000;

    parametersPolygonalFaceset.points[2].x = 47.859639;
    parametersPolygonalFaceset.points[2].y = 12.423785;
    parametersPolygonalFaceset.points[2].z = 15.000000;

    parametersPolygonalFaceset.points[3].x = 47.859639;
    parametersPolygonalFaceset.points[3].y = 12.423785;
    parametersPolygonalFaceset.points[3].z = 0.000000;

    parametersPolygonalFaceset.points[4].x = 57.859639;
    parametersPolygonalFaceset.points[4].y = 0.973380;
    parametersPolygonalFaceset.points[4].z = 0.000000;

    parametersPolygonalFaceset.points[5].x = 57.859639;
    parametersPolygonalFaceset.points[5].y = 0.973380;
    parametersPolygonalFaceset.points[5].z = 15.000000;

    parametersPolygonalFaceset.points[6].x = 57.859639;
    parametersPolygonalFaceset.points[6].y = 12.423785;
    parametersPolygonalFaceset.points[6].z = 15.000000;

    parametersPolygonalFaceset.points[7].x = 57.859639;
    parametersPolygonalFaceset.points[7].y = 12.423785;
    parametersPolygonalFaceset.points[7].z = 0.000000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 2;
    parametersPolygonalFaceset.indices[5] = 1;
    parametersPolygonalFaceset.indices[6] = 5;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 7;
    parametersPolygonalFaceset.indices[9] = 3;
    parametersPolygonalFaceset.indices[10] = 2;
    parametersPolygonalFaceset.indices[11] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 4;
    parametersPolygonalFaceset.indices[17] = 8;
    parametersPolygonalFaceset.indices[18] = 5;
    parametersPolygonalFaceset.indices[19] = 1;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 8;
    parametersPolygonalFaceset.indices[21] = 7;
    parametersPolygonalFaceset.indices[22] = 6;
    parametersPolygonalFaceset.indices[23] = 5;

    webifc::Geometry geometry6 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry6);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;

    parametersPolygonalFaceset.numPoints = 8;
    parametersPolygonalFaceset.points = new glm::dvec3[parametersPolygonalFaceset.numPoints];
    parametersPolygonalFaceset.points[0].x = 62.0000;
    parametersPolygonalFaceset.points[0].y = -11.4504;
    parametersPolygonalFaceset.points[0].z = 0.0000;

    parametersPolygonalFaceset.points[1].x = 62.0000;
    parametersPolygonalFaceset.points[1].y = -11.4504;
    parametersPolygonalFaceset.points[1].z = 15.0000;

    parametersPolygonalFaceset.points[2].x = 62.0000;
    parametersPolygonalFaceset.points[2].y = 0.0000;
    parametersPolygonalFaceset.points[2].z = 15.0000;

    parametersPolygonalFaceset.points[3].x = 62.0000;
    parametersPolygonalFaceset.points[3].y = 0.0000;
    parametersPolygonalFaceset.points[3].z = 0.0000;

    parametersPolygonalFaceset.points[4].x = 72.0000;
    parametersPolygonalFaceset.points[4].y = -11.4504;
    parametersPolygonalFaceset.points[4].z = 15.0000;

    parametersPolygonalFaceset.points[5].x = 72.0000;
    parametersPolygonalFaceset.points[5].y = -11.4504;
    parametersPolygonalFaceset.points[5].z = 0.0000;

    parametersPolygonalFaceset.points[6].x = 72.0000;
    parametersPolygonalFaceset.points[6].y = 0.0000;
    parametersPolygonalFaceset.points[6].z = 15.0000;

    parametersPolygonalFaceset.points[7].x = 72.0000;
    parametersPolygonalFaceset.points[7].y = 0.0000;
    parametersPolygonalFaceset.points[7].z = 0.0000;

    parametersPolygonalFaceset.indicesPerFace = 4;
    parametersPolygonalFaceset.numIndices = 6 * parametersPolygonalFaceset.indicesPerFace;
    parametersPolygonalFaceset.indices = new uint32_t[parametersPolygonalFaceset.numIndices];
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[0] = 1;
    parametersPolygonalFaceset.indices[1] = 2;
    parametersPolygonalFaceset.indices[2] = 3;
    parametersPolygonalFaceset.indices[3] = 4;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[4] = 5;
    parametersPolygonalFaceset.indices[5] = 2;
    parametersPolygonalFaceset.indices[6] = 1;
    parametersPolygonalFaceset.indices[7] = 6;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[8] = 2;
    parametersPolygonalFaceset.indices[9] = 5;
    parametersPolygonalFaceset.indices[10] = 7;
    parametersPolygonalFaceset.indices[11] = 3;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[12] = 4;
    parametersPolygonalFaceset.indices[13] = 3;
    parametersPolygonalFaceset.indices[14] = 7;
    parametersPolygonalFaceset.indices[15] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[16] = 6;
    parametersPolygonalFaceset.indices[17] = 1;
    parametersPolygonalFaceset.indices[18] = 4;
    parametersPolygonalFaceset.indices[19] = 8;
    //IFCINDEXEDPOLYGONALFACE
    parametersPolygonalFaceset.indices[20] = 5;
    parametersPolygonalFaceset.indices[21] = 6;
    parametersPolygonalFaceset.indices[22] = 8;
    parametersPolygonalFaceset.indices[23] = 7;

    webifc::Geometry geometry7 = conwayGeometryProcessor.getPolygonalFaceSetGeometry(parametersPolygonalFaceset);

    geometryVec.push_back(geometry7);

    //free memory 
    delete parametersPolygonalFaceset.points;
    delete parametersPolygonalFaceset.indices;

    auto time = webifc::ms() - start;

    std::cout << "Processing geometry t " << time << "ms" << std::endl;

    std::cout << "Testing individual obj export..." << std::endl;

    for ( int geometryIndex = 0; geometryIndex < geometryVec.size(); geometryIndex++ )
    {
        size_t offset = 0;
        std::string singleObj = conwayGeometryProcessor.GeometryToObj(geometryVec[geometryIndex], offset, NormalizeMat);

        //filthy but just testing quick
        std::string fileName = "./";
        fileName += std::to_string(geometryIndex);
        fileName += "_conway.obj";


        std::wstring wsTmp(fileName.begin(), fileName.end());

        webifc::writeFile(wsTmp, singleObj);
    }

    std::cout << "Testing complete obj export..." << std::endl;

    std::string completeObj = "";
    size_t offset = 0;

    for ( int geometryIndex = 0; geometryIndex < geometryVec.size(); geometryIndex++ )
    {
        completeObj += conwayGeometryProcessor.GeometryToObj(geometryVec[geometryIndex], offset, NormalizeMat);
    }

    //filthy but just testing quick
    std::string fileName = "./index_ifc_full_conway.obj";

    std::wstring wsTmp(fileName.begin(), fileName.end());

    webifc::writeFile(wsTmp, completeObj);

    geometryVec.clear();

    std::cout << "Done" << std::endl;
}
