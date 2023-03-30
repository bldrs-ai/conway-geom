/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stack>
#include <memory>
#include <map>

#include <emscripten/bind.h>

#include "include/web-ifc.h"
#include "include/conway-geometry.h"

#include "version.h"

std::unique_ptr<conway::ConwayGeometryProcessor> conwayProcessor;// = std::make_unique<conway::ConwayGeometryProcessor>();

uint32_t GLOBAL_MODEL_ID_COUNTER = 0;

#ifdef __EMSCRIPTEN_PTHREADS__
    constexpr bool MT_ENABLED = true;
#else
    constexpr bool MT_ENABLED = false;
#endif

bool shown_version_header = false;

// use to construct API placeholders
int main() 
{
    conwayProcessor = std::make_unique<conway::ConwayGeometryProcessor>();

    return 0;
}

conway::IfcGeometry GetGeometry(conway::ConwayGeometryProcessor::ParamsPolygonalFaceSet parameters)
{
    return conwayProcessor->getPolygonalFaceSetGeometry(parameters);
}

EMSCRIPTEN_BINDINGS(my_module) {

    emscripten::class_<conway::IfcGeometry>("IfcGeometry")
        .constructor<>()
        .function("GetVertexData", &conway::IfcGeometry::GetVertexData)
        .function("GetVertexDataSize", &conway::IfcGeometry::GetVertexDataSize)
        .function("GetIndexData", &conway::IfcGeometry::GetIndexData)
        .function("GetIndexDataSize", &conway::IfcGeometry::GetIndexDataSize)
        ;


    emscripten::value_object<glm::dvec4>("dvec4")
        .field("x", &glm::dvec4::x)
        .field("y", &glm::dvec4::y)
        .field("z", &glm::dvec4::z)
        .field("w", &glm::dvec4::w)
        ;

    emscripten::value_array<std::array<double, 16>>("array_double_16")
            .element(emscripten::index<0>())
            .element(emscripten::index<1>())
            .element(emscripten::index<2>())
            .element(emscripten::index<3>())
            .element(emscripten::index<4>())
            .element(emscripten::index<5>())
            .element(emscripten::index<6>())
            .element(emscripten::index<7>())
            .element(emscripten::index<8>())
            .element(emscripten::index<9>())
            .element(emscripten::index<10>())
            .element(emscripten::index<11>())
            .element(emscripten::index<12>())
            .element(emscripten::index<13>())
            .element(emscripten::index<14>())
            .element(emscripten::index<15>())
            ;

    emscripten::register_vector<std::string>("stringVector");
    emscripten::register_vector<uint32_t>("UintVector");
    emscripten::function("GetGeometry", &GetGeometry);
}
