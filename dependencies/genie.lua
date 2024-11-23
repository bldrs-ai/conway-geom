solution "dependencies"
configurations {"Debug", "Release"}
includedirs {"include"}
location(_ACTION)

platforms {"x64", "Emscripten"}

configuration {"vs*", "Debug"}
buildoptions {"/bigobj"}

configuration {"Debug*"}
targetdir "bin/debug"

configuration {"Release*"}
targetdir "bin/release"

configuration {}
project "draco"
language "C++"
kind "StaticLib"solution "dependencies"
configurations {"Debug", "Release"}
includedirs {"include"}
location(_ACTION)

platforms {"x64", "Emscripten"}

configuration {"vs*", "Debug"}
buildoptions {"/bigobj"}

configuration {"Debug*"}
targetdir "bin/debug"

configuration {"Release*"}
targetdir "bin/release"

configuration {}
project "draco"
language "C++"
kind "StaticLib"
files {}

DracoSourceFiles = {
    "../external/draco/src/draco/animation/*.cc",
    "../external/draco/src/draco/attributes/*.cc",
    "../external/draco/src/draco/compression/*.cc",
    "../external/draco/src/draco/compression/attributes/*.cc",
    "../external/draco/src/draco/compression/attributes/prediction_schemes/*.cc",
    "../external/draco/src/draco/compression/bit_coders/*.cc",
    "../external/draco/src/draco/compression/config/*.cc",
    "../external/draco/src/draco/compression/entropy/*.cc",
    "../external/draco/src/draco/compression/mesh/*.cc",
    "../external/draco/src/draco/compression/mesh/traverser/*.cc",
    "../external/draco/src/draco/compression/point_cloud/*.cc",
    "../external/draco/src/draco/compression/point_cloud/algorithms/*.cc",
    "../external/draco/src/draco/core/*.cc",
    "../external/draco/src/draco/io/*.cc",
    "../external/draco/src/draco/material/*.cc",
    "../external/draco/src/draco/mesh/*.cc",
    "../external/draco/src/draco/meshdata/*.cc",
    "../external/draco/src/draco/metadata/*.cc",
    "../external/draco/src/draco/point_cloud/*.cc",
    "../external/draco/src/draco/scene/*.cc",
    "../external/draco/src/draco/texture/*.cc"
}

configuration {"windows or linux or macosx or ios or gmake"}
buildoptions_cpp {
    "-O3",
    "-DNDEBUG",
    "-Wall",
    "-fexceptions",
    "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
    "-std=c++17",
    "-pthread"
}

configuration {"windows or macosx or linux"}
files {
    DracoSourceFiles
}

configuration {"gmake and not macosx and not windows"}
linkoptions {
    "--bind",
    "-03",
    "-flto",
    '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=Draco -s MODULARIZE=1 -s SIDE_MODULE=1'
}
configuration {}
libdirs {}
links {}
flags {"Symbols", "FullSymbols"}

includedirs {
    "../external/draco/src",
}

excludes {
    --Draco Source Files
    "../external/draco/src/draco/javascript/**.*",
    "../external/draco/src/draco/maya/**.*",
    "../external/draco/src/draco/tools/**.*",
    "../external/draco/src/draco/unity/**.*",
    "../external/draco/src/draco/animation/**.*",
    "../external/draco/src/draco/io/**.*",
    --Draco Test Files
    "../external/draco/src/draco/**/*test*cc"
}

configuration {"Debug"}

configuration {"Release", "gmake"}

configuration {"Emscripten", "Release"}
postbuildcommands {"cp ../bin/release/draco.bc ../wasm/libdraco.a"}

configuration {"macosx", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"cp ../bin/64/release/libdraco.a ../macOS-arm64/libdraco.a"}
flags {"EnableAVX2"}

configuration {"macosx", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"cp ../bin/64/release/libdraco.a ../macOS-arm64/libdraco.a"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libdraco.a ..\\win\\libdraco.a* /Y"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libdraco.a ..\\win\\libdraco.a* /Y"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libdraco.a ..\\win\\libdraco.a* /Y"}
flags {"EnableAVX2"}

project "manifold"
language "C++"
kind "StaticLib"
files {}

ManifoldSrcFiles = {
    "../external/manifold/src/**.*",
    "../external/manifold/src/collider/include/*.h",
    "../external/manifold/src/utilities/include/*.h"
}

configuration {"windows or linux or macosx or ios or gmake"}
buildoptions_cpp {
    "-O3",
    "-DNDEBUG",
    "-Wall",
    "-fexceptions",
    "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
    "-std=c++17",
    "-pthread"
}

configuration {"windows or macosx or linux"}
files {
    ManifoldSrcFiles
}

configuration {"gmake and not macosx and not windows"}
linkoptions {
    "--bind",
    "-03",
    "-flto",
    '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=Manifold -s MODULARIZE=1 -s SIDE_MODULE=1'
}
configuration {}
libdirs {}
links {}
flags {"Symbols", "FullSymbols"}

includedirs {
    "../external/manifold/src",
    "../external/manifold/src/utilities/include",
    "../external/manifold/src/collider/include",
    "../external/manifold/src/utilities/include",
    "../external/manifold/src/third_party/thrust",
    "../external/manifold/src/manifold/include",
    "../external/manifold/src/polygon/include",
    "../external/manifold/src/sdf/include",
    "../external/manifold/src/third_party/graphlite/include",
    "../external/manifold/src/third_party/glm",
}

excludes {
    --Manifold Test files
    "../external/manifold/src/third_party/glm/test/**.*",
    "../external/manifold/src/third_party/thrust/examples/**.*",
    "../external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
    "../external/manifold/src/third_party/glm/test/gtc/**.*",
}

configuration {"Debug"}

configuration {"Release", "gmake"}

configuration {"Emscripten", "Release"}
postbuildcommands {"cp ../bin/release/manifold.bc ../wasm/libmanifold.a"}

configuration "Release*"
flags {"OptimizeSpeed", "NoIncrementalLink"}

configuration {"macosx", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"cp ../bin/64/release/libmanifold.a ../macOS-arm64/libmanifold.a"}
flags {"EnableAVX2"}

configuration {"macosx", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"cp ../bin/64/release/libmanifold.a ../macOS-arm64/libmanifold.a"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libmanifold.a ..\\win\\libmanifold.a* /Y"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libmanifold.a ..\\win\\libmanifold.a* /Y"}
flags {"EnableAVX2"}

project "gltfsdk"
language "C++"
kind "StaticLib"
files {}

glTFSDKSrcFiles = {"../external/gltf-sdk/GLTFSDK/source/**.*"}

configuration {"linux or macosx or ios or gmake"}
buildoptions_cpp {
    "-O3",
    "-DNDEBUG",
    "-Wall",
    "-fexceptions",
    "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
    "-std=c++17",
    "-pthread"
}

configuration {"windows or macosx or linux"}
files {
    glTFSDKSrcFiles
}

configuration {"gmake", "Emscripten"}
linkoptions {
    "-O3",
    "--bind",
    "--dts",
    "-03",
    "-flto",
    '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=gltfsdk -s ENVIRONMENT=web -s SINGLE_FILE=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=["FS, WORKERFS"] -lworkerfs.js'
}

configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

configuration {}
libdirs {}
links {}
flags {"Symbols", "FullSymbols","UseObjectResponseFile"}

includedirs {
    "../external/gltf-sdk/GLTFSDK/Inc",
    "../external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
}

excludes {
    --glTF-SDK Source Files
    "../external/gltf-sdk/GLTFSDK/source/Version.cpp",
}

configuration {"Debug"}

configuration {"Release", "gmake"}

configuration {"Emscripten", "Release"}
postbuildcommands {"cp ../bin/release/gltfsdk.bc ../wasm/libgltfsdk.a"}

configuration "Release*"
flags {"OptimizeSpeed", "NoIncrementalLink"}

configuration {"macosx", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"cp ../bin/64/release/libgltfsdk.a ../macOS-arm64/libgltfsdk.a"}
flags {"EnableAVX2"}

configuration {"macosx", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"cp ../bin/64/release/libgltfsdk.a ../macOS-arm64/libgltfsdk.a"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Debug"}
targetdir(path.join("bin", "64", "debug"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libgltfsdk.a ..\\win\\libgltfsdk.a* /Y"}
flags {"EnableAVX2"}

configuration {"windows", "x64", "Release"}
targetdir(path.join("bin", "64", "release"))
postbuildcommands {"xcopy ..\\bin\\64\\release\\libgltfsdk.a ..\\win\\libgltfsdk.a* /Y"}
flags {"EnableAVX2"}