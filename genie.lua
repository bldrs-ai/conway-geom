solution "conway_geom"
    configurations {
        "Debug",
        "Release"
    }

    includedirs {"include"}
    location(_ACTION)

    platforms {
        "x64",
        "Emscripten"
    }

    configuration {"vs*", "Debug"}
        buildoptions { "/bigobj" }


    configuration {"Debug*"}
        targetdir "bin/debug"


    configuration {"Release*"}
        targetdir "bin/release"

    configuration {}

project "conway_geom_native"
    language "C++"
    kind "ConsoleApp"
    files {}

    ConwayCoreFiles = {
        "conway_geometry/*.h",
        "conway_geometry/*.cpp",
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geomeetry/legacy/**.*"
    }
    WebIfcSourceFiles = {"web-ifc-api.cpp"}
    WebIfcTestSourceFiles = {"test/*.cpp"}
    WebIfcTestingMain = {"web-ifc-test.cpp"}
    ConwayNativeMain = {"conway-native.cpp"}

    configuration {"windows or linux or macosx or ios or gmake"}
        buildoptions_cpp {
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-fexceptions",
            "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
            "-std=c++20"
        }

    configuration {"windows or macosx or linux"}
        files {
            ConwayCoreFiles,
            ConwayNativeMain
        }

    configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

    configuration {"gmake and not macosx and not windows"}
        linkoptions {
            "--bind",
            "-03",
            "-flto",
            '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=["FS, WORKERFS"] -lworkerfs.js'
        }

    configuration {}
        libdirs {}
        links {}
        flags {
            "Symbols",
            "FullSymbols",
            "UseObjectResponseFile"
        }

        includedirs {
            "external/tinynurbs/include",
            "external/manifold/src",
            "external/manifold/src/utilities/include",
            "external/glm",
            "external/earcut.hpp/include",
            "external/TinyCppTest/Sources",
            "external/manifold/src/collider/include",
            "external/manifold/src/utilities/include",
            "external/manifold/src/third_party/thrust",
            "external/manifold/src/manifold/include",
            "external/manifold/src/polygon/include",
            "external/manifold/src/sdf/include",
            "external/manifold/src/third_party/graphlite/include",
            "external/manifold/src/third_party/glm",
            "external/gltf-sdk/GLTFSDK/Inc",
            "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
            "external/draco/src",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt",
            "external/csgjs-cpp"
        }

        excludes {
            -- Manifold Test files
            "external/**/**cc",
            "external/manifold/src/third_party/glm/test/**.*",
            "external/manifold/src/third_party/thrust/examples/**.*",
            "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
            "external/manifold/src/third_party/glm/test/gtc/**.*",
            -- Draco Source Files
            "external/draco/**/*cc",
            -- glTF-SDK Source Files
            "external/gltf-sdk/**/**cpp",
            "external/fuzzy-bools/fuzzy/main.cpp"
        }

    configuration {"Debug"}

    configuration {"Release", "gmake"}

    configuration "Release*"
        flags {
            "OptimizeSpeed",
            "NoIncrementalLink"
        }

    configuration {"Emscripten", "Debug"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"Emscripten", "Release"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

project "conway_geom_native_tests"
    language "C++"
    kind "ConsoleApp"
    files {}

    ConwayCoreFiles = {
        "conway_geometry/*.h",
        "conway_geometry/*.cpp",
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geomeetry/legacy/**.*"
    }
    ConwayTestSourceFiles = {"test/*.cpp"}

    configuration {"windows or linux or macosx or ios or gmake"}
        buildoptions_cpp {
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-fexceptions",
            "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
            "-std=c++20"
        }

    configuration {"windows or macosx or linux"}
        files {
            ConwayCoreFiles,
            ConwayTestSourceFiles
        }

    configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

    configuration {"gmake and not macosx and not windows"}
        linkoptions {
            "--bind",
            "-03",
            "-flto",
            '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -sSTACK_SIZE=5MB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native_tests -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=["FS, WORKERFS"] -lworkerfs.js'
        }

    configuration {}
        libdirs {}
        links {}
        flags {
            "Symbols",
            "FullSymbols",
            "UseObjectResponseFile"
        }

        includedirs {
            "external/tinynurbs/include",
            "external/manifold/src",
            "external/manifold/src/utilities/include",
            "external/glm",
            "external/earcut.hpp/include",
            "external/TinyCppTest/Sources",
            "external/manifold/src/collider/include",
            "external/manifold/src/utilities/include",
            "external/manifold/src/third_party/thrust",
            "external/manifold/src/manifold/include",
            "external/manifold/src/polygon/include",
            "external/manifold/src/sdf/include",
            "external/manifold/src/third_party/graphlite/include",
            "external/manifold/src/third_party/glm",
            "external/gltf-sdk/GLTFSDK/Inc",
            "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
            "external/draco/src",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt",
            "external/csgjs-cpp"
        }

        excludes {
            -- Manifold Test files
            "external/manifold/src/third_party/glm/test/**.*",
            "external/manifold/src/third_party/thrust/examples/**.*",
            "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
            "external/manifold/src/third_party/glm/test/gtc/**.*",
            -- Draco Source Files
            "external/draco/src/draco/javascript/**.*",
            "external/draco/src/draco/maya/**.*",
            "external/draco/src/draco/tools/**.*",
            "external/draco/src/draco/unity/**.*",
            "external/draco/src/draco/animation/**.*",
            "external/draco/src/draco/io/**.*",
            -- Draco Test Files
            "external/draco/src/draco/**/*test*cc",
            -- glTF-SDK Source Files
            "external/gltf-sdk/GLTFSDK/source/Version.cpp",
            "external/fuzzy-bools/fuzzy/main.cpp"
        }

    configuration {"Debug"}

    configuration {"Release", "gmake"}

    configuration "Release*"
        flags {"OptimizeSpeed", "NoIncrementalLink"}

    configuration {"Emscripten", "Debug"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"Emscripten", "Release"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

project "webifc_native"
    language "C++"
    kind "ConsoleApp"
    files {}

    WebIfcCoreFiles = {
        "geometry/**.*",
        "parsing/**.*",
        "utility/**.*",
        "schema/**.*",
        "test/io_helpers.cpp",
        "test/io_helpers.h"
    }
    WebIfcTestingMain = {"web-ifc-test.cpp"}

    configuration {"windows or linux or macosx or ios or gmake"}
        buildoptions_cpp {
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-fexceptions",
            "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
            "-std=c++20"
        }

    configuration {"windows or macosx or linux"}
        files {WebIfcCoreFiles, WebIfcTestingMain}

    configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

    configuration {"gmake and not macosx and not windows"}
        linkoptions {
            "--bind",
            "-03",
            "-flto",
            '--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=webifc_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=["FS, WORKERFS"] -lworkerfs.js'
        }

    configuration {}
        libdirs {}
        links {}
        flags {
            "Symbols",
            "FullSymbols",
            "UseObjectResponseFile"
        }

        includedirs {
            "external/tinynurbs/include",
            "external/manifold/src",
            "external/manifold/src/utilities/include",
            "external/glm",
            "external/earcut.hpp/include",
            "external/TinyCppTest/Sources",
            "external/manifold/src/collider/include",
            "external/manifold/src/utilities/include",
            "external/manifold/src/third_party/thrust",
            "external/manifold/src/manifold/include",
            "external/manifold/src/polygon/include",
            "external/manifold/src/sdf/include",
            "external/manifold/src/third_party/graphlite/include",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt",
            "external/fast_float/include"
        }

        excludes {
            "external/manifold/src/third_party/glm/test/**.*",
            "external/manifold/src/third_party/thrust/examples/**.*",
            "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
            "external/manifold/src/third_party/glm/test/gtc/**.*",
            "external/fuzzy-bools/fuzzy/main.cpp"
        }

    configuration {"Debug"}

    configuration {"Release", "gmake"}

    configuration {"Emscripten", "Debug"}
        libdirs {"./dependencies/wasm"}
        links {"manifold"}

    configuration {"Emscripten", "Release"}
        libdirs {"./dependencies/wasm"}
        links {"manifold"}

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/macOS-arm64"}
        links {"manifold"}
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/macOS-arm64"}
        links {"manifold"}
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {"manifold"}
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/win"}
        links {"manifold"}
        flags {"EnableAVX2"}

project "ConwayGeomWasmNode"
    language "C++"
    kind "ConsoleApp"
    files {}

    targetname "ConwayGeomWasm"

    targetextension ".js"

    ConwayCoreFiles = {
        "conway_geometry/*.h",
        "conway_geometry/*.cpp",
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geomeetry/legacy/**.*"
    }
    ConwaySourceFiles = {"conway-api.cpp"}

    configuration {"linux or macosx or ios or gmake"}
        buildoptions_cpp {
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-fexceptions",
            "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
            "-std=c++20"
        }

    configuration {"windows or macosx or linux"}
        files {
            ConwayCoreFiles,
            ConwaySourceFiles
        }

    configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

if _ARGS[1] == "profile" and _ARGS[2] ~= nil then
    configuration {"gmake"}
    linkoptions {
        "-g -O0",
        "-gdwarf-5",
        "-gpubnames",
        "--bind",
        "--dts",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s ENVIRONMENT=node",
        "-s ALLOW_MEMORY_GROWTH=1",
        "-s MAXIMUM_MEMORY=4GB",
        "-s STACK_SIZE=5MB",
        "-s FORCE_FILESYSTEM=1",
        "-gsource-map",
        "--source-map-base " .. _ARGS[2],
        "-sASSERTIONS",
        "-s SAFE_HEAP=1",
        "-s EXPORT_NAME=ConwayGeomWasm",
        --"-s USE_ES6_IMPORT_META=0",
       "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
        "-s EXPORT_ES6=1",
        "-s MODULARIZE=1",
        "-sNO_DISABLE_EXCEPTION_CATCHING"
    }
else 
    configuration {"gmake"}
    linkoptions {
        "-O3",
        "--bind",
        "--dts",
        "-03",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s ALLOW_MEMORY_GROWTH=1",
        "-s MAXIMUM_MEMORY=4GB",
        "-s STACK_SIZE=5MB",
        "-s FORCE_FILESYSTEM=1",
        "-s EXPORT_NAME=ConwayGeomWasm",
        "-s ENVIRONMENT=node",
        "-s SINGLE_FILE=1",
        --"-s USE_ES6_IMPORT_META=0",
        "-s EXPORT_ES6=1",
        "-s MODULARIZE=1",
        "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
        "-lworkerfs.js"
    }
end

    configuration {}
        libdirs {}
        links {}
        flags {
            "Symbols",
            "FullSymbols",
            "UseObjectResponseFile"
        }

        includedirs {
            "external/tinynurbs/include",
            "external/manifold/src",
            "external/manifold/src/utilities/include",
            "external/glm",
            "external/earcut.hpp/include",
            "external/TinyCppTest/Sources",
            "external/manifold/src/collider/include",
            "external/manifold/src/utilities/include",
            "external/manifold/src/third_party/thrust",
            "external/manifold/src/manifold/include",
            "external/manifold/src/polygon/include",
            "external/manifold/src/sdf/include",
            "external/manifold/src/third_party/graphlite/include",
            "external/manifold/src/third_party/glm",
            "external/gltf-sdk/GLTFSDK/Inc",
            "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
            "external/draco/src",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt",
            "external/csgjs-cpp"
        }

        excludes {
            -- Manifold Test files
            "external/**/**cc",
            "external/manifold/src/third_party/glm/test/**.*",
            "external/manifold/src/third_party/thrust/examples/**.*",
            "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
            "external/manifold/src/third_party/glm/test/gtc/**.*",
            -- Draco Source Files
            "external/draco/**/*cc",
            -- glTF-SDK Source Files
            "external/gltf-sdk/**/**cpp",
            "external/fuzzy-bools/fuzzy/main.cpp"
        }

    configuration {"Debug"}

    configuration "Release*"
        flags {
            "OptimizeSpeed",
            "NoIncrementalLink"
        }

    configuration {"Emscripten", "Debug"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"Emscripten", "Release"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

project "ConwayGeomWasmWeb"
    language "C++"
    kind "ConsoleApp"
    files {}

    targetname "ConwayGeomWasm"

    targetextension ".js"

    ConwayCoreFiles = {
        "conway_geometry/*.h",
        "conway_geometry/*.cpp",
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geomeetry/legacy/**.*"
    }
    ConwaySourceFiles = {"conway-api.cpp"}

    configuration {"linux or macosx or ios or gmake"}
        buildoptions_cpp {
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-fexceptions",
            "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
            "-std=c++20"
        }

    configuration {"windows or macosx or linux"}
        files {
            ConwayCoreFiles,
            ConwaySourceFiles
        }

    configuration {"windows"}
        prelinkcommands {
            "$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
            "$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))",
            "$(foreach string,$(OBJECTS),\
                $(file >> $(NEWOBJRESP),$(string) )\
                )"
        }

if _ARGS[1] == "profile" and _ARGS[2] ~= nil then
    configuration {"gmake"}
    linkoptions {
        "-g -O0",
        "-gdwarf-5",
        "-gpubnames",
        "--bind",
        "--dts",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s ENVIRONMENT=web",
        "-s ALLOW_MEMORY_GROWTH=1",
        "-s MAXIMUM_MEMORY=4GB",
        "-s STACK_SIZE=5MB",
        "-s FORCE_FILESYSTEM=1",
        "-gsource-map",
        "--source-map-base " .. _ARGS[2],
        "-sASSERTIONS",
        "-s SAFE_HEAP=1",
        "-s EXPORT_NAME=ConwayGeomWasm",
        "-s USE_ES6_IMPORT_META=0",
        "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
        "-s EXPORT_ES6=1",
        "-s MODULARIZE=1",
        "-sNO_DISABLE_EXCEPTION_CATCHING"
    }
else 
    configuration {"gmake"}
    linkoptions {
        "-O3",
        "--bind",
        "--dts",
        "-03",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s ALLOW_MEMORY_GROWTH=1",
        "-s MAXIMUM_MEMORY=4GB",
        "-s STACK_SIZE=5MB",
        "-s FORCE_FILESYSTEM=1",
        "-s EXPORT_NAME=ConwayGeomWasm",
        "-s ENVIRONMENT=web",
        "-s SINGLE_FILE=1",
        "-s USE_ES6_IMPORT_META=0",
        "-s EXPORT_ES6=1",
        "-s MODULARIZE=1",
        "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
        "-lworkerfs.js"
    }
end

    configuration {}
        libdirs {}
        links {}
        flags {
            "Symbols",
            "FullSymbols",
            "UseObjectResponseFile"
        }

        includedirs {
            "external/tinynurbs/include",
            "external/manifold/src",
            "external/manifold/src/utilities/include",
            "external/glm",
            "external/earcut.hpp/include",
            "external/TinyCppTest/Sources",
            "external/manifold/src/collider/include",
            "external/manifold/src/utilities/include",
            "external/manifold/src/third_party/thrust",
            "external/manifold/src/manifold/include",
            "external/manifold/src/polygon/include",
            "external/manifold/src/sdf/include",
            "external/manifold/src/third_party/graphlite/include",
            "external/manifold/src/third_party/glm",
            "external/gltf-sdk/GLTFSDK/Inc",
            "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
            "external/draco/src",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt",
            "external/csgjs-cpp"
        }

        excludes {
            -- Manifold Test files
            "external/**/**cc",
            "external/manifold/src/third_party/glm/test/**.*",
            "external/manifold/src/third_party/thrust/examples/**.*",
            "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
            "external/manifold/src/third_party/glm/test/gtc/**.*",
            -- Draco Source Files
            "external/draco/**/*cc",
            -- glTF-SDK Source Files
            "external/gltf-sdk/**/**cpp",
            "external/fuzzy-bools/fuzzy/main.cpp"
        }

    configuration {"Debug"}

    configuration "Release*"
        flags {
            "OptimizeSpeed",
            "NoIncrementalLink"
        }

    configuration {"Emscripten", "Debug"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"Emscripten", "Release"}
        libdirs {"./dependencies/wasm"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/macOS-arm64"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}
