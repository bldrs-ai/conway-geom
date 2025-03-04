local nodeCores = "Math.max(require('os').cpus().length - 1, 1)"
local webCores  = "Math.max(navigator.hardwareConcurrency - 1, 1)"

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


project "conway_csg_native"
    language "C++"
    kind "ConsoleApp"
    files {}
    
    flags { "NoPCH", "NoRTTI", "ExtraWarnings" }

    configuration { "vs* or windows" }
      defines	"_CRT_SECURE_NO_WARNINGS"
        
    configuration "Release*"
      flags { "OptimizeSpeed", "NoIncrementalLink" }

    configuration { "vs*" }
        buildoptions { "/std:c++latest" }
        
    configuration { "Release and vs*" }
        buildoptions { "/Zi" }

    ConwayCoreFiles = {
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geometry/structures/**.*",
        "conway_geometry/csg/**.*",
        "logging/**.*"
      --  "conway_geometry/legacy/**.*"
    }
    ConwayNativeMain = {"conway_geometry/utilities/csg_command_line.cpp"}

    configuration {"windows or macosx or linux"}
        files {
            ConwayCoreFiles,
            ConwayNativeMain
        }

    configuration {"gmake"}
        buildoptions_cpp {
          "-Wall",
          "-fexceptions",
          "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
          "-std=c++20",
          "-pthread"
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
          "conway_geometry",
          "./",
      --    "external/fuzzy-bools",
          "external/tinynurbs/include",
          "external/glm",
          "external/earcut.hpp/include",
          "external/TinyCppTest/Sources",
          "external/CDT/CDT/include",
          "external/tinyobjloader"
      }

      excludes {
          -- Manifold Test files
          "external/**/**cc",
      }

    configuration {"Debug"}

    configuration {"Release", "gmake"}

    configuration "Release*"
        flags {
            "OptimizeSpeed",
            "NoIncrementalLink"
        }

    configuration {"Emscripten", "Debug"}

    configuration {"Emscripten", "Release"}

    configuration {"macosx", "x64", "Debug"}
        targetdir(path.join("bin", "64", "debug"))
        flags {"EnableAVX2"}

    configuration {"macosx", "x64", "Release"}
        targetdir(path.join("bin", "64", "release"))
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Debug"}
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
        targetdir(path.join("bin", "64", "debug"))
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
        targetdir(path.join("bin", "64", "release"))
        flags {"EnableAVX2"}

project "conway_geom_native"
    language "C++"
    kind "ConsoleApp"
    files {}

    ConwayCoreFiles = {
        "conway_geometry/*.h",
        "conway_geometry/*.cpp",
        "conway_geometry/operations/**.*",
        "conway_geometry/representation/**.*",
        "conway_geometry/structures/**.*",
        "conway_geometry/csg/**.*",
        "logging/**.*"
      --  "conway_geometry/legacy/**.*"
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
            "conway_geometry",
            "logging",
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
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
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
        "logging/**.*"
      --  "conway_geometry/legacy/**.*"
    }
    ConwayTestSourceFiles = {"test/ConwayGeometryProcessor_test.cpp", "test/main.cpp"}

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
            "conway_geometry",
            "logging",
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
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
        targetdir(path.join("bin", "64", "debug"))
        libdirs {"./dependencies/win"}
        links {
            "draco",
            "manifold",
            "gltfsdk"
        }
        flags {"EnableAVX2"}

    configuration {"windows", "x64", "Release"}
        defines "_USE_MATH_DEFINES" --required by legacy boolean library csgjs
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
    
        targetname "ConwayGeomWasmNode"
    
        targetextension ".js"
    
        ConwayCoreFiles = {
            "conway_geometry/*.h",
            "conway_geometry/*.cpp",
            "conway_geometry/operations/**.*",
            "conway_geometry/representation/**.*",
            "conway_geometry/structures/**.*",
            "conway_geometry/csg/**.*",
            "logging/**.*"
          --  "conway_geometry/legacy/**.*"
        }
        ConwaySourceFiles = {"conway-api.cpp"}
    
        configuration {"linux or macosx or ios or gmake"}
            buildoptions_cpp {
                "-Wall",
                "-fexceptions",
                "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
                -- TODO(Conor): I want threads for performance reasons.
                -- "-pthread",
                "-std=c++20",
                "-fexperimental-library",
                -- TODO(pablo): https://github.com/bldrs-ai/conway/wiki/Performance#simd
                -- "-msimd128",
                -- "-DGLM_FORCE_INTRINSICS=1"
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
            "-s PRECISE_F32=1",
            "--define-macro=REAL_T_IS_DOUBLE",
            "-s ENVIRONMENT=node",
            "-s ALLOW_MEMORY_GROWTH=1",
            "-s MAXIMUM_MEMORY=4GB",
            "-s STACK_SIZE=10MB",
            "-s FORCE_FILESYSTEM=1",
            "-gsource-map",
            "--source-map-base " .. _ARGS[2],
            "-s NODERAWFS=1",
            "-s EXPORT_NAME=ConwayGeomWasm",
            "-s ABORTING_MALLOC=0",
            --"-s USE_ES6_IMPORT_META=0",
            "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
            "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
            "-s EXPORT_ES6=1",
            "-s MODULARIZE=1",
            "-sNO_DISABLE_EXCEPTION_CATCHING",
            "-lworkerfs.js",
      --      "-pthread"
        }
    else 
        configuration {"gmake"}
        buildoptions_cpp {
          "-O3",
          "-DNDEBUG"
        }
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
            "-s PRECISE_F32=1",
            "-s NODERAWFS=1",
            "-s EXPORT_NAME=ConwayGeomWasm",
            "-s ENVIRONMENT=node",
            "-s SINGLE_FILE=1",
            --"-s USE_ES6_IMPORT_META=0",
            "-s EXPORT_ES6=1",
            "-s MODULARIZE=1",
            "-s ABORTING_MALLOC=0",
            "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
            "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
            "-lworkerfs.js",
            "-sNO_DISABLE_EXCEPTION_CATCHING",
        --    "-pthread"
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
                ".",
                "utility",
                "conway_geometry",
                "logging",
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
           --     "external/fuzzy-bools",
          --      "external/csgjs-cpp",
                "external/CDT/CDT/include",
                "external/fast_float/include"--,
               -- "external/tinyobjloader"
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

project "ConwayGeomWasmNodeMT"
language "C++"
kind "ConsoleApp"
files {}

targetname "ConwayGeomWasmNodeMT"

targetextension ".js"

ConwayCoreFiles = {
    "conway_geometry/*.h",
    "conway_geometry/*.cpp",
    "conway_geometry/operations/**.*",
    "conway_geometry/representation/**.*",
    "conway_geometry/structures/**.*",
    "conway_geometry/csg/**.*",
    "logging/**.*"
  --  "conway_geometry/legacy/**.*"
}
ConwaySourceFiles = {"conway-api.cpp"}

configuration {"linux or macosx or ios or gmake"}
    buildoptions_cpp {
        "-Wall",
        "-fexceptions",
        "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
        -- TODO(Conor): I want threads for performance reasons.

        "-s USE_PTHREADS=1",
        "-std=c++20",
        "-fexperimental-library",
        -- TODO(pablo): https://github.com/bldrs-ai/conway/wiki/Performance#simd
        -- "-msimd128",
        -- "-DGLM_FORCE_INTRINSICS=1"
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
    "-s PRECISE_F32=1",
    "--define-macro=REAL_T_IS_DOUBLE",
    "-s ENVIRONMENT=node",
    "-s ALLOW_MEMORY_GROWTH=1",
    "-s MAXIMUM_MEMORY=4GB",
    "-s STACK_SIZE=10MB",
    "-s FORCE_FILESYSTEM=1",
    "-gsource-map",
    "--source-map-base \"" .. _ARGS[2] .. "\"",
    "-s NODERAWFS=1",
    --"-sASSERTIONS",
--     "-s SAFE_HEAP=1",
    "-s EXPORT_NAME=ConwayGeomWasm",
    "-s ABORTING_MALLOC=0",
    --"-s USE_ES6_IMPORT_META=0",
    "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
    "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
    "-s EXPORT_ES6=1",
    "-s MODULARIZE=1",
    "-sNO_DISABLE_EXCEPTION_CATCHING",
    "-lworkerfs.js",
--      "-pthread"
}
else 
configuration {"gmake"}
buildoptions_cpp {
  "-O3",
  "-DNDEBUG",
  "-pthread"
}
linkoptions {
    "-O3",
    "--bind",
    "--dts",
    "-03",
    "-flto",
    "-pthread",
    "-sSHARED_MEMORY",
    "-s USE_PTHREADS=1",
    "-sPTHREAD_POOL_SIZE=\"" .. nodeCores .. "\"",
    "--define-macro=REAL_T_IS_DOUBLE",
    "-s ALLOW_MEMORY_GROWTH=1",
    "-s MAXIMUM_MEMORY=4GB",
    "-s STACK_SIZE=5MB",
    "-s FORCE_FILESYSTEM=1",
    "-s PRECISE_F32=1",
    "-s NODERAWFS=1",
    "-s EXPORT_NAME=ConwayGeomWasm",
    "-s ENVIRONMENT=node",
    "-s SINGLE_FILE=1",
    --"-s USE_ES6_IMPORT_META=0",
    "-s EXPORT_ES6=1",
    "-s MODULARIZE=1",
    "-s ABORTING_MALLOC=0",
    "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
    "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
    "-lworkerfs.js",
    "-sNO_DISABLE_EXCEPTION_CATCHING",
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
        ".",
        "utility",
        "conway_geometry",
        "logging",
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
    --     "external/fuzzy-bools",
  --      "external/csgjs-cpp",
        "external/CDT/CDT/include",
         "external/fast_float/include"--,
        -- "external/tinyobjloader"
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
      "draco_MT",
      "manifold_MT",
      "gltfsdk_MT"
    }

configuration {"Emscripten", "Release"}
    libdirs {"./dependencies/wasm"}
    links {
        "draco_MT",
        "manifold_MT",
        "gltfsdk_MT"
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

    targetname "ConwayGeomWasmWeb"

    targetextension ".js"

    ConwayCoreFiles = {
      "conway_geometry/*.h",
      "conway_geometry/*.cpp",
      "conway_geometry/operations/**.*",
      "conway_geometry/representation/**.*",
      "conway_geometry/structures/**.*",
      "conway_geometry/csg/**.*",
      "logging/**.*"
    --  "conway_geometry/legacy/**.*"
  }
  ConwaySourceFiles = {"conway-api.cpp"}

    configuration {"linux or macosx or ios or gmake"}
      buildoptions_cpp {
        "-Wall",
        "-fexceptions",
        "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
        -- TODO(Conor): I want threads for performance reasons.
        -- "-pthread",
        "-std=c++20",
        "-fexperimental-library",
        -- TODO(pablo): https://github.com/bldrs-ai/conway/wiki/Performance#simd
        -- "-msimd128",
        -- "-DGLM_FORCE_INTRINSICS=1"
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
    buildoptions_cpp {
      "-fsanitize=address"
    }
    linkoptions {
        "-g -O0",
        "-gdwarf-5",
        "-gpubnames",
        "--bind",
        "--dts",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s PRECISE_F32=1",
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
        "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
        "-s EXPORT_ES6=1",
        "-s MODULARIZE=1",
        "-sNO_DISABLE_EXCEPTION_CATCHING",
        "-s DISABLE_EXCEPTION_CATCHING=0",
    }
else 
    configuration {"gmake"}
    buildoptions_cpp {
      "-O3",
      "-DNDEBUG"
    }
    linkoptions {
        "-O3",
        "--bind",
        "--dts",
        "-03",
        "-flto",
        "--define-macro=REAL_T_IS_DOUBLE",
        "-s PRECISE_F32=1",
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
        "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
        "-lworkerfs.js",
        "-sNO_DISABLE_EXCEPTION_CATCHING",
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
            ".",
            "utility",
            "conway_geometry",
            "logging",
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
       --     "external/fuzzy-bools",
            "external/csgjs-cpp",
            "external/CDT/CDT/include",
            "external/fast_float/include"--,
           -- "external/tinyobjloader"
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

project "ConwayGeomWasmWebMT"
  language "C++"
  kind "ConsoleApp"
  files {}

targetname "ConwayGeomWasmWebMT"

targetextension ".js"

ConwayCoreFiles = {
  "conway_geometry/*.h",
  "conway_geometry/*.cpp",
  "conway_geometry/operations/**.*",
  "conway_geometry/representation/**.*",
  "conway_geometry/structures/**.*",
  "conway_geometry/csg/**.*",
  "logging/**.*"
--  "conway_geometry/legacy/**.*"
}
ConwaySourceFiles = {"conway-api.cpp"}

configuration {"linux or macosx or ios or gmake"}
    buildoptions_cpp {
        "-O3",
        "-Wall",
        "-fexceptions",
        "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP",
        "-std=c++20",
        "-fexperimental-library",
        "-matomics",
        "-mbulk-memory",
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
    "-pthread",
    "-s PRECISE_F32=1",
    "-sSHARED_MEMORY",
    "-s USE_PTHREADS=1",
    "--define-macro=REAL_T_IS_DOUBLE",
    "-sENVIRONMENT=web,worker",
    "-s ALLOW_MEMORY_GROWTH=1",
    "-s MAXIMUM_MEMORY=4GB",
    "-s STACK_SIZE=5MB",
    "-s FORCE_FILESYSTEM=1",
    "-gsource-map",
    "--source-map-base " .. _ARGS[2],
    --"-sASSERTIONS",
    --"-s SAFE_HEAP=1",
    "-s EXPORT_NAME=ConwayGeomWasm",
    "-s USE_ES6_IMPORT_META=0",
    "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
    "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
    "-s EXPORT_ES6=1",
    "-s MODULARIZE=1",
    "-sNO_DISABLE_EXCEPTION_CATCHING",
    "-s DISABLE_EXCEPTION_CATCHING=0",
}
else 
configuration {"gmake"}
buildoptions_cpp {
  "-O3",
  "-DNDEBUG",
  "-pthread"
}
linkoptions {
    "-O3",
    "--bind",
    "--dts",
    "-flto",
    "-pthread",
    "-s PRECISE_F32=1",
    "-sSHARED_MEMORY",
    "-s USE_PTHREADS=1",
    "-sPTHREAD_POOL_SIZE=\"" .. webCores .. "\"",
    "--define-macro=REAL_T_IS_DOUBLE",
    "-s ALLOW_MEMORY_GROWTH=1",
    "-s MAXIMUM_MEMORY=4GB",
    "-s STACK_SIZE=5MB",
    "-s FORCE_FILESYSTEM=1",
    "-s EXPORT_NAME=ConwayGeomWasm",
    "-s ENVIRONMENT=web,worker",
    "-s USE_ES6_IMPORT_META=0",
    "-s EXPORT_ES6=1",
    "-s MODULARIZE=1",
    "-s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"]",
    "-s EXPORTED_FUNCTIONS=[\"_malloc, _free\"]",
    "-lworkerfs.js",
    "-sNO_DISABLE_EXCEPTION_CATCHING",
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
        ".",
        "utility",
        "conway_geometry",
        "logging",
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
    --     "external/fuzzy-bools",
        "external/csgjs-cpp",
        "external/CDT/CDT/include",
        "external/fast_float/include"--,
        -- "external/tinyobjloader"
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
        "draco_MT",
        "manifold_MT",
        "gltfsdk_MT"
    }

configuration {"Emscripten", "Release"}
    libdirs {"./dependencies/wasm"}
    links {
        "draco_MT",
        "manifold_MT",
        "gltfsdk_MT"
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
