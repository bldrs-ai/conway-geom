
solution "conway_geom"
    configurations { "Debug", "Release"} --"WasmDebug", "WasmRelease" }
    includedirs    { "include" }
    location       ( _ACTION )

    platforms { "x64", "Emscripten" }


    configuration { "vs*", "Debug" }
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

        --targetextension ".js"

        WebIfcCoreFiles       = { "include/*.h" }
        WebIfcMathFiles       = { "include/math/*.h" }
        WebIfcParsingFiles    = { "include/parsing/*.h" }
        WebIfcSourceFiles     = { "web-ifc-api.cpp" }
        WebIfcTestSourceFiles = { "test/*.cpp" }
        WebIfcTestingMain     = { "web-ifc-test.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}

        configuration { "linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { WebIfcCoreFiles,
                    WebIfcMathFiles,
                    WebIfcParsingFiles,
                    --WebIfcSourceFiles,
                    ManifoldSrcFiles,
                    WebIfcTestSourceFiles,
                     WebIfcTestingMain
                    }

        configuration {"gmake and not macosx"}
        linkoptions { "--bind", "-03", "-flto", "--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"] -lworkerfs.js" }
        configuration {}
        libdirs {  }
        links {  }
        flags { "Symbols", "FullSymbols" }

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
            "external/manifold/src/third_party/glm"
            --"/Users/soar/Documents/GitHub/emsdk/upstream/emscripten/system/include"
            
            --"$(EMSDK)/upstream/emscripten/system/include"
        }    

        excludes { "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*" }     

        configuration { "Debug" }

        configuration { "Release", "gmake" }
                        
        configuration "Release*"
            flags { "OptimizeSpeed", "NoIncrementalLink" }
                             
        configuration { "x64", "Debug" }
            targetdir ( path.join( "bin", "64", "debug" ) )            
            flags { "EnableAVX2" }

        configuration { "x64", "Release" }
            targetdir ( path.join( "bin", "64", "release" ) )
            flags { "EnableAVX2" }


    project "conway_geom_wasm"
    language "C++"
        kind "ConsoleApp"
        files {}

        targetextension ".js"

        WebIfcCoreFiles       = { "include/*.h" }
        WebIfcMathFiles       = { "include/math/*.h" }
        WebIfcParsingFiles    = { "include/parsing/*.h" }
        WebIfcSourceFiles     = { "web-ifc-api.cpp" }
        WebIfcTestSourceFiles = { "test/*.cpp" }
        WebIfcTestingMain     = { "web-ifc-test.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}

        configuration { "linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { WebIfcCoreFiles,
                    WebIfcMathFiles,
                    WebIfcParsingFiles,
                    WebIfcSourceFiles,
                    ManifoldSrcFiles
                    --WebIfcTestSourceFiles,
                    -- WebIfcTestingMain
                    }

        configuration {"gmake"}
        linkoptions { "--bind", "-03", "-flto", "--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"] -lworkerfs.js" }
        configuration {}
        libdirs {  }
        links {  }
        flags { "Symbols", "FullSymbols" }

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
            "external/manifold/src/third_party/glm"
            --"/Users/soar/Documents/GitHub/emsdk/upstream/emscripten/system/include"
            
            --"$(EMSDK)/upstream/emscripten/system/include"
        }    

        excludes { "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*" }     

        configuration { "Debug" }

        configuration { "Release", "gmake" }
                        
        configuration "Release*"
            flags { "OptimizeSpeed", "NoIncrementalLink" }
                             
        configuration { "x64", "Debug" }
            targetdir ( path.join( "bin", "64", "debug" ) )            
            flags { "EnableAVX2" }

        configuration { "x64", "Release" }
            targetdir ( path.join( "bin", "64", "release" ) )
            flags { "EnableAVX2" }


    project "conway_geom_wasm_mt"
    language "C++"
        kind "ConsoleApp"
        files {}

        targetextension ".js"


        WebIfcCoreFiles       = { "include/*.h" }
        WebIfcMathFiles       = { "include/math/*.h" }
        WebIfcParsingFiles    = { "include/parsing/*.h" }
        WebIfcSourceFiles     = { "web-ifc-api.cpp" }
        WebIfcTestSourceFiles = { "test/*.cpp" }
        WebIfcTestingMain     = { "web-ifc-test.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}

        configuration { "linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-pthread", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { WebIfcCoreFiles,
                    WebIfcMathFiles,
                    WebIfcParsingFiles,
                    WebIfcSourceFiles,
                    ManifoldSrcFiles
                    --WebIfcTestSourceFiles,
                    --WebIfcTestingMain
                    }

        configuration { "gmake" }
        linkoptions { "-pthread", "-s PTHREAD_POOL_SIZE=navigator.hardwareConcurrency", "--bind", "-03", "-flto", "--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"] -lworkerfs.js" }
        configuration {}
        libdirs {  }
        links {  }
        flags { "Symbols", "FullSymbols" }

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
            "external/manifold/src/third_party/glm"
            --"/Users/soar/Documents/GitHub/emsdk/upstream/emscripten/system/include"
            
            --"$(EMSDK)/upstream/emscripten/system/include"
        } 

        excludes { "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*" }       

        configuration { "Debug" }

        configuration { "Release", "gmake" }
                        
        configuration "Release*"
            flags { "OptimizeSpeed", "NoIncrementalLink" }
                             
        configuration { "x64", "Debug" }
            targetdir ( path.join( "bin", "64", "debug" ) )            
            flags { "EnableAVX2" }

        configuration { "x64", "Release" }
            targetdir ( path.join( "bin", "64", "release" ) )
            flags { "EnableAVX2" }

    project "genie"
        kind "StaticLib"
        language "C"
        files { "**.lua" }

        configuration { "macos or ios"}
            postbuildcommands { "$(SolutionDir)../macos_genie/./genie gmake" }

        configuration { "windows" }
            postbuildcommands { "$(SolutionDir)../windows_genie/genie.exe gmake" }

        configuration { "Debug"}
            removeplatforms { "x64", "Native", "Universal64", "ARM64" }
        configuration { "Release"}
            removeplatforms { "x64", "Native", "Universal64", "ARM64" }
