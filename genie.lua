
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

        ConwayCoreFiles       = { "conway_geometry/*.h", "conway_geometry/*.cpp", "conway_geometry/operations/**.*", "conway_geometry/representation/**.*"}
        WebIfcSourceFiles     = { "web-ifc-api.cpp" }
        WebIfcTestSourceFiles = { "test/*.cpp" }
        WebIfcTestingMain     = { "web-ifc-test.cpp" }
        ConwayNativeMain     = { "conway-native.cpp"}
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}
        glTFSDKSrcFiles       = {"external/gltf-sdk/GLTFSDK/source/**.*"}
        DracoSourceFiles      = {"external/draco/src/draco/animation/*.cc",
                                "external/draco/src/draco/attributes/*.cc",
                                "external/draco/src/draco/compression/*.cc",
                                "external/draco/src/draco/compression/attributes/*.cc",
                                "external/draco/src/draco/compression/attributes/prediction_schemes/*.cc",
                                "external/draco/src/draco/compression/bit_coders/*.cc",
                                "external/draco/src/draco/compression/config/*.cc",
                                "external/draco/src/draco/compression/entropy/*.cc",
                                "external/draco/src/draco/compression/mesh/*.cc",
                                "external/draco/src/draco/compression/mesh/traverser/*.cc",
                                "external/draco/src/draco/compression/point_cloud/*.cc",
                                "external/draco/src/draco/compression/point_cloud/algorithms/*.cc",
                                "external/draco/src/draco/core/*.cc",
                                "external/draco/src/draco/io/*.cc",
                                "external/draco/src/draco/material/*.cc",
                                "external/draco/src/draco/mesh/*.cc",
                                "external/draco/src/draco/meshdata/*.cc",
                                "external/draco/src/draco/metadata/*.cc",
                                "external/draco/src/draco/point_cloud/*.cc",
                                "external/draco/src/draco/scene/*.cc",
                                "external/draco/src/draco/texture/*.cc"}

        configuration { "windows or linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { ConwayCoreFiles,
                    ManifoldSrcFiles,
                    glTFSDKSrcFiles,
                    DracoSourceFiles,
                    ConwayNativeMain
                }

        configuration {"gmake and not macosx and not windows"}
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
            "external/manifold/src/third_party/glm",
            "external/gltf-sdk/GLTFSDK/Inc",
            "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
            "external/draco/src",
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt"
            --"/Users/soar/Documents/GitHub/emsdk/upstream/emscripten/system/include"
            
            --"$(EMSDK)/upstream/emscripten/system/include"
        }    

        excludes 
        { 
        --Manifold Test files 
        "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*",
        --Draco Source Files
        "external/draco/src/draco/javascript/**.*",
        "external/draco/src/draco/maya/**.*",
        "external/draco/src/draco/tools/**.*",
        "external/draco/src/draco/unity/**.*",
		"external/draco/src/draco/animation/**.*",
		"external/draco/src/draco/io/**.*",
        --Draco Test Files
        "external/draco/src/draco/animation/*test*cc",
        "external/draco/src/draco/attributes/*test*cc",
        "external/draco/src/draco/core/*test*cc",
        "external/draco/src/draco/io/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/mesh/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/point_cloud/*test*cc",
        "external/draco/src/draco/scene/*test*cc",
        "external/draco/src/draco/texture/*test*cc",
        "external/draco/src/draco/metadata/*test*cc",
        "external/draco/src/draco/compression/*test*cc",
        "external/draco/src/draco/compression/attributes/*test*cc",
        "external/draco/src/draco/compression/attributes/prediction_schemes/*test*cc",
        "external/draco/src/draco/compression/bit_coders/*test*cc",
        "external/draco/src/draco/compression/config/*test*cc",
        "external/draco/src/draco/compression/entropy/*test*cc",
        "external/draco/src/draco/compression/mesh/*test*cc",
        "external/draco/src/draco/compression/mesh/traverser/*test*cc",
        "external/draco/src/draco/compression/point_cloud/*test*cc",
        "external/draco/src/draco/compression/point_cloud/algorithms/*test*cc",
		--glTF-SDK Source Files
		"external/gltf-sdk/GLTFSDK/source/Version.cpp",
        "external/fuzzy-bools/fuzzy/main.cpp"}     

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

            project "conway_geom_native_tests"
            language "C++"
                kind "ConsoleApp"
                files {}
        
                ConwayCoreFiles       = { "conway_geometry/*.h", "conway_geometry/*.cpp", "conway_geometry/operations/**.*", "conway_geometry/representation/**.*"}
                ConwayTestSourceFiles = { "test/*.cpp" }
                --ConwayTestingMain     = { "conway-test.cpp"}
                ManifoldSrcFiles      = { "external/manifold/src/**.*", 
                "external/manifold/src/collider/include/*.h", 
                "external/manifold/src/utilities/include/*.h"}
                glTFSDKSrcFiles       = {"external/gltf-sdk/GLTFSDK/source/**.*"}
                DracoSourceFiles      = {"external/draco/src/draco/animation/*.cc",
                                        "external/draco/src/draco/attributes/*.cc",
                                        "external/draco/src/draco/compression/*.cc",
                                        "external/draco/src/draco/compression/attributes/*.cc",
                                        "external/draco/src/draco/compression/attributes/prediction_schemes/*.cc",
                                        "external/draco/src/draco/compression/bit_coders/*.cc",
                                        "external/draco/src/draco/compression/config/*.cc",
                                        "external/draco/src/draco/compression/entropy/*.cc",
                                        "external/draco/src/draco/compression/mesh/*.cc",
                                        "external/draco/src/draco/compression/mesh/traverser/*.cc",
                                        "external/draco/src/draco/compression/point_cloud/*.cc",
                                        "external/draco/src/draco/compression/point_cloud/algorithms/*.cc",
                                        "external/draco/src/draco/core/*.cc",
                                        "external/draco/src/draco/io/*.cc",
                                        "external/draco/src/draco/material/*.cc",
                                        "external/draco/src/draco/mesh/*.cc",
                                        "external/draco/src/draco/meshdata/*.cc",
                                        "external/draco/src/draco/metadata/*.cc",
                                        "external/draco/src/draco/point_cloud/*.cc",
                                        "external/draco/src/draco/scene/*.cc",
                                        "external/draco/src/draco/texture/*.cc"}
        
                configuration { "windows or linux or macosx or ios or gmake" }
                buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }
        
                configuration { "windows or macosx or linux"}
                    files { ConwayCoreFiles,
                            ManifoldSrcFiles,
                            glTFSDKSrcFiles,
                            DracoSourceFiles,
                            ConwayTestSourceFiles
                        }
        
                configuration {"gmake and not macosx and not windows"}
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
                    "external/manifold/src/third_party/glm",
                    "external/gltf-sdk/GLTFSDK/Inc",
                    "external/gltf-sdk/External/RapidJSON/232389d4f1012dddec4ef84861face2d2ba85709/include",
                    "external/draco/src",
                    "external/fuzzy-bools",
                    "external/fuzzy-bools/deps/cdt"
                    --"/Users/soar/Documents/GitHub/emsdk/upstream/emscripten/system/include"
                    
                    --"$(EMSDK)/upstream/emscripten/system/include"
                }    
        
                excludes 
                { 
                --Manifold Test files 
                "external/manifold/src/third_party/glm/test/**.*",
                "external/manifold/src/third_party/thrust/examples/**.*",
                "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
                "external/manifold/src/third_party/glm/test/gtc/**.*",
                --Draco Source Files
                "external/draco/src/draco/javascript/**.*",
                "external/draco/src/draco/maya/**.*",
                "external/draco/src/draco/tools/**.*",
                "external/draco/src/draco/unity/**.*",
                "external/draco/src/draco/animation/**.*",
                "external/draco/src/draco/io/**.*",
                --Draco Test Files
                "external/draco/src/draco/animation/*test*cc",
                "external/draco/src/draco/attributes/*test*cc",
                "external/draco/src/draco/core/*test*cc",
                "external/draco/src/draco/io/*test*cc",
                "external/draco/src/draco/material/*test*cc",
                "external/draco/src/draco/mesh/*test*cc",
                "external/draco/src/draco/material/*test*cc",
                "external/draco/src/draco/point_cloud/*test*cc",
                "external/draco/src/draco/scene/*test*cc",
                "external/draco/src/draco/texture/*test*cc",
                "external/draco/src/draco/metadata/*test*cc",
                "external/draco/src/draco/compression/*test*cc",
                "external/draco/src/draco/compression/attributes/*test*cc",
                "external/draco/src/draco/compression/attributes/prediction_schemes/*test*cc",
                "external/draco/src/draco/compression/bit_coders/*test*cc",
                "external/draco/src/draco/compression/config/*test*cc",
                "external/draco/src/draco/compression/entropy/*test*cc",
                "external/draco/src/draco/compression/mesh/*test*cc",
                "external/draco/src/draco/compression/mesh/traverser/*test*cc",
                "external/draco/src/draco/compression/point_cloud/*test*cc",
                "external/draco/src/draco/compression/point_cloud/algorithms/*test*cc",
                --glTF-SDK Source Files
                "external/gltf-sdk/GLTFSDK/source/Version.cpp",
                "external/fuzzy-bools/fuzzy/main.cpp"}     
        
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

    project "webifc_native"
    language "C++"
        kind "ConsoleApp"
        files {}

        WebIfcCoreFiles       = { "geometry/**.*", "parsing/**.*", "utility/**.*", "schema/**.*" }
        WebIfcSourceFiles     = { "web-ifc-api.cpp" }
        --WebIfcTestSourceFiles = { "test/*.cpp" }
        WebIfcTestingMain     = { "web-ifc-test.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", "external/manifold/src/collider/include/*.h", "external/manifold/src/utilities/include/*.h"}

        configuration { "windows or linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { WebIfcCoreFiles, ManifoldSrcFiles, WebIfcTestingMain}

        configuration {"gmake and not macosx and not windows"}
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
            "external/fuzzy-bools",
            "external/fuzzy-bools/deps/cdt"
        }    

        excludes { "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*",
        "external/fuzzy-bools/fuzzy/main.cpp" }     

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

        ConwayCoreFiles       = { "conway_geometry/*.h", "conway_geometry/*.cpp", "conway_geometry/operations/**.*", "conway_geometry/representation/**.*"}
        ConwaySourceFiles     = { "conway-api.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}
        glTFSDKSrcFiles       = {"external/gltf-sdk/GLTFSDK/source/**.*"}
        DracoSourceFiles      = {"external/draco/src/draco/animation/*.cc",
                                "external/draco/src/draco/attributes/*.cc",
                                "external/draco/src/draco/compression/*.cc",
                                "external/draco/src/draco/compression/attributes/*.cc",
                                "external/draco/src/draco/compression/attributes/prediction_schemes/*.cc",
                                "external/draco/src/draco/compression/bit_coders/*.cc",
                                "external/draco/src/draco/compression/config/*.cc",
                                "external/draco/src/draco/compression/entropy/*.cc",
                                "external/draco/src/draco/compression/mesh/*.cc",
                                "external/draco/src/draco/compression/mesh/traverser/*.cc",
                                "external/draco/src/draco/compression/point_cloud/*.cc",
                                "external/draco/src/draco/compression/point_cloud/algorithms/*.cc",
                                "external/draco/src/draco/core/*.cc",
                                "external/draco/src/draco/io/*.cc",
                                "external/draco/src/draco/material/*.cc",
                                "external/draco/src/draco/mesh/*.cc",
                                "external/draco/src/draco/meshdata/*.cc",
                                "external/draco/src/draco/metadata/*.cc",
                                "external/draco/src/draco/point_cloud/*.cc",
                                "external/draco/src/draco/scene/*.cc",
                                "external/draco/src/draco/texture/*.cc"}

        configuration { "linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { ConwayCoreFiles,
                    ConwaySourceFiles,
                    ManifoldSrcFiles,
                    glTFSDKSrcFiles,
                    DracoSourceFiles,
                    }
					
		configuration {"windows"}
		prelinkcommands {"$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
		"$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))" ,
		"$(foreach string,$(OBJECTS),\
		$(file >> $(NEWOBJRESP),$(string) )\
		)"}

        configuration {"gmake"}
        linkoptions {"-O3", "--bind", "--dts", "-03", "-flto", "--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_wasm -s ENVIRONMENT=web -s SINGLE_FILE=1 -s EXPORT_ES6=1 -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"] -lworkerfs.js" }
		
        configuration {}
        libdirs {  }
        links {  }
        flags { "Symbols", "FullSymbols", "UseObjectResponseFile" }
		

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
            "external/fuzzy-bools/deps/cdt"
        }    

        excludes 
        { 
        --Manifold Test files 
        "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*",
        --Draco Source Files
        "external/draco/src/draco/javascript/**.*",
        "external/draco/src/draco/maya/**.*",
        "external/draco/src/draco/tools/**.*",
        "external/draco/src/draco/unity/**.*",
		"external/draco/src/draco/animation/**.*",
		"external/draco/src/draco/io/**.*",
        --Draco Test Files
        "external/draco/src/draco/animation/*test*cc",
        "external/draco/src/draco/attributes/*test*cc",
        "external/draco/src/draco/core/*test*cc",
        "external/draco/src/draco/io/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/mesh/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/point_cloud/*test*cc",
        "external/draco/src/draco/scene/*test*cc",
        "external/draco/src/draco/texture/*test*cc",
        "external/draco/src/draco/metadata/*test*cc",
        "external/draco/src/draco/compression/*test*cc",
        "external/draco/src/draco/compression/attributes/*test*cc",
        "external/draco/src/draco/compression/attributes/prediction_schemes/*test*cc",
        "external/draco/src/draco/compression/bit_coders/*test*cc",
        "external/draco/src/draco/compression/config/*test*cc",
        "external/draco/src/draco/compression/entropy/*test*cc",
        "external/draco/src/draco/compression/mesh/*test*cc",
        "external/draco/src/draco/compression/mesh/traverser/*test*cc",
        "external/draco/src/draco/compression/point_cloud/*test*cc",
        "external/draco/src/draco/compression/point_cloud/algorithms/*test*cc",
		--glTF-SDK Source Files
		"external/gltf-sdk/GLTFSDK/source/Version.cpp",
        "external/fuzzy-bools/fuzzy/main.cpp"}         

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


        ConwayCoreFiles       = { "conway_geometry/*.h", "conway_geometry/*.cpp", "conway_geometry/operations/**.*", "conway_geometry/representation/**.*"}
        ConwaySourceFiles     = { "conway-api.cpp" }
        ManifoldSrcFiles      = { "external/manifold/src/**.*", 
        "external/manifold/src/collider/include/*.h", 
        "external/manifold/src/utilities/include/*.h"}
        glTFSDKSrcFiles       = {"external/gltf-sdk/GLTFSDK/source/**.*"}
        DracoSourceFiles      = {"external/draco/src/draco/animation/*.cc",
                                "external/draco/src/draco/attributes/*.cc",
                                "external/draco/src/draco/compression/*.cc",
                                "external/draco/src/draco/compression/attributes/*.cc",
                                "external/draco/src/draco/compression/attributes/prediction_schemes/*.cc",
                                "external/draco/src/draco/compression/bit_coders/*.cc",
                                "external/draco/src/draco/compression/config/*.cc",
                                "external/draco/src/draco/compression/entropy/*.cc",
                                "external/draco/src/draco/compression/mesh/*.cc",
                                "external/draco/src/draco/compression/mesh/traverser/*.cc",
                                "external/draco/src/draco/compression/point_cloud/*.cc",
                                "external/draco/src/draco/compression/point_cloud/algorithms/*.cc",
                                "external/draco/src/draco/core/*.cc",
                                "external/draco/src/draco/io/*.cc",
                                "external/draco/src/draco/material/*.cc",
                                "external/draco/src/draco/mesh/*.cc",
                                "external/draco/src/draco/meshdata/*.cc",
                                "external/draco/src/draco/metadata/*.cc",
                                "external/draco/src/draco/point_cloud/*.cc",
                                "external/draco/src/draco/scene/*.cc",
                                "external/draco/src/draco/texture/*.cc"}

        configuration { "linux or macosx or ios or gmake" }
        buildoptions_cpp { "-O3", "-DNDEBUG", "-pthread", "-Wall", "-fexceptions", "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP", "-std=c++17" }

        configuration { "windows or macosx or linux"}
            files { ConwayCoreFiles,
                    ConwaySourceFiles,
                    ManifoldSrcFiles,
                    glTFSDKSrcFiles,
                    DracoSourceFiles,
                    }
					
		configuration {"windows"}
		prelinkcommands {"$(eval NEWLINKOBJS=$(LINKOBJS)_) $(eval NEWOBJRESP=$(OBJRESP)_) $(eval LINKCMD=$(CXX) -o $(TARGET) $(NEWLINKOBJS) $(RESOURCES) $(ARCH) $(ALL_LDFLAGS) $(LIBS))",
		"$(if $(wildcard $(NEWOBJRESP)), $(shell del $(subst /,\\,$(NEWOBJRESP))))" ,
		"$(foreach string,$(OBJECTS),\
		$(file >> $(NEWOBJRESP),$(string) )\
		)"}

        configuration { "gmake" }
        linkoptions { "-pthread", "-s PTHREAD_POOL_SIZE=navigator.hardwareConcurrency", "--bind", "-03", "-flto", "--define-macro=REAL_T_IS_DOUBLE -s ALLOW_MEMORY_GROWTH=1 -s MAXIMUM_MEMORY=4GB -s FORCE_FILESYSTEM=1 -s EXPORT_NAME=conway_geom_native -s MODULARIZE=1 -s EXPORTED_RUNTIME_METHODS=[\"FS, WORKERFS\"] -lworkerfs.js" }
        configuration {}
        libdirs {  }
        links {  }
        flags { "Symbols", "FullSymbols", "UseObjectResponseFile" }

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
            "external/fuzzy-bools/deps/cdt"
        } 

        excludes 
        { 
        --Manifold Test files 
        "external/manifold/src/third_party/glm/test/**.*",
        "external/manifold/src/third_party/thrust/examples/**.*",
        "external/manifold/src/third_party/thrust/dependencies/cub/test/**.*",
        "external/manifold/src/third_party/glm/test/gtc/**.*",
        --Draco Source Files
        "external/draco/src/draco/javascript/**.*",
        "external/draco/src/draco/maya/**.*",
        "external/draco/src/draco/tools/**.*",
        "external/draco/src/draco/unity/**.*",
		"external/draco/src/draco/animation/**.*",
		"external/draco/src/draco/io/**.*",
        --Draco Test Files
        "external/draco/src/draco/animation/*test*cc",
        "external/draco/src/draco/attributes/*test*cc",
        "external/draco/src/draco/core/*test*cc",
        "external/draco/src/draco/io/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/mesh/*test*cc",
        "external/draco/src/draco/material/*test*cc",
        "external/draco/src/draco/point_cloud/*test*cc",
        "external/draco/src/draco/scene/*test*cc",
        "external/draco/src/draco/texture/*test*cc",
        "external/draco/src/draco/metadata/*test*cc",
        "external/draco/src/draco/compression/*test*cc",
        "external/draco/src/draco/compression/attributes/*test*cc",
        "external/draco/src/draco/compression/attributes/prediction_schemes/*test*cc",
        "external/draco/src/draco/compression/bit_coders/*test*cc",
        "external/draco/src/draco/compression/config/*test*cc",
        "external/draco/src/draco/compression/entropy/*test*cc",
        "external/draco/src/draco/compression/mesh/*test*cc",
        "external/draco/src/draco/compression/mesh/traverser/*test*cc",
        "external/draco/src/draco/compression/point_cloud/*test*cc",
        "external/draco/src/draco/compression/point_cloud/algorithms/*test*cc",
		--glTF-SDK Source Files
		"external/gltf-sdk/GLTFSDK/source/Version.cpp",
        "external/fuzzy-bools/fuzzy/main.cpp"}         

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
