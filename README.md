# conway-geom

### Build Steps (Mac OS)

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Make sure gmake is installed. 
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Run ```build_osx.sh debug``` or ```build_osx.sh release``` which will run genie and use gmake to build wasm modules (single and multithreaded) + native test executables in debug or release mode.

### Build Steps (Windows)
1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Make sure gmake is installed. 
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Add EMSCRIPTEN environment variable, example: ```EMSCRIPTEN=D:\emsdk\upstream\emscripten```
7. [Install MinGW-64](https://github.com/msys2/msys2-installer/releases/download/2022-06-03/msys2-x86_64-20220603.exe) and add ```g++.exe``` location to your PATH variable. 
8. Run ```build_win.bat debug``` or ```build_win.bat release``` which will run genie and use gmake to build wasm modules (single and multithreaded) + native test executables in debug or release mode.

## Current IFC Schema Coverage (For Ifc-JS Feature Parity)

ifc::IFCPLANE
ifc::IFCBSPLINESURFACE
ifc::IFCBSPLINESURFACEWITHKNOTS
ifc::IFCRATIONALBSPLINESURFACEWITHKNOTS
ifc::IFCAXIS2PLACEMENT2D
ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2D
ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM
ifc::IFCPOLYLINE
ifc::IFCLINE
ifc::IFCINDEXEDPOLYCURVE
ifc::IFCCIRCLE
ifc::IFCELLIPSE
ifc::IFCBSPLINECURVE
ifc::IFCBSPLINECURVEWITHKNOTS
ifc::IFCRATIONALBSPLINECURVEWITHKNOTS
ifc::IFCAXIS1PLACEMENT
ifc::IFCAXIS2PLACEMENT3D
ifc::IFCLOCALPLACEMENT
ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3D
ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM
ifc::IFCCONNECTEDFACESET
ifc::IFCCLOSEDSHELL
ifc::IFCOPENSHELL
ifc::IFCFACE
ifc::IFCADVANCEDFACE
ifc::IFCPOLYLOOP
ifc::IFCINDEXEDPOLYGONALFACE
ifc::IFCPOLYGONALFACESET
~~ifc::IFCMAPPEDITEM~~
~~ifc::IFCBOOLEANCLIPPINGRESULT~~
~~ifc::IFCBOOLEANRESULT~~
~~ifc::IFCHALFSPACESOLID~~
~~ifc::IFCPOLYGONALBOUNDEDHALFSPACE~~
~~ifc::IFCREPRESENTATIONMAP~~
~~ifc::IFCFACEBASEDSURFACEMODEL~~
~~ifc::IFCSHELLBASEDSURFACEMODEL~~
~~ifc::IFCADVANCEDBREP~~
~~ifc::IFCFACETEDBREP~~
~~ifc::IFCPRODUCTREPRESENTATION~~
~~ifc::IFCPRODUCTDEFINITIONSHAPE~~
~~ifc::IFCSHAPEREPRESENTATION~~
~~ifc::IFCTRIANGULATEDFACESET~~
~~ifc::IFCSURFACECURVESWEPTAREASOLID~~
~~ifc::IFCSWEPTDISKSOLID~~
~~ifc::IFCREVOLVEDAREASOLID~~
~~ifc::IFCEXTRUDEDAREASOLID~~
~~ifc::IFCINDEXEDPOLYGONALFACEWITHVOIDS~~
~~ifc::IFCPRESENTATIONSTYLEASSIGNMENT~~
~~ifc::IFCSURFACESTYLE~~
~~ifc::IFCSURFACESTYLERENDERING~~
~~ifc::IFCSURFACESTYLESHADING~~
~~ifc::IFCSTYLEDREPRESENTATION~~
~~ifc::IFCSTYLEDITEM~~
~~ifc::IFCCOLOURRGB~~
~~ifc::IFCMATERIALLAYERSETUSAGE~~
~~ifc::IFCMATERIALLAYERSET~~
~~ifc::IFCMATERIALLAYER~~
~~ifc::IFCMATERIAL~~
~~ifc::IFCFILLAREASTYLE~~
~~ifc::IFCMATERIALLIST~~
~~ifc::IFCMATERIALCONSTITUENTSET~~
~~ifc::IFCMATERIALCONSTITUENT~~
~~ifc::IFCMATERIALPROFILESETUSAGE~~
~~ifc::IFCMATERIALPROFILE~~
~~ifc::IFCMATERIALPROFILESET~~
~~ifc::IFCFACEOUTERBOUND~~
~~ifc::IFCFACEBOUND~~
~~ifc::IFCEDGELOOP~~
~~ifc::IFCEDGECURVE~~
~~ifc::IFCARBITRARYOPENPROFILEDEF~~
~~ifc::IFCARBITRARYCLOSEDPROFILEDEF~~
~~ifc::IFCARBITRARYPROFILEDEFWITHVOIDS~~
~~ifc::IFCRECTANGLEPROFILEDEF~~
~~ifc::IFCROUNDEDRECTANGLEPROFILEDEF~~
~~ifc::IFCRECTANGLEHOLLOWPROFILEDEF~~
~~ifc::IFCCIRCLEPROFILEDEF~~
~~ifc::IFCELLIPSEPROFILEDEF~~
~~ifc::IFCCIRCLEHOLLOWPROFILEDEF~~
~~ifc::IFCISHAPEPROFILEDEF~~
~~ifc::IFCLSHAPEPROFILEDEF~~
~~ifc::IFCUSHAPEPROFILEDEF~~
~~ifc::IFCDERIVEDPROFILEDEF~~
~~ifc::IFCCOMPOSITEPROFILEDEF~~
~~ifc::IFCARBITRARYOPENPROFILEDEF~~
~~ifc::IFCCYLINDRICALSURFACE~~
~~ifc::IFCSURFACEOFREVOLUTION~~
~~ifc::IFCSURFACEOFLINEAREXTRUSION~~
~~ifc::IFCCOMPOSITECURVE~~
~~ifc::IFCCOMPOSITECURVESEGMENT~~
~~ifc::IFCTRIMMEDCURVE~~