# conway-geom

### Build Steps (Mac OS)

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Install the `gmake` and `node` dependencies via Homebrew (```brew install gmake node```).
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update --init --recursive```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Run ```build_osx.sh debug``` or ```build_osx.sh release``` which will run genie and use gmake to build wasm modules (single and multithreaded) + native test executables in debug or release mode.

### Build Steps (Windows)
1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Make sure gmake is installed. 
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update --init --recursive```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Add EMSCRIPTEN environment variable, example: ```EMSCRIPTEN=D:\emsdk\upstream\emscripten```
7. [Install MinGW-64](https://github.com/msys2/msys2-installer/releases/download/2022-06-03/msys2-x86_64-20220603.exe) and add ```g++.exe``` location to your PATH variable. 
8. Run ```build_win.bat debug``` or ```build_win.bat release``` which will run genie and use gmake to build wasm modules (single and multithreaded) + native test executables in debug or release mode.

## conway_native Usage
1. Running application currently parses geometry from index.ifc and outputs individual obj files + a single complete obj file. 

## webifc_native Usage
1. Launching executable with no arguments loads index.ifc from repository root and parses ifc type information.
2. -i /Path/To/IFC.ifc - Loads external ifc files and parses ifc type information.
3. -objs - outputs geometry from ifc to individudal obj files
4. -stats - outputs currently supported ifc type coverage for a given ifc file
5. -statsv - Outputs currently supported ifc type coverage + prints type name info for a given ifc file
6. -t - Traces and prints IFC type information gathered from parsing a given ifc file
7. -c - Prints generated code for conway geometry processor (Internal only - used for debugging / development - WIP)
8. -h - Displays help.

## Current IFC Schema Coverage (For Ifc-JS Feature Parity)
### 38.88% Complete (51/90)

1.  IFCPLANE
2.  ~~IFCBSPLINESURFACE~~
3.  ~~IFCBSPLINESURFACEWITHKNOTS~~
4.  ~~IFCRATIONALBSPLINESURFACEWITHKNOTS~~
5.  IFCAXIS2PLACEMENT2D
6.  IFCCARTESIANTRANSFORMATIONOPERATOR2D
7.  IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM
8.  IFCPOLYLINE
9.  IFCLINE
10. IFCINDEXEDPOLYCURVE
11. IFCCIRCLE
12. IFCELLIPSE
13. ~~IFCBSPLINECURVE~~
14. ~~IFCBSPLINECURVEWITHKNOTS~~
15. ~~IFCRATIONALBSPLINECURVEWITHKNOTS~~
16. IFCAXIS1PLACEMENT
17. IFCAXIS2PLACEMENT3D
18. IFCLOCALPLACEMENT
19. IFCCARTESIANTRANSFORMATIONOPERATOR3D
20. IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM
21. IFCCONNECTEDFACESET
22. IFCCLOSEDSHELL
23. IFCOPENSHELL
24. IFCFACE
25. IFCADVANCEDFACE
26. IFCPOLYLOOP
27. IFCINDEXEDPOLYGONALFACE
28. IFCPOLYGONALFACESET
29. IFCMAPPEDITEM
30. IFCBOOLEANCLIPPINGRESULT
31. IFCBOOLEANRESULT
32. IFCHALFSPACESOLID
33. IFCPOLYGONALBOUNDEDHALFSPACE
34. ~~IFCREPRESENTATIONMAP~~
35. IFCFACEBASEDSURFACEMODEL
36. IFCSHELLBASEDSURFACEMODEL
37. ~~IFCADVANCEDBREP~~
38. IFCFACETEDBREP
39. ~~IFCPRODUCTREPRESENTATION~~
40. ~~IFCPRODUCTDEFINITIONSHAPE~~
41. IFCSHAPEREPRESENTATION
42. ~~IFCTRIANGULATEDFACESET~~
43. ~~IFCSURFACECURVESWEPTAREASOLID~~
44. ~~IFCSWEPTDISKSOLID~~
45. ~~IFCREVOLVEDAREASOLID~~
46. ~~IFCEXTRUDEDAREASOLID~~
47. ~~IFCINDEXEDPOLYGONALFACEWITHVOIDS~~
48. ~~IFCPRESENTATIONSTYLEASSIGNMENT~~
49. IFCSURFACESTYLE
50. ~~IFCSURFACESTYLERENDERING~~
51. ~~IFCSURFACESTYLESHADING~~
52. ~~IFCSTYLEDREPRESENTATION~~
53. IFCSTYLEDITEM
54. IFCCOLOURRGB
55. IFCMATERIALLAYERSETUSAGE
56. IFCMATERIALLAYERSET
57. IFCMATERIALLAYER
58. IFCMATERIAL
59. ~~IFCFILLAREASTYLE~~
60. IFCMATERIALLIST
61. IFCMATERIALCONSTITUENTSET
62. IFCMATERIALCONSTITUENT
63. ~~IFCMATERIALPROFILESETUSAGE~~
64. IFCMATERIALPROFILE
65. ~~IFCMATERIALPROFILESET~~
66. IFCFACEOUTERBOUND
67. IFCFACEBOUND
68. ~~IFCEDGELOOP~~
69. ~~IFCEDGECURVE~~
70. ~~IFCARBITRARYOPENPROFILEDEF~~
71. IFCARBITRARYCLOSEDPROFILEDEF
72. IFCARBITRARYPROFILEDEFWITHVOIDS
73. IFCRECTANGLEPROFILEDEF
74. IFCROUNDEDRECTANGLEPROFILEDEF
75. ~~IFCRECTANGLEHOLLOWPROFILEDEF~~
76. IFCCIRCLEPROFILEDEF
77. ~~IFCELLIPSEPROFILEDEF~~
78. ~~IFCCIRCLEHOLLOWPROFILEDEF~~
79. ~~IFCISHAPEPROFILEDEF~~
80. ~~IFCLSHAPEPROFILEDEF~~
81. ~~IFCUSHAPEPROFILEDEF~~
82. ~~IFCDERIVEDPROFILEDEF~~
83. ~~IFCCOMPOSITEPROFILEDEF~~
84. ~~IFCARBITRARYOPENPROFILEDEF~~
85. ~~IFCCYLINDRICALSURFACE~~
86. ~~IFCSURFACEOFREVOLUTION~~
87. ~~IFCSURFACEOFLINEAREXTRUSION~~
88. IFCCOMPOSITECURVE
89. ~~IFCCOMPOSITECURVESEGMENT~~
90. IFCTRIMMEDCURVE
