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

1.  ifc::IFCPLANE
2.  ifc::IFCBSPLINESURFACE
3.  ifc::IFCBSPLINESURFACEWITHKNOTS
4.  ifc::IFCRATIONALBSPLINESURFACEWITHKNOTS
5.  ifc::IFCAXIS2PLACEMENT2D
6.  ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2D
7.  ifc::IFCCARTESIANTRANSFORMATIONOPERATOR2DNONUNIFORM
8.  ifc::IFCPOLYLINE
9.  ifc::IFCLINE
10. ifc::IFCINDEXEDPOLYCURVE
11. ifc::IFCCIRCLE
12. ifc::IFCELLIPSE
13. ifc::IFCBSPLINECURVE
14. ifc::IFCBSPLINECURVEWITHKNOTS
15. ifc::IFCRATIONALBSPLINECURVEWITHKNOTS
16. ifc::IFCAXIS1PLACEMENT
17. ifc::IFCAXIS2PLACEMENT3D
18. ifc::IFCLOCALPLACEMENT
19. ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3D
20. ifc::IFCCARTESIANTRANSFORMATIONOPERATOR3DNONUNIFORM
21. ifc::IFCCONNECTEDFACESET
22. ifc::IFCCLOSEDSHELL
23. ifc::IFCOPENSHELL
24. ifc::IFCFACE
25. ifc::IFCADVANCEDFACE
26. ifc::IFCPOLYLOOP
27. ifc::IFCINDEXEDPOLYGONALFACE
28. ifc::IFCPOLYGONALFACESET
29. ~~ifc::IFCMAPPEDITEM~~
30. ~~ifc::IFCBOOLEANCLIPPINGRESULT~~
31. ~~ifc::IFCBOOLEANRESULT~~
32. ~~ifc::IFCHALFSPACESOLID~~
33. ~~ifc::IFCPOLYGONALBOUNDEDHALFSPACE~~
34. ~~ifc::IFCREPRESENTATIONMAP~~
35. ~~ifc::IFCFACEBASEDSURFACEMODEL~~
36. ~~ifc::IFCSHELLBASEDSURFACEMODEL~~
37. ~~ifc::IFCADVANCEDBREP~~
38. ~~ifc::IFCFACETEDBREP~~
39. ~~ifc::IFCPRODUCTREPRESENTATION~~
40. ~~ifc::IFCPRODUCTDEFINITIONSHAPE~~
41. ~~ifc::IFCSHAPEREPRESENTATION~~
42. ~~ifc::IFCTRIANGULATEDFACESET~~
43. ~~ifc::IFCSURFACECURVESWEPTAREASOLID~~
44. ~~ifc::IFCSWEPTDISKSOLID~~
45. ~~ifc::IFCREVOLVEDAREASOLID~~
46. ~~ifc::IFCEXTRUDEDAREASOLID~~
47. ~~ifc::IFCINDEXEDPOLYGONALFACEWITHVOIDS~~
48. ~~ifc::IFCPRESENTATIONSTYLEASSIGNMENT~~
49. ~~ifc::IFCSURFACESTYLE~~
50. ~~ifc::IFCSURFACESTYLERENDERING~~
51. ~~ifc::IFCSURFACESTYLESHADING~~
52. ~~ifc::IFCSTYLEDREPRESENTATION~~
53. ~~ifc::IFCSTYLEDITEM~~
54. ~~ifc::IFCCOLOURRGB~~
55. ~~ifc::IFCMATERIALLAYERSETUSAGE~~
56. ~~ifc::IFCMATERIALLAYERSET~~
57. ~~ifc::IFCMATERIALLAYER~~
58. ~~ifc::IFCMATERIAL~~
59. ~~ifc::IFCFILLAREASTYLE~~
60. ~~ifc::IFCMATERIALLIST~~
61. ~~ifc::IFCMATERIALCONSTITUENTSET~~
62. ~~ifc::IFCMATERIALCONSTITUENT~~
63. ~~ifc::IFCMATERIALPROFILESETUSAGE~~
64. ~~ifc::IFCMATERIALPROFILE~~
65. ~~ifc::IFCMATERIALPROFILESET~~
66. ~~ifc::IFCFACEOUTERBOUND~~
67. ~~ifc::IFCFACEBOUND~~
68. ~~ifc::IFCEDGELOOP~~
69. ~~ifc::IFCEDGECURVE~~
70. ~~ifc::IFCARBITRARYOPENPROFILEDEF~~
71. ~~ifc::IFCARBITRARYCLOSEDPROFILEDEF~~
72. ~~ifc::IFCARBITRARYPROFILEDEFWITHVOIDS~~
73. ~~ifc::IFCRECTANGLEPROFILEDEF~~
74. ~~ifc::IFCROUNDEDRECTANGLEPROFILEDEF~~
75. ~~ifc::IFCRECTANGLEHOLLOWPROFILEDEF~~
76. ~~ifc::IFCCIRCLEPROFILEDEF~~
77. ~~ifc::IFCELLIPSEPROFILEDEF~~
78. ~~ifc::IFCCIRCLEHOLLOWPROFILEDEF~~
79. ~~ifc::IFCISHAPEPROFILEDEF~~
80. ~~ifc::IFCLSHAPEPROFILEDEF~~
81. ~~ifc::IFCUSHAPEPROFILEDEF~~
82. ~~ifc::IFCDERIVEDPROFILEDEF~~
83. ~~ifc::IFCCOMPOSITEPROFILEDEF~~
84. ~~ifc::IFCARBITRARYOPENPROFILEDEF~~
85. ~~ifc::IFCCYLINDRICALSURFACE~~
86. ~~ifc::IFCSURFACEOFREVOLUTION~~
87. ~~ifc::IFCSURFACEOFLINEAREXTRUSION~~
88. ~~ifc::IFCCOMPOSITECURVE~~
89. ~~ifc::IFCCOMPOSITECURVESEGMENT~~
90. ~~ifc::IFCTRIMMEDCURVE~~
