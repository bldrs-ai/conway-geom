@echo off

rem add path to emsdk environment below
call D:\emsdk\emsdk_env.bat
call windows_genie\genie gmake

if "%1" NEQ "" (
if "%1" == "debug" (
cd gmake
make config=debugemscripten conway_geom_native
make config=debugemscripten conway_geom_native_mt
)

if "%1" == "release" (
	cd gmake
	make config=releaseemscripten conway_geom_native
	make config=releaseemscripten conway_geom_native_mt
))

echo Finished.