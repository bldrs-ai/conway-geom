@echo off

set EMSDKPATH=%EMSDK%
set FILENAME="\emsdk_env.bat"

rem Combine two paths
set COMBINED="%EMSDK:"=%\%FILENAME:"=%"
call %COMBINED%
call windows_genie\genie gmake

if "%1" NEQ "" (
if "%1" == "debug" (
cd gmake
	make config=debug64 conway_geom_native
	make config=debug64 webifc_native
	make config=debugemscripten conway_geom_wasm
	
	cd ..
)

if "%1" == "release" (
	cd gmake
	make config=release64 conway_geom_native
	make config=release64 webifc_native
	make config=releaseemscripten conway_geom_wasm
	
	cd ..
))

echo Finished.