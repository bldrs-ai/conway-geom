#!/bin/sh

#add path to emsdk environment below
envPath=$EMSDK

envPath+='/emsdk_env.sh'

echo ${envPath}

source ${envPath}
./macos_genie/genie gmake

if [ "$1" != "" ];
then
	if [ "$1" = "debug" ] 
	then
		cd gmake
		make config=debug64 conway_geom_native
		make config=debug64 webifc_native
		make config=debugemscripten conway_geom_wasm
		make config=debugemscripten conway_geom_wasm_mt
	fi

	if [ "$1" = "release" ]
	then
		cd gmake
		make config=release64 conway_geom_native
		make config=release64 webifc_native
		make config=releaseemscripten conway_geom_wasm
		make config=releaseemscripten conway_geom_wasm_mt
	fi
fi

echo "Finished."