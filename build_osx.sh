#!/bin/sh

#add path to emsdk environment below
source /Users/soar/Documents/GitHub/emsdk/emsdk_env.sh
./macos_genie/genie gmake

if [ "$1" != "" ];
then
	if [ "$1" = "debug" ] 
	then
		cd gmake
		gmake config=debugemscripten conway_geom_native
		gmake config=debugemscripten conway_geom_native_mt
	fi

	if [ "$1" = "release" ]
	then
		cd gmake
		gmake config=releaseemscripten conway_geom_native
		gmake config=releaseemscripten conway_geom_native_mt
	fi
fi

echo "Finished."