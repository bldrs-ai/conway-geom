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
		make config=debugemscripten conway_geom_native
		make config=debugemscripten conway_geom_native_mt
	fi

	if [ "$1" = "release" ]
	then
		cd gmake
		make config=releaseemscripten conway_geom_native
		make config=releaseemscripten conway_geom_native_mt
	fi
fi

echo "Finished."