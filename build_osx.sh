#!/bin/sh

case "${1}" in
    "test")
        native_config=debug64
        wasm_config=debugemscripten
        ;;
    "clean")
        ( cd gmake && \
            gmake config=debug64 clean && \
            gmake config=debugemscripten clean && \
            gmake config=release64 clean && \
            gmake config=releaseemscripten clean
        ) && exit 0 || exit 1
        ;;
    "debug")
        native_config=debug64
        wasm_config=debugemscripten
        ;;
    "release")
        native_config=release64
        wasm_config=releaseemscripten
        ;;
    *)
        echo "$0 <clean | debug | release | test> <native | wasm>" 1>&2
        exit 1
        ;;
esac

if [ -z "${EMSDK}" ]; then
    echo "! EMSDK environment variable not defined" 1>&2
    echo "! Did you forget to activate the emsdk environment?" 1>&2
    echo "!! See https://emscripten.org/docs/getting_started/downloads.html" 1>&2
    exit 1
fi

source "${EMSDK}/emsdk_env.sh"
export EMSCRIPTEN="${EMSDK}/upstream/emscripten"

if [ "$(/usr/bin/arch)" == "arm64" ]; then
    sed -I '' -e "s@^NODE_JS = .*@NODE_JS = '$(brew --prefix node)/bin/node'@" "${EMSDK}/.emscripten"
fi

if [ ! -d ".git/modules" ]; then
    git submodule update --init
    if [ $? -ne 0 ]; then
        echo "! Could not initialize Git submodules" 1>&2
        exit 1
    fi
fi

./macos_genie/genie-$(uname -m) gmake
if [ $? -ne 0 ]; then
    echo "! Could not generate makefiles" 1>&2
    exit 1
fi

if [ "$1" = "test" ]; then
    # Your code for the 'test' case goes here
    ( cd gmake &&
    make config=${native_config} conway_geom_native_tests &&
    ../bin/64/debug/conway_geom_native_tests
    )
    if [ $? -ne 0 ]; then
        echo "! Build failed" 1>&2
        exit 1
    fi
else
    if [ -z "$2" ]; then
        echo "No platform specified, building for native + wasm"
        # Your code for other cases goes here
        ( cd gmake &&
        make config=${native_config} conway_geom_native webifc_native &&
        make config=${wasm_config} ConwayGeomWasm )
    else
        echo $2
        if [ "$2" = "native" ]; then
            ( cd gmake &&
            make config=${native_config} conway_geom_native webifc_native )
        elif [ "$2" = "wasm" ]; then
            ( cd gmake &&
            make config=${wasm_config} ConwayGeomWasm )
        else
            echo "Platform invalid, building for native + wasm"
            # Your code for other cases goes here
            ( cd gmake &&
            make config=${native_config} conway_geom_native webifc_native &&
            make config=${wasm_config} ConwayGeomWasm )
        fi
    fi
    if [ $? -ne 0 ]; then
        echo "! Build failed" 1>&2
        exit 1
    fi
fi

echo "Finished."
