#!/bin/sh

case "${1}" in
    "debug")
        native_config=debug64
        wasm_config=debugemscripten
        ;;
    "release")
        native_config=release64
        wasm_config=releaseemscripten
        ;;
    *)
        echo "$0 <debug | release>" 1>&2
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

( cd gmake && \
    make config=${native_config} conway_geom_native webifc_native && \
    make config=${wasm_config} conway_geom_wasm conway_geom_wasm_mt
)
if [ $? -ne 0 ]; then
    echo "! Build failed" 1>&2
    exit 1
fi

echo "Finished."
