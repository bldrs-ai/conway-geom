#!/bin/sh

case "${1}" in
    "test")
        native_config=debug64
        wasm_config=debugemscripten
        ;;
    "clean")
        if [ -d "gmake" ]; then
        ( cd gmake && \
            make config=debug64 clean && \
            make config=debugemscripten clean && \
            make config=release64 clean && \
            make config=releaseemscripten clean
        ) && exit 0 || exit 1
        else
            echo "Directory 'gmake' does not exist, nothing to clean."
            exit 0
        fi
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

if [ -z "$3" ]; then
    ./macos_genie/genie-$(uname -m) gmake
elif [ "$3" = "profile" ]; then
    work_dir="$(pwd)"
    compiled_path="file://$work_dir/../../compiled/dependencies/conway-geom/Dist/"
    ./macos_genie/genie-$(uname -m) gmake profile "$compiled_path"
fi 
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
        ( cd gmake &&
        make config=${native_config} conway_geom_native webifc_native &&
        make config=${wasm_config} ConwayGeomWasm )
    else
        if [ "$2" = "native" ]; then
            ( cd gmake &&
            make config=${native_config} webifc_native )
        elif [ "$2" = "wasmNode" ]; then
            ( cd gmake &&
            make config=${wasm_config} ConwayGeomWasmNode )
        elif [ "$2" = "wasmWeb" ]; then
            ( cd gmake &&
            make config=${wasm_config} ConwayGeomWasmWeb )
        elif [ "$2" = "wasmNodeMT" ]; then
            ( cd gmake &&
            make config=${wasm_config} ConwayGeomWasmNodeMT )
        elif [ "$2" = "wasmWebMT" ]; then
            ( cd gmake &&
            make config=${wasm_config} ConwayGeomWasmWebMT )
        else
            echo "Platform invalid!" 1>&2
            exit 1
        fi
    fi
    if [ $? -ne 0 ]; then
        echo "! Build failed" 1>&2
        exit 1
    fi
fi

if [ -n "$3" ] && [ "$3" = "profile" ]; then
    echo "Remapping source map..."
    
    if [ "$2" = "wasmNode" ]; then
        SOURCE_MAP="./bin/release/ConwayGeomWasmNode.wasm.map"
    elif [ "$2" = "wasmWeb" ]; then
        SOURCE_MAP="./bin/release/ConwayGeomWasmWeb.wasm.map"
    else
        SOURCE_MAP="./bin/release/ConwayGeomWasm.wasm.map"
    fi

    # Use sed to replace the entire path up to and including /emsdk/ with $EMSDK
    sed -i.bak -e 's|[^"]*\/emsdk\/|'"$EMSDK"'\/|g' "$SOURCE_MAP"

    # Use sed to replace the entire path up to and including /system/lib with $EMSDK/upstream/emscripten/system/lib
    sed -i.bak -e 's|[^"]*/system/lib|'"$EMSDK"'/upstream/emscripten/system/lib|g' "$SOURCE_MAP"

    # Prepend current working directory to paths that start with ../../../
    sed -i.bak -e 's|"../../../|"'"$(pwd)"'/|g' "$SOURCE_MAP"

    # Prepend current working directory to paths that start with ../../
    sed -i.bak -e 's|"../../|"'"$(pwd)"'/|g' "$SOURCE_MAP" 

    # Cleanup backup file
    rm "$SOURCE_MAP.bak"
fi

echo "Finished."
