# conway-geom

## Init

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Install the `gmake` and `node` dependencies via Homebrew (```brew install gmake node```).
3. Clone this repo and navigate to the root of the repository.

### Init Win
1. [Install MinGW-64](https://github.com/msys2/msys2-installer/releases/download/2022-06-03/msys2-x86_64-20220603.exe) and add ```g++.exe``` location to your PATH variable. 

## Build

These will genie and use gmake to build wasm modules (single and multithreaded) + native test executables in debug or release mode.

### Build (Mac OS)

1. Add EMSDK environment variable to your terminal path. Example: ```. $EMSDK_DIR/emsdk_env.sh```
1. Run ```git submodule update --init --recursive```
1. `cd dependencies/wasm ; unzip dependencies.zip ; cd -`
1. Run ```build_osx.sh release (wasmNode|wasmWeb) profile``` or ```build_osx.sh release (wasmNode|wasmWeb)```

### Build (Windows)
1. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
1. Add EMSCRIPTEN environment variable, example: ```EMSCRIPTEN=D:\emsdk\upstream\emscripten```
1. Run ```git submodule update --init --recursive```
1. Run ```build_win.bat release (wasmNode|wasmWeb) profile``` or ```build_win.bat release (wasmNode|wasmWeb)```

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

## IFC Schema Coverage

See the wiki page for coverage https://github.com/bldrs-ai/conway/wiki/IFC-Coverage
