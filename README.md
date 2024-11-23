# conway-geom

`conway-geom` is Bldrs’ rewrite of [web-ifc](https://github.com/ThatOpen/engine_web-ifc), engineered for high-performance and to support the full suite of open CAD standards within the IFC and STEP families.

Features:
- *Fast near-exact CSG*: Our cutting-edge implementation of Constructive Solid Geometry (CSG) partly based on “Exact predicates, exact constructions and combinatorics for mesh CSG” Lévy et al. 2024 [https://arxiv.org/pdf/2405.12949], with robust handling of real-world model data, tailored for the web environment.
- *Seamless NURBS and Advanced BREP*: Conway ensures smooth handling of Non-Uniform Rational B-Splines (NURBS) and Advanced Boundary Representation (BREP), with adaptive tesselation preserving accuracy and quality across complex geometries to streamline workflows for intricate designs.
- *Advanced Model Introspection*: More than pretty pictures, Conway gives type-safe access and powerful querying including geometric relations.
- *Enhanced Memory Handling*: Enables smooth, lag-free interaction with large, highly detailed 3D models across devices, even with complex assets or high-poly designs.
- *Web Geom*: Common libraries shared with including Earcut, Manifold, csgjs-cpp, fast_float, GLM, and tinynurbs for robust, high-speed processing.
- *STEP AP2xx*: Initial support for [STEP AP2xx](https://en.wikipedia.org/wiki/ISO_10303-21) expands Conway’s capabilities to include Automotive and 3D-printing applications.

## IFC Schema Coverage

[Conway](https://github.com/bldrs-ai/conway) generates code bindings for the entire IFC 4.x specs, of 1000+ types. Around 10% of those types have geometric representations.  `conway-geom` covers most of these (80% Nov 2024) and is progressing towards full coverage.  See [IFC_Coverage](https://github.com/bldrs-ai/conway/wiki/IFC-Coverage) for a detailed coverage report.

## Getting Started

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Install the `gmake` and `node` dependencies via Homebrew (```brew install gmake node```).
3. Clone this repo and navigate to the root of the repository.

We've noticed some sentitivity in versions, so here's the most recent working ver:
```
emsdk % ./emsdk activate latest
Resolving SDK alias 'latest' to '3.1.48'
Resolving SDK version '3.1.48' to 'sdk-releases-694434b6d47c5f6eff2c8fbd9eeb016c977ae9dc-64bit'
Setting the following tools as active:
   node-16.20.0-64bit
   python-3.9.2-64bit
   releases-694434b6d47c5f6eff2c8fbd9eeb016c977ae9dc-64bit
emscripten % ./emcc -v
emcc (Emscripten gcc/clang-like replacement + linker emulating GNU ld) 3.1.48 (e967e20b4727956a30592165a3c1cde5c67fa0a8)
clang version 18.0.0 (https://github.com/llvm/llvm-project a54545ba6514802178cf7cf1c1dd9f7efbf3cde7)
Target: wasm32-unknown-emscripten
Thread model: posix
```

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
