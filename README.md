# conway-geom

### Build Steps (Mac OS)

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Make sure gmake is installed. 
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Run ```build_osx.sh debug``` or ```build_osx.sh release``` which will run genie and use gmake to build wasm modules (single and multithreaded) in debug or release mode.

### Build Steps (Windows)
1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Make sure gmake is installed. 
3. Clone this repo and navigate to the root of the repository.
4. Run ```git submodule update```
5. Add EMSDK environment variable to your terminal path. Example: ```EMSDK=D:\emsdk``` (no trailing slash)
6. Add EMSCRIPTEN environment variable, example: ```EMSCRIPTEN=D:\emsdk\upstream\emscripten```
7. Run ```build_win.bat debug``` or ```build_win.bat release``` which will run genie and use gmake to build wasm modules (single and multithreaded) in debug or release mode.