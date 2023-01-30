# conway-geom

### Build Steps

1. [Install EMSDK](https://github.com/emscripten-core/emsdk) and follow the instructions for your platform. 
2. Clone this repo and navigate to the root of the repository. 
3. Edit build_osx.sh script to use your ```emsdk/emsdk_env.sh``` path. (TODO: make this dynamic based on environment variables). 
4. Run ```build_osx.sh debug``` or ```build_osx.sh release``` which will run genie and use gmake to build wasm modules (single and multithreaded) in debug or release mode.
