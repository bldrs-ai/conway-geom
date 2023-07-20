@echo off
setlocal enabledelayedexpansion

if "%1"=="" (
    echo Usage: %0 ^<clean ^| debug ^| release ^| test^>
    exit /b 1
)

if /i "%1"=="clean" (
    cd gmake && (
        make config=debug64 clean && make config=debugemscripten clean && make config=release64 clean && make config=releaseemscripten clean
    )
    if errorlevel 1 (
        exit /b 1
    ) else (
        exit /b 0
    )
)

if not defined EMSDK (
    echo ! EMSDK environment variable not defined
    echo ! Did you forget to activate the emsdk environment?
    echo !! See https://emscripten.org/docs/getting_started/downloads.html
    exit /b 1
)

call "%EMSDK%\emsdk_env.bat"
set EMSCRIPTEN=%EMSDK%\upstream\emscripten

if not exist ".git\modules" (
    git submodule update --init
    if errorlevel 1 (
        echo ! Could not initialize Git submodules
        exit /b 1
    )
)

.\windows_genie\genie gmake
if errorlevel 1 (
    echo ! Could not generate makefiles
    exit /b 1
)

if /i "%1"=="debug" (
    set native_config=debug64
    set wasm_config=debugemscripten
) else if /i "%1"=="release" (
    set native_config=release64
    set wasm_config=releaseemscripten
) else if /i "%1"=="test" (
    set native_config=debug64
    set wasm_config=debugemscripten
) else (
    echo %0 ^<clean ^| debug ^| release^> 1>&2
    exit /b 1
)

IF "%1"=="test" (
    cd gmake
    make config=%native_config% conway_geom_native_tests
    ..\bin\64\debug\conway_geom_native_tests
    if %errorlevel% neq 0 (
        echo "! Build failed" >&2
        exit /b 1
    )
) ELSE (
    IF "%2"=="" (
        echo "No platform specified, building for native + wasm"
        cd gmake
        make config=%native_config% conway_geom_native webifc_native
        make config=%wasm_config% ConwayGeomWasm
    ) ELSE (
        IF "%2"=="native" (
            cd gmake
            make config=%native_config% conway_geom_native webifc_native
        ) ELSE IF "%2"=="wasm" (
            cd gmake
            make config=%wasm_config% ConwayGeomWasm
        ) ELSE (
            echo "Platform invalid"
        )
    )
    if %errorlevel% neq 0 (
        echo "! Build failed" >&2
        exit /b 1
    )
)

cd ..

if errorlevel 1 (
    echo ! Build failed
    exit /b 1
) else (
    echo Finished.
    exit /b 0
)