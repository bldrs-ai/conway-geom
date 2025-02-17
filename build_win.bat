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

REM if not exist ".git\modules" (
REM     git submodule update --init
REM     if errorlevel 1 (
REM         echo ! Could not initialize Git submodules
REM         exit /b 1
REM     )
REM )

REM Check if the third argument is empty
if "%~3"=="" (
    .\windows_genie\genie gmake
) else (
    REM Check if the third argument is "profile"
    if "%~3"=="profile" (
       SET compiled_path=file:///%CD%\..\..\compiled\dependencies\conway-geom\Dist\
        .\windows_genie\genie gmake profile !compiled_path!
    )
)
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

if "%1"=="test" (
    cd gmake
    make config=%native_config% conway_geom_native_tests
    ..\bin\64\debug\conway_geom_native_tests
    if %errorlevel% neq 0 (
        echo "! Build failed" >&2
        exit /b 1
    )
) else (
    if "%2"=="" (
        echo "No platform specified, building for native + wasm"
        cd gmake
        make config=%native_config% conway_geom_native webifc_native
        make config=%wasm_config% ConwayGeomWasm
    ) else (
        if "%2"=="native" (
            cd gmake
            make config=%native_config% conway_geom_native webifc_native
        ) else if "%2"=="wasmNode" (
            cd gmake
            make config=%wasm_config% ConwayGeomWasmNode
        ) else if "%2"=="wasmWeb" (
            cd gmake
            make config=%wasm_config% ConwayGeomWasmWeb
        ) else if "%2"=="wasmNodeMT" (
            cd gmake
            make config=%wasm_config% ConwayGeomWasmNodeMT
        ) else if "%2"=="wasmWebMT" (
            cd gmake
            make config=%wasm_config% ConwayGeomWasmWebMT
        ) else (
            echo "Platform invalid" >&2
            exit /b 1
        )
    )
    if %errorlevel% neq 0 (
        echo "! Build failed" >&2
        exit /b 1
    )
)

REM Check if the third argument is "profile"
IF "%~3" NEQ "" IF "%~3" == "profile" (
    echo Remapping source map paths...

    REM Pass the correct source map directly to PowerShell based on the second argument
    IF "%~2" == "wasmWeb" (
        powershell -ExecutionPolicy Bypass -NoLogo -NoProfile -File "%CD%\..\remap_source_map_win.ps1" -SOURCE_MAP "..\bin\release\ConwayGeomWasmWeb.wasm.map"
    ) ELSE IF "%~2" == "wasmNode" (
        powershell -ExecutionPolicy Bypass -NoLogo -NoProfile -File "%CD%\..\remap_source_map_win.ps1" -SOURCE_MAP "..\bin\release\ConwayGeomWasmNode.wasm.map"
    ) ELSE (
        powershell -ExecutionPolicy Bypass -NoLogo -NoProfile -File "%CD%\..\remap_source_map_win.ps1" -SOURCE_MAP "..\bin\release\ConwayGeomWasm.wasm.map"
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

endlocal