param (
    [string]$SOURCE_MAP
)

# Ensure the source map file path is passed as an argument
if (-not $SOURCE_MAP) {
    Write-Host "Error: SOURCE_MAP parameter is not provided."
    exit 1
}

if (-not (Test-Path $SOURCE_MAP)) {
    Write-Host "Source map file not found: $SOURCE_MAP"
    exit 1
}

# Replace escaped JSON \ with /
$EMSDK = $env:EMSDK -replace '\\', '/'

# Replace escaped JSON \ with /
$CWD = (Get-Location).Path -replace '\\', '/'

# Read the content once
$content = Get-Content $SOURCE_MAP -Raw

# Extract the section between the first [ and the first ]
$match = [regex]::Match($content, '\[(.*?)\]')
$section = $match.Groups[1].Value

# Process the section (find beginning of header to /emsdk/)
$processedSection = $section -replace '[^"]*\/emsdk\/', "$EMSDK/"

# NEW: Process paths leading up to /system/lib and replace them with $EMSDK/upstream/emscripten/system/lib
$processedSection = $processedSection -replace '[^"]*/system/lib', "$EMSDK/upstream/emscripten/system/lib"

# Process second section (../../../)
$processedSection = $processedSection -replace '"../../../', "`"$CWD/../"

# Process final section (../../)
$processedSection = $processedSection -replace '"../../', "`"$CWD/../"

# Replace the old section with the processed one in the original content
$content = $content.Replace($section, $processedSection)

# Write the content back
$content | Set-Content $SOURCE_MAP
