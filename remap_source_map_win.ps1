$SOURCE_MAP = "..\bin\release\ConwayGeomWasm.wasm.map"

#replace escaped JSON \ with /
$EMSDK = $env:EMSDK -replace '\\', '/'

#replace escaped JSON \ with /
$CWD = (Get-Location).Path -replace '\\', '/'

# Read the content once
$content = Get-Content $SOURCE_MAP -Raw

# Extract the section between the first [ and the first ]
$match = [regex]::Match($content, '\[(.*?)\]')
$section = $match.Groups[1].Value

# Process the section (find beginning of header to /emsdk/
# TODO: will /emsdk/ folder always be emsdk? Maybe target /upstream/?
$processedSection = $section -replace '[^"]*\/emsdk\/', "$EMSDK/"

# Process second section (../../../)
$processedSection = $processedSection -replace '"../../../', "`"$CWD/../"

# Process final section (../../)
$processedSection = $processedSection -replace '"../../', "`"$CWD/../"

# Replace the old section with the processed one in the original content
$content = $content.Replace($section, $processedSection)

# Write the content back
$content | Set-Content $SOURCE_MAP