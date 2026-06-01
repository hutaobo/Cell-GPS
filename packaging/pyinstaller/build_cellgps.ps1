param(
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $scriptDir "..\..")
Push-Location $repoRoot

try {
    if (Test-Path ".\build\cellgps") {
        Remove-Item ".\build\cellgps" -Recurse -Force
    }

    if (Test-Path ".\dist\cellgps.exe") {
        Remove-Item ".\dist\cellgps.exe" -Force
    }

    if (Test-Path ".\dist\error.log") {
        Remove-Item ".\dist\error.log" -Force
    }

    & $PythonExe -m PyInstaller --clean --noconfirm ".\packaging\pyinstaller\cellgps.spec"
}
finally {
    Pop-Location
}
