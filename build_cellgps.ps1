param(
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
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

    & $PythonExe -m PyInstaller --clean --noconfirm ".\cellgps.spec"
}
finally {
    Pop-Location
}
