[CmdletBinding()]
param(
    [ValidateSet('dev', 'built')]
    [string]$Mode = 'dev',
    [string]$RunsDirectory = 'M:\ChatGPT\LinkerForge 2\runs\folds',
    [string]$PymolWslProject = 'M:\ChatGPT\LinkerForgeTools\pymol-open-source',
    [string]$VsDevCmd = 'M:\MicrosoftBuildTools\Common7\Tools\VsDevCmd.bat'
)

$ErrorActionPreference = 'Stop'
$repository = Split-Path -Parent $PSScriptRoot
$frontend = Join-Path $repository 'frontend'
$builtExecutable = Join-Path $frontend 'src-tauri\target\release\linkerforge-desktop.exe'

if (-not (Test-Path -LiteralPath $RunsDirectory)) {
    New-Item -ItemType Directory -Path $RunsDirectory | Out-Null
}

$env:LINKERFORGE_ENABLE_WSL = '1'
$env:LINKERFORGE_RUNS_DIR = $RunsDirectory
$env:CARGO_INCREMENTAL = '0'
$env:CARGO_PROFILE_DEV_DEBUG = '0'

if (Test-Path -LiteralPath (Join-Path $PymolWslProject '.pixi\envs\default\bin\pymol')) {
    $env:LINKERFORGE_PYMOL_WSL_PROJECT = $PymolWslProject
}

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    throw 'WSL2 is unavailable; the audited LocalColabFold backend cannot start.'
}
& wsl.exe --exec which colabfold_batch | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw 'colabfold_batch is not available on the WSL2 PATH.'
}
if (-not (Get-Command nvidia-smi.exe -ErrorAction SilentlyContinue)) {
    throw 'nvidia-smi is unavailable; verify the NVIDIA driver before local folding.'
}

Write-Host "Runs: $RunsDirectory"
Write-Host "ColabFold: WSL2 PATH"
Write-Host $(if ($env:LINKERFORGE_PYMOL_WSL_PROJECT) { "PyMOL: WSL2 open-source project at $PymolWslProject" } else { 'PyMOL: Windows detection/fallback' })

if ($Mode -eq 'built') {
    if (-not (Test-Path -LiteralPath $builtExecutable)) {
        throw "Built executable not found at $builtExecutable. Run with -Mode dev or build it first."
    }
    Start-Process -FilePath $builtExecutable
    exit 0
}

if (-not (Test-Path -LiteralPath $VsDevCmd)) {
    throw "Visual Studio developer environment was not found at $VsDevCmd."
}
if (-not (Get-Command npm.cmd -ErrorAction SilentlyContinue)) {
    throw 'npm.cmd is unavailable on PATH.'
}

$environmentAssignments = @(
    "set LINKERFORGE_ENABLE_WSL=$env:LINKERFORGE_ENABLE_WSL",
    "set LINKERFORGE_RUNS_DIR=$env:LINKERFORGE_RUNS_DIR",
    "set LINKERFORGE_PYMOL_WSL_PROJECT=$env:LINKERFORGE_PYMOL_WSL_PROJECT",
    "set CARGO_INCREMENTAL=$env:CARGO_INCREMENTAL",
    "set CARGO_PROFILE_DEV_DEBUG=$env:CARGO_PROFILE_DEV_DEBUG"
) -join ' && '
$command = "call $VsDevCmd -arch=x64 -host_arch=x64 && $environmentAssignments && cd /d `"$frontend`" && npm run desktop"
& cmd.exe /d /s /c $command
exit $LASTEXITCODE
