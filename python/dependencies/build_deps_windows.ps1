$ErrorActionPreference = "Stop"

$VcpkgRoot = $env:VCPKG_ROOT
if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) { $VcpkgRoot = "C:\vcpkg" }

# Enable vcpkg binary cache so rebuilt packages are reused across runs
# Put it OUTSIDE the vcpkg repo dir so a reclone does not wipe it
$BinaryCache = "C:\vcpkg_binary_cache"
New-Item -ItemType Directory -Force -Path $BinaryCache | Out-Null
$env:VCPKG_DEFAULT_BINARY_CACHE = $BinaryCache
Write-Host "VCPKG_DEFAULT_BINARY_CACHE=$env:VCPKG_DEFAULT_BINARY_CACHE"

# Triplet
$Triplet = $env:VCPKG_TARGET_TRIPLET
if ([string]::IsNullOrWhiteSpace($Triplet)) { $Triplet = "x64-windows-static" }

# Repo root where vcpkg.json and vcpkg-configuration.json live
$RepoRoot = (Resolve-Path ".").Path
$vcpkgJson = Join-Path $RepoRoot "vcpkg.json"
$vcpkgCfg  = Join-Path $RepoRoot "vcpkg-configuration.json"

Write-Host "=== vcpkg manifest setup ==="
Write-Host "RepoRoot: $RepoRoot"
Write-Host "Triplet:  $Triplet"
Write-Host "vcpkg.json: $(Test-Path $vcpkgJson)  ($vcpkgJson)"
Write-Host "vcpkg-configuration.json: $(Test-Path $vcpkgCfg)  ($vcpkgCfg)"
if (-not (Test-Path $vcpkgJson)) { throw "Missing vcpkg.json at $vcpkgJson" }
if (-not (Test-Path $vcpkgCfg))  { throw "Missing vcpkg-configuration.json at $vcpkgCfg" }

Write-Host ""
Write-Host "=== vcpkg bootstrap ==="

$BootstrapBat = Join-Path $VcpkgRoot "bootstrap-vcpkg.bat"
$VcpkgExe     = Join-Path $VcpkgRoot "vcpkg.exe"

function Ensure-VcpkgCheckout {
  param([string]$Root)

  $bootstrap = Join-Path $Root "bootstrap-vcpkg.bat"
  $exe = Join-Path $Root "vcpkg.exe"

  if ((Test-Path $bootstrap) -or (Test-Path $exe)) {
    Write-Host "Using existing vcpkg at $Root"
    return
  }

  Write-Host "Directory exists but no bootstrap/vcpkg.exe found. Recreating $Root"
  if (Test-Path $Root) {
    Remove-Item -Recurse -Force $Root
  }

  Write-Host "Cloning vcpkg into $Root"
  git clone https://github.com/microsoft/vcpkg $Root
  if ($LASTEXITCODE -ne 0) { throw "git clone vcpkg failed with exit code $LASTEXITCODE" }
}

Ensure-VcpkgCheckout -Root $VcpkgRoot

if (-not (Test-Path $VcpkgExe)) {
  & $BootstrapBat -disableMetrics
  if ($LASTEXITCODE -ne 0) { throw "bootstrap-vcpkg failed with exit code $LASTEXITCODE" }
}

& $VcpkgExe version
if ($LASTEXITCODE -ne 0) { throw "vcpkg version failed with exit code $LASTEXITCODE" }

Write-Host ""
Write-Host "=== vcpkg install (manifest mode) ==="
Push-Location $RepoRoot
& $VcpkgExe install --triplet $Triplet --clean-after-build
if ($LASTEXITCODE -ne 0) { throw "vcpkg install failed with exit code $LASTEXITCODE" }
Pop-Location

Write-Host ""
Write-Host "=== locate installed tree ==="
$installedRoot = Join-Path $RepoRoot "vcpkg_installed\$Triplet"
Write-Host "Installed root: $installedRoot"
Write-Host "Exists? " (Test-Path $installedRoot)
if (-not (Test-Path $installedRoot)) { throw "Could not find vcpkg_installed for triplet at $installedRoot" }

$binDir     = Join-Path $installedRoot "bin"
$libDir     = Join-Path $installedRoot "lib"
$includeDir = Join-Path $installedRoot "include"

Write-Host "bin:     $binDir (exists: $(Test-Path $binDir))"
Write-Host "lib:     $libDir (exists: $(Test-Path $libDir))"
Write-Host "include: $includeDir (exists: $(Test-Path $includeDir))"

Write-Host ""
Write-Host "=== quick HDF5 verification ==="
Write-Host "hdf5.h:"
Get-ChildItem -Recurse -ErrorAction SilentlyContinue $includeDir -Filter "hdf5.h" |
  Select-Object -First 5 FullName | ForEach-Object { Write-Host "  $_" }

Write-Host "libs containing hdf5:"
Get-ChildItem -ErrorAction SilentlyContinue $libDir -Filter "*hdf5*" |
  Select-Object -First 20 Name,Length | ForEach-Object { Write-Host ("  {0} ({1} bytes)" -f $_.Name,$_.Length) }

Write-Host "DLLs named hdf5*.dll in bin (dynamic triplet):"
if (Test-Path $binDir) {
  Get-ChildItem -ErrorAction SilentlyContinue $binDir -Filter "hdf5*.dll" |
    Select-Object -First 20 Name,Length | ForEach-Object { Write-Host ("  {0} ({1} bytes)" -f $_.Name,$_.Length) }
} else {
  Write-Host "  (no bin dir)"
}
