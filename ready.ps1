# Ready Script for Parabola Containment Branch Development on Windows
# This script automates downloading Git, GCC toolchain, configuring zlib, cloning and building the parabola repository.

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"

$CurrentDir = $PSScriptRoot
$ToolsDir = "$CurrentDir\tools"
$GitDir = "$ToolsDir\git"
$MinGWDir = "$ToolsDir\mingw64"
$RepoDir = $CurrentDir

Write-Host "=== 1. Creating directories ==="
New-Item -ItemType Directory -Force -Path $ToolsDir | Out-Null

Write-Host "=== 2. Downloading & Installing MinGit ==="
if (-not (Test-Path "$GitDir\cmd\git.exe")) {
    Write-Host "Downloading MinGit..."
    $GitZip = "$ToolsDir\git.zip"
    Invoke-WebRequest -Uri "https://github.com/git-for-windows/git/releases/download/v2.55.0.windows.3/MinGit-2.55.0.3-64-bit.zip" -OutFile $GitZip
    Write-Host "Extracting MinGit..."
    New-Item -ItemType Directory -Force -Path $GitDir | Out-Null
    Expand-Archive -Path $GitZip -DestinationPath $GitDir -Force
    Remove-Item -Path $GitZip -Force
    Write-Host "MinGit installed successfully."
} else {
    Write-Host "MinGit is already installed."
}

# Configure Git to ignore certificate revocation checks to prevent CRYPT_E_NO_REVOCATION_CHECK errors on Windows
& "$GitDir\cmd\git.exe" config --global http.schannelCheckRevoke false


Write-Host "=== 3. Downloading & Installing MinGW-w64 ==="
if (-not (Test-Path "$MinGWDir\bin\gcc.exe")) {
    Write-Host "Downloading MinGW-w64 (WinLibs)..."
    $MinGWZip = "$ToolsDir\mingw.zip"
    Invoke-WebRequest -Uri "https://github.com/brechtsanders/winlibs_mingw/releases/download/16.1.0posix-14.0.0-ucrt-r3/winlibs-x86_64-posix-seh-gcc-16.1.0-mingw-w64ucrt-14.0.0-r3.zip" -OutFile $MinGWZip
    Write-Host "Extracting MinGW-w64..."
    tar.exe -xf $MinGWZip -C $ToolsDir
    Remove-Item -Path $MinGWZip -Force
    # Copy mingw32-make.exe as make.exe
    Copy-Item "$MinGWDir\bin\mingw32-make.exe" "$MinGWDir\bin\make.exe" -Force
    Write-Host "MinGW-w64 installed successfully."
} else {
    Write-Host "MinGW-w64 is already installed."
}

Write-Host "=== 4. Setting up PATH environment variables ==="
$GitBinPath = "$GitDir\cmd"
$MinGWBinPath = "$MinGWDir\bin"
$UserPath = [System.Environment]::GetEnvironmentVariable("Path", "User")
$PathUpdated = $false

if ($UserPath -notlike "*$GitBinPath*") {
    $UserPath = "$GitBinPath;$UserPath"
    $PathUpdated = $true
}
if ($UserPath -notlike "*$MinGWBinPath*") {
    $UserPath = "$MinGWBinPath;$UserPath"
    $PathUpdated = $true
}

if ($PathUpdated) {
    [System.Environment]::SetEnvironmentVariable("Path", $UserPath, "User")
    Write-Host "PATH updated successfully. You may need to restart your terminal to apply changes."
} else {
    Write-Host "PATH is already configured."
}

# Update current process PATH for immediate execution
$env:PATH = "$MinGWBinPath;$GitBinPath;" + $env:PATH

Write-Host "=== 5. Compiling and Installing zlib ==="
if (-not (Test-Path "$MinGWDir\include\zlib.h") -or -not (Test-Path "$MinGWDir\lib\libz.a") -or -not (Test-Path "$MinGWDir\bin\zlib1.dll")) {
    Write-Host "Downloading zlib source..."
    $ZlibTempDir = "$ToolsDir\zlib"
    if (Test-Path $ZlibTempDir) { Remove-Item -Path $ZlibTempDir -Recurse -Force }
    git clone https://github.com/madler/zlib.git $ZlibTempDir
    
    Write-Host "Configuring and compiling zlib..."
    Start-Process -FilePath "make.exe" -ArgumentList "-f win32/Makefile.gcc" -WorkingDirectory $ZlibTempDir -NoNewWindow -Wait
    
    Write-Host "Installing zlib headers and libs..."
    Copy-Item "$ZlibTempDir\zlib.h" "$MinGWDir\include\" -Force
    Copy-Item "$ZlibTempDir\zconf.h" "$MinGWDir\include\" -Force
    Copy-Item "$ZlibTempDir\libz.a" "$MinGWDir\lib\" -Force
    Copy-Item "$ZlibTempDir\zlib1.dll" "$MinGWDir\bin\" -Force
    Copy-Item "$ZlibTempDir\libz.dll.a" "$MinGWDir\lib\" -Force
    
    Remove-Item -Path $ZlibTempDir -Recurse -Force
    Write-Host "zlib configured successfully."
} else {
    Write-Host "zlib is already configured."
}

Write-Host "=== 6. Cloning Parabola Repository ==="
if (-not (Test-Path $RepoDir)) {
    Write-Host "Cloning repository..."
    git clone https://github.com/hsl2001/parabola.git $RepoDir
    git -C $RepoDir checkout containment
    Write-Host "Cloned and switched to containment branch."
} else {
    Write-Host "Parabola repository directory already exists."
    # Ensure on containment branch
    git -C $RepoDir checkout containment
}

Write-Host "=== 7. Applying Windows fixes ==="
$KthreadFile = "$RepoDir\klib\kthread.c"
if (Test-Path $KthreadFile) {
    $Content = Get-Content -Path $KthreadFile -Raw
    if ($Content -notlike "*#include <stdint.h>*") {
        $Content = $Content.Replace("#include <limits.h>", "#include <limits.h>`n#include <stdint.h>")
        Set-Content -Path $KthreadFile -Value $Content
        Write-Host "Patched klib/kthread.c to include <stdint.h>"
    } else {
        Write-Host "klib/kthread.c is already patched."
    }
}

Write-Host "=== 8. Compiling Parabola ==="
if (Test-Path "$RepoDir\Makefile") {
    Write-Host "Running make..."
    # Run make inside RepoDir
    Start-Process -FilePath "make.exe" -WorkingDirectory $RepoDir -NoNewWindow -Wait
    if (Test-Path "$RepoDir\reverb.exe") {
        Write-Host "Parabola compiled successfully! reverb.exe is ready."
    } else {
        Write-Warning "Compilation finished, but reverb.exe was not found."
    }
} else {
    Write-Error "Makefile not found in parabola repository."
}

Write-Host "=== 9. Configuring IDE and .gitignore ==="
# Add tools/ to .gitignore
$GitignoreFile = "$CurrentDir\.gitignore"
if (Test-Path $GitignoreFile) {
    $GitignoreContent = Get-Content -Path $GitignoreFile
    if ($GitignoreContent -notcontains "tools/") {
        Add-Content -Path $GitignoreFile -Value "`ntools/"
        Write-Host "Added 'tools/' to .gitignore"
    } else {
        Write-Host "'tools/' already in .gitignore"
    }
} else {
    Set-Content -Path $GitignoreFile -Value "tools/`n"
    Write-Host "Created .gitignore and added 'tools/'"
}

# Update .vscode/settings.json
$VscodeDir = "$CurrentDir\.vscode"
$SettingsFile = "$VscodeDir\settings.json"
if (-not (Test-Path $VscodeDir)) {
    New-Item -ItemType Directory -Force -Path $VscodeDir | Out-Null
}
$SettingsObj = @{}
if (Test-Path $SettingsFile) {
    try {
        $SettingsObj = Get-Content -Path $SettingsFile -Raw | ConvertFrom-Json -AsHashtable
    } catch {
        Write-Warning "Failed to parse existing settings.json, overwriting."
    }
}
$GitPathNormalized = "$GitDir\cmd\git.exe".Replace("\", "/")
$SettingsObj["git.path"] = $GitPathNormalized
$SettingsObj | ConvertTo-Json -Depth 100 | Set-Content -Path $SettingsFile -Force
Write-Host "Updated .vscode/settings.json with git.path: $GitPathNormalized"

Write-Host "=== Setup Complete! ==="
