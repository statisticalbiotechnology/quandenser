:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::                                                                           :::
::: tested on a Windows 10 64-bit installation with VS2015 Community          :::
:::                                                                           :::
::: N.B.: VS2015 Express cannot compile x64 natively and will break the build :::
:::                                                                           :::
::: To support reading Thermo RAW Files you need the MSFileReader executable, :::
::: which can be downloaded from https://thermo.flexnetoperations.com/.       :::
::: Search for "MSFileReader" and download and install the 3.0 SP3 version    :::
:::                                                                           :::
::: To support reading Waters RAW Files, Visual C++ runtime is necessary to   :::
::: run http://www.microsoft.com/en-us/download/details.aspx?id=30679         :::
:::                                                                           :::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off

set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build\win64
set RELEASE_DIR=%SRC_DIR%\release\win64
set BUILD_TYPE=Release
set VENDOR=false

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
IF "%~1"=="-e" (set VENDOR=true)
SHIFT
GOTO parse
:endparse

call %SRC_DIR%\quandenser\ext\maracluster\admin\builders\_init_msvc_.bat 64bit

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal

call %SRC_DIR%\quandenser\ext\maracluster\admin\builders\_urls_and_file_names_.bat
call %SRC_DIR%\quandenser\admin\builders\_urls_and_file_names_.bat

set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set ZIP_DIR=%INSTALL_DIR%\%ZIP_BASE%
if not exist "%ZIP_DIR%" (
  echo Downloading and installing 7-Zip
  call :downloadfile %ZIP_URL% %INSTALL_DIR%\7zip.exe
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%ZIP_DIR%\7z.exe"

set CMAKE_DIR=%INSTALL_DIR%\%CMAKE_BASE%
if not exist "%CMAKE_DIR%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%CMAKE_DIR%\bin\cmake.exe"

call %SRC_DIR%\quandenser\ext\maracluster\admin\builders\install_proteowizard.bat 64bit

set MVN_DIR=%INSTALL_DIR%\%MVN_BASE%
if not exist "%MVN_DIR" (
  echo Downloading and installing Maven
  call :downloadfile %MVN_URL% %INSTALL_DIR%\mvn.zip
  %ZIP_EXE% x "%INSTALL_DIR%\mvn.zip" -o"%INSTALL_DIR%" > NUL
)
setlocal
set PATH=%PATH%;%MVN_DIR%\bin

set JAVA_DIR=%INSTALL_DIR%\%JAVA_BASE%
if not exist "%JAVA_DIR" (
  echo Downloading and installing Java JDK   
  call :downloadfile %JAVA_URL% %INSTALL_DIR%\java.zip
  %ZIP_EXE% x "%INSTALL_DIR%\java.zip" -o"%INSTALL_DIR%" > NUL
)
set JAVA_HOME=%JAVA_DIR%

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  call :downloadfile "%NSIS_URL%" %INSTALL_DIR%\nsis.exe
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)
set PATH=%PATH%;%INSTALL_DIR%\nsis

::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: END INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::
:::::::::::: START BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::

if not exist "%BUILD_DIR%" (md "%BUILD_DIR%")

::::::: Building quandenser :::::::
if not exist "%BUILD_DIR%\quandenser" (md "%BUILD_DIR%\quandenser")
cd /D "%BUILD_DIR%\quandenser"
echo cmake quandenser.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -A x64 -DBOOST_ROOT="%PWIZ_DIR%\libraries\%BOOST_BASE%" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\%ZLIB_BASE%" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%\quandenser"

echo build quandenser (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m

if "%VENDOR%" == "true" (
  ::::::: Building quandenser with vendor support :::::::
  if not exist "%BUILD_DIR%\quandenser-vendor-support" (md "%BUILD_DIR%\quandenser-vendor-support")
  cd /D "%BUILD_DIR%\quandenser-vendor-support"
  echo cmake quandenser with vendor support.....
  %CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -A x64 -DBOOST_ROOT="%PWIZ_DIR%\libraries\%BOOST_BASE%" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\%ZLIB_BASE%" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" -DVENDOR_SUPPORT=ON "%SRC_DIR%\quandenser"

  echo build quandenser with vendor support.....
  msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

  ::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
  ::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m
)

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
set /A exit_code=0
call :copytorelease "%BUILD_DIR%\quandenser\quan*.exe"

if "%VENDOR%" == "true" (
  call :copytorelease "%BUILD_DIR%\quandenser-vendor-support\quan*.exe"
)

echo Finished buildscript execution in build directory %BUILD_DIR%

cd /D "%SRC_DIR%"

EXIT /B %exit_code%

:downloadfile
echo Downloading %1 to %2
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B

:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
