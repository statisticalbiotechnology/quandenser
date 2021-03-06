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
set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set ZIP_URL=http://downloads.sourceforge.net/sevenzip/7z920.exe
if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  call :downloadfile %ZIP_URL% %INSTALL_DIR%\7zip.exe
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

set CMAKE_VERSION=3.16.6
set CMAKE_BASE=cmake-%CMAKE_VERSION%-win32-x86
set CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v%CMAKE_VERSION%/%CMAKE_BASE%.zip
if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"

::: https://teamcity.labkey.org/viewType.html?buildTypeId=bt81 :::
::: without-t = without tests :::
set PWIZ_VERSION_URL=https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/VERSION
call :downloadfile %PWIZ_VERSION_URL% %INSTALL_DIR%\VERSION
set /p PWIZ_VERSION_STRING=<%INSTALL_DIR%\VERSION
set PWIZ_BASE=pwiz-src-without-t-%PWIZ_VERSION_STRING: =_%
set PWIZ_URL=https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/%PWIZ_BASE%.tar.bz2
set PWIZ_DIR=%INSTALL_DIR%\proteowizard
::: Boost asio library
set BOOST_ASIO_BASE=boost_asio_1_10_4
set BOOST_ASIO_URL=https://sourceforge.net/projects/asio/files/asio/1.10.4 (Stable)/%BOOST_ASIO_BASE%.zip/download
if not exist "%PWIZ_DIR%\lib" (
  echo Downloading and installing ProteoWizard
  if not exist "%PWIZ_DIR%" (
    call :downloadfile %PWIZ_URL% %INSTALL_DIR%\pwiz.tar.bz2
    if not exist "%INSTALL_DIR%\pwiz.tar" (
      %ZIP_EXE% x "%INSTALL_DIR%\pwiz.tar.bz2" -o"%INSTALL_DIR%" -aoa > NUL
    )
    %ZIP_EXE% x "%INSTALL_DIR%\pwiz.tar" -o"%PWIZ_DIR%" -aoa > NUL
  )
  
  cd /D "%PWIZ_DIR%"
  echo Starting Proteowizard and Boost build, this can take a while...
  call quickbuild.bat address-model=64 -j4 --toolset=msvc-%MSVC_VER%.0 --i-agree-to-the-vendor-licenses ^
                pwiz/data/common//pwiz_data_common ^
                pwiz/data/identdata//pwiz_data_identdata ^
                pwiz/data/identdata//pwiz_data_identdata_version ^
                pwiz/data/msdata//pwiz_data_msdata ^
                pwiz/data/msdata//pwiz_data_msdata_version ^
                pwiz/data/msdata/mz5//pwiz_data_msdata_mz5 ^
                pwiz/data/proteome//pwiz_data_proteome ^
                pwiz/utility/chemistry//pwiz_utility_chemistry ^
                pwiz/utility/minimxml//pwiz_utility_minimxml ^
                pwiz/utility/misc//SHA1 ^
                pwiz/utility/misc//pwiz_utility_misc ^
                pwiz/data/vendor_readers//pwiz_data_vendor_readers ^
                pwiz/data/vendor_readers/Thermo//pwiz_reader_thermo ^
                pwiz_aux/msrc/utility/vendor_api/thermo//pwiz_vendor_api_thermo ^
                pwiz/data/vendor_readers/Shimadzu//pwiz_reader_shimadzu ^
                pwiz_aux/msrc/utility/vendor_api/Shimadzu//pwiz_vendor_api_shimadzu ^
                pwiz/data/vendor_readers/UIMF//pwiz_reader_uimf ^
                pwiz_aux/msrc/utility/vendor_api/UIMF//pwiz_vendor_api_uimf ^
                pwiz/data/vendor_readers/UNIFI//pwiz_reader_unifi ^
                pwiz_aux/msrc/utility/vendor_api/UNIFI//pwiz_vendor_api_unifi ^
                pwiz/data/vendor_readers/Agilent//pwiz_reader_agilent ^
                pwiz_aux/msrc/utility/vendor_api/Agilent//pwiz_vendor_api_agilent ^
                pwiz/data/vendor_readers/Waters//pwiz_reader_waters ^
                pwiz/data/vendor_readers/Bruker//pwiz_reader_bruker ^
                pwiz_aux/msrc/utility/vendor_api/Bruker//pwiz_vendor_api_bruker ^
                pwiz/data/vendor_readers/ABI/T2D//pwiz_reader_abi_t2d ^
                pwiz_aux/msrc/utility/vendor_api/ABI/T2D//pwiz_vendor_api_abi_t2d ^
                pwiz/data/vendor_readers/ABI//pwiz_reader_abi ^
                pwiz_aux/msrc/utility/vendor_api/ABI//pwiz_vendor_api_abi ^
                /ext/zlib//z ^
                /ext/hdf5//hdf5 ^
                /ext/hdf5//hdf5pp ^
                libraries/SQLite//sqlite3 ^
                libraries/SQLite//sqlite3pp ^
                /ext/boost//system ^
                /ext/boost//thread ^
                /ext/boost//chrono ^
                /ext/boost//regex ^
                /ext/boost//filesystem ^
                /ext/boost//iostreams ^
                /ext/boost//program_options ^
                /ext/boost//nowide ^
                /ext/boost//serialization > pwiz_installation.log 2>&1
    
  mkdir lib
  for /r build-nt-x86_64 %%x in (*.lib) do copy "%%x" lib\ /Y > NUL
  cd lib
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace '-vc%MSVC_VER%0-mt','' })"
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace 'libpwiz_','pwiz_' })"
  Ren libzlib.lib zlib.lib
  Ren libhdf5pp.lib hdf5pp.lib
  Ren libhdf5.lib hdf5.lib
  Ren libSHA1.lib SHA1.lib
  Ren libsqlite3pp.lib sqlite3pp.lib
  Ren libsqlite3.lib sqlite3.lib
  ::: these DLLs might not work, as they are for VS2013 :::
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Waters\vc12_x64\* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x64\baf2sql_c.* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x64\timsdata.* . > NUL
  
  ::: Generate lib from dll for cdt.dll
  setlocal enableDelayedExpansion
  set DLL_BASE=cdt
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  lib /DEF:"!DEF_FILE!" /MACHINE:X64
  endlocal
  
  cd ..

  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y > NUL
  
  ::: copy the boost::asio library, which is not included by the ProteoWizard boost tar but is needed for maracluster
  call :downloadfile "%BOOST_ASIO_URL%" %INSTALL_DIR%\boost_asio.zip
  %ZIP_EXE% x "%INSTALL_DIR%\boost_asio.zip" -o"%INSTALL_DIR%" -aoa > NUL
  PowerShell "Copy-Item -Path '%INSTALL_DIR%\%BOOST_ASIO_BASE%\boost' -Destination '%PWIZ_DIR%\libraries\boost_1_67_0' -Recurse -Force"
)

set MVN_BASE=apache-maven-3.6.3
set MVN_URL=http://ftp.unicamp.br/pub/apache/maven/maven-3/3.6.3/binaries/%MVN_BASE%-bin.zip
if not exist "%INSTALL_DIR%\%MVN_BASE%" (
  echo Downloading and installing Maven
  call :downloadfile %MVN_URL% %INSTALL_DIR%\mvn.zip
  %ZIP_EXE% x "%INSTALL_DIR%\mvn.zip" -o"%INSTALL_DIR%" > NUL
)
setlocal
set PATH=%PATH%;%INSTALL_DIR%\%MVN_BASE%\bin

set JAVA_BASE=openjdk-8u212-b04
set JAVA_URL=https://github.com/AdoptOpenJDK/openjdk8-upstream-binaries/releases/download/jdk8u212-b04/OpenJDK8u-x64_windows_8u212b04.zip
if not exist "%INSTALL_DIR%\%JAVA_BASE%" (
  echo Downloading and installing Java JDK   
  call :downloadfile %JAVA_URL% %INSTALL_DIR%\java.zip
  %ZIP_EXE% x "%INSTALL_DIR%\java.zip" -o"%INSTALL_DIR%" > NUL
)
set JAVA_HOME=%INSTALL_DIR%\%JAVA_BASE%

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS 3/3.04/nsis-3.04-setup.exe/download
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
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -A x64 -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%\quandenser"

echo build quandenser (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m

if "%VENDOR%" == "true" (
  ::::::: Building quandenser with vendor support :::::::
  if not exist "%BUILD_DIR%\quandenser-vendor-support" (md "%BUILD_DIR%\quandenser-vendor-support")
  cd /D "%BUILD_DIR%\quandenser-vendor-support"
  echo cmake quandenser with vendor support.....
  %CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -A x64 -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" -DVENDOR_SUPPORT=ON "%SRC_DIR%\quandenser"

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
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B

:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
