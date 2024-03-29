#############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required(VERSION 2.6)

if(APPLE AND ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.0)
  # Identify AppleClang's COMPILER_ID as AppleClang instead of just Clang 
  # (https://gitlab.kitware.com/cmake/cmake/issues/19067)
  cmake_policy(SET CMP0025 NEW)
endif()

project(QUANDENSER)
set(QUANDENSER_SOURCE_DIR ${CMAKE_SOURCE_DIR})
set(MARACLUSTER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/maracluster)
set(PERCOLATOR_SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/percolator)
include("CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Debug" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

option(TARBALL_BUILD "Create a tarball instead of an installer package." OFF)
if(TARBALL_BUILD)
  add_definitions(-DTARBALL_BUILD) # triggers if statements in CPack.txt
endif(TARBALL_BUILD)

option(VENDOR_SUPPORT "Add support for reading vendor RAW files using Proteowizard." OFF)
if(VENDOR_SUPPORT)
  add_definitions(-DVENDOR_SUPPORT) # triggers ifdef statements in the C++ code
endif(VENDOR_SUPPORT)

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building QUANDENSER:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}, ${CMAKE_MODULE_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)

MESSAGE( STATUS "QUANDENSER options:" )
MESSAGE( STATUS "VENDOR_SUPPORT = ${VENDOR_SUPPORT}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

###############################################################################
# COMPILING MARACLUSTER
###############################################################################

include("${MARACLUSTER_SOURCE_DIR}/CommonCMake.txt") # get the "local" version number
load_config_files_with_output_dir("${MARACLUSTER_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/ext/maracluster")

add_subdirectory(ext/maracluster/src ./maracluster EXCLUDE_FROM_ALL)

foreach (TARGET maracluster batchlibrary maraclusterlibrary)
  target_include_directories(${TARGET} PRIVATE ${MARACLUSTER_SOURCE_DIR})
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ext/maracluster/src) # needed for CMake generated files
endforeach()

###############################################################################
# COMPILING PERCOLATOR
###############################################################################

include("${PERCOLATOR_SOURCE_DIR}/CommonCMake.txt") # get the "local" version number
load_config_files_with_output_dir("${PERCOLATOR_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/ext/percolator")

add_subdirectory(ext/percolator/src ./percolator EXCLUDE_FROM_ALL)

foreach (TARGET perclibrary blas fido picked_protein percolator qvality)
  target_include_directories(${TARGET} PRIVATE ${PERCOLATOR_SOURCE_DIR})
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/ext/percolator/src) # needed for CMake generated files
endforeach()

target_compile_options(perclibrary PRIVATE -UBOOST_ERROR_CODE_HEADER_ONLY)

#############################################################################
# COMPILING DINOSAUR
#############################################################################

if (MSVC)
  set(MVN_CMD mvn.cmd)
else (MSVC)
  set(MVN_CMD mvn)
endif (MSVC)

MESSAGE( STATUS "Building Dinosaur with Maven" )

execute_process(COMMAND ${MVN_CMD} package -Pconf -DskipTests -Ddir=${CMAKE_CURRENT_BINARY_DIR}/dinosaur WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/ext/dinosaur OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/dinosaur_maven_stdout.txt ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/dinosaur_maven_stderr.txt RESULT_VARIABLE MVN_RESULT)

if(NOT "${MVN_RESULT}" STREQUAL "0")
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/dinosaur_maven_stdout.txt DINOSAUR_STDOUT)
  file(READ ${CMAKE_CURRENT_BINARY_DIR}/dinosaur_maven_stderr.txt DINOSAUR_STDERR)
  
  message(STATUS "DINOSAUR_STDOUT = ${DINOSAUR_STDOUT}")
  message(STATUS "DINOSAUR_STDERR = ${DINOSAUR_STDERR}")
  
  message(FATAL_ERROR "Cannot create module Dinosaur, mvn execution failed: ${MVN_RESULT}")
endif()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dinosaur/Dinosaur-1.2.1.free.jar ${CMAKE_SOURCE_DIR}/ext/advParams_dinosaur.txt ${CMAKE_SOURCE_DIR}/ext/advParams_dinosaur_targeted.txt DESTINATION ${JAR_PATH})

###############################################################################
# COMPILING QUANDENSER
###############################################################################

include("CommonCMake.txt")
# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files_skip_ext_dir("${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(src)

###############################################################################
# TESTING
###############################################################################

# Enabling system level tests (Ctest)
#enable_testing()
# Scheduling system level tests
#add_subdirectory(data/system_tests/percolator)
# Scheduling unit level tests
#if(GOOGLE_TEST)
#  add_subdirectory(data/unit_tests/percolator)
#endif()

###############################################################################
# INSTALLING
###############################################################################


if( MINGW )
  message( STATUS "    Installing system-libraries: MinGW DLLs." )
  #set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libxerces-c-3-1.dll)
  include( InstallRequiredSystemLibraries )
endif( MINGW )

if( MSVC )
  message( STATUS "    Installing system-libraries: Win32 DLLs." )
  set( CMAKE_INSTALL_OPENMP_LIBRARIES TRUE )
  if (VENDOR_SUPPORT)
    set( DLL_LIBRARIES "${PWIZ_DATA_VENDOR_BRUKER_BAF2SQL_API_LIBRARY};${PWIZ_DATA_VENDOR_WATERS_API_LIBRARY};${PWIZ_DATA_VENDOR_WATERS_IMS_API_LIBRARY}" )
    STRING(REPLACE ".lib" ".dll" DLL_LIBRARIES "${DLL_LIBRARIES}")
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${DLL_LIBRARIES} )
  endif(VENDOR_SUPPORT)
  include( InstallRequiredSystemLibraries )
endif( MSVC )


# ###############################################################################
# # PACKAGING
# ###############################################################################
# # PACKAGING OPTIONS: GENERAL

if(VENDOR_SUPPORT)
  set(NICKNAME_SUFFIX "-vendor-support")
endif(VENDOR_SUPPORT)

set(PACKAGE_NICKNAME "quandenser${NICKNAME_SUFFIX}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Clustering and quantifying mass spectrometry data")
include("CPack.txt")

set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")


include(CPack)
