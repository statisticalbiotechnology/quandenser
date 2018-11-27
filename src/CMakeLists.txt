

###############################################################################
# FIND EXTERNAL LIBRARIES REQUIRED BY QUANDENSER
###############################################################################

# INCLUDE HEADER FILES IN SRC DIRECTORY
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(../ext)

# the percolator header files point across subdirectories...
include_directories(../ext/percolator/src)
include_directories(../ext/percolator/src/fido)
include_directories(../ext/percolator/src/picked_protein)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../ext)

MESSAGE( STATUS "Checking packages:" )

if(MSVC)
  if(VENDOR_SUPPORT)
    add_definitions(-DPWIZ_READER_THERMO)
  endif(VENDOR_SUPPORT)
  add_definitions(-DBOOST_ALL_NO_LIB)
  set(BOOST_NOWIDE "nowide")
else(MSVC)
  set(BOOST_NOWIDE "")
endif(MSVC)

find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS filesystem iostreams regex thread serialization system chrono ${BOOST_NOWIDE} REQUIRED)
if(Boost_FOUND)
  message(STATUS "Package Boost found at: ${Boost_INCLUDE_DIRS}")
else(Boost_FOUND)
  message(FATAL_ERROR "Boost has not been found")
endif(Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  message(STATUS "Zlib found at: ${ZLIB_INCLUDE_DIR}; ${ZLIB_LIBRARY}")
else(ZLIB_FOUND)
  message(FATAL_ERROR "The package Zlib has not been found")
endif(ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(ProteoWizard REQUIRED)
if(PWIZ_FOUND)
  message(STATUS "ProteoWizard found at: ${PWIZ_INCLUDE_DIR}")
else(PWIZ_FOUND)
  message(FATAL_ERROR "The package ProteoWizard has not been found")
endif(PWIZ_FOUND)
include_directories(${PWIZ_INCLUDE_DIR})

find_package(Threads REQUIRED)
if(MSVC)
  if(THREADS_FOUND)
    message(STATUS  "Threads found")
  else(THREADS_FOUND)
    message(FATAL_ERROR "The package Threads has not been not found")
  endif()
  include_directories(${Threads_INCLUDE_DIR})
  set(CMAKE_THREAD_LIBS_INIT ${CMAKE_THREAD_LIBS_INIT} ${THREADS_LIBRARY})
endif()

#############################################################################
# COMPILE QUANDENSER
#############################################################################

add_library(quandenserlibrary STATIC Quandenser.cpp MaRaClusterAdapter.cpp MaRaClusterIO.cpp PercolatorAdapter.cpp FeatureAlignment.cpp FeatureGroups.cpp AlignRetention.cpp SpectrumToPrecursorMap.cpp SpectrumFiles.cpp SplineRegression.cpp DinosaurIO.cpp ConsensusMerger.cpp Globals.cpp MyException.cpp Option.cpp)

add_executable(quandenser main.cpp)

if(APPLE)
  set_property(TARGET quandenser PROPERTY LINK_SEARCH_START_STATIC FALSE)
  set_property(TARGET quandenser PROPERTY LINK_SEARCH_END_STATIC TRUE)
#  set_property(TARGET mergespec PROPERTY MACOSX_BUNDLE ON)
elseif(MINGW)
#  set_target_properties(quandenser PROPERTIES LINK_FLAGS "-Wl,-Bstatic -static-libgcc -lcurl -liconv")
  set_target_properties(quandenser PROPERTIES LINK_FLAGS "-static-libgcc")
endif()

if(PROFILING AND UNIX)
  set_target_properties(quandenser PROPERTIES COMPILE_FLAGS -pg LINK_FLAGS -pg)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DNDEBUG -Wno-unused-result -Wno-overflow -Wall -g -pg -fprofile-arcs -ftest-coverage")
endif() 

set(COMMON_LIBRARIES ${PWIZ_LIBRARIES} ${Boost_LIBRARIES} ${ZLIB_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

# SET LIBRARIES FOR QUANDENSER

target_link_libraries(quandenser quandenserlibrary batchlibrary maraclusterlibrary perclibrary fido picked_protein ${COMMON_LIBRARIES})
#if(MSVC)
#  target_link_libraries(extractraw batchlibrary quandenserlibrary ${COMMON_LIBRARIES})
#endif(MSVC)

# INSTALL QUANDENSER
install(TARGETS quandenser EXPORT QUANDENSER DESTINATION bin) # Important to use relative path here (used by CPack)!
#install(TARGETS extractspec msgffixmzml EXPORT QUANDENSER DESTINATION bin) # Important to use relative path here (used by CPack)!
#if(MSVC)
#  install(TARGETS extractraw EXPORT QUANDENSER DESTINATION bin) # Important to use relative path here (used by CPack)!
#endif(MSVC)
