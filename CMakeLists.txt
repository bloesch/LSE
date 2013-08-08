cmake_minimum_required (VERSION 2.6)
project(LSE)

##################### Set some variables #####################
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
SET(CMAKE_BUILD_TYPE "Release")
set (TINYXML "tinyxml")


##################### Find, include, and compile library #####################
find_package(Eigen REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})
#include_directories($ENV{EIGEN_ROOT})
include_directories(./include)
include_directories(${TINYXML})

add_library(LSE src/Manager.cpp src/FilterOCEKF.cpp src/FilterSync.cpp src/DelayCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)

##################### Handle FindLSE.cmake file #####################
file(WRITE cmake/FindLSE.cmake
"set(LSE_INCLUDE_DIRS
	${CMAKE_INSTALL_PREFIX}/include
)

set(LSE_LIBRARIES
  LSE
)

set(LSE_FOUND TRUE)")

##################### Install stuff #####################
install(TARGETS LSE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/LSE)")
install(FILES include/Manager.hpp include/FilterOCEKF.hpp include/FilterSync.hpp include/Rotations.hpp include/DelayCalibration.hpp include/Common.hpp DESTINATION ${CMAKE_INSTALL_PREFIX}/include/LSE)
install(FILES cmake/FindLSE.cmake DESTINATION ${CMAKE_ROOT}/Modules)