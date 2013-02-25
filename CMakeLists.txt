cmake_minimum_required (VERSION 2.6)
project(LSE)

##################### Set some variables #####################
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
SET(CMAKE_BUILD_TYPE "Release")
set (TINYXML "tinyxml")


##################### Find, include, and compile library #####################
find_package(Eigen REQUIRED)
include_directories(./include)
include_directories(${Eigen_INCLUDE_DIRS})
include_directories(${TINYXML})

# Build a library to be imported as a python module.
set(WRAP_PYTHON TRUE CACHE BOOL "Build Python Wrapper")
if(WRAP_PYTHON)
        LINK_LIBRARIES(boost_python)
        INCLUDE_DIRECTORIES("/usr/include/python2.7")
        LINK_DIRECTORIES("/usr/lib/python2.7")
        ADD_LIBRARY(_PythonManager SHARED src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/DelayCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)
        SET_TARGET_PROPERTIES(_PythonManager PROPERTIES PREFIX "")   
endif(WRAP_PYTHON)

add_library(LSE src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/DelayCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)

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
install(FILES include/Manager.hpp include/PythonManager.hpp include/FilterBase.hpp include/FilterOCEKF.hpp include/FilterVUKF.hpp include/Rotations.hpp include/DelayCalibration.hpp include/Common.hpp DESTINATION ${CMAKE_INSTALL_PREFIX}/include/LSE)
install(FILES cmake/FindLSE.cmake DESTINATION ${CMAKE_ROOT}/Modules)