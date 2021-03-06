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

set(USE_CERES TRUE CACHE BOOL "Use Ceres")
if(USE_CERES)
	add_definitions( -DUSE_CERES=ON )
	find_package(Ceres REQUIRED)
	include_directories(${CERES_INCLUDES})
endif(USE_CERES)

# Build a library to be imported as a python module.
set(WRAP_PYTHON TRUE CACHE BOOL "Build Python Wrapper")
if(WRAP_PYTHON)
	add_definitions( -DWRAP_PYTHON=ON )
	find_package(PythonLibs REQUIRED) 
	LINK_LIBRARIES(boost_python)
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
	LINK_DIRECTORIES(${PYTHON_INCLUDE_DIR})
	
	if(USE_CERES)
		ADD_LIBRARY(_PythonManager SHARED src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/FilterVUKF2.cpp src/FilterInertialOF.cpp src/FilterFLS.cpp src/DelayCalibration.cpp src/RobotCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)
		TARGET_LINK_LIBRARIES(_PythonManager ${CERES_LIBRARIES_SHARED})
	else(USE_CERES)
		ADD_LIBRARY(_PythonManager SHARED src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/FilterVUKF2.cpp src/DelayCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)
	endif(USE_CERES)
	SET_TARGET_PROPERTIES(_PythonManager PROPERTIES PREFIX "")
endif(WRAP_PYTHON)

if(USE_CERES)
	add_library(LSE src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/FilterVUKF2.cpp src/FilterInertialOF.cpp src/FilterFLS.cpp src/DelayCalibration.cpp src/RobotCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)
	TARGET_LINK_LIBRARIES(LSE ${CERES_LIBRARIES_SHARED})
else(USE_CERES)
	add_library(LSE src/Manager.cpp src/PythonManager.cpp src/FilterOCEKF.cpp src/FilterVUKF.cpp src/FilterVUKF2.cpp src/DelayCalibration.cpp ${TINYXML}/tinyxml.cpp ${TINYXML}/tinyxmlerror.cpp ${TINYXML}/tinyxmlparser.cpp)
endif(USE_CERES)

##################### Handle FindLSE.cmake file #####################
if(USE_CERES)
	file(WRITE cmake/FindLSE.cmake
	"find_package(Ceres REQUIRED)
	
	set(LSE_INCLUDE_DIRS
		${CMAKE_INSTALL_PREFIX}/include
		${Eigen_INCLUDE_DIRS}
		${CERES_INCLUDES}
	)
	
	set(LSE_LIBRARIES
	  LSE
	  ${CERES_LIBRARIES_SHARED}
	)
	
	set(LSE_FOUND TRUE)")
else(USE_CERES)
	file(WRITE cmake/FindLSE.cmake
	"set(LSE_INCLUDE_DIRS
		${CMAKE_INSTALL_PREFIX}/include
		${Eigen_INCLUDE_DIRS}
	)
	
	set(LSE_LIBRARIES
	  LSE
	)
	
	set(LSE_FOUND TRUE)")
endif(USE_CERES)

##################### Install stuff #####################
install(TARGETS LSE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/include/LSE)")
install(FILES include/Manager.hpp include/PythonManager.hpp include/FilterBase.hpp include/FilterOCEKF.hpp include/FilterVUKF.hpp include/FilterVUKF2.hpp include/FilterInertialOF.hpp include/FilterFLS.hpp include/Rotations.hpp include/OptimizationFramework.hpp include/DelayCalibration.hpp include/RobotCalibration.hpp include/Common.hpp DESTINATION ${CMAKE_INSTALL_PREFIX}/include/LSE)
install(FILES cmake/FindLSE.cmake DESTINATION ${CMAKE_ROOT}/Modules)
message(${CMAKE_ROOT})