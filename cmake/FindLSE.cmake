find_package(Ceres REQUIRED)

set(LSE_INCLUDE_DIRS
	/usr/local/include
	/usr/include/eigen3
	/usr/local/include
)

set(LSE_LIBRARIES
  LSE
  ceres_shared
)

set(LSE_FOUND TRUE)