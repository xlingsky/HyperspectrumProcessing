cmake_minimum_required(VERSION 3.1)
project(HSP)

SET(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS  filesystem REQUIRED)
find_package(glog REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(src/)

SET(OpenMP_LIBS "")
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	ADD_DEFINITIONS(-D_USE_OPENMP)
	SET(_USE_OPENMP TRUE)
	if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang" AND (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.0"))
		SET(OpenMP_LIBS ${OpenMP_libomp_LIBRARY})
	endif()
endif()

	ADD_DEFINITIONS(-D_LOGGING)

set(hsp
  src/main.cc
  src/decode/decode.cc
  src/bigfileio.cc
  src/pos.cc
  src/lx_geometry_rpc.cc
  )

add_executable(hsp ${hsp})

target_link_libraries(hsp PRIVATE ${GDAL_LIBRARY} gflags glog::glog ${Boost_LIBRARIES} opencv_photo opencv_core ${OpenMP_LIBS})
