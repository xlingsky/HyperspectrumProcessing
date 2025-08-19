cmake_minimum_required(VERSION 3.10)
project(HSP)

# SET(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL CONFIG REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS  filesystem CONFIG REQUIRED)
find_package(glog REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(src/)

message(STATUS "BOOST: ${Boost_INCLUDE_DIRS}")

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

file(GLOB sofa src/SOFA/*.c)

set(raster
  src/hsp.cc
  src/raster/detail/lsd.c
  )

set(j2cgcs
  src/j2cgcs.cc
  src/geometry/xsofa.cc
  )

set(tracker
  src/tracking.cc
  )

add_executable(hsp ${raster})
target_link_libraries(hsp PRIVATE GDAL::GDAL gflags glog::glog Boost::filesystem opencv_photo opencv_core opencv_imgcodecs ${OpenMP_LIBS} )

add_executable(j2cgcs ${j2cgcs} ${sofa})
target_link_libraries(j2cgcs PRIVATE GDAL::GDAL gflags glog::glog ${Boost_LIBRARIES} )

add_executable(tracker ${tracker} )
target_link_libraries(tracker PRIVATE gflags glog::glog Boost::filesystem opencv_video )