cmake_minimum_required(VERSION 3.1)
project(HSP)

SET(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GDAL CONFIG REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS  filesystem REQUIRED)
find_package(glog REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
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

set(raster
  src/hsp.cc
  src/raster/detail/lsd.c
  )

set(geo
  src/geo.cc
  src/geometry/pos.cc
  src/geometry/lx_geometry_rpc.cc
  src/util/bigfileio.cc
  )

set(extractor
  src/extractor.cc
  src/decode/decode.cc
  src/util/bigfileio.cc
  src/geometry/pos.cc
  )

add_executable(hsp_raster ${raster})
target_link_libraries(hsp_raster PRIVATE GDAL::GDAL gflags glog::glog ${Boost_LIBRARIES} opencv_photo opencv_core opencv_imgcodecs ${OpenMP_LIBS} )

add_executable(hsp_geo ${geo})
target_link_libraries(hsp_geo PRIVATE GDAL::GDAL ${Boost_LIBRARIES})

add_executable(hsp_extractor ${extractor})
target_link_libraries(hsp_extractor PRIVATE GDAL::GDAL gflags glog::glog ${Boost_LIBRARIES})
