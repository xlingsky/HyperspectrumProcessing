include(ExternalProject)

ExternalProject_Add(
    ext_sofa
    PREFIX sofa
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        ${ExternalProject_CMAKE_ARGS}
    BUILD_BYPRODUCTS
        <INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}${lib_name}${lib_suffix}${CMAKE_STATIC_LIBRARY_SUFFIX}
}

ExternalProject_Get_property(ext_sofa INSTALL_DIR)
set(SOFA_INCLUDE_DIR ${INSTALL_DIR}/include/sofa )
set(SOFA_LIB_DIR ${INSTALL_DIR}/lib)
set(SOFA_LIBRARIES ${lib_name})

message(STATUS "SOFA_INCLUDE_DIR: ${SOFA_INCLUDE_DIR}")
message(STATUS "SOFA_LIB_DIR ${SOFA_LIB_DIR}")
message(STATUS "SOFA_LIBRARIES: ${SOFA_LIBRARIES}")