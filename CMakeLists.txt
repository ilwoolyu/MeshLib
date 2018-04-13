cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME Mesh)
project(${PROJECT_NAME})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

SET(VERSION_MAJOR "1")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")

# recursive call: header directories
MACRO(HEADER_DIRECTORIES return_list)
	FILE(GLOB_RECURSE new_list include/*.h)
	SET(dir_list "")
	FOREACH(file_path ${new_list})
		GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
		SET(dir_list ${dir_list} ${dir_path})
	ENDFOREACH()
	LIST(REMOVE_DUPLICATES dir_list)
	SET(${return_list} ${dir_list})
ENDMACRO()
HEADER_DIRECTORIES(HEADER_DIR)
include_directories(${HEADER_DIR})

set(SRC_FILES src/Mesh.cpp src/MeshIO.cpp src/MNIObjIO.cpp src/VtkIO.cpp)
set(SRC_UTIL_FILES src/Util/Geom.cpp src/Util/AABB.cpp src/Util/AABB_Sphere.cpp src/Util/Slicer.cpp src/Util/SphericalHarmonics.cpp src/Util/SurfaceUtil.cpp)

# External libraries
# geodesic distance
set(SRC_GEOD_FILES src/Geodesic/Geodesic.cpp src/Geodesic/GeodesicA.cpp
	src/Geodesic/GeodesicPath.cpp
	src/Geodesic/gw/gw_core/GW_Config.cpp
	src/Geodesic/gw/gw_core/GW_FaceIterator.cpp
	src/Geodesic/gw/gw_core/GW_SmartCounter.cpp
	src/Geodesic/gw/gw_core/GW_VertexIterator.cpp
	src/Geodesic/gw/gw_core/GW_Face.cpp
	src/Geodesic/gw/gw_core/GW_Mesh.cpp
	src/Geodesic/gw/gw_core/GW_Vertex.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicFace.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicMesh.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicMeshA.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicPath.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicPoint.cpp
	src/Geodesic/gw/gw_geodesic/GW_TriangularInterpolation_Cubic.cpp
	src/Geodesic/gw/gw_geodesic/GW_GeodesicVertex.cpp
	src/Geodesic/gw/gw_geodesic/GW_TriangularInterpolation_Linear.cpp
	src/Geodesic/gw/gw_geodesic/GW_TriangularInterpolation_Quadratic.cpp)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

if(BUILD_SHARED_LIBS)
	set(BUILD_SHARED SHARED)
endif(BUILD_SHARED_LIBS)

add_library(Mesh ${BUILD_SHARED} ${SRC_FILES} ${SRC_UTIL_FILES} ${SRC_GEOD_FILES})
set_target_properties(Mesh PROPERTIES SOVERSION 1.0)

if(APPLE)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
       SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
endif(APPLE)

FILE(GLOB_RECURSE MESHLIB_HEADERS include/*.h)

install(FILES
    ${MESHLIB_HEADERS}
    DESTINATION include
    COMPONENT dev)

install(TARGETS Mesh
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)


configure_file(MeshLibConfig.cmake.in
 "${PROJECT_BINARY_DIR}/MeshLibConfig.cmake" @ONLY)
install(FILES
    ${PROJECT_BINARY_DIR}/MeshLibConfig.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}-${VERSION_MAJOR}.${VERSION_MINOR}"
    COMPONENT dev)
