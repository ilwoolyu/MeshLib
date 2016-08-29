cmake_minimum_required(VERSION 2.8)

project(Mesh)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)

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
set(SRC_UTIL_FILES src/Util/Geom.cpp src/Util/AABB.cpp src/Util/Slicer.cpp src/Util/SphericalHarmonics.cpp src/Util/SurfaceUtil.cpp)

# External libraries
# geodesic distance
set(SRC_GEOD_FILES src/Geodesic/Geodesic.cpp src/Geodesic/GeodesicA.cpp
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

add_library(Mesh SHARED ${SRC_FILES} ${SRC_UTIL_FILES} ${SRC_GEOD_FILES})
set_target_properties(Mesh PROPERTIES SOVERSION 1.0)
