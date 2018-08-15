cmake_minimum_required(VERSION 3.5)

project(MeshLib VERSION "1.0.0")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(MeshLib_INSTALL_CONFIG_DIR "lib/cmake/Mesh-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

set(SRC_FILES
  src/Mesh.cpp
  src/MeshIO.cpp
  src/MNIObjIO.cpp
  src/VtkIO.cpp
  )
set(SRC_UTIL_FILES 
  src/Util/Geom.cpp
  src/Util/AABB.cpp
  src/Util/Slicer.cpp
  src/Util/SphericalHarmonics.cpp
  src/Util/SurfaceUtil.cpp
  )
set(SRC_GEOD_FILES
	src/Geodesic/Geodesic.cpp
	src/Geodesic/GeodesicA.cpp
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
	src/Geodesic/gw/gw_geodesic/GW_TriangularInterpolation_Quadratic.cpp
  )

add_library(Mesh ${SRC_FILES} ${SRC_UTIL_FILES} ${SRC_GEOD_FILES})
set_target_properties(Mesh PROPERTIES SOVERSION 1.0)
target_include_directories(Mesh PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Geodesic>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Geodesic/gw>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Geodesic/gw/gw_core>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Geodesic/gw/gw_geodesic>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Geodesic/gw/gw_maths>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/Util>
  $<INSTALL_INTERFACE:include>
  )

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

file(GLOB_RECURSE MESHLIB_HEADERS include/*.h)
install(
  FILES ${MESHLIB_HEADERS}
  DESTINATION include COMPONENT Development
  )

install(TARGETS Mesh EXPORT MeshLibTargets
  LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION lib COMPONENT Development
  )

# Configure 'MeshLibTargets.cmake'
export(EXPORT MeshLibTargets FILE ${PROJECT_BINARY_DIR}/MeshLibTargets.cmake)

# Configure 'MeshLibConfig.cmake' for a build tree
include(CMakePackageConfigHelpers)
set(build_config ${PROJECT_BINARY_DIR}/MeshLibConfig.cmake)
configure_package_config_file(
  MeshLibConfig.cmake.in
  ${build_config}
  INSTALL_DESTINATION ${PROJECT_BINARY_DIR}
  INSTALL_PREFIX ${PROJECT_BINARY_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

# Configure 'MeshLibConfig.cmake' for an install tree
set(install_config ${PROJECT_BINARY_DIR}/CMakeFiles/MeshLibConfig.cmake)
configure_package_config_file(
  MeshLibConfig.cmake.in
  ${install_config}
  INSTALL_DESTINATION ${MeshLib_INSTALL_CONFIG_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

# Install 'MeshLibConfig.cmake'
install(
  FILES ${install_config}
  DESTINATION ${MeshLib_INSTALL_CONFIG_DIR}
  COMPONENT Development
  )

# Install 'MeshLibTargets.cmake'
install(EXPORT MeshLibTargets
  FILE MeshLibTargets.cmake
  DESTINATION ${MeshLib_INSTALL_CONFIG_DIR}
  COMPONENT Development
  )

