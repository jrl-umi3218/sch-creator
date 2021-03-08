# -*- cmake -*-
cmake_minimum_required(VERSION 3.1.4)
if(POLICY CMP0063)
  cmake_policy(SET CMP0063 NEW)
endif()

SET(PROJECT_NAME sch_creator)
SET(PROJECT_DESCRIPTION "Creation of Strictly Convex Hulls (Sphere-Torus Patches Bounding Volumes)")
set(PROJECT_URL "https://github.com/jrl-umi3218/sch-core")
set(PROJECT_USE_CMAKE_EXPORT TRUE)

set(PROJECT_VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
set(CXX_DISABLE_WERROR TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set_property(GLOBAL PROPERTY USE_FOLDERS OFF)


include(cmake/base.cmake)
#INCLUDE(cmake/eigen.cmake)
include(cmake/msvc-specific.cmake)
list(APPEND CMAKE_MODULE_PATH "C:/devel/install/share/eigen3/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

project(
  ${PROJECT_NAME}
  LANGUAGES CXX
  VERSION ${PROJECT_VERSION})

set(DOXYGEN_USE_MATHJAX "YES")
set(DOXYGEN_EXAMPLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/doc")

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_project_dependency(Boost REQUIRED COMPONENTS program_options)
add_project_dependency(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  add_project_dependency(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
else()
  add_project_dependency(Eigen3 MODULE REQUIRED)
endif()

ADD_REQUIRED_DEPENDENCY("sch-core")

file(
  GENERATE
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/include/sch-creator/config.h"
  INPUT "${CMAKE_CURRENT_SOURCE_DIR}/include/sch-creator/config.h.in")
install(FILES "${CMAKE_BINARY_DIR}/$<CONFIG>/include/sch-creator/config.h"
        DESTINATION ${INCLUDE_INSTALL_DESTINATION})

add_subdirectory(src)

ADD_SUBDIRECTORY(script)
