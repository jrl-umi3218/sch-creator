# -*- cmake -*-
cmake_minimum_required(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/qhull.cmake)
INCLUDE(cmake/cpack.cmake)

set(CXX_DISABLE_WERROR True)
set(CMAKE_VERBOSE_MAKEFILE ON)


SET(PROJECT_NAME sch_creator)
SET(PROJECT_DESCRIPTION "Creation of Sphere-Torus Patches Bounding Volumes")
SET(PROJECT_URL "")

SETUP_PROJECT()

####################################
# DEPENDENCIES
find_package(Boost COMPONENTS program_options REQUIRED)

ADD_REQUIRED_DEPENDENCY("sch-core")

####################################
## section: include directory

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

## section: source files
# Add your source files here (one file per line), please SORT in alphabetical order for future maintenance
SET (${PROJECT_NAME}_SOURCE_FILES
./Main.cpp
./SmoothHullGeneratorVVR.cpp
)

## section: header files
# Add your header files here(one file per line), please SORT in alphabetical order for future maintenance!
SET(${PROJECT_NAME}_HEADER_FILES
./SmoothHullGeneratorVVR.h
	./vector3.h
)

SET_SOURCE_FILES_PROPERTIES(${PROJECT_NAME}_HEADER_FILES
                            PROPERTIES HEADER_FILE_ONLY TRUE)
LIST(APPEND ${PROJECT_NAME}_SOURCE_FILES ${${PROJECT_NAME}_HEADER_FILES})


## section: add target
ADD_EXECUTABLE (${PROJECT_NAME} ${${PROJECT_NAME}_SOURCE_FILES} )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} )

PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} sch-core)

ADD_SUBDIRECTORY(script)

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION bin)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()

