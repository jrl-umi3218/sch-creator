# -*- cmake -*-
cmake_minimum_required(VERSION 2.6)
# written by Darren Ha(darren.ha@samsung.com)

# replace SmoothCollisionDetection with binary name. e.g. mbase or FBase
SET (this_target STPBV_generator)
PROJECT(${this_target})

find_package(Boost COMPONENTS program_options REQUIRED)

## section: include directory

INCLUDE_DIRECTORIES(
  ../STPBV/include/SCD/Matrix
  ${Boost_INCLUDE_DIRS}
)

## section: source files
# Add your source files here (one file per line), please SORT in alphabetical order for future maintenance
SET (${this_target}_SOURCE_FILES
./Main.cpp
./SmoothHullGeneratorVVR.cpp
)

## section: header files
# Add your header files here(one file per line), please SORT in alphabetical order for future maintenance!
SET(${this_target}_HEADER_FILES
./SmoothHullGeneratorVVR.h
	./vector3.h
)

SET_SOURCE_FILES_PROPERTIES(${this_target}_HEADER_FILES
                            PROPERTIES HEADER_FILE_ONLY TRUE)
LIST(APPEND ${this_target}_SOURCE_FILES ${${this_target}_HEADER_FILES})


## section: add target
ADD_EXECUTABLE (${this_target} ${${this_target}_SOURCE_FILES} )

TARGET_LINK_LIBRARIES(${this_target} ${Boost_LIBRARIES} )

INSTALL(TARGETS ${this_target} DESTINATION bin)

