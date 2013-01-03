# -*- cmake -*-
  cmake_minimum_required(VERSION 2.6)

  
SET (this_target STP_Builder)
#PROJECT(${this_target})


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

## section: add definitions
# add prefix -D. example> -DSHP
# - DO NOT add the following definitions(already defined in ${OSP_DEFINITIONS}:
# -DSHP, -DWIN32, -D_WINDOWS, -D_DEBUG, -D_USRDLL, -D_CRT_SECURE_NO_DEPRECATE
ADD_DEFINITIONS(

)

## section: add target
ADD_EXECUTABLE (${this_target} ${${this_target}_SOURCE_FILES} )




