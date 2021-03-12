set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME})

set(HEADERS
    ${PROJECT_INCLUDE_DIR}/SmoothHullGeneratorVVR.h
    ${PROJECT_INCLUDE_DIR}/vector3.h)

set(COMMON_SRC
    SmoothHullGeneratorVVR.cpp)

macro(ADD_SCH_CREATOR_EXEC name)
  add_executable(${name} ${name}.cpp ${COMMON_SRC} ${HEADERS} )
  target_link_libraries(${name} PUBLIC Boost::program_options)
  INSTALL(TARGETS ${name} DESTINATION bin)
endmacro()

ADD_SCH_CREATOR_EXEC(sch-creator)
