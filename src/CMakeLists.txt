set(HEADERS
    ${PROJECT_INCLUDE_DIR}/${PROJECT_NAME}/SmoothHullGeneratorVVR.h
    ${PROJECT_INCLUDE_DIR}/${PROJECT_NAME}/vector3.h)

set(COMMON_SRC
    SmoothHullGeneratorVVR.cpp)

macro(ADD_SCH_CREATOR_EXEC name)
  add_executable(${name} ${name}.cpp ${COMMON_SRC} ${HEADERS} )
  target_link_libraries(${name} Boost::program_options)
  target_include_directories(${name} PUBLIC "${EIGEN3_INCLUDE_DIR}")
  INSTALL(TARGETS ${name} DESTINATION bin)

  target_link_libraries(
    ${name}
    yaml-cpp
  )
  PKG_CONFIG_USE_DEPENDENCY(${name} sch-core)

endmacro()

ADD_SCH_CREATOR_EXEC(sch-creator)
ADD_SCH_CREATOR_EXEC(sch-creator-2d)
ADD_SCH_CREATOR_EXEC(sch-creator-3d)
ADD_SCH_CREATOR_EXEC(supportCheck)
