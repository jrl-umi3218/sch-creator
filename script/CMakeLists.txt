MACRO(INSTALL_SCRIPT SCRIPT)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT}.sh 
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        RENAME ${SCRIPT}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ
       )
ENDMACRO()

CONFIG_FILES(cloud2sch.sh)
INSTALL_SCRIPT(cloud2sch)

CONFIG_FILES(vrml2sch.sh)
INSTALL_SCRIPT(vrml2sch)

# Install packages depending on qhull
SEARCH_FOR_QHULL()
IF(QHULL_FOUND)
  find_program(QCONVEX_EXEC
    NAMES qconvex
    HINTS "${QHULL_ROOT}" "$ENV{QHULL_ROOT}"
    PATH_SUFFIXES project build bin lib)

  CONFIG_FILES(cloud2qhull.sh)
  INSTALL_SCRIPT(cloud2qhull)

  CONFIG_FILES(vrml2qhull.sh)
  INSTALL_SCRIPT(cloud2qhull)
ENDIF(QHULL_FOUND)

# Install packages depending on blender
SET(BLENDER_EXECUTABLE  blender)
EXEC_PROGRAM("${BLENDER_EXECUTABLE} --version >> /dev/null" RETURN_VALUE BlenderExecResult)
IF(${BlenderExecResult} EQUAL 0)
  MESSAGE("Blender found")
  SET(BLENDER_FOUND TRUE)
ELSE()
  MESSAGE("Blender not found")
  SET(BLENDER_FOUND FALSE)
ENDIF(${BlenderExecResult} EQUAL 0)


