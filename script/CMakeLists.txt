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

  IF(QCONVEX_EXEC)
    CONFIG_FILES(cloud2qhull.sh)
    INSTALL_SCRIPT(cloud2qhull)

    CONFIG_FILES(vrml2qhull.sh)
    INSTALL_SCRIPT(cloud2qhull)
  ENDIF(QCONVEX_EXEC)
ENDIF(QHULL_FOUND)

# Install packages depending on blender
IF(UNIX)
  SET(BLENDER_EXECUTABLE  blender)
  EXECUTE_PROCESS(COMMAND "${BLENDER_EXECUTABLE}" --version
                OUTPUT_VARIABLE BlenderExecOutput
                RESULT_VARIABLE BlenderExecResult)
  IF(${BlenderExecResult} EQUAL 0)
    STRING(REGEX MATCH "[0-9.]+" BLENDER_VERSION "${BlenderExecOutput}")
    MESSAGE("Blender found ${BLENDER_VERSION}")
    SET(BLENDER_FOUND TRUE)

    INSTALL(FILES blender/io_qconvex/export_qc.py
                  blender/io_qconvex/__init__.py
            DESTINATION ~/.blender/${BLENDER_VERSION}/scripts/addons/io_qconvex/)
  ELSE()
    MESSAGE("Blender not found")
    SET(BLENDER_FOUND FALSE)
  ENDIF(${BlenderExecResult} EQUAL 0)
ENDIF(UNIX)


### Add some tests.
# Unit testing:
# Create volumes for simple surfaces
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tmp/)
SET(FILE_LIST cube2x2 s2x2)
FOREACH(f ${FILE_LIST})
  #convert the file from vrml to cloud
  ADD_TEST(test_vrml_parsing_${f} ${PYTHON_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/script/vrml2cloud.py
    ${CMAKE_SOURCE_DIR}/examples/wrl/${f}.wrl
    ${CMAKE_BINARY_DIR}/tmp/${f}.cloud)
  ADD_TEST(test_vrml_parsing_${f}_cmp ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_BINARY_DIR}/tmp/${f}.cloud
    ${CMAKE_SOURCE_DIR}/examples/cloud/${f}.cloud
  )

  #convert the file from cloud to sch
  ADD_TEST(NAME test_build_from_cloud_${f}
    COMMAND ${PROJECT_NAME}
    ${CMAKE_SOURCE_DIR}/examples/cloud/${f}.cloud
    ${CMAKE_BINARY_DIR}/tmp/${f}.txt)
  ADD_TEST(test_build_from_cloud_${f}_cmp ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_BINARY_DIR}/tmp/${f}.txt
    ${CMAKE_SOURCE_DIR}/examples/sch/${f}.txt
  )
ENDFOREACH()

