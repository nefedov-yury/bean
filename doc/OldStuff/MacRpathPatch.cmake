######################################################################
#
# This "patch" is only for dyld of MacOS
#
# MAC_INSTALL_NAME_TOOL is "install_name_tool" to change dynamic
# shared library install names.
#
######################################################################

IF( CMAKE_SYSTEM_NAME MATCHES Darwin ) # MacOS
  # see cmake --help-property MACOSX_RPATH
  # https://blog.kitware.com/upcoming-in-cmake-2-8-12-osx-rpath-support/
  # MESSAGE( STATUS " +++ CMAKE_MACOSX_RPATH= ${CMAKE_MACOSX_RPATH}" )
  IF ( CMAKE_MACOSX_RPATH )
    # unset so that the macros in this module are not executed
    UNSET( MAC_INSTALL_NAME_TOOL )
  ELSE()
    FIND_PROGRAM( MAC_INSTALL_NAME_TOOL NAMES install_name_tool )
    IF( NOT MAC_INSTALL_NAME_TOOL )
      MESSAGE( STATUS
        "I can not find install_name_tool so @rpath will not work."
        " Modify the DYLD_LIBRARY_PATH to include the ROOT, CLHEP and"
        " your \"workdir/BeanLib_x.x.x\" directories."
        " You may want to use workdir/setup.sh script."
        )
    ENDIF()
  ENDIF()
ENDIF()

######################################################################
#
# MACRO for "patching" RPATH name of libraries
#  Parameters:
#       tag_libname - the library target name (see ADD_LIBRARY)
#
# The  @rpath  is the path prefix of a library install name which
# is substituted (by dyld) with each path in the run path list until
# a loadable dylib if found.
#
######################################################################
MACRO( MAC_LIB_RPATH taglibname )
  IF( MAC_INSTALL_NAME_TOOL )
    GET_TARGET_PROPERTY( LIB_PATH_NAME ${taglibname} LOCATION )

    # find actual library name (like libBeanCore.dylib)
    GET_FILENAME_COMPONENT( LIB_NAME ${LIB_PATH_NAME} NAME )
    # print for debugging:
    # MESSAGE( " LIB_PATH_NAME= ${LIB_PATH_NAME}" )
    # MESSAGE( " LIB_NAME= ${LIB_NAME}" )

    # after build lib
    ADD_CUSTOM_COMMAND( TARGET ${libname}
      POST_BUILD
      COMMAND ${MAC_INSTALL_NAME_TOOL}
      ARGS -id @rpath/${install_lib_dir}/${LIB_NAME} ${LIB_PATH_NAME}
      )

  ENDIF()
ENDMACRO( MAC_LIB_RPATH )

######################################################################
#
# MACRO for "patching" RPATH of prog
#  Parameter:
#       prog - the executable target <prog> (see ADD_EXECUTABLE)
#
######################################################################
MACRO( MAC_PROG_RPATH prog )
  IF( MAC_INSTALL_NAME_TOOL )
    GET_TARGET_PROPERTY( PROG_PATH_NAME ${prog} LOCATION )

    # after build prog
    # add both directories to LC_RPATH because there is no
    # good solution for calling commands after installation.
    ADD_CUSTOM_COMMAND( TARGET ${prog}
      POST_BUILD
      COMMAND ${MAC_INSTALL_NAME_TOOL}
      ARGS -add_rpath ${PROJECT_BINARY_DIR} ${PROG_PATH_NAME}
      COMMAND ${MAC_INSTALL_NAME_TOOL}
      ARGS -add_rpath ${CMAKE_INSTALL_PREFIX} ${PROG_PATH_NAME}
      )

  ENDIF()
ENDMACRO( MAC_PROG_RPATH )
