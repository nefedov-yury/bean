######################################################################
# ** Script for searching for the ROOT program **
# -----------------------------------------------
#
# For UNIX-like system this module tries to find the root-config
# script. Then, using this script extracts the meta-information
# about ROOT installation.
# The root-config script requires bash-shell and a number of standard
# UNIX tools like test, sed, awk and so on.
#
# In Windows we rely on environment variable ROOTSYS. We check that
# the "${ROOTSYS}" directory exists and assume the position of the
# "\bin", "\lib" and "\include" directories in relation to it.
#
# Input:
#       The ROOT_CONFIG_SEARCHPATH variable, if defined, must contain
#       the exact path to the root-config script.
#       Otherwise, we check directory $ROOTSYS/bin using the
#       environment variable ROOTSYS.
#
# Output variables:
#       ROOTSYS - We found ROOT, and ROOTSYS is the path to
#                 the directory in which it is installed.
#
#       ROOT_CONFIG_EXECUTABLE - root-config with full path to it
#       ROOT_CINT_EXECUTABLE   - rootcint (or rootcling for ROOT-6)
#
#       ROOT_LIBRARY_DIR - the library directory
#       ROOT_BINARY_DIR  - the executable directory
#       ROOT_INCLUDE_DIR - the header directory
#       ROOT_LIBRARIES   - regular ROOT libraries
#       ROOT_GLIBS       - regular libs and GUI ROOT libraries
#       ROOT_CFLAGS      - extra compiler flags
#
#       ROOTVERSION      - the version of ROOT (like 6.30.06) and
#       ROOT_MAJOR_VERS, ROOT_MINOR_VERS, ROOT_PATCH_VERS   - are
#               the corresponding components (like 6 30 6)
#
#       ROOT5_GENERATE_DICTIONARY and
#       ROOT6_GENERATE_DICTIONARY   -  cmake functions to create
#               dictionaries for ROOT-5 and ROOT-6 respectively
#
# Nefedov dot Yury at jinr dot ru
# Based on a script written by the F.Uhlig@gsi.de (fairroot.gsi.de)
######################################################################

# if we already found ROOT do nothing
IF( NOT ROOTSYS )
  MESSAGE( STATUS "Looking for Root..." )

  # function to check the existence of important files or directories
  FUNCTION( _root_check  File )
    IF( NOT EXISTS "${File}" )
      SET( ROOTSYS, NOTFOUND )
      MESSAGE( FATAL_ERROR
        "${File} does not exist. "
        "ROOT must not be installed correctly."
      )
    ENDIF()
  ENDFUNCTION( _root_check )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS

    IF( NOT ROOT_CONFIG_SEARCHPATH )
      SET( ROOT_CONFIG_SEARCHPATH $ENV{ROOTSYS}/bin )
    ENDIF()

    FIND_PROGRAM( ROOT_CONFIG_EXECUTABLE
      NAMES root-config
      PATHS ${ROOT_CONFIG_SEARCHPATH}
      NO_DEFAULT_PATH
    )

    IF( ROOT_CONFIG_EXECUTABLE )
      MESSAGE( STATUS
        "Looking for Root... found root-config:\n\t"
        " ${ROOT_CONFIG_EXECUTABLE}"
      )

      # Set ROOTSYS
      STRING( REGEX REPLACE "(^.*)/bin/root-config" "\\1"
        ROOTSYS ${ROOT_CONFIG_EXECUTABLE}
      )

      # function to call 'root-config' and remove the trailing
      # newline from output
      FUNCTION( _root_config ARG OUT )
        EXECUTE_PROCESS( COMMAND ${ROOT_CONFIG_EXECUTABLE} ${ARG}
          OUTPUT_VARIABLE tmp
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        STRING(REPLACE "\n" "" tmp "${tmp}")      # UNIX
        #  STRING(REPLACE "\r\n" "" tmp "${tmp}")    # Win
        SET( ${OUT} ${tmp} PARENT_SCOPE )
        MESSAGE( DEBUG
          "DEBUG: _root_config ${ARG} => ${OUT} = ${tmp}" )
      ENDFUNCTION( _root_config )

      # Set ROOTVERSION
      _root_config( --version ROOTVERSION )

      # Set ROOT_LIBRARY_DIR
      _root_config( --libdir ROOT_LIBRARY_DIR )
      _root_check( ${ROOT_LIBRARY_DIR} )

      # Set ROOT_BINARY_DIR
      _root_config( --bindir ROOT_BINARY_DIR )
      _root_check( ${ROOT_BINARY_DIR} )

      # Set ROOT_INCLUDE_DIR
      _root_config( --incdir ROOT_INCLUDE_DIR )
      _root_check( ${ROOT_INCLUDE_DIR} )

      # Set ROOT_LIBRARIES and ROOT_GLIBS
      _root_config( --libs ROOT_LIBRARIES )
      _root_config( --glibs ROOT_GLIBS )

      # Set ROOT_CFLAGS
      _root_config( --cflags ROOT_CFLAGS )

    ELSE ( ROOT_CONFIG_EXECUTABLE )
      SET( ROOTSYS, NOTFOUND )
      MESSAGE( FATAL_ERROR
        "Could not find the directory where ROOT is installed. "
        "Please set the ROOTSYS environment variable or "
        "set the ROOT_CONFIG_SEARCHPATH cmake variable to "
        "the canonical path to the root-config program."
      )
    ENDIF( ROOT_CONFIG_EXECUTABLE )

  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows ) # Windows

    IF( MSVC )
      MESSAGE( VERBOSE "Compiler: MSVC, version: " ${MSVC_VERSION} )
    ELSE()
      MESSAGE( FATAL_ERROR "Only MSVC is supported on Windows" )
    ENDIF()

    # Set ROOTSYS
    IF( ROOT_CONFIG_SEARCHPATH )
      STRING( REGEX REPLACE "(^.*)[/\]bin" "\\1"
        ROOTSYS ${ROOT_CONFIG_SEARCHPATH}
      )
    ELSE()
      SET( ROOTSYS $ENV{ROOTSYS} )
    ENDIF()

    # converts a native path into a cmake-style path with (/)
    FILE( TO_CMAKE_PATH "${ROOTSYS}" ROOTSYS )
    MESSAGE( STATUS "Looking for Root... ROOTSYS=\n\t${ROOTSYS}" )

    IF( ROOTSYS )
      _root_check( ${ROOTSYS} )
    ELSE()
      MESSAGE( FATAL_ERROR
        "The environment variable ROOTSYS _must be_ defined" )
    ENDIF()

    # Set ROOT_BINARY_DIR
    SET( ROOT_BINARY_DIR "${ROOTSYS}/bin" )
    _root_check( ${ROOT_BINARY_DIR} )

    # Set ROOT_LIBRARY_DIR
    SET( ROOT_LIBRARY_DIR "${ROOTSYS}/lib" )
    _root_check( ${ROOT_LIBRARY_DIR} )

    # Set ROOT_INCLUDE_DIR
    SET( ROOT_INCLUDE_DIR "${ROOTSYS}/include" )
    _root_check( ${ROOT_INCLUDE_DIR} )

    # Check files ROOT/RVersion.hxx (root-6) or RVersion.h (root-5)
    # and set ROOTVERSION
    IF( EXISTS "${ROOT_INCLUDE_DIR}/ROOT/RVersion.hxx" )
      FILE( READ "${ROOT_INCLUDE_DIR}/ROOT/RVersion.hxx" contents )
      MESSAGE( DEBUG "DEBUG: RVersion.hxx=\n${contents}\n" )
      STRING(
        REGEX REPLACE ".*define[ ]+ROOT_VERSION_MAJOR[ ]+([0-9]*).*"
        "\\1" ROOT_MAJOR_VERS "${contents}"
      )
      STRING(
        REGEX REPLACE ".*define[ ]+ROOT_VERSION_MINOR[ ]+([0-9]*).*"
        "\\1" ROOT_MINOR_VERS "${contents}"
      )
      STRING(
        REGEX REPLACE ".*define[ ]+ROOT_VERSION_PATCH[ ]+([0-9]*).*"
        "\\1" ROOT_PATCH_VERS "${contents}"
      )
      SET( ROOTVERSION
        "${ROOT_MAJOR_VERS}.${ROOT_MINOR_VERS}.${ROOT_PATCH_VERS}" )
    ELSEIF( EXISTS "${ROOT_INCLUDE_DIR}/RVersion.h" )
      FILE( READ "${ROOT_INCLUDE_DIR}/RVersion.h" contents )
      MESSAGE( DEBUG "DEBUG: RVersion.h=\n${contents}\n" )
      STRING(
        REGEX REPLACE ".*define[ ]+ROOT_RELEASE[ ]+[\"]([^\"]*).*"
        "\\1" ROOTVERSION  "${contents}"
      )
    ELSE()
      MESSAGE( FATAL_ERROR
        "Can not find RVersion.h(xx) file."
        "ROOT must not be installed correctly."
      )
    ENDIF()

    # Set ROOT_LIBRARIES and ROOT_GLIBS == ROOT_LIBRARIES
    FILE( GLOB ROOT_LIBRARIES ${ROOT_LIBRARY_DIR}/*.lib )
    SET( ROOT_GLIBS ${ROOT_LIBRARIES} )

    # Set ROOT_CFLAGS
    STRING( CONCAT ROOT_CFLAGS
      "-Zc:__cplusplus -std:c++17 -MD -GR -EHs"
      " -D_WIN32 -I${ROOT_INCLUDE_DIR}" )

  ENDIF() # Linux, Darwin, Windows

  MESSAGE( STATUS
    "Looking for Root... found version: ${ROOTVERSION}" )

  IF( NOT DEFINED ROOT_MAJOR_VERS )
    # parse the ROOTVERSION string into variables
    STRING( REGEX REPLACE "^([0-9]+)[.][0-9]+[./][0-9]+.*$" "\\1"
      ROOT_MAJOR_VERS     "${ROOTVERSION}"
    )
    STRING( REGEX REPLACE "^[0-9]+[.]([0-9]+)[./][0-9]+.*$" "\\1"
      ROOT_MINOR_VERS     "${ROOTVERSION}"
    )
    STRING( REGEX REPLACE "^[0-9]+[.][0-9]+[./]([0-9]+).*$" "\\1"
      ROOT_PATCH_VERS     "${ROOTVERSION}"
    )
  ENDIF()

  # check the executables of rootcint and set ROOT_CINT_EXECUTABLE
  SET( _rootcint rootcint )
  IF( ROOT_MAJOR_VERS GREATER 5 )
    SET( _rootcint rootcling )
  ENDIF()
  FIND_PROGRAM( ROOT_CINT_EXECUTABLE
    NAMES ${_rootcint}
    PATHS ${ROOT_BINARY_DIR}
    NO_DEFAULT_PATH
  )
  MESSAGE( STATUS "Looking for Root... found ${_rootcint}:\n\t "
    "${ROOT_CINT_EXECUTABLE}"
  )

  MESSAGE( VERBOSE 
    "Summary of ROOT search results: \n"
    "\t version: ${ROOTVERSION} : <${ROOT_MAJOR_VERS}>, "
    "<${ROOT_MINOR_VERS}>, <${ROOT_PATCH_VERS}>\n"
    "\t ROOT_LIBRARY_DIR= ${ROOT_LIBRARY_DIR}\n"
    "\t ROOT_BINARY_DIR= ${ROOT_BINARY_DIR}\n"
    "\t ROOT_INCLUDE_DIR= ${ROOT_INCLUDE_DIR}\n"
    "\t ROOT_LIBRARIES= ${ROOT_LIBRARIES}\n"
    "\t ROOT_GLIBS= ${ROOT_GLIBS}\n"
    "\t ROOT_CFLAGS= ${ROOT_CFLAGS}\n"
  )

ENDIF()


######################################################################
#
# function for building ROOT dictionaries
#  Parameters:
#       DictNameCxx  - the C++ name of CINT dictionary file
#       IncludeFiles - the list of include files
#       IncludeDirs  - the include file directories to be searched
#       LinkDefH     - the LinkDef.h file (optional)
#
#  Output files:  two dictionary files ${DictNameCxx} and
#                 complimentary header file with ".h" suffix
#
######################################################################

FUNCTION( ROOT5_GENERATE_DICTIONARY
    DictNameCxx IncludeFiles IncludeDirs LinkDefH )

  # add "-I" before each include file directory
  SET( INCLUDE_DIRS )
  FOREACH( dir ${IncludeDirs} )
    SET( INCLUDE_DIRS ${INCLUDE_DIRS} -I${dir} )
  ENDFOREACH()

  STRING( REGEX REPLACE "^(.*)[.](.*)$" "\\1.h"
    DictNameH "${DictNameCxx}" )
  SET( OUTFILES ${DictNameCxx} ${DictNameH} )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
      ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
      -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Darwin )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
      ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
      -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
      -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ENDIF()

ENDFUNCTION( ROOT5_GENERATE_DICTIONARY )


######################################################################
#
# function for building ROOT6 dictionaries
#  Parameters:
#       TagLibName   - the target library name
#       IncludeFiles - the list of include files
#       IncludeDirs  - the include file directories to be searched
#       LinkDefH     - the LinkDef.h file
#
#  Output files:  two dictionary files with names:
#                 TagLibName_rdict.cxx & TagLibName_rdict.pcm
#
######################################################################

FUNCTION( ROOT6_GENERATE_DICTIONARY
    TagLibName IncludeFiles IncludeDirs LinkDefH )

  # add "-I" before each include file directory
  SET( INCLUDE_DIRS )
  FOREACH( dir ${IncludeDirs} )
    SET( INCLUDE_DIRS ${INCLUDE_DIRS} -I${dir} )
  ENDFOREACH()

  SET( DictNameCxx ${TagLibName}_rdict.cxx )
  SET( DictNamePcm ${TagLibName}_rdict.pcm )
  SET( OUTFILES ${DictNameCxx} ${DictNamePcm} )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
      ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
      ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Darwin )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
      ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
      ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
      ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ENDIF()

ENDFUNCTION( ROOT6_GENERATE_DICTIONARY )
