######################################################################
# ** Script for searching for the CLHEP library **
# --------------------------------------------_---
#
# For UNIX-like system this module tries to find the clhep-config
# script in the directory ${CLHEP_SEARCHPATH}/bin.
# Then, using this script, extracts the meta-information about CLHEP
# installation.
# Note that some version of clhep-config returns quoted strings.
# If the information returned by the script is not correct, then
# print out an informational message and try to specify the position
# of 'include' and 'lib' directories relative to the
# ${CLHEP_SEARCHPATH} directory.
#
# In Windows we rely on value CLHEP_SEARCHPATH only and
# set the directories 'include' and 'lib' with respect to
# the ${CLHEP_SEARCHPATH) directory.
#
# Input:
#   CLHEP_SEARCHPATH    Set this variable to the canonical path to
#                       the CLHEP installation directory (MUST BE)
#
# Output variables:
#   CLHEP_FOUND                 - we found the CLHEP library
#   CLHEP_VERSION               - version (like CLHEP 2.0.4.5)
#   CLHEP_CONFIG_EXECUTABLE     - clhep-config with full path to it
#   CLHEP_INCLUDE_DIR           - the header directory
#   CLHEP_LIBRARY_DIR           - the library directory
#   CLHEP_INCLUDES              - the include path (with -I)
#   CLHEP_LIBRARIES             - the libraries needed to use CLHEP
#
# Nefedov dot Yury at jinr dot ru
# Based on the RootFind script
######################################################################

# if we already found CLHEP do nothing
IF( NOT CLHEP_FOUND )
  MESSAGE( STATUS "Looking for CLHEP..." )

  # function to check the existence of important files or directories
  FUNCTION( _clhep_check  File )
    IF( NOT EXISTS "${File}" )
      SET( CLHEP_FOUND, NOTFOUND )
      MESSAGE( FATAL_ERROR
        "${File} does not exist. "
        "CLHEP must not be installed correctly."
      )
    ENDIF()
  ENDFUNCTION( _clhep_check )

  # macro to set the position of 'lab' and 'include' directories
  # relative to the ${CLHEP_SEARCHPATH}
  MACRO( _set_relative )
    # Set CLHEP_INCLUDE_DIR
    SET( CLHEP_INCLUDE_DIR "${CLHEP_SEARCHPATH}/include" )
    _clhep_check( ${CLHEP_INCLUDE_DIR} )

    # Set CLHEP_VERSION (parse file CLHEP/ClhepVersion.h)
    _clhep_check( "${CLHEP_INCLUDE_DIR}/CLHEP/ClhepVersion.h" )
    FILE( READ "${CLHEP_INCLUDE_DIR}/CLHEP/ClhepVersion.h" contents )
    MESSAGE( DEBUG "DEBUG: contents=\n${contents}\n" )
    STRING( REGEX REPLACE ".*[ ]String[^\"]+[\"]([^\"]*).*" "\\1"
      CLHEP_VERSION  "${contents}"
    )
    MESSAGE( STATUS
      "Looking for CLHEP... - found version: ${CLHEP_VERSION}" )

    # Set CLHEP_LIBRARY_DIR
    SET( CLHEP_LIBRARY_DIR "${CLHEP_SEARCHPATH}/lib" )
    _clhep_check( ${CLHEP_LIBRARY_DIR} )

  ENDMACRO( _set_relative )

  IF( NOT CLHEP_SEARCHPATH OR NOT EXISTS "${CLHEP_SEARCHPATH}" )
    MESSAGE( FATAL_ERROR
      "Could not find the directory where CLHEP is installed. "
      "Please specify the CLHEP_SEARCHPATH variable "
      "to the canonical path to CLHEP installation directory."
    )
  ENDIF()

  # Set CLHEP_FOUND
  SET( CLHEP_FOUND ${CLHEP_SEARCHPATH} )
  MESSAGE( DEBUG
    "DEBUG: Looking for CLHEP... CLHEP_FOUND= ${CLHEP_FOUND}" )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS

    FIND_PROGRAM( CLHEP_CONFIG_EXECUTABLE
      NAMES clhep-config
      PATHS "${CLHEP_FOUND}/bin"
      NO_DEFAULT_PATH
    )

    # function to call 'clhep-config' and remove the trailing
    # newline from output
    FUNCTION( _clhep_config ARG OUT )
      EXECUTE_PROCESS( COMMAND ${CLHEP_CONFIG_EXECUTABLE} ${ARG}
        OUTPUT_VARIABLE tmp
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      STRING(REPLACE "\n" "" tmp "${tmp}")      # UNIX
      #  STRING(REPLACE "\r\n" "" tmp "${tmp}")    # Win
      SET( ${OUT} ${tmp} PARENT_SCOPE)
      MESSAGE( DEBUG
        "DEBUG: _clhep_config ${ARG} => ${OUT} = ${tmp}" )
    ENDFUNCTION( _clhep_config )

    # check that clhep-config return the credible information
    SET( is_credible "YES" )
    IF( CLHEP_CONFIG_EXECUTABLE )
      _clhep_config( --prefix _prefix_dir )
      STRING( REGEX REPLACE "^\"(.*)\"$" "\\1"
        _prefix_dir "${_prefix_dir}" )
      IF( NOT EXISTS "${_prefix_dir}" )
        MESSAGE( WARNING
          "clhep-config --prefix gives directory: ${_prefix_dir} "
          "which does not exist."
        )
        MESSAGE( WARNING "CLHEP must not be installed correctly" )
        SET( is_credible "NO" )
      ENDIF()
    ELSE( CLHEP_CONFIG_EXECUTABLE )
      MESSAGE( WARNING
        "clhep-config is not found in the ${CLHEP_FOUND}/bin "
        "Let's try to find the 'include' and 'lib' directories "
        "in ${CLHEP_FOUND} ..."
      )
      SET( is_credible "NO" )
    ENDIF( CLHEP_CONFIG_EXECUTABLE )

    IF( is_credible )
      MESSAGE( STATUS "Looking for CLHEP... found clhep-config:\n\t"
        " ${CLHEP_CONFIG_EXECUTABLE}"
      )

      # Set CLHEP_VERSION
      _clhep_config( --version CLHEP_VERSION )
      STRING( REGEX REPLACE "^[^0-9]*([0-9.]+).*$" "\\1"
        CLHEP_VERSION "${CLHEP_VERSION}" )
      MESSAGE( STATUS
        "Looking for CLHEP... found version: ${CLHEP_VERSION}" )

      # Set CLHEP_INCLUDES and CLHEP_INCLUDE_DIR
      _clhep_config( --include CLHEP_INCLUDES )
      STRING( REGEX REPLACE "^-I(.*)$" "\\1"
        CLHEP_INCLUDE_DIR "${CLHEP_INCLUDES}" )
      STRING( REGEX REPLACE "^\"(.*)\"$" "\\1"
        CLHEP_INCLUDE_DIR "${CLHEP_INCLUDE_DIR}" )
      _clhep_check( ${CLHEP_INCLUDE_DIR} )

      # Set CLHEP_LIBRARIES
      _clhep_config( --libs CLHEP_LIBRARIES )

      # Set CLHEP_LIBRARY_DIR
      STRING(REGEX REPLACE "^-L(.*) -l.*$" "\\1"
        CLHEP_LIBRARY_DIR "${CLHEP_LIBRARIES}" )
      STRING( REGEX REPLACE "^\"(.*)\"$" "\\1"
        CLHEP_LIBRARY_DIR "${CLHEP_LIBRARY_DIR}" )
      _clhep_check( ${CLHEP_LIBRARY_DIR} )

    ELSE() # is not credible
      _set_relative()

      # Set CLHEP_INCLUDES
      SET( CLHEP_INCLUDES "-I${CLHEP_INCLUDE_DIR}" )

      # Set CLHEP_LIBRARIES
      SET( CLHEP_LIBRARIES
        "-L${CLHEP_LIBRARY_DIR} -lCLHEP-${CLHEP_VERSION}" )

    ENDIF( is_credible )

  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows ) # Windows
    IF( MSVC )
      MESSAGE( VERBOSE "Compiler: MSVC, version: " ${MSVC_VERSION} )
    ELSE()
      MESSAGE( FATAL_ERROR "Only MSVC is supported on Windows" )
    ENDIF()

    _set_relative()

    # Set CLHEP_INCLUDES
    SET( CLHEP_INCLUDES "/I ${CLHEP_INCLUDE_DIR}" )

    # Set CLHEP_LIBRARIES
    SET( CLHEP_LIBRARIES
      "${CLHEP_LIBRARY_DIR}/CLHEP-${CLHEP_VERSION}.lib" )
    _clhep_check( ${CLHEP_LIBRARIES} )

  ENDIF() # Linux, Darwin, Windows

  MESSAGE( VERBOSE
    "Summary of CLHEP search results: \n"
    "\t version: ${CLHEP_VERSION}\n"
    "\t CLHEP_LIBRARY_DIR= ${CLHEP_LIBRARY_DIR}\n"
    "\t CLHEP_INCLUDE_DIR= ${CLHEP_INCLUDE_DIR}\n"
    "\t CLHEP_LIBRARIES= ${CLHEP_LIBRARIES}\n"
  )

ENDIF( NOT CLHEP_FOUND )
