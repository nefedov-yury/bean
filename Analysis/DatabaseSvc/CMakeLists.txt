#---------------------------------------------------------------------
# DatabaseSvc: interface library to sqlite database
#              <=> boss/Database/DatabaseSvc/
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# this part is kept for debugging
#---------------------------------------------------------------------
# PROJECT( db )
# CMAKE_MINIMUM_REQUIRED( VERSION 3.15...3.31 )
# OPTION( BUILTIN_SQLITE
#   "Use included sqlite, or system sqlite library" ON )
# ADD_COMPILE_DEFINITIONS( BEAN )
#---------------------------------------------------------------------

SET( libname DatabaseSvc )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#---------------------------------------------------------------------
# sqlite library

UNSET( SQLIB )
IF( BUILTIN_SQLITE )

  IF( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.c OR
      NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.h )
    MESSAGE( FATAL_ERROR
      "There are no local files \"sqlite3.c\" and \"sqlite3.h\". "
      "Please download the SQLite-amalgamation version of them "
      "from the site:  http://www.sqlite.org  and put in "
      "\"${CMAKE_CURRENT_SOURCE_DIR}\" directory.\n"
      "To do this, you can use the \"get_sqlite3.sh\" script "
      "from the BEAN/scripts folder. "
      "Alternatively, you can try using the SQLITE system library "
      "by setting -DBUILTIN_SQLITE=OFF when call the cmake."
    )
  ENDIF()

  # add sqlite3.c to the list of source files
  SET( CXX_FILES ${CXX_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.c )
  MESSAGE( DEBUG "DEBUG: DatabaseSvc CXX_FILES= ${CXX_FILES}" )

ELSE()

  # try to find SQLITE3 library
  FIND_LIBRARY( SQLIB sqlite3 )
  IF( NOT SQLIB )
    MESSAGE( FATAL_ERROR
      "The system sqlite library was not found. "
      "Try to set BUILTIN_SQLITE to ON and rerun cmake."
    )
  ENDIF()

ENDIF( BUILTIN_SQLITE )

#---------------------------------------------------------------------
# updatedb
IF( ${CMAKE_VERSION} VERSION_LESS "3.12.0" )
  find_package( PythonInterp )
  IF( PYTHONINTERP_FOUND )
    ADD_CUSTOM_TARGET( updatedb
      COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
      COMMAND "${PYTHON_EXECUTABLE}" "./dbupdate2.py" "../dat"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
    )
  ENDIF()
ELSE()
  find_package (Python3 COMPONENTS Interpreter)
  IF( Python3_Interpreter_FOUND )
    ADD_CUSTOM_TARGET( updatedb
      COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
      COMMAND "${Python3_EXECUTABLE}" "./dbupdate.py" "../dat"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
    )
  ELSE()
    IF( Python2_Interpreter_FOUND )
      ADD_CUSTOM_TARGET( updatedb
        COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
        COMMAND "${Python2_EXECUTABLE}" "./dbupdate2.py" "../dat"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
      )
    ENDIF()
  ENDIF()
ENDIF()

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )  # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ENDIF()

# We can not assume that pthread library will be included
# automatically with "root-config"
FIND_PACKAGE( Threads )
IF( SQLIB )
  TARGET_LINK_LIBRARIES( ${libname} ${SQLIB}
    ${CMAKE_THREAD_LIBS_INIT} )
ELSE()
  TARGET_LINK_LIBRARIES( ${libname} ${CMAKE_THREAD_LIBS_INIT} )
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

