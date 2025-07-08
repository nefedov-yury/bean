#---------------------------------------------------------------------
# CMake configuration file to build Bean: see README.md how to use it
#---------------------------------------------------------------------

# message for users of IHEP cluster
IF( CMAKE_VERSION VERSION_LESS "3.15.0" )
  CMAKE_HOST_SYSTEM_INFORMATION(RESULT _host_name QUERY HOSTNAME)
  if( _host_name MATCHES "lxslc7[0-9]+.ihep.ac.cn" )
    MESSAGE("\n+++ Please use cmake3 on IHEP cluster +++\n")
  ENDIF()
ENDIF()

# the project has been tested on versions 3.15...3.29 and there is no
# need to warn about compatibility issues for them
CMAKE_MINIMUM_REQUIRED( VERSION 3.15...3.31 )

# CMAKE_HOST_SYSTEM_NAME is about the system where cmake is being run
# CMAKE_SYSTEM_NAME is about the platform the script is building for
# and is only defined after the project directive
MESSAGE( STATUS
  "+++ BEAN: You are working on ${CMAKE_HOST_SYSTEM_NAME} system\n"
  "   +++ The program was tested under Linux, Darwin and Windows\n"
)

PROJECT( BEAN )

#---------------------------------------------------------------------
# Configuration parameters
#---------------------------------------------------------------------
# installation directory:
IF( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  SET( CMAKE_INSTALL_PREFIX
    ${CMAKE_SOURCE_DIR}/workdir
    CACHE PATH
    "Choose the directory to install program"
    FORCE
  )
ENDIF()
MESSAGE( STATUS
  "+++ BEAN: CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}" )

# Set build type for single-config generators (makefiles)
SET( default_build_type "RelWithDebInfo" )
SET( CMAKE_BUILD_TYPE
  ${default_build_type}
  CACHE STRING
  "Choose the type of build: RelWithDebInfo Release Debug"
  FORCE
)
MESSAGE( STATUS "+++ BEAN: CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}" )

# BOSS_VERSION is the version of the BOSS program used to create the
# DST files you want to analyze. This value defines the following:
# 1) the ROOTEVENTDATA_VERSION
# 2) the preprocessor variable BOSS_VER, which determines
#    the versions of the algorithms in the Analysis folder
# Tested values: 6.6.3-6.6.5;
#                7.0.3;7.0.4;7.0.5-7.0.9
#                7.1.0-2
#
# The string BOSS_VERSION should contain three digits separated by
# dots. Please DO NOT SPECIFY THE PATCH VERSION of BOSS production.
SET( BOSS_VERSION
  6.6.4
  CACHE STRING
  "Version of BOSS in dst production"
)

OPTION( USE_PROOF
  "Building BEAN with the ability to run in PROOF environment" OFF
)
IF( USE_PROOF AND (CMAKE_SYSTEM_NAME MATCHES Windows) )  # Windows
  MESSAGE( VERBOSE "+++ BEAN: PROOF is not available on Windows" )
  SET ( USE_PROOF OFF )
ENDIF()

OPTION( BUILTIN_SQLITE
  "Use sqlite included in BEAN (default) or system library" ON
)

IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  OPTION( ENABLE_CLHEP_XLINK
    "Enable RPATH-based link for CLHEP library" ON )
  OPTION( ENABLE_ROOT_XLINK
    "Enable RPATH-based link for ROOT libraries" ON )
ENDIF()

#---------------------------------------------------------------------
# Set ROOTEVENTDATA_VERSION and BOSS_VER
#---------------------------------------------------------------------
STRING( REGEX REPLACE "^([0-9]+)[.][0-9]+[.][0-9]+.*$" "\\1"
  BOSS_VER_MAJOR "${BOSS_VERSION}"
)
STRING( REGEX REPLACE "^[0-9][.]([0-9]+)[.][0-9]+.*$" "\\1"
  BOSS_VER_MINOR "${BOSS_VERSION}"
)
STRING( REGEX REPLACE "^[0-9]+[.][0-9]+[.]([0-9]+).*$" "\\1"
  BOSS_VER_MICRO "${BOSS_VERSION}"
)
SET( BOSS_VER "${BOSS_VER_MAJOR}${BOSS_VER_MINOR}${BOSS_VER_MICRO}" )
MESSAGE( STATUS "+++ BEAN: BOSS_VER is ${BOSS_VER}" )

# ROOTEVENTDATA_VERSION defines format of dst-file:
#   BOSS_VERSION                        ROOTEVENTDATA_VERSION
#   6.5.2; 6.5.3; 6.5.5; 6.6.1          = BOSS_VERSION
#   6.6.2                               6.6.1
#   >= 6.6.3 < 7.0.5                    6.6.3
#   >= 7.0.5 <= 7.0.9                   7.0.5
#   7.1.0, 7.1.1, 7.1.2                 7.1.0
#
SET( ROOTEVENTDATA_VERSION ${BOSS_VERSION} )
IF( BOSS_VERSION  VERSION_EQUAL  "6.6.2" )
  SET( ROOTEVENTDATA_VERSION "6.6.1" )
ELSEIF( BOSS_VERSION  VERSION_GREATER  "6.6.2" AND
    BOSS_VERSION  VERSION_LESS  "7.0.5" )
  SET( ROOTEVENTDATA_VERSION "6.6.3" )
ELSEIF( BOSS_VERSION  VERSION_LESS  "7.0.10" )
  SET( ROOTEVENTDATA_VERSION "7.0.5" )
ELSE()
  SET( ROOTEVENTDATA_VERSION "7.1.0" )
ENDIF()
MESSAGE( DEBUG
  "DEBUG: ROOTEVENTDATA_VERSION is ${ROOTEVENTDATA_VERSION}" )

SET( ROOTEVENTDATA RootEventData_${ROOTEVENTDATA_VERSION} )
IF( NOT EXISTS "${CMAKE_SOURCE_DIR}/${ROOTEVENTDATA}" )
  MESSAGE( FATAL_ERROR
    "'${ROOTEVENTDATA}' does not exist in folder:\n"
    "'${CMAKE_SOURCE_DIR}'\n"
    "The BOSS_VERSION=${BOSS_VERSION} is not supported or "
    "is specified incorrectly.\n"
  )
ENDIF()
MESSAGE( STATUS "+++ BEAN: ROOTEVENTDATA is '${ROOTEVENTDATA}'" )

#---------------------------------------------------------------------
# Searching for ROOT
#---------------------------------------------------------------------
SET( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake )

SET( ROOT_CONFIG_SEARCHPATH
  ""
  CACHE PATH
  "Path to directory with root-config program"
)
INCLUDE( FindROOT )

IF( ROOTVERSION VERSION_GREATER 6.31 )
  IF( USE_PROOF )
    SET ( USE_PROOF OFF )
    MESSAGE( STATUS
      "+++ BEAN: The root with version ${ROOTVERSION} > 6.31 "
      "no longer supports the PROOF framework" )
    MESSAGE( STATUS
      "+++ BEAN: USE_PROOF option is now disabled" )
  ENDIF()
ENDIF()

#---------------------------------------------------------------------
# Searching for CLHEP
#---------------------------------------------------------------------
SET( CLHEP_SEARCHPATH
  "$ENV{CLHEP_DIR}"
  CACHE PATH
  "Path to directory where clhep is installed"
)
INCLUDE( FindCLHEP )

#---------------------------------------------------------------------
# Compilation options
# We need a compiler that supports at least C++11, so we retrieve the
# corresponding settings from ROOT_CFLAGS. In most cases this is
# enough.
#---------------------------------------------------------------------
# CMAKE_COMPILER_IS_GNUCC is depricated
IF( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR # g++
  (CMAKE_CXX_COMPILER_ID MATCHES "Clang") )   # Clang or AppleClang

  IF( ROOT_CFLAGS MATCHES "-std=[^ ]*" )
    # cut out a section with -std=...
    STRING( REGEX REPLACE "^.*(-std=[^ ]*).*$" "\\1"
      ROOT_STDCXX "${ROOT_CFLAGS}" )
    MESSAGE( DEBUG "DEBUG: ROOT_STDCXX= ${ROOT_STDCXX}" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_STDCXX}" )
  ELSE()
    # 'old' root versions
    IF( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
      MESSAGE( DEBUG
        "DEBUG: g++ version is ${CMAKE_CXX_COMPILER_VERSION}" )
      # c++14 mode is the default in GCC 6.1 up 10 (including)
      # c++17 mode is the default in GCC >= 11
      IF( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1 )
        IF( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.7 )
          SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
        ELSEIF( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.3 )
          SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )
        ELSE()
          MESSAGE( FATAL_ERROR
            "Insufficient g++ version ${CMAKE_CXX_COMPILER_VERSION}\n"
            "BEAN uses C++11 features"
          )
        ENDIF()
      ENDIF()
    ELSEIF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
    ENDIF()
  ENDIF()

  # to see compilation warnings
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )
ENDIF()

IF( CMAKE_SYSTEM_NAME MATCHES Windows )         # Windows
  IF( ROOT_CFLAGS MATCHES "std:[^ ]*" )
    STRING( REGEX REPLACE "^.*std:([^ ]*).*$" "\\1"
      ROOT_STDCXX "${ROOT_CFLAGS}" )
    SET( ROOT_STDCXX "-Zc:__cplusplus -std:${ROOT_STDCXX}" )
    MESSAGE( DEBUG "DEBUG: ROOT_STDCXX= ${ROOT_STDCXX}" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_STDCXX}" )
    IF( ROOT_STDCXX  MATCHES "17" )
      # std::byte is enabled by /std:c++17, conflict with SDK
      ADD_COMPILE_DEFINITIONS( _HAS_STD_BYTE=0 )
    ENDIF()
  ELSE()
    SET( CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Zc:__cplusplus -std:c++11" )
  ENDIF()

  # use run-time library (MD or MDd)
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -MD -GR -EHs" )
  # to see compilation warnings
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W3" )

  # Math Constants like M_PI...
  ADD_COMPILE_DEFINITIONS( _USE_MATH_DEFINES )
  # To eliminate warnings for the older, `less` secure functions
  # like ctime, fopen, sscanf, strcpy and so on
  ADD_COMPILE_DEFINITIONS( _CRT_SECURE_NO_WARNINGS )
ENDIF()

# Define system-specific compiler macros:
# __unix__  for Linux
# __APPLE__ for Darwin
# _WIN32    for Windows (32 or 64 bit)
IF( CMAKE_SYSTEM_NAME MATCHES Linux )           # Linux
  ADD_COMPILE_DEFINITIONS( __unix__ )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Darwin )      # MacOS
  ADD_COMPILE_DEFINITIONS( __APPLE__ )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )     # Windows
  ADD_COMPILE_DEFINITIONS( _WIN32 )
ENDIF()

# Preprocessor definitions specific for the BEAN
ADD_COMPILE_DEFINITIONS( BEAN )
ADD_COMPILE_DEFINITIONS( BOSS_VER=${BOSS_VER} )
IF( USE_PROOF )
  ADD_COMPILE_DEFINITIONS( USE_PROOF=1 )
ELSEIF()
  # this is equivalent to leaving the variable undefined
  ADD_COMPILE_DEFINITIONS( USE_PROOF=0 )
ENDIF()

# print C++ compiler flags and preprocessor definitions
MESSAGE( STATUS "+++ BEAN: CMAKE_CXX_FLAGS are '${CMAKE_CXX_FLAGS}'" )
GET_DIRECTORY_PROPERTY( PREDEFS COMPILE_DEFINITIONS )
IF ( PREDEFS  )
  MESSAGE( STATUS
    "+++ BEAN: Preprocessor definitions are '${PREDEFS}'" )
ENDIF()

#---------------------------------------------------------------------
# Libraries options
#---------------------------------------------------------------------
# Set the common place where all libraries to be built.
# Please do not to be confused it with the place where the libraries
# will be installed later.
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

# Set actual place where the libraries will be installed
SET( install_lib_dir "BeanLib_${BOSS_VERSION}" )

#---------------------------------------------------------------------
# Shared library settings: RPATH
#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Darwin ) # MacOS
  # initialize the MACOSX_RPATH property on all targets
  SET ( CMAKE_MACOSX_RPATH ON )
ENDIF()

IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS

  IF( ENABLE_ROOT_XLINK )
    IF( ROOT_LIBRARIES MATCHES ",-rpath," )
      MESSAGE( VERBOSE
        "+++ BEAN: the ROOT_LIBRARIES variable already contains "
        "the rpath setting, so the ENABLE_ROOT_XLINK parameter "
        "is redundant"
      )
    ELSE()
      SET( ROOTXLINK "-Wl,-rpath,${ROOT_LIBRARY_DIR}" )
      SET( CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} ${ROOTXLINK}" )
    ENDIF()
  ENDIF()

  IF( ENABLE_CLHEP_XLINK )
    SET( CLHEPXLINK "-Wl,-rpath,${CLHEP_LIBRARY_DIR}" )
    SET( CMAKE_SHARED_LINKER_FLAGS
      "${CMAKE_SHARED_LINKER_FLAGS} ${CLHEPXLINK}" )
  ENDIF()

  SET( CMAKE_SKIP_BUILD_RPATH FALSE )

  # append directories in the linker search path
  SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

  SET( CMAKE_INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${install_lib_dir}" )
ENDIF()

#---------------------------------------------------------------------
# build libraries
#---------------------------------------------------------------------
ADD_SUBDIRECTORY( ${ROOTEVENTDATA} )
SET( ROOTEVENTDATA_DIR  ${CMAKE_SOURCE_DIR}/${ROOTEVENTDATA} )

ADD_SUBDIRECTORY( BeanCore )
SET( BeanCore_DIR  ${CMAKE_SOURCE_DIR}/BeanCore )

# -- Analysis --
IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  ADD_SUBDIRECTORY( Analysis/Identifier )
  SET( Identifier_DIR  ${CMAKE_SOURCE_DIR}/Analysis/Identifier )
ENDIF()

ADD_SUBDIRECTORY( Analysis/AbsCor )
SET( AbsCor_DIR  ${CMAKE_SOURCE_DIR}/Analysis/AbsCor )

ADD_SUBDIRECTORY( Analysis/DatabaseSvc )
SET( DatabaseSvc_DIR  ${CMAKE_SOURCE_DIR}/Analysis/DatabaseSvc )

ADD_SUBDIRECTORY( Analysis/MagneticField )
SET( MagneticField_DIR  ${CMAKE_SOURCE_DIR}/Analysis/MagneticField )

ADD_SUBDIRECTORY( Analysis/EventTag )
SET( EventTag_DIR  ${CMAKE_SOURCE_DIR}/Analysis/EventTag )

IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  ADD_SUBDIRECTORY( Analysis/ParticleID )
  SET( ParticleID_DIR  ${CMAKE_SOURCE_DIR}/Analysis/ParticleID )
ELSE()
  ADD_SUBDIRECTORY( Analysis/ParticleID_boss6 )
  SET( ParticleID_DIR  ${CMAKE_SOURCE_DIR}/Analysis/ParticleID_boss6 )
ENDIF()

ADD_SUBDIRECTORY( Analysis/VertexFit )
SET( VertexFit_DIR  ${CMAKE_SOURCE_DIR}/Analysis/VertexFit )

ADD_SUBDIRECTORY( Analysis/RscanDQ )
SET( RscanDQ_DIR  ${CMAKE_SOURCE_DIR}/Analysis/RscanDQ )

ADD_SUBDIRECTORY( Analysis/TrackCorrection )
SET(TrackCorrection_DIR ${CMAKE_SOURCE_DIR}/Analysis/TrackCorrection )

# -- User --
ADD_SUBDIRECTORY( BeanUser )

#---------------------------------------------------------------------
# BEAN executable
#---------------------------------------------------------------------
# program name
SET( prog bean_${BOSS_VERSION}.exe )
SET( uninstall_prog "${prog}"  )
IF( CMAKE_SYSTEM_NAME MATCHES Windows )
  # in windows the `exe` extention is added by default
  SET( prog bean_${BOSS_VERSION} )
ENDIF()

INCLUDE_DIRECTORIES(
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CLHEP_INCLUDE_DIR}
  ${DatabaseSvc_DIR}
)

SET( MAIN_SRCS "BeanCore/main.cxx" )
SET_SOURCE_FILES_PROPERTIES( "BeanCore/main.cxx"
  PROPERTIES COMPILE_DEFINITIONS BEANBASE="${CMAKE_SOURCE_DIR}" )

# getopt() port for MSVC: https://www.ludvikjerabek.com/software/
IF( CMAKE_SYSTEM_NAME MATCHES Windows )
  LIST( APPEND MAIN_SRCS "BeanCore/win_getopt/getopt.c" )
  ADD_COMPILE_DEFINITIONS( STATIC_GETOPT )
ENDIF()

ADD_EXECUTABLE( ${prog} ${MAIN_SRCS} )

IF( CMAKE_SYSTEM_NAME MATCHES Linux )
  # this linker option only affects ELF files: man ld /--no-as-needed
  SET( CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed" )
ENDIF()

TARGET_LINK_LIBRARIES( ${prog}
  RootEventData
  BeanCore
  BeanUser
  AbsCor
  DatabaseSvc
  MagneticField
  EventTag
  ParticleID
  VertexFit
)

IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  SET_TARGET_PROPERTIES( ${prog}
    PROPERTIES  LINK_FLAGS  "${ROOTXLINK} ${CLHEPXLINK}"
  )
ENDIF()

# final message (must be simple!) after compilation is done
ADD_CUSTOM_COMMAND( TARGET ${prog} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red --bold
  "Use   make install   to install the program!"
  COMMENT "Compilation is over!"
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( DIRECTORY DESTINATION ${install_lib_dir} )
INSTALL( TARGETS ${prog}  DESTINATION ${CMAKE_INSTALL_PREFIX} )

ADD_CUSTOM_TARGET( install_bean
  "${CMAKE_COMMAND}" -P
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
)

# generate 'setup.sh' and install it if asked
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  SET( LD_LIBRARY_PATH "LD_LIBRARY_PATH" )
  SET( DYLD_LIBRARY_PATH "DYLD_LIBRARY_PATH" )

  CONFIGURE_FILE(
    ${CMAKE_SOURCE_DIR}/cmake/setup.sh.cmake
    ${PROJECT_BINARY_DIR}/setup_${BOSS_VERSION}.sh
    @ONLY
  )
  ADD_CUSTOM_TARGET( setup_file
    COMMAND "${CMAKE_COMMAND}" -E copy
    ${PROJECT_BINARY_DIR}/setup_${BOSS_VERSION}.sh
    ${CMAKE_INSTALL_PREFIX}
  )
ENDIF()

#---------------------------------------------------------------------
# PROOF
#---------------------------------------------------------------------
IF ( USE_PROOF )
  # Note: the name of par-file MUST match the name of the directory

  # PROOF_A_FILES is the list of data files that could be used
  # in the BEAN from the directory "Analysis"
  FILE( GLOB PROOF_A_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/PROOF-INF/*.*" )
  FILE( GLOB TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/DatabaseSvc/dat/*.db" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  FILE( GLOB TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/MagneticField/dat/*.dat" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  FILE( GLOB TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/MagneticField/dat/*.txt" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  FILE( GLOB_RECURSE TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/ParticleID/share/*.txt" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  FILE( GLOB TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/AbsCor/dat/*.txt" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  FILE( GLOB TMP_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "Analysis/EventTag/share/*.codes" )
  LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
  LIST( APPEND PROOF_A_FILES "Analysis/EventTag/share/pdt.table" )

  MESSAGE( DEBUG "DEBUG: PROOF_A_FILES are\n${PROOF_A_FILES}\n" )

  SET( ANALYSIS_PAR "${CMAKE_INSTALL_PREFIX}/par/Analysis.par" )
  ADD_CUSTOM_COMMAND(
    OUTPUT   ${ANALYSIS_PAR}
    COMMAND "${CMAKE_COMMAND}" -E
    make_directory ${CMAKE_INSTALL_PREFIX}/par
    COMMAND "${CMAKE_COMMAND}" -E
    tar -czf ${ANALYSIS_PAR} ${PROOF_A_FILES}
    DEPENDS  ${PROOF_A_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Creating ${ANALYSIS_PAR}"
  )

  # PROOF_B_FILES is the list of data files that could be used
  # in the BEAN from the directory "BeanUser"
  FILE( GLOB PROOF_B_FILES RELATIVE "${CMAKE_SOURCE_DIR}"
    "BeanUser/PROOF-INF/*.*" )
  LIST( APPEND PROOF_B_FILES "BeanUser/mydecay.codes" )

  SET( BEANUSER_PAR "${CMAKE_INSTALL_PREFIX}/par/BeanUser.par" )
  ADD_CUSTOM_COMMAND(
    OUTPUT   ${BEANUSER_PAR}
    COMMAND "${CMAKE_COMMAND}" -E
    make_directory ${CMAKE_INSTALL_PREFIX}/par
    COMMAND "${CMAKE_COMMAND}" -E
    tar -czf ${BEANUSER_PAR} ${PROOF_B_FILES}
    DEPENDS  ${PROOF_B_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Creating ${BEANUSER_PAR}"
  )

  SET( LIB_PAR "${CMAKE_INSTALL_PREFIX}/par/${install_lib_dir}.par" )
  SET( PROOF_INF_LIB
    "${CMAKE_INSTALL_PREFIX}/${install_lib_dir}/PROOF-INF" )
  ADD_CUSTOM_COMMAND(
    OUTPUT   ${PROOF_INF_LIB}/BUILD.sh ${PROOF_INF_LIB}/SETUP.C
    COMMAND "${CMAKE_COMMAND}" -E make_directory ${PROOF_INF_LIB}
    COMMAND "${CMAKE_COMMAND}" -E
    copy BUILD_lib.sh ${PROOF_INF_LIB}/BUILD.sh
    COMMAND "${CMAKE_COMMAND}" -E
    copy SETUP_lib.C ${PROOF_INF_LIB}/SETUP.C
    DEPENDS ${CMAKE_SOURCE_DIR}/cmake/BUILD_lib.sh
    ${CMAKE_SOURCE_DIR}/cmake/SETUP_lib.C
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/cmake"
    COMMENT "Creating ${PROOF_INF_LIB}"
  )

  ADD_CUSTOM_COMMAND(
    OUTPUT   ${LIB_PAR}
    COMMAND "${CMAKE_COMMAND}" -E
    make_directory ${CMAKE_INSTALL_PREFIX}/par
    COMMAND "${CMAKE_COMMAND}" -E
    tar -czf ${LIB_PAR} ${install_lib_dir}
    DEPENDS ${PROOF_INF_LIB}/BUILD.sh ${PROOF_INF_LIB}/SETUP.C
    WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
    COMMENT "Creating ${LIB_PAR}"
  )

  ADD_CUSTOM_TARGET( proofbean
    DEPENDS ${ANALYSIS_PAR} ${BEANUSER_PAR} ${LIB_PAR}
  )
  ADD_DEPENDENCIES( proofbean install_bean )
ELSE( )
  MESSAGE( DEBUG "DEBUG: USE_PROOF is ${USE_PROOF}" )
ENDIF( USE_PROOF )

#---------------------------------------------------------------------
# uninstall or delete "by hand"
#---------------------------------------------------------------------
ADD_CUSTOM_TARGET( uninstall
  COMMAND "${CMAKE_COMMAND}" -E remove -f ${uninstall_prog}
  COMMAND "${CMAKE_COMMAND}" -E remove_directory ${install_lib_dir}
  COMMAND "${CMAKE_COMMAND}" -E remove -f setup_${BOSS_VERSION}.sh
  COMMAND "${CMAKE_COMMAND}" -E remove_directory par
  WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
)
