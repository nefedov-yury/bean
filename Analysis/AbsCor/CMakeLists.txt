#---------------------------------------------------------------------
# This is rewrote version of AbsCor from Boss:
#                               <- boss/Analysis/PhotonCor/AbsCor
#---------------------------------------------------------------------

SET( libname AbsCor )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${CLHEP_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
)
IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  INCLUDE_DIRECTORIES( ${Identifier_DIR} )
ENDIF()

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )  # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
)

IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  TARGET_LINK_LIBRARIES( ${libname} Identifier)
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

