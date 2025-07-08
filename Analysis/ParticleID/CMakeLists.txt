#---------------------------------------------------------------------
# ParticleID: library for particle identification
#              <=> boss/Analysis/ParticleID for BOSS >= 7
#---------------------------------------------------------------------

SET( libname ParticleID )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${CLHEP_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
)

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  #  SET( ParticleID_ROOT_LIBRARIES "${ROOT_LIBRARIES} -lMLP" )
  SET( ParticleID_ROOT_LIBRARIES "-L ${ROOT_LIBRARY_DIR} -lMLP" )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )  # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
  SET( ParticleID_ROOT_LIBRARIES "${ROOT_LIBRARIES}" )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
  ${ParticleID_ROOT_LIBRARIES}
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

