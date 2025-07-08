#---------------------------------------------------------------------
# TrackCorrection:
#       this class based on "TrackCorrection" package modified for use
#       in a BEAN (see ~guoyp/public/TrackCorrection)
#
#---------------------------------------------------------------------

SET( libname TrackCorrection )

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
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )  # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

