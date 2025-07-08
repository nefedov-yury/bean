#---------------------------------------------------------------------
# MagneticField: interface library to magnetic field database
#              <=> boss/MagneticField/
#---------------------------------------------------------------------

SET( libname MagneticField )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${DatabaseSvc_DIR}
  ${CLHEP_INCLUDE_DIR}
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
  DatabaseSvc
  ${CLHEP_LIBRARIES}
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

