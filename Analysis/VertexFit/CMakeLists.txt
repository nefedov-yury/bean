#---------------------------------------------------------------------
# VertexFit: <=> boss/Analysis/VertexFit
#---------------------------------------------------------------------

SET( libname VertexFit )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${MagneticField_DIR}
  ${DatabaseSvc_DIR}
  ${CLHEP_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
)

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

# the file VertexParameter.cxx MUST be in the end of list
# to prevent MSVC Linker Tools Warning LNK4221
FILE( GLOB LAST_FILE src/VertexParameter.cxx )
LIST( REMOVE_ITEM CXX_FILES ${LAST_FILE} )
LIST( APPEND CXX_FILES ${LAST_FILE} )
MESSAGE( DEBUG "DEBUG: VertexFit src files are ${CXX_FILES}" )

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )  # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
  MagneticField
  DatabaseSvc
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

