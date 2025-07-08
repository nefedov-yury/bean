#---------------------------------------------------------------------
# RootEventData: library to read DST-files (root-files)
#---------------------------------------------------------------------

SET( libname RootEventData )

SET( INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOT_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} )
MESSAGE( DEBUG
  "DEBUG: RootEventData INCLUDES are ${INCLUDE_DIRECTORIES}" )

#---------------------------------------------------------------------
# the list of headers are needed to generate the root dictionary
FILE( GLOB H_FILES "RootEventData/*.h" )
FOREACH( file ${H_FILES} )
  IF( NOT file MATCHES LinkDef )
    SET( LINKDEFHDRS ${LINKDEFHDRS} ${file} )
  ELSE()
    MESSAGE( DEBUG "DEBUG: RootEventData filter out file ${file}" )
  ENDIF()
ENDFOREACH()

#---------------------------------------------------------------------
IF( ROOT_MAJOR_VERS LESS 6 )
  # the name of the dictionary
  SET( DICTIONARY ${libname}_rootcint.cxx )
  ROOT5_GENERATE_DICTIONARY(
    ${DICTIONARY}
    "${LINKDEFHDRS}"
    "${INCLUDE_DIRECTORIES}"
    RootEventData/LinkDef.h
  )
ELSE()
  ROOT6_GENERATE_DICTIONARY(
    ${libname}
    "${LINKDEFHDRS}"
    "${INCLUDE_DIRECTORIES}"
    RootEventData/LinkDef.h
  )
  SET( DICTIONARY ${libname}_rdict.cxx )
  SET( DictNamePcm ${libname}_rdict.pcm )
ENDIF()

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

# add the dictionary to the list of source files
LIST( APPEND CXX_FILES ${DICTIONARY} )

MESSAGE( DEBUG "DEBUG: RootEventData src files are ${CXX_FILES}\n" )

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  TARGET_LINK_LIBRARIES( ${libname}
    ${ROOT_LIBRARIES}
  )
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows ) # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
  TARGET_LINK_LIBRARIES( ${libname}
    ${ROOT_LIBRARIES}
  )
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

IF( ROOT_MAJOR_VERS GREATER 5 )
  INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/${DictNamePcm}
    DESTINATION ${install_lib_dir} )
ENDIF()

