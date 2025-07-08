#---------------------------------------------------------------------
# BeanCore: library provides the basic functionality of the bean
#           and contains transitional classes Root->Analysis
#---------------------------------------------------------------------

SET( libname BeanCore )

SET( INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CLHEP_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} )

MESSAGE( DEBUG
  "DEBUG: BeanCore INCLUDES are ${INCLUDE_DIRECTORIES}" )

#---------------------------------------------------------------------
# the list of headers are needed to generate the root dictionary
SET( LINKDEFHDRS  DstFormat.h ReadDst.h Bean.h TMergeableMap.h )
#  IF( USE_PROOF )
  #  SET( LINKDEFHDRS  ${LINKDEFHDRS} MasterInfoHack.cxx )
#  ENDIF()

IF( ROOT_MAJOR_VERS LESS 6 )
  # the name of the dictionary
  SET( DICTIONARY ${libname}_rootcint.cxx )
  ROOT5_GENERATE_DICTIONARY(
    ${DICTIONARY}
    "${LINKDEFHDRS}"
    "${INCLUDE_DIRECTORIES}"
    ""
  )
ELSE()
  ROOT6_GENERATE_DICTIONARY(
    ${libname}
    "${LINKDEFHDRS}"
    "${INCLUDE_DIRECTORIES}"
    ""
  )
  SET( DICTIONARY ${libname}_rdict.cxx )
  SET( DictNamePcm ${libname}_rdict.pcm )
ENDIF()

#---------------------------------------------------------------------
# the list of source files
SET( CXX_FILES  DstFormat.cxx ReadDst.cxx Bean.cxx
  TMergeableMap.cxx DstEvtRecTracks.cxx TofHitStatus.cxx
  BMdcTrack.cxx BMdcKalTrack.cxx BEmcTrack.cxx BExtTrack.cxx )

IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # Windows
  SET_SOURCE_FILES_PROPERTIES( "Bean.cxx"
    PROPERTIES COMPILE_DEFINITIONS BEANLIBDIR="${install_lib_dir}" )
ENDIF()

IF( USE_PROOF )
  #  SET( CXX_FILES  ${CXX_FILES} MasterInfoHack.cxx )
  SET( CXX_FILES  ${CXX_FILES} unix_cp.cxx )
ENDIF()

# add the dictionary to the list of source files
LIST( APPEND CXX_FILES ${DICTIONARY} )

#---------------------------------------------------------------------
IF( CMAKE_SYSTEM_NAME MATCHES Linux|Darwin ) # Linux or MacOS
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  IF( USE_PROOF )
    SET( BeanCore_ROOT_LIBRARIES
      "-L ${ROOT_LIBRARY_DIR} -lEG -lProofPlayer -lProof" )
  ELSE()
    SET( BeanCore_ROOT_LIBRARIES
      "-L ${ROOT_LIBRARY_DIR} -lEG" )
  ENDIF()
ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows ) # Windows
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
  SET( BeanCore_ROOT_LIBRARIES "${ROOT_LIBRARIES}" )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  ${BeanCore_ROOT_LIBRARIES}
  ${CLHEP_LIBRARIES}
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

IF( ROOT_MAJOR_VERS GREATER 5 )
  INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/${DictNamePcm}
    DESTINATION ${install_lib_dir} )
ENDIF()
