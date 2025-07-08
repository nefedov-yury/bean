#---------------------------------------------------------------------
# BeanUser is the place for user analysis and examples
#---------------------------------------------------------------------

# Explicit list of files that will be compiled
# --------------------------------------------
SET( CXX_FILES
  # ++ Several functions with examples and tests
  User1.cxx
  UserTest.cxx
  TestAbsCor.cxx
  TestPID.cxx
  MagField.cxx
  TestEventTag.cxx
  TestDb.cxx
  JobInfo.cxx
  EntryList.cxx
  #
  # ++ Rhopi.cxx is the 'Jpsi -> rho pi' analysis rewritten for BEAN
  #    using the Boss cobe 'Analysis/Physics/RhopiAlg' as a template
  Rhopi.cxx
  #
  # ++ Here are possibly useful old files
  #     Etaetagamma.cxx
  #     Semileptonic.cxx
  #     DoubleDTag.cxx
  #     Bhabha.cxx
  #     RadBhabha.cxx
  #
  #     SelectOmegaEta.cxx
  #     SelectOmegapi0.cxx
  #     Eta_TwoPi.cxx
  #     PipPimPi0Eta.cxx
  #     SelectJpsiX.cxx
  #
  # ++ User programs are here:
  SelectKKgg.cxx
  PsipJpsiPhiEta.cxx
  PsipPiPiKK.cxx
  PsipJpsiGammaEta.cxx
  PsipHelixParCor.cxx
  DelClonedTrks.cxx
)

#---------------------------------------------------------------------
SET( libname BeanUser )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${ParticleID_DIR}
  ${MagneticField_DIR}
  ${VertexFit_DIR}
  ${AbsCor_DIR}
  ${DatabaseSvc_DIR}
  ${EventTag_DIR}
  ${RscanDQ_DIR}
  ${TrackCorrection_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CLHEP_INCLUDE_DIR}
)

ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
  ParticleID
  MagneticField
  VertexFit
  AbsCor
  DatabaseSvc
  EventTag
  RscanDQ
  TrackCorrection
)

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

