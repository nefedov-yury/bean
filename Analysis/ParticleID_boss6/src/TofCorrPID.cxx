#include <cmath>

#include "ParticleID/TofCorrPID.h"

#ifndef BEAN
#include "EvtRecEvent/EvtRecTrack.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "ExtEvent/RecExtTrack.h"
#include "TofRecEvent/RecTofTrack.h"
#include "DstEvent/TofHitStatus.h"
#else
#include "TofHitStatus.h"
#endif
#include <fstream>


TofCorrPID * TofCorrPID::m_pointer = 0;
TofCorrPID * TofCorrPID::instance() {
   if(!m_pointer) m_pointer = new TofCorrPID();
   return m_pointer;
}

TofCorrPID::TofCorrPID():ParticleIDBase() {
   m_runBegin = 0;
   m_runEnd   = 0;
}

void TofCorrPID::init() {

   for(int i = 0; i < 5; i++) {
      m_chi[i]    = -100.0;
      m_prob[i]   = -1.0;
      m_sigma[i]  = 10.0;
      m_offset[i] = -1000.0;
   }
   m_chimin = 99.;
   m_pdfmin = 99.;
   m_ndof = 0;

   for( unsigned int i=0; i<5; i++ ) {
     for( unsigned int j=0; j<7; j++ ) {
       m_deltaT[j][i] = -1000.0;
     }
   }

   int run = getRunNo();
   if( !( abs(run)>=m_runBegin && abs(run)<=m_runEnd ) ) {
     inputParameter( getRunNo() );
   }

   return;
}


void TofCorrPID::calculate() {
   if(particleIDCalculation() == 0) m_ndof=1;
}


int TofCorrPID::particleIDCalculation() {
   int irc=-1;

   EvtRecTrack* recTrk = PidTrk();
   if(!(recTrk->isMdcTrackValid())) return irc;
   if(!(recTrk->isMdcKalTrackValid())) return irc;
   if(!(recTrk->isExtTrackValid())) return irc;
   if(!(recTrk->isTofTrackValid())) return irc;

#ifndef BEAN
   SmartRefVector<RecTofTrack> tofTrk = recTrk->tofTrack();
   SmartRefVector<RecTofTrack>::iterator it;
#else
   const std::vector<RecTofTrack* >& tofTrk = recTrk->tofTrack();
   std::vector<RecTofTrack* >::const_iterator it;
#endif

   RecMdcTrack* mdcTrk = recTrk->mdcTrack();
   int charge = mdcTrk->charge();

   double p[5], betaGamma[5];
   RecMdcKalTrack* kalTrk = recTrk->mdcKalTrack();
   for( int i=0; i<5; i++ ) {
     if( i==0 ) {
       kalTrk->setPidType(RecMdcKalTrack::electron);
     }
     else if( i==1 ) {
       kalTrk->setPidType(RecMdcKalTrack::muon);
     }
     else if( i==2 ) {
       kalTrk->setPidType(RecMdcKalTrack::pion);
     }
     else if( i==3 ) {
       kalTrk->setPidType(RecMdcKalTrack::kaon);
     }
     else if( i==4 ) {
       kalTrk->setPidType(RecMdcKalTrack::proton);
     }
#ifndef BEAN
     HepLorentzVector ptrk = kalTrk->p4();
#else
     HepLorentzVector ptrk = kalTrk->p4( xmass(i) );
#endif
     p[i] = ptrk.rho();
     betaGamma[i] = p[i]/xmass(i);
   }

   double zrhit[2];
   RecExtTrack* extTrk = recTrk->extTrack();
   zrhit[0] = extTrk->tof1Position().z();
   zrhit[1] = extTrk->tof2Position().z();

   int tofid[2] = { -9, -9 };

   bool readFile = false;
   bool signal[2] = { false, false };
   TofHitStatus *hitst = new TofHitStatus;
   for( it = tofTrk.begin(); it!=tofTrk.end(); it++ ) {
     unsigned int st = (*it)->status();
     hitst->setStatus(st);
     if( hitst->is_raw() ) return irc;  // TOF no hit
     bool barrel  = hitst->is_barrel();
     if( !barrel ) { zrhit[0] = extTrk->tof1Position().rho(); }
     bool readout = hitst->is_readout();
     bool counter = hitst->is_counter();
     bool cluster = hitst->is_cluster();
     int  layer   = hitst->layer();
     tofid[layer-1] = (*it)->tofID();
     bool east      = hitst->is_east();
     double tof     = (*it)->tof();
     if( readout ) {
       // default 0: end cap
       unsigned int ipmt = 0;
       // barrel: 0: inner east, 1:inner west, 2:outer east, 3: outer west
       // endcap: 0 
       if( barrel ) { ipmt = ( ( st & 0xC0 ) >> 5 ) + ( ( ( st ^ 0x20 ) & 0x20 ) >> 5 ) - 2; }
       for( unsigned int i=0; i<5; i++ ) {
	 double offset = (*it)->toffset(i);
	 double texp   = (*it)->texp(i);
	 if( texp<0.0 ) continue;
	 double dt     = tof - offset - texp;
	 m_deltaT[ipmt][i] = offsetTof( i, barrel, ipmt, betaGamma[i], charge, zrhit[layer-1], dt );
       }
       if( counter && cluster ) {
	 for( unsigned int i=0; i<5; i++ ) {
	   if( ((*it)->texp(i))>0.0 ) {
	     m_offset[i] = m_deltaT[ipmt][i];
	     m_sigma[i]  = sigmaTof( i, charge, barrel, ipmt, &tofid[0], &zrhit[0], betaGamma[i] );
	   }
	 }
       }
     }
     else {
       if( counter ) {
	 for( unsigned int i=0; i<5; i++ ) {
	   double offset = (*it)->toffset(i);
	   double texp   = (*it)->texp(i);
	   if( texp<0.0 ) continue;
	   double dt     = tof - offset - texp;
	   m_deltaT[layer+3][i] = offsetTof( i, tofid[layer-1], zrhit[layer-1], dt );
	   m_offset[i] = m_deltaT[layer+3][i];
	   m_sigma[i]  = sigmaTof( i, charge, barrel, layer+3, &tofid[0], &zrhit[0], betaGamma[i] );
	 }
	 signal[layer-1] = correlationCheck();
       }
       else {
	 if( cluster ) {
	   if( signal[0] && signal[1] ) {
	     for( unsigned int i=0; i<5; i++ ) {
	       double offset = (*it)->toffset(i);
	       double texp   = (*it)->texp(i);
	       if( texp<0.0 ) continue;
	       double dt     = tof - offset - texp;
	       m_deltaT[6][i] = offsetTof( i, tofid[0], tofid[1], zrhit[0], zrhit[1], dt );
	       m_offset[i] = m_deltaT[6][i];
	       m_sigma[i]  = sigmaTof( i, charge, barrel, 6, &tofid[0], &zrhit[0], betaGamma[i] );
	     }
	   }
	   else if( signal[0] && !signal[1] ) {
	     for( unsigned int i=0; i<5; i++ ) {
	       m_offset[i] = m_deltaT[4][i];
	       m_sigma[i]  = sigmaTof( i, charge, barrel, 4, &tofid[0], &zrhit[0], betaGamma[i] );
	     }
	   }
	   else if( !signal[0] && signal[1] ) {
	     for( unsigned int i=0; i<5; i++ ) {
	       m_offset[i] = m_deltaT[5][i];
	       m_sigma[i]  = sigmaTof( i, charge, barrel, 5, &tofid[1], &zrhit[1], betaGamma[i] );
	     }
	   }
	   else return irc;
	 }
       }
     }
   }

   bool good = false;
   double pdftemp = 0;
   for( unsigned int i=0; i<5; i++ ) {
     m_chi[i] = m_offset[i]/m_sigma[i];
     if( ( m_chi[i]>-4.0 && m_chi[i]<6.0 ) && !good ) { good = true; }
     double ppp = pdfCalculate(m_chi[i],1);
     if( fabs(ppp) > pdftemp) { pdftemp = fabs(ppp); }
   }
   m_pdfmin = pdftemp;
   if( pdftemp < pdfCalculate(pdfMinSigmaCut(),1) ) return irc;
   if( !good ) return irc;
   for(int i = 0; i < 5; i++) {
      m_prob[i] = probCalculate(m_chi[i]*m_chi[i], 1);
   }
   irc = 0;
   return irc;
}


void TofCorrPID::inputParameter( int run ) {

  std::string filePath = path + "/share/TofCorrPID/";

  if( abs(run)>=8093 && abs(run)<=9025 ) {
    filePath = filePath + "psip2009";
    m_runBegin = 8093;
    m_runEnd   = 9025;
  }
  else if( abs(run)>=9947 && abs(run)<=10878 ) {
    filePath = filePath + "jpsi2009";
    m_runBegin = 9947;
    m_runEnd   = 10878;
  }

  if( run>0 ) {
    filePath = filePath + "/data/";
  }
  else {
    filePath = filePath + "/mc/";
  }


  // weight from tof calibration
  std::string fileWeight = filePath + "calib_barrel_sigma.txt";
  ifstream inputWeight( fileWeight.c_str(), std::ios_base::in );
  if( !inputWeight ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileWeight << endl;
    exit(1);
  }

  for( unsigned int tofid=0; tofid<176; tofid++ ) {
    for( unsigned int readout=0; readout<3; readout++ ) {
      for( unsigned int p_i=0; p_i<5; p_i++ ) {
	inputWeight >> m_p_weight[tofid][readout][p_i];
      }
    }
  }
  //  cout << "finish read " << fileWeight << endl;

  // common item, from bunch size and bunch time
  std::string fileCommon = filePath + "calib_barrel_common.txt";
  ifstream inputCommon( fileCommon.c_str(), std::ios_base::in );
  if( !inputCommon ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileCommon << endl;
    exit(1);
  }
  inputCommon >> m_p_common;
  //  cout << "finish read " << fileCommon << endl;

  // endcap sigma
  std::string fileEcSigma = filePath + "calib_endcap_sigma.txt";
  ifstream inputEcSigma( fileEcSigma.c_str(), std::ios_base::in );
  if( !inputEcSigma ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileEcSigma << endl;
    exit(1);
  }

  for( unsigned int tofid=0; tofid<96; tofid++ ) {
    for( unsigned int p_i=0; p_i<3; p_i++ ) {
      inputEcSigma >> m_ec_sigma[tofid][p_i];
    }
  }
  //  cout << "finish read " << fileEcSigma << endl;

  // curve of betaGamma versus Q0
  std::string fileQ0BetaGamma = filePath + "curve_Q0_BetaGamma.txt";
  ifstream inputQ0BetaGamma( fileQ0BetaGamma.c_str(), std::ios_base::in );
  if( !inputQ0BetaGamma ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileQ0BetaGamma << endl;
    exit(1);
  }
  // barrel layer1 layer2 and endcap
  for( unsigned int layer=0; layer<3; layer++ ) {
    for( unsigned int ipar=0; ipar<5; ipar++ ) {
      inputQ0BetaGamma >> m_q0_bg[layer][ipar];
    }
  }
  //  cout << "finish read " << fileQ0BetaGamma << endl;

  // paramter of A and B
  std::string fileParAB = filePath + "parameter_A_B.txt";
  ifstream inputParAB( fileParAB.c_str(), std::ios_base::in );
  if( !inputParAB ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileParAB << endl;
    exit(1);
  }
  // 0: inner east, 1: inner west, 2: outer east, 3: outer west, 4: west endcap 
  // 0: parameter A, 1: parameter B
  for( unsigned int ipmt=0; ipmt<5; ipmt++ ) {
    for( unsigned int iab=0; iab<2; iab++ ) {
      for( unsigned int ipar=0; ipar<5; ipar++ ) {
	inputParAB >> m_par_ab[ipmt][iab][ipar];
      }
    }
  }
  for( unsigned int ipmt=0; ipmt<5; ipmt++ ) {
    for( unsigned int iab=0; iab<2; iab++ ) {
      for( unsigned int ipar=0; ipar<5; ipar++ ) {
	inputParAB >> m_par_pbar_ab[ipmt][iab][ipar];
      }
    }
  }
  //  cout << "finish read " << fileParAB << endl;

  // sigma for pion, kaon and proton
  std::string fileSigma = filePath + "parameter_sigma.txt";
  ifstream inputSigma( fileSigma.c_str(), std::ios_base::in );
  if( !inputSigma ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileSigma << endl;
    exit(1);
  }
  // 0: pion,  1: kaon,  2: proton,  3: anti-proton
  // 0: inner east, 1: inner west, 2: outer east, 3: outer west
  // 4: inner layer, 5: outer layer, 6: double layer
  // 7: west endcap
  for( unsigned int ispecies=0; ispecies<4; ispecies++ ) {
    for( unsigned int ipmt=0; ipmt<8; ipmt++ ) {
      for( unsigned int ipar=0; ipar<9; ipar++ ) {
	inputSigma >> m_par_sigma[ispecies][ipmt][ipar];
      }
    }
  }
  //  cout << "finish read " << fileSigma << endl;

  // chi for pion, kaon and proton
  /*
  std::string fileChi = filePath + "parameter_sigma_momentum.txt";
  ifstream inputChi( fileChi.c_str(), std::ios_base::in );
  if( !inputChi ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileChi << endl;
    exit(1);
  }
  // 0: inner east, 1: inner west, 2: outer east, 3: outer west
  // 4: inner layer, 5: outer layer, 6: double layer
  for( unsigned int ispecies=0; ispecies<3; ispecies++ ) {
    for( unsigned int ipmt=0; ipmt<7; ipmt++ ) {
      for( unsigned int ipar=0; ipar<3; ipar++ ) {
	inputChi >> m_par_sig_mom[ispecies][ipmt][ipar];
      }
    }
  }
  */
  //  cout << "finish read " << fileChi << endl;


  // offset for low momentum proton and anti-proton
  std::string fileProtonOffset = filePath + "parameter_offset_proton.txt";
  ifstream inputProtonOffset( fileProtonOffset.c_str(), std::ios_base::in );
  if( !inputProtonOffset ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileProtonOffset << endl;
    exit(1);
  }
  // 0: proton,  1: anti-proton
  // 0: inner east, 1: inner west, 2: outer east, 3: outer west
  for( unsigned int ispecies=0; ispecies<2; ispecies++ ) {
    for( unsigned int ipmt=0; ipmt<4; ipmt++ ) {
      for( unsigned int ipar=0; ipar<10; ipar++ ) {
	for( unsigned int jpar=0; jpar<20; jpar++ ) {
	  inputProtonOffset >> m_p_offset[ispecies][ipmt][ipar][jpar];
	}
      }
    }
  }
  //  cout << "finish read " << fileProtonOffset << endl;

  // sigma for low momentum proton and anti-proton
  std::string fileProtonSigma = filePath + "parameter_sigma_proton.txt";
  ifstream inputProtonSigma( fileProtonSigma.c_str(), std::ios_base::in );
  if( !inputProtonSigma ) {
    cout << "ParticleID::TofCorrPID: Can NOT open file: " << fileProtonSigma << endl;
    exit(1);
  }
  // 0: proton,  1: anti-proton
  // 0: inner east, 1: inner west, 2: outer east, 3: outer west
  // 4: inner layer, 5: outer layer, 6: double layer
  // 7: west endcap
  for( unsigned int ispecies=0; ispecies<2; ispecies++ ) {
    for( unsigned int ipmt=0; ipmt<7; ipmt++ ) {
      for( unsigned int ipar=0; ipar<10; ipar++ ) {
	for( unsigned int jpar=0; jpar<20; jpar++ ) {
	  inputProtonSigma >> m_p_sigma[ispecies][ipmt][ipar][jpar];
	}
      }
    }
  }
  //  cout << "finish read " << fileProtonSigma << endl;


  return;
}


double TofCorrPID::offsetTof( unsigned int ispecies, bool barrel, unsigned int ipmt, double betaGamma, int charge, double zrhit, double dt ) {
  if( ispecies==0 || ispecies==1 ) { return dt; }

  double deltaT = -1000.0;
  if( ( ipmt>= 4 && barrel ) || ( ipmt!=0 && !barrel ) || betaGamma<0.0 || abs(charge)!=1 || fabs(zrhit)>120.0 ) {
    cout << "Particle::TofCorrPID: offsetTof for single end: the input parameter are NOT correct! Please check them!" << endl;
    return deltaT;
  }
  unsigned int layer=0;
  if( barrel ) {
    if( ipmt==0 || ipmt==1 ) { layer=0; }
    else if( ipmt==2 || ipmt==3 ) { layer=1; }
  }
  else { layer=2; }
  double q0 = qCurveFunc( layer, betaGamma );

  unsigned int inumber=ipmt;
  if( !barrel ) { inumber=4; }

  double func[5];
  func[0] = 1.0;
  func[1] = zrhit;
  func[2] = zrhit*zrhit;
  func[3] = zrhit*zrhit*zrhit;
  func[4] = zrhit*zrhit*zrhit*zrhit;

  double parA = 0.0;
  double parB = 0.0;
  // anti-proton
  if( ispecies==4 && charge==-1 ) {
    for( unsigned int i=0; i<5; i++ ) {
      parA += m_par_pbar_ab[inumber][0][i]*func[i];
      parB += m_par_pbar_ab[inumber][1][i]*func[i];
    }
  }
  // proton
  else {
    for( unsigned int i=0; i<5; i++ ) {
      parA += m_par_ab[inumber][0][i]*func[i];
      parB += m_par_ab[inumber][1][i]*func[i];
    }
  }

  double tcorr = parA + parB/sqrt(q0);

  // barrel TOF low momentum proton and anti-proton
  if( barrel && ispecies==4 && betaGamma<0.8 ) {
    int    nzbin  = 10;
    double zbegin = -115.0;
    double zend   = 115.0;
    double zstep  = (zend-zbegin)/nzbin;
    
    double nbgbin = 20.0;
    double bgbegin[2] = { 0.4, 0.55 };
    double bgend[2]   = { 0.8, 0.8  };
    double bgstep[2];
    bgstep[0] = (bgend[0]-bgbegin[0])/nbgbin;
    bgstep[1] = (bgend[1]-bgbegin[1])/nbgbin;

    int izbin = static_cast<int>((zrhit-zbegin)/zstep);
    if( izbin<0 )           { izbin=0;       }
    else if( izbin>=nzbin ) { izbin=nzbin-1; }
    unsigned int layer=0;
    if( ipmt==2 || ipmt==3 ) { layer=1; }
    int ibgbin = static_cast<int>((betaGamma-bgbegin[layer])/bgstep[layer]);
    if( ibgbin<0 )            { ibgbin=0;        }
    else if( ibgbin>=nbgbin ) { ibgbin=nbgbin-1; }

    if( charge==1 ) {
      tcorr += m_p_offset[0][ipmt][izbin][ibgbin];
    }
    else {
      tcorr += m_p_offset[1][ipmt][izbin][ibgbin];
    }
  }

  deltaT = dt - tcorr;

  return deltaT;
}


double TofCorrPID::offsetTof( unsigned int ispecies, int tofid, double zrhit, double dt ) {
  if( ispecies==0 || ispecies==1 ) { return dt; }

  double deltaT = -1000.0;
  if( tofid<0 || tofid>=176 ) {
    cout << "Particle::TofCorrPID: offsetTof for single layer: the input parameter are NOT correct! Please check them!" << endl;
    exit(1);
  }
  int pmt[3] = { -9, -9, -9 };
  if( tofid>=0 && tofid<=87 ) {
    pmt[0] = 0;
    pmt[1] = 1;
    pmt[2] = 4;
  }
  else {
    pmt[0] = 2;
    pmt[1] = 3;
    pmt[2] = 5;
  }

  double sigmaCorr2  = m_p_common*m_p_common;
  double sigmaLeft   = bSigma( 0, tofid, zrhit );
  double sigmaLeft2  = sigmaLeft*sigmaLeft;
  double sigmaRight  = bSigma( 1, tofid, zrhit );
  double sigmaRight2 = sigmaRight*sigmaRight;

  double fraction = ( sigmaRight2 - sigmaCorr2 )/( sigmaLeft2 + sigmaRight2 - 2.0*sigmaCorr2);
  deltaT = fraction*m_deltaT[pmt[0]][ispecies] + (1.0-fraction)*m_deltaT[pmt[1]][ispecies];

  return deltaT;
}


double TofCorrPID::offsetTof( unsigned int ispecies, int tofid1, int tofid2, double zrhit1, double zrhit2, double dt ) {
  if( ispecies==0 || ispecies==1 ) { return dt; }

  double deltaT = -1000.0;
  if( tofid1<0 || tofid1>=88 || tofid2<88 || tofid2>=176 ) {
    cout << "Particle::TofCorrPID: offsetTof for double layer: the input parameter are NOT correct! Please check them!" << endl;
    exit(1);
  }

  double sigmaCorr2  = m_p_common*m_p_common;
  double sigmaInner  = bSigma( 2, tofid1, zrhit1 );
  double sigmaInner2 = sigmaInner*sigmaInner;
  double sigmaOuter  = bSigma( 2, tofid2, zrhit2 );
  double sigmaOuter2 = sigmaOuter*sigmaOuter;
  double sigma = sqrt( (sigmaInner2*sigmaOuter2-sigmaCorr2*sigmaCorr2)/(sigmaInner2+sigmaOuter2-2.0*sigmaCorr2) );

  m_sigma[0] = sigma;
  m_sigma[1] = sigma;

  double fraction = ( sigmaOuter2 - sigmaCorr2 )/( sigmaInner2 + sigmaOuter2 - 2.0*sigmaCorr2);
  deltaT = fraction*m_deltaT[4][ispecies] + (1.0-fraction)*m_deltaT[5][ispecies];

  return deltaT;
}


double TofCorrPID::sigmaTof( unsigned int ispecies, int charge, bool barrel, unsigned int ipmt, int tofid[2], double zrhit[2], double betaGamma ) {

  double sigma = 1.0e-6;

  if( ispecies==0 || ispecies==1 ) {
    if( barrel ) {
      if( ipmt==0 ) {
	sigma = bSigma( 0, tofid[0], zrhit[0] );
      }
      else if( ipmt==1 ) {
	sigma = bSigma( 1, tofid[0], zrhit[0] );
      }
      else if( ipmt==2 ) {
	sigma = bSigma( 0, tofid[1], zrhit[1] );
      }
      else if( ipmt==3 ) {
	sigma = bSigma( 1, tofid[1], zrhit[1] );
      }
      else if( ipmt==4 ) {
	sigma = bSigma( 2, tofid[0], zrhit[0] );
      }
      else if( ipmt==5 ) {
	sigma = bSigma( 2, tofid[1], zrhit[1] );
      }
      else if( ipmt==6 ) {
	sigma = bSigma( &tofid[0], &zrhit[0] );
      }
    }
    else {
      sigma = eSigma( tofid[0], zrhit[0] );
    }
  }
  else {
    int ibgbin = -1;
    unsigned int iz = 0;
    if( ipmt==2 || ipmt==3 || ipmt==5 ) { iz=1; }
    // low momentum proton and anti-proton
    if( barrel && ispecies==4 && betaGamma<0.8 ) {
      double nbgbin = 20.0;
      double bgbegin[2] = { 0.4, 0.55 };
      double bgend[2]   = { 0.8, 0.8  };
      double bgstep[2];
      bgstep[0] = (bgend[0]-bgbegin[0])/nbgbin;
      bgstep[1] = (bgend[1]-bgbegin[1])/nbgbin;

      ibgbin = static_cast<int>((betaGamma-bgbegin[iz])/bgstep[iz]);
      if( ibgbin<0 )            { ibgbin=0;        }
      else if( ibgbin>=nbgbin ) { ibgbin=nbgbin-1; }
    }

    sigma = sigmaTof( ispecies, charge, barrel, ipmt, zrhit[iz], ibgbin );
  }

  return sigma;
}


double TofCorrPID::sigmaTof( unsigned int ispecies, int charge, bool barrel, unsigned int ipmt, double zrhit, int ibgbin ) {

  int izbin=0;
  if( barrel && ispecies==4 && ibgbin!=-1 ) {
    int    nzbin  = 10;
    double zbegin = -115.0;
    double zend   = 115.0;
    double zstep  = (zend-zbegin)/nzbin;

    izbin = static_cast<int>((zrhit-zbegin)/zstep);
    if( izbin<0 )           { izbin=0;       }
    else if( izbin>=nzbin ) { izbin=nzbin-1; }
  }

  double func[10];
  for( unsigned int i=0; i<9; i++ ) {
    if( i==0 ) { func[i] = 1.0; }
    else {
      func[i] = func[i-1]*zrhit;
    }
  }
  func[9] = -1.0/(zrhit*zrhit-115.0*115.0);

  unsigned int inumber = ipmt;
  if( !barrel ) { inumber=7; }

  double sigma = 0.0;
  // barrel
  if( barrel ) {
    // barrel inner layer east/west
    if( ipmt==0 || ipmt==1 ) {
      if( ispecies==2 || ispecies==3 ) { // pion / kaon
	for( unsigned int i=0; i<5; i++ ) {
	  if( i!=4 ) {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
	  }
	  else {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[9];
	  }
	}
      }
      else if( ispecies==4 ) {
	if( ibgbin!=-1 ) {
	  if( charge==1 ) {
	    sigma = m_p_sigma[0][inumber][izbin][ibgbin];
	  }
	  else {
	    sigma = m_p_sigma[1][inumber][izbin][ibgbin];
	  }
	}
	else {
	  for( unsigned int i=0; i<7; i++ ) {
	    if( charge==1 ) {
	      sigma += m_par_sigma[2][inumber][i]*func[i];
	    }
	    else {
	      sigma += m_par_sigma[3][inumber][i]*func[i];
	    }
	  }
	}
      }
    }
    // barrel outer layer east/west
    else if( ipmt==2 || ipmt==3 ) {
      if( ispecies==2 || ispecies==3 ) { // pion / kaon
	for( unsigned int i=0; i<4; i++ ) {
	  if( i!=3 ) {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
	  }
	  else {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[9];
	  }
	}
      }
      else if( ispecies==4 ) {
	if( ibgbin!=-1 ) {
	  if( charge==1 ) {
	    sigma = m_p_sigma[0][inumber][izbin][ibgbin];
	  }
	  else {
	    sigma = m_p_sigma[1][inumber][izbin][ibgbin];
	  }
	}
	else {
	  for( unsigned int i=0; i<5; i++ ) {
	    if( charge==1 ) {
	      sigma += m_par_sigma[2][inumber][i]*func[i];
	    }
	    else {
	      sigma += m_par_sigma[3][inumber][i]*func[i];
	    }
	  }
	}
      }
    }
    // barrel inner layer and outer layer
    else if( ipmt==4 || ipmt==5 ) {
      if( ispecies==2 ) { // pion
	for( unsigned int i=0; i<5; i++ ) {
	  if( i!=4 ) {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
	  }
	  else {
	    sigma += m_par_sigma[ispecies-2][inumber][i]*func[9];
	  }
	}
      }
      else if( ispecies==3 ) { // kaon
	for( unsigned int i=0; i<4; i++ ) {
	  sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
	}
      }
      else if( ispecies==4 ) {
	if( ibgbin!=-1 ) {
	  if( charge==1 ) {
	    sigma = m_p_sigma[0][inumber][izbin][ibgbin];
	  }
	  else {
	    sigma = m_p_sigma[1][inumber][izbin][ibgbin];
	  }
	}
	else {
	  for( unsigned int i=0; i<9; i++ ) {
	    if( charge==1 ) {
	      sigma += m_par_sigma[2][inumber][i]*func[i];
	    }
	    else {
	      sigma += m_par_sigma[3][inumber][i]*func[i];
	    }
	  }
	}
      }
    }
    // barrel double layer
    else if( ipmt==6 ) {
      if( ispecies==2 || ispecies==3 ) { // pion / kaon
	for( unsigned int i=0; i<5; i++ ) {
	  sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
	}
      }
      else if( ispecies==4 ) {
	if( ibgbin!=-1 ) {
	  if( charge==1 ) {
	    sigma = m_p_sigma[0][inumber][izbin][ibgbin];
	  }
	  else {
	    sigma = m_p_sigma[1][inumber][izbin][ibgbin];
	  }
	}
	else {
	  for( unsigned int i=0; i<9; i++ ) {
	    if( charge==1 ) {
	      sigma += m_par_sigma[2][inumber][i]*func[i];
	    }
	    else {
	      sigma += m_par_sigma[3][inumber][i]*func[i];	      
	    }
	  }
	}
      }
    }
  }
  // endcap
  else {
    if( ispecies==2 || ispecies==3 ) { // pion / kaon
      for( unsigned int i=0; i<3; i++ ) {
	sigma += m_par_sigma[ispecies-2][inumber][i]*func[i];
      }
    }
    else if( ispecies==4 ) {
      for( unsigned int i=0; i<4; i++ ) {
	int iparticle=4;
	if( charge==-1 ) { iparticle=5; }
	if( i!=3 ) {
	  sigma += m_par_sigma[iparticle-2][inumber][i]*func[i];
	}
	else {
	  sigma += m_par_sigma[iparticle-2][inumber][i]*func[9];
	}
      }
    }
  }

  /*
  double chi = 1.0;
  if( barrel ) {
    if( ( ispecies==3 || ispecies==4 ) && ipmt<4 ) {
      chi = m_par_sig_mom[ispecies-2][inumber][0];
    }
    else {
      chi = m_par_sig_mom[ispecies-2][inumber][0]*exp(0.0-m_par_sig_mom[ispecies-2][inumber][1]*p)+m_par_sig_mom[ispecies-2][inumber][2];
    }
  }

  sigma = sigma*chi;
  */

  return sigma;
}


double TofCorrPID::qCurveFunc( unsigned int layer, double betaGamma ) {
  double q0 = -100.0;
  if( layer>=3 || betaGamma<0.0 ) {
    cout << "Particle::TofCorrPID::qCurveFunc: the input parameter are NOT correct! Please check them!" << endl;
    return q0;
  }

  double beta = betaGamma/sqrt(1.0+betaGamma*betaGamma);
  double logterm = log( m_q0_bg[layer][2]+pow((1.0/betaGamma), m_q0_bg[layer][4] ) );
  q0 = m_q0_bg[layer][0]*( m_q0_bg[layer][1]-pow( beta, m_q0_bg[layer][3] ) - logterm )/pow( beta, m_q0_bg[layer][3] );

  return q0;
}


double TofCorrPID::bSigma( unsigned int end, int tofid, double zrhit ) {

  double func[5];
  func[0] = 1.0;
  func[1] = zrhit;
  func[2] = zrhit*zrhit;
  func[3] = zrhit*zrhit*zrhit;
  func[4] = zrhit*zrhit*zrhit*zrhit;

  double sigma = 0.0;
  for( unsigned int i=0; i<5; i++ ) {
    sigma += m_p_weight[tofid][end][i]*func[i];
  }

  return sigma;
}


double TofCorrPID::bSigma( int tofid[2], double zrhit[2] ) {

  double sigma1 = bSigma( 2, tofid[0], zrhit[0] );
  double sigma2 = bSigma( 2, tofid[1], zrhit[1] );
  double sigmaCorr2  = m_p_common*m_p_common;
  double sigma = ( sigma1*sigma1*sigma2*sigma2 - sigmaCorr2*sigmaCorr2 )/( sigma1*sigma1 + sigma2*sigma2 - 2.0*sigmaCorr2 );
  sigma = sqrt(fabs(sigma));

  return sigma;
}


double TofCorrPID::eSigma( int tofid, double zrhit ) {

  double func[5];
  func[0] = 1.0;
  func[1] = zrhit;
  func[2] = zrhit*zrhit;

  double sigma = 0.0;
  for( unsigned int i=0; i<3; i++ ) {
    sigma += m_ec_sigma[tofid][i]*func[i];
  }

  return sigma;
}


bool TofCorrPID::correlationCheck() {
  bool chiCut = false;
  bool good = false;
  double pdftemp = 0;
  for( unsigned int i=0; i<5; i++ ) {
    m_chi[i] = m_offset[i]/m_sigma[i];
    if( ( m_chi[i]>-4.0 && m_chi[i]<6.0 ) && !good ) { good=true; }
    double ppp = pdfCalculate(m_chi[i],1);
    if( fabs(ppp) > pdftemp) { pdftemp = fabs(ppp); }
  }
  m_pdfmin = pdftemp;
  if( pdftemp < pdfCalculate(pdfMinSigmaCut(),1) ) return chiCut;
  if( !good ) return chiCut;
  chiCut = true;
  return chiCut;
}
