//==============================================================================================//
//
// SelectJsiX - search for e+ e- -> J/psi X-> (mu+ mu-) X
//
//==============================================================================================//
//
// 3-09-2021  (+BOSS 7.0.7)
//
//
//
//N psi'X(ISR)  integral info     18-01-2021
//16-02-2021  + 3 sigma m_Jpsi, nGood MDC
//31-03-2021   sigma_chic1,2 = 0.008 GeV
//5-04-2021    N J/psi X(ISR)  integral info
//8-04-2021    N chic integral info, sigma = 7.5 MeV
//3-06-2021    EventNavigator + EMC  ->cellId()
//3-09-2021    J/psi: 3.0 - 3.2
//3-09-2021    chic1: 3.47-3.53343, chic2: 3.53343-3.6
//8-09-2021   psi': 3.65-3.72,  psi' ISR: 3.59-3.77
//14-12-2021  iMC for bkgr study
//14-12-2021   J/psi Data vs MC
//==============================================================================================//





#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <fstream>              // to save momentums in file
#include <cmath>
#include <vector>
#include <map>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/TwoVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/LorentzRotation.h>
#include <CLHEP/Geometry/Point3D.h>

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using CLHEP::pi;
using CLHEP::twopi;
using CLHEP::HepLorentzRotation;


#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TMcParticle.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

//=====================================
#include "RootEventData/TEvtNavigator.h"
#include "DstFormat.h"
//=====================================

#include "AbsCor/AbsCor.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "EventTag/EventTagSvc.h"

#include "DstEvtRecTracks.h"
#include "TofHitStatus.h"
#include "ReadDst.h"

using namespace std;




//for eMC chi_c1,2
//#define MC_chi_c1 OK
//#define MC_chi_c2 OK
//#define MC_psi2S_2pi OK



struct Select_JpsiX {
  // charge tracks:
  vector<int> qtrksPos;  //+trk# in EvtRecTrkCol (Kalman)
  vector<int> qtrksNeg;  //-trk# in EvtRecTrkCol (Kalman)
  vector<int> qtrks;     // trk# in EvtRecTrkCol (MDC)
  vector<int> qmups;
  vector<int> qmums;
  vector<int> qeps;
  vector<int> qems;
  vector<int> qmuEPp;
  vector<int> qmuEPm;
  vector<int> qeEPp;
  vector<int> qeEPm;
  vector<int> qpps;
  vector<int> qpms;
  vector<int> qpips;
  vector<int> qpims;
  vector<int> qpipsFromPsi2S;
  vector<int> qpimsFromPsi2S;
  vector<int> qpipPsi2StrkID;
  vector<int> qpimPsi2StrkID;
  vector<int> qmupJpsiTrkID;
  vector<int> qmumJpsiTrkID;
  vector<int> qgChiCshwrID;
  vector<int> qpipENsort;
  vector<int> qpimENsort;
  vector<int> qgENsort;
  vector<int> qmupENsort;
  vector<int> qmumENsort;
  vector<int> qkps;
  vector<int> qkms;
  vector<int> qgchic;
  vector<double> qRxy0;
  vector<double> qRxy0Pos;
  vector<double> qRxy0Neg;
  vector<double> qkpRxy0;
  vector<double> qkmRxy0;
  vector<double> qpmRxy0;
  vector<double> qpipRxy0;
  vector<double> qpimRxy0;
  vector<double> qRz0;
  vector<double> qRz0Pos;
  vector<double> qRz0Neg;
  vector<double> qkpRz0;
  vector<double> qkmRz0;
  vector<double> qpipRz0;
  vector<double> qpimRz0;
  vector<double> qpmRz0;
  vector<double> qCosTheta; //all MDC trks
  vector<double> qMom; //all MDC trks

  int ipip, ipim, iKp, iKm, imup, imum, iep, iem, ilepp, ilepm;               // trk# for pi+ & pi-
  HepLorentzVector ppip, ppim, pmup, pmum, pmup_1C, pmum_1C, pep, pem, plepp, plepm, pKp, pKm, pKPos[2], pKNeg[2], ppiPos[2], ppiNeg[2];  // 4-momentum of pi+/-
  int ipip_mc, ipim_mc;         // correspondent mcParticle
  int pip_pdg, pim_pdg;         // PDG code of these mc-particles (or 0)

  int id, iad, ip[2], iap[2];
  double Eg_chic1, Eg_chic2;

  // neutral (gammas) tracks
  vector<int> gammas;           // trk# in EvtRecTrkCol
  vector<HepLorentzVector> Pg;  // 4-momentum of gammas
  vector<HepLorentzVector> Pmups;
  vector<HepLorentzVector> Pmums;
  vector<HepLorentzVector> Peps;
  vector<HepLorentzVector> Pems;
  vector<HepLorentzVector> PmuEPp;
  vector<HepLorentzVector> PmuEPm;
  vector<HepLorentzVector> PeEPp;
  vector<HepLorentzVector> PeEPm;
  vector<HepLorentzVector> Ppps;
  vector<HepLorentzVector> Ppms;
  vector<HepLorentzVector> Ppips;
  vector<HepLorentzVector> Ppims;
  vector<HepLorentzVector> PpipsFromPsi2S;
  vector<HepLorentzVector> PpimsFromPsi2S;
  vector<HepLorentzVector> PgFromChiC;
  vector<HepLorentzVector> Pkps;
  vector<HepLorentzVector> Pkms;

  int ig1, ig2, ig3, ig4, ig5, ig6;                 // the best gamma pair for pi0
  HepLorentzVector ppi0;        // 4-momentum of pi0 = g+g

  HepLorentzVector Ptot;        // 4-momentum of the overall system

  // momentums of pi(+,-,0) after kinematic fit
  HepLorentzVector Ppip, Ppim, Pmum, Pmup, PISR_missing, Pem, Pep, Pgg, PKp, PKm, Pphi, Pf1;

  Select_JpsiX() {
    ipip = ipim = -1;
    imup = imum = -1;
    ipip_mc = ipim_mc = -1;
    pip_pdg = pim_pdg = 0;
    ig1 = ig2 = ig3 = ig4 = ig5 = ig6 = -1;
    id = iad = -1;
  }

  int nGood() { return qtrks.size(); }
  int nPos() { return qtrksPos.size(); }
  int nNeg() { return qtrksNeg.size(); }
  int mc_trk(int i) { return (qtrks[i] == ipip) ? ipip_mc : ipim_mc; }
  int isPiPlus(int i) { return qtrks[i] == ipip; }
  int isPiMinus(int i) { return qtrks[i] == ipim; }
  int nGam() { return gammas.size(); }
  int gamma1() { return (ig1 >= 0) ? gammas[ig1] : -1; }
  int gamma2() { return (ig2 >= 0) ? gammas[ig2] : -1; }
  int gamma3() { return (ig3 >= 0) ? gammas[ig3] : -1; }
  int gamma4() { return (ig4 >= 0) ? gammas[ig4] : -1; }
  int gamma5() { return (ig5 >= 0) ? gammas[ig5] : -1; }
  int gamma6() { return (ig6 >= 0) ? gammas[ig6] : -1; }
  bool GoodGammas() { return (ig1 != -1 && ig2 != -1); }
  bool GoodGammas6() { return (ig1 != -1 && ig2 != -1  && ig3 != -1  && ig4 != -1  && ig5 != -1  && ig6 != -1); }
};

typedef Select_JpsiX Select;

//-----------------------------------------------------------------------------
// Global file variables
//-----------------------------------------------------------------------------
const static bool init_selection=false;


//for ISR cut
const static double cosThetaISR = 1.95; // 0.95
const static double pFracISR = 1.85; //0.85


const static double mk = 0.493677;
const static double mpi = 0.13957;
const static double mpi0 = 0.13498;
const static double meta = 0.547862;
const static double meta958 = 0.95778;
const static double mphi = 1.019461;
const static double mLambda = 1.115683;
const static double mKs = 0.497614;
const static double mf1 = 1.2819;
const static double mmu = 0.1056583715;
const static double me = 0.00051099828;
const static double mJpsi = 3.096916;
const static double mpsi2S = 3.686097;
const static double mchic1 = 3.51066;
const static double mchic2 = 3.5562;
const static double mp = 0.938272046;
const static double md = 1.875612928;
const static double mDpm = 1.86961;
const static double mD0 = 1.86484;
const static double mDspm = 1.9683;


const static double beam_angle = 0.011; // 11 mrad


//const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
//const double velc = 299.792458;   // tof path unit in mm


static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;

static int ERROR_WARNING = 0;
static bool isMC = false;
static bool iMCJpsiISR = false;
static bool iMCJpsi = false;
static bool iMCpsi2SToJpsi2pi = false;
static bool iMCchiC1ToJpsiGamma = false;
static bool iMCchiC2ToJpsiGamma = false;
static bool MCJpsiTo2mu = false;
static bool JpsiVertexFit = false;
static bool Jpsi1CFit = false;
static bool leptonISmuon=false;
static bool Jpsi2piFromPsi2S=false;
static bool JpsiGammaFromChiC1=false;
static bool JpsiGammaFromChiC2=false;
static bool MuonsFromJpsi=false;
static bool MuonsFromJpsi3sigma=false;
static bool MuonsFromJpsi3032=false;
static bool ElectronsFromJpsi=false;
static bool MuonsFromPsi2S=false;
static bool ElectronsFromPsi2S=false;
static bool isAntiproton = false;
static bool isAntideuteron = false;
static bool isDeuteron = false;
static bool isLambda = false;
static bool isAntiLambda = false;
static double pJpsi_iMC = 0.;
static bool iMCpsi2SInEvent = false;
static bool iMCchiC1InEvent = false;
static bool iMCchiC2InEvent = false;
static bool JpsiWithoutPcut = false;
static bool JpsiWithPcut = false;
static bool SBJpsi=false;
static bool SBpsi2S=false;
static bool ISRcutRegion=false;
static bool isISRevent=false;
static bool isTrk4psi=false;
static bool isGoodMCJpsi24pi = true;
static bool isJpsiX=false;
static bool isJpsiISR=false;
static bool ispsiX=false;
static bool ispsiISR=false;
static bool ischic1X=false;
static bool ischic2X=false;
static bool JpsiSelectedForDataMCcompar=false;
static double Ecms = 0.;


// Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
inline Double_t RtoD(Double_t ang){return ang*180./M_PI;}
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------------- Initialization -------------------------------------------------------------
void SelectJpsiXStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << " SelectJpsiXStartJob() " << endl;
  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }

  his1.resize(1000,(TH1D*)0);
  his2.resize(1000,(TH2D*)0);

  // init Absorption Correction
  m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

  //initialize EventTag
  m_EventTagSvc = EventTagSvc::instance();

  // set paths to pdg & decayCodes files:
  m_EventTagSvc->setPdtFile(
        selector->AbsPath("Analysis/EventTag/share/pdt.table") );
  m_EventTagSvc->setDecayTabsFile(
        selector->AbsPath("Analysis/EventTag/share/decay.codes") );

  if ( selector->Verbose() ) m_EventTagSvc->setVerbose(1);

  m_EventTagSvc->initialize();

//-------------------------------------------------------------------------------- We have to initialize DatabaseSvc-------------------------------------------

  DatabaseSvc* dbs = DatabaseSvc::instance();
  if ( (dbs->GetDBFilePath()).empty() ) {
    // set path to directory with databases:
    dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
  }

  //------------------------------------------------------------------------------ We have to initialize Magnetic field----------------------------------------

  MagneticFieldSvc* mf = MagneticFieldSvc::instance();
  if ( !(mf->GetPath()).empty() ) {
    cerr << " RhopiStartJob::WARNING:"
         << " MagneticFieldSvc has already been initialized" << endl
         << "                         path = " << mf->GetPath() << endl;
  }
  // set path to directory with magnetic fields tables
  mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
  mf->UseDBFlag(false); // like in the boss program
  mf->RunMode(3); // like in the boss program
  //------------------------------------------------------------------------------------------ Book histograms--------------------------------------------------

  his1[0] = new TH1D("Runs", "run numbers", 1000, 0.,50000);
  his1[1] = new TH1D("PidInfo", "0 - no pid, 1 - pid", 2, -0.5, 1.5);
  his1[2] = new TH1D("Rxy","R_{xy}", 100, -5., 5.);
  his1[3] = new TH1D("Rz","R_{z}", 100, -50., 50.);
  his1[4] = new TH1D("cos_theta","cos(#theta)", 100, 0., 1.);
  his1[5] = new TH1D("Ntrks","N_{charged tracks} in event", 11, -0.5, 10.5);

  his1[6] = new TH1D("Ncharge","N_{charge} in event", 11, -5.5, 5.5);
  his1[7] = new TH1D("AngleTrack","The angle between the photon & the nearest charged track, degrees", 180, 0., 180.);
  his1[8] = new TH1D("LogCLpi","lg(CL_{#pi})", 100, -5., 0.);
  his1[28] = new TH1D("LogCLmu","lg(CL_{#mu})", 100, -5., 0.);
  his1[9] = new TH1D("QPtrack","Charged momentum (Q*P)", 100, -2.5, 2.5);
  his1[10] = new TH1D("Ngamma","N_{gamma} in event", 11, -0.5, 10.5);

  his1[11] = new TH1D("E_p","#frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);
  his1[12] = new TH1D("E_p_new","NEW  #frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);
  his1[13] = new TH1D("E_p_fin","FINAL  #frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);

  his1[14] = new TH1D("vtx_fit", "0 - no fit, 1 - good fit", 2, -0.5, 1.5);
  his1[15] = new TH1D("vtx_chi2", "vertex fit #chi^2 for J/#psi-> mu^{+}mu^{-}", 100, 0., 100.);
  his1[25] = new TH1D("vtx_chi2_eta958", "vertex fit #chi^2 for  #eta'", 200, 0., 200.);
  his1[16] = new TH1D("Mgg_eta_4C", "4C-fit: M_{#gamma#gamma}(for #eta)", 500, 0.3, 0.8);
  his1[17] = new TH1D("Etot_4C","4C-fit: E_{tot}", 100, 2.5, 3.5);
  his1[317] = new TH1D("Etot_4C_4","4C-fit: E_{tot}, 4 charged tracks", 100, 2.5, 3.5);
  his1[318] = new TH1D("Etot_4C_2","4C-fit: E_{tot}, 2 charged tracks", 100, 2.5, 3.5);
  his1[18] = new TH1D("Ptot_4C","4C-fit: P_{tot}", 200, 0, 1.0);
  his1[318] = new TH1D("Ptot_4C_4","4C-fit: P_{tot},  4 charged tracks", 200, 0, 1.0);
  his1[319] = new TH1D("Ptot_4C_2","4C-fit: P_{tot},  2 charged tracks", 200, 0, 1.0);

  his1[21] = new TH1D("ini_P","Charged momentum (Q*P)", 100, -2., 2.);
  his1[22] = new TH1D("ini_Psum","P(K^{+})+P(K^{-})", 100, 0., 2.);
  his1[23] = new TH1D("ini_Pgmax","energy of maximal gammas", 100, 0., 2.);
  his1[24] = new TH1D("ini_Pgsum","energy of sum mom. of gammas", 100, 0., 2.);

  his2[1] = new TH2D("depth", "Depth_{-} vs Depth_{+}", 100, -10., 60., 100, -10., 60.);
  his2[2] = new TH2D("pi_charge","K+ vs K-",10,-0.5,9.5,10,-0.5,9.5);

  his1[30] = new TH1D("Ecms","E_{cms}", 200, 3.6, 5.6);
  his1[31] = new TH1D("chi2_1C", "#chi^{2} for 1C-Fit J/#psi -> mu^{+}mu^{-}", 300, 0, 300 );
  his1[331] = new TH1D("chi2_4C_4", "min #chi^{2} for 4C-Fit,  4 charged tracks", 300, 0, 300 );
  his1[332] = new TH1D("chi2_4C_2", "min #chi^{2} for 4C-Fit,  2 charged tracks", 300, 0, 300 );
  his1[32] = new TH1D("Etot","E_{KK#pi0} rec.", 100, 0., 3.5);
  his1[33] = new TH1D("Ptot","P_{KK#pi0} rec.", 100, 0., 2.);
  his1[34] = new TH1D("Etot_chicut","E_{KK#pi0} rec. after chi cut", 100, 0., 3.5);
  his1[35] = new TH1D("Ptot_chicut","P_{KK#pi0} rec. after chi cut", 100, 0., 2.);

  his2[51] = new TH2D("Mpm_M0_fromChi2", "M^{2}(K^{+}K^{-}) - y, M^{2}(#pi^{0}) - x from chi2", 100,0.,1.0, 100,0.,9.);
  his2[52] = new TH2D("Mp0_Mm0_fromChi2", "M^{2}(K^{+}#pi^{0}) - y, M^{2}(K^{-}#pi^{0}) - x from chi2", 100,0.,9., 100,0.,9.);
  his2[53] = new TH2D("Mpm_M0_fromChi2Mgg", "M^{2}(K^{+}K^{-}) - y, M^{2}(#pi^{0}) - xfrom chi2 & Mgg", 100,0.,1.0, 100,0.,9.);
  his2[54] = new TH2D("Mp0_Mm0_fromChi2Mgg", "M^{2}(K^{+}#pi^{0}) - y, M^{2}(K^{-}#pi^{0}) - x from chi2 & Mgg", 100,0.,9., 100,0.,9.);
  his2[55] = new TH2D("Mpm_M0_SB", "M^{2}(K^{+}K^{-}) - y, M^{2}(#pi^{0}) - x SB", 100,0.,1.0, 100,0.,9.);
  his2[56] = new TH2D("Mp0_Mm0_SB", "M^{2}(K^{+}#pi^{0}) - y, M^{2}(K^{-}#pi^{0}) - x SB", 100,0.,9., 100,0.,9.);


  // Monte Carlo histograms
  his1[101] = new TH1D("mc_type", "0 - bkg, 1 - 3pi, 2 - pi+rho", 3, -0.5, 2.5);
  his1[102] = new TH1D("mc_deccode", "decCode",258,-1.5,256.5);
  his1[103] = new TH1D("mc_deccode_final", "decCode after select",258,-1.5,256.5);
  his1[106] = new TH1D("mc_deccode_P", "decCode for Ptot_rec > 0.2",258,-1.5,256.5);

  his1[104] = new TH1D("mc_Jpsi", "MC: 0 - J/#psi, 1 - #mu from J/#psi, 2 - e from J/#psi, 3 - J/#psi + #gamma_{ISR}", 7, -3.5, 3.5);
  his1[105] = new TH1D("mc_Psi2S", "MC: 0 - #psi(2S), 1 - #mu from #psi(2S), 2 - e from #psi(2S), 3 - J/#psi#pi^{+}#pi^{-}", 7, -3.5, 3.5);
  his1[109] = new TH1D("mc_muQP_Jpsi","MC: Charged momentum (Q*P) for #mu from J/#psi", 1000, -5.0, 5.0);
  his1[110] = new TH1D("mc_muQP_Psi2S","MC: Charged momentum (Q*P) for #mu from #psi(2S)", 1000, -5.0, 5.0);
  his1[60]  = new TH1D("mc_cosTheta_Jpsi","MC: cos(#theta) for #mu^{+/-} from J/#psi", 200, -1., 1.);
  his1[61]  = new TH1D("mc_cosTheta_Psi2S","MC: cos(#theta) for #mu^{+/-} from #psi(2S)", 200, -1., 1.);
  his1[62]  = new TH1D("mc_Phi_Jpsi","MC: #phi for #mu^{+/-} from J/#psi; degrees", 720, -360., 360.);
  his1[63]  = new TH1D("mc_Phi_Psi2S","MC: #phi for #mu^{+/-} from #psi(2S); degrees", 720, -360., 360.);

  his2[60]  = new TH2D("mc_phi_vs_cosTheta_mu_ini","MC: phi vs cos(#theta) for #mu all init", 720, -360., 360., 200, -1., 1.);
  his2[61]  = new TH2D("mc_phi_vs_cosTheta_mu_fin","MC: phi vs cos(#theta) for #mu all fin", 720, -360., 360., 200, -1., 1.);
  his2[62]  = new TH2D("mc_phi_vs_cosTheta_mu_nonJpsi_ini","MC: phi vs cos(#theta) for #mu non J/#psi init", 720, -360., 360., 200, -1., 1.);
  his2[63]  = new TH2D("mc_phi_vs_cosTheta_mu_nonJpsi_fin","MC: phi vs cos(#theta) for #mu non J/#psi fin", 720, -360., 360., 200, -1., 1.);
  his2[64]  = new TH2D("mc_phi_vs_cosTheta_mu_Jpsi_ini","MC: phi vs cos(#theta) for #mu from J/#psi init", 720, -360., 360., 200, -1., 1.);
  his2[65]  = new TH2D("mc_phi_vs_cosTheta_mu_Jpsi_fin","MC: phi vs cos(#theta) for #mu from J/#psi fin", 720, -360., 360., 200, -1., 1.);
  his2[66]  = new TH2D("mc_phi_vs_cosTheta_ee_ini","MC: phi vs cos(#theta) for e from J/#psi init", 720, -360., 360., 200, -1., 1.);
  his2[67]  = new TH2D("mc_phi_vs_cosTheta_ee_fin","MC: phi vs cos(#theta) for e from J/#psi fin", 720, -360., 360., 200, -1., 1.);
  his2[68]  = new TH2D("mc_phi_vs_cosTheta_Jpsi_ini","MC: phi vs cos(#theta) for J/#psi init", 720, -360., 360., 200, -1., 1.);
  his2[69]  = new TH2D("mc_phi_vs_cosTheta_Jpsi_fin","MC: phi vs cos(#theta) for J/#psi fin", 720, -360., 360., 200, -1., 1.);
  his2[70]  = new TH2D("mc_P_vs_cosTheta_mu_fin","MC: p vs cos(#theta) for #mu all fin", 80, 0., 4., 200, -1., 1.);
  his2[71]  = new TH2D("mc_pJpsi_vs_phimu","MC: p J/#psi vs #phi #mu^{+/-} from J/#psi; #phi #mu^{+/-}; p J/#psi", 720, -360., 360., 80, 0., 4.);
  his2[72]  = new TH2D("mc_P_vs_cosTheta_mu_init","MC: p vs cos(#theta) for #mu all init", 80, 0., 4., 200, -1., 1.);
  his2[73]  = new TH2D("mc_P_vs_cosTheta_muJpsi_init","MC: p vs cos(#theta) for #mu from J/#psi init; p #mu^{+/-}; cos(#theta) #mu^{+/-}", 80, 0., 4., 200, -1., 1.);
  his2[74]  = new TH2D("mc_P_vs_cosTheta_muJpsi_fin","MC: p vs cos(#theta) for #mu from J/#psi fin; p #mu^{+/-}; cos(#theta) #mu^{+/-}", 80, 0., 4., 200, -1., 1.);
  his2[75]  = new TH2D("P_vs_cosTheta_muJpsi_rec","p vs cos(#theta) for #mu from J/#psi rec.; p #mu^{+/-}; cos(#theta) #mu^{+/-}", 80, 0., 4., 200, -1., 1.);
  his2[76]  = new TH2D("P_vs_M2mu_Jpsi_rec","p vs M_{#mu^{+}#mu^{-}} from J/#psi rec.; p J/#psi; M_{#mu^{+}#mu^{-}}", 3000, 0., 3., 200, 3.0, 3.2);
  his2[77]  = new TH2D("PJpsi_vs_Ntrk","p J/#psi rec. vs. N_{charged tracks}; p J/#psi; N_{charged tracks}", 3000, 0., 3., 11, -0.5, 10.5);
  his2[78]  = new TH2D("M2mu_vs_Ntrk","M_{#mu^{+}#mu^{-}} vs. N_{charged tracks}; M_{#mu^{+}#mu^{-}}; N_{charged tracks}", 200, 3.0, 3.2, 11, -0.5, 10.5);
  his2[79]  = new TH2D("P_vs_M2mu_Jpsi_pISRcut","p vs M_{#mu^{+}#mu^{-}} from J/#psi rec. (0.85*p_{J/#psi} cut); p J/#psi; M_{#mu^{+}#mu^{-}}", 3000, 0., 3., 200, 3.0, 3.2);
  his2[80]  = new TH2D("P_vs_CosTheta_Jpsi","p vs cos(#Theta) J/#psi; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[81]  = new TH2D("P_vs_Mpsi2S_rec","p vs M_{#psi(2S)} rec.; p #psi(2S); M_{#psi(2S)} rec.", 3000, 0., 3., 200, 3.59, 3.79);
  his2[82]  = new TH2D("Ppsi2S_vs_Ntrk","p #psi(2S) rec. vs. N_{charged tracks}; p #psi(2S); N_{charged tracks}", 3000, 0., 3., 11, -0.5, 10.5);
  his2[83]  = new TH2D("Mpsi2S_vs_Ntrk","M_{#psi(2S)} rec. vs. N_{charged tracks}; M_{#psi(2S)} rec.; N_{charged tracks}", 200, 3.59, 3.79, 11, -0.5, 10.5);
  his2[84]  = new TH2D("P_vs_Mpsi2S_pISRcut","p vs M_{#psi(2S)} rec. (0.85*p_{#psi(2S)} cut); p #psi(2S); M_{#psi(2S)} rec.", 3000, 0., 3., 200, 3.59, 3.79);
  his2[85]  = new TH2D("P_vs_CosTheta_psi2S","p vs cos(#Theta) #psi(2S); p #psi(2S); cos(#Theta) #psi(2S)", 3000, 0., 3., 200, -1., 1.);
  his2[86]  = new TH2D("mc_phi_vs_cosTheta_psi2S_ini","MC: phi vs cos(#theta) for #psi(2S) init", 720, -360., 360., 200, -1., 1.);
  his2[87]  = new TH2D("mc_phi_vs_cosTheta_psi2S_fin","MC: phi vs cos(#theta) for #psi(2S) fin", 720, -360., 360., 200, -1., 1.);
  his2[88]  = new TH2D("mc_phi_vs_cosTheta_chic1_ini","MC: phi vs cos(#theta) for #chi_{c1} init", 720, -360., 360., 200, -1., 1.);
  his2[89]  = new TH2D("mc_phi_vs_cosTheta_chic1_fin","MC: phi vs cos(#theta) for #chi_{c1} fin", 720, -360., 360., 200, -1., 1.);
  his2[90]  = new TH2D("mc_phi_vs_cosTheta_chic2_ini","MC: phi vs cos(#theta) for #chi_{c2} init", 720, -360., 360., 200, -1., 1.);
  his2[91]  = new TH2D("mc_phi_vs_cosTheta_chic2_fin","MC: phi vs cos(#theta) for #chi_{c2} fin", 720, -360., 360., 200, -1., 1.);
  his2[92]  = new TH2D("P_vs_Mchic1_rec_1C","1C: p vs M_{#chi_{c1}} rec.; p #chi_{c1}; M_{#chi_{c1}} rec.", 3000, 0., 3., 200, 3.45, 3.65);
  his2[93]  = new TH2D("Pchic1_vs_Ntrk_1C","1C: p #chi_{c1} rec. vs. N_{charged tracks}; p #chi_{c1}; N_{charged tracks}", 3000, 0., 3., 11, -0.5, 10.5);
  his2[94]  = new TH2D("Mchic1_vs_Ntrk_1C","1C: M_{#chi_{c1}} rec. vs. N_{charged tracks}; M_{#chi_{c1}} rec.; N_{charged tracks}", 200, 3.45, 3.65, 11, -0.5, 10.5);
  his2[95]  = new TH2D("P_vs_CosTheta_chic1_1C","1C: p vs cos(#Theta) #chi_{c1}; p #chi_{c1}; cos(#Theta) #chi_{c1}", 3000, 0., 3., 200, -1., 1.);
  his2[96]  = new TH2D("P_vs_Mchic2_rec_1C","1C: p vs M_{#chi_{c2}} rec.; p #chi_{c2}; M_{#chi_{c2}} rec.", 3000, 0., 3., 200, 3.45, 3.65);
  his2[97]  = new TH2D("Pchic2_vs_Ntrk_1C","1C: p #chi_{c2} rec. vs. N_{charged tracks}; p #chi_{c2}; N_{charged tracks}", 3000, 0., 3., 11, -0.5, 10.5);
  his2[98]  = new TH2D("Mchic2_vs_Ntrk_1C","1C: M_{#chi_{c2}} rec. vs. N_{charged tracks}; M_{#chi_{c2}} rec.; N_{charged tracks}", 200, 3.45, 3.65, 11, -0.5, 10.5);
  his2[99]  = new TH2D("P_vs_CosTheta_chic2_1C","1C: p vs cos(#Theta) #chi_{c2}; p #chi_{c2}; cos(#Theta) #chi_{c2}", 3000, 0., 3., 200, -1., 1.);
  his2[100]  = new TH2D("MJpsi2pi_vs_Ntrk_X","M_{J/#psi#pi^{+}#pi^{-}} rec. vs. N_{charged tracks} (M_{J/#psi#pi^{+}#pi^{-}}: 3.72-3.9 GeV); M_{J/#psi#pi^{+}#pi^{-}} rec.; N_{charged tracks}", 200, 3.7, 3.9, 11, -0.5, 10.5);
  his2[101]  = new TH2D("mc_P_vs_M2mu_Jpsi_rec","MC(J/#psi + ISR): p vs M_{#mu^{+}#mu^{-}} from J/#psi rec.; p J/#psi; M_{#mu^{+}#mu^{-}}", 3000, 0., 3., 200, 3.0, 3.2);
  his2[102]  = new TH2D("mc__PJpsi_vs_Ntrk","MC(J/#psi + ISR): p J/#psi rec. vs. N_{charged tracks}; p J/#psi; N_{charged tracks}", 3000, 0., 3., 11, -0.5, 10.5);
  his2[103]  = new TH2D("mc_P_vs_CosTheta_Jpsi_rec","MC(J/#psi + ISR): p vs cos(#Theta) J/#psi rec.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[104]  = new TH2D("mc_P_vs_CosTheta_nonJpsi","MC(J/#psi + ISR): p vs cos(#Theta) non J/#psi; p #mu^{+}#mu^{-}; cos(#Theta) #mu^{+}#mu^{-}", 3000, 0., 3., 200, -1., 1.);
  his2[105]  = new TH2D("mc_phi_vs_cosTheta_Jpsi_rec","MC: phi vs cos(#theta) for J/#psi rec.", 720, -360., 360., 200, -1., 1.);
  his2[106]  = new TH2D("mc_P_vs_CosTheta_Jpsi_ini","MC(J/#psi + ISR): p vs cos(#Theta) J/#psi ini.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[107]  = new TH2D("mc_P_vs_CosTheta_Jpsi_fin","MC(J/#psi + ISR): p vs cos(#Theta) J/#psi fin.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[108]  = new TH2D("mc_P_vs_CosTheta_Jpsi_mu_ini","MC(J/#psi + ISR): p vs cos(#Theta) J/#psi->#mu^{+}#mu^{-} ini.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[109]  = new TH2D("mc_P_vs_CosTheta_Jpsi_mu_fin","MC(J/#psi + ISR): p vs cos(#Theta) J/#psi->#mu^{+}#mu^{-} fin.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
  his2[110]  = new TH2D("mc_Prec_vs_Pmc_Jpsi","MC(J/#psi + ISR): p rec. vs p MC  J/#psi; p rec. J/#psi; p MC J/#psi", 3000, 0., 3., 3000, 0., 3.);

  his1[107] = new TH1D("cuts","cuts: 0-all, 1-good charged, 2-J/#psi->2#mu, 3 - >3 trk, 4-2-#sigma J/#psi, 5- 2-#sigma #psi, 10 - 1C", 20,-0.5,19.5);
  his1[108] = new TH1D("cuts70","cuts #phi#eta' (decCode = 70): 0-all, 1-good charged, 2-good gamma, 3 - PID, 4-fit, 5-#chi^{2}, 6-M_{gg} (eta), 7 - SBc, 8- SBxy, 9 - SBd", 10,-0.5,9.5);
  his1[207] = new TH1D("cuts_EtaPhi","cuts: 0-all,1-good charged,2-good gamma,3-muon, 4-fit, 5-#chi^{2}, 6 - SBc, 7- SBxy, 8 - SBd ", 10,-0.5,9.5);


  his1[111] = new TH1D("mc_Eg", "energy of one gamma ISR (MeV)",5000,0.,5000.);
  his1[112] = new TH1D("mc_Egs", "E(ISR) (MeV)",1000,0.,1000.);
  his1[113] = new TH1D("mc_E3pi", "E(3#pi) (MeV)",1000,0.,3100.);
  his1[114] = new TH1D("mc_Etot", "E(3#pi)+E(ISR) (MeV)",1000,0.,3500.);
  his2[111] = new TH2D("mc_Mpm_Mp0", "M^{2}(#pi^{+}#pi^{-}) vs M^{2}(#pi^{+}#pi^{0})", 100,0.,9., 100,0.,9.);
  his2[112] = new TH2D("mc_Mpm_Mm0", "M^{2}(#pi^{+}#pi^{-}) vs M^{2}(#pi^{-}#pi^{0})", 100,0.,9., 100,0.,9.);
  his2[113] = new TH2D("mc_Mp0_Mm0", "M^{2}(#pi^{+}#pi^{0}) vs M^{2}(#pi^{-}#pi^{0})", 100,0.,9., 100,0.,9.);

  his1[121] = new TH1D("mc_m_ang","angle #pi^{+/-} rec-MC track",100,-0.1,0.1);
  his1[122] = new TH1D("mc_m_dp", "delta momentum #pi^{+/-} rec-MC track",100,-0.2,0.2);
  his1[123] = new TH1D("mc_m_pdg","good/bad (1/2) match for #pi",5,-2.5,2.5);
  his1[125] = new TH1D("mc_m_pip","mc pdg for rec. #pi^{+}",999,-499.5,499.5);
  his1[126] = new TH1D("mc_m_pim","mc pdg for rec. #pi^{-}",999,-499.5,499.5);

  his1[131] = new TH1D("mc_E_p_pi","E/p for pions(mc)", 100, 0, 1.1);
  his1[132] = new TH1D("mc_E_p_e","E/p for electrons(mc)", 100, 0, 1.1);
  his1[133] = new TH1D("mc_E_p_oth","E/p for others(mc)", 100, 0, 1.1);
  his1[134] = new TH1D("mc_E_p_pi_new","NEW  E/p for pions(mc)", 100, 0, 1.1);
  his1[135] = new TH1D("mc_E_p_e_new","NEW  E/p for electrons(mc)", 100, 0, 1.1);
  his1[136] = new TH1D("mc_E_p_oth_new","NEW  E/p for others(mc)", 100, 0, 1.1);
  his1[137] = new TH1D("mc_E_p_pi_fin","FINAL  E/p for pions(mc)", 100, 0, 1.1);
  his1[138] = new TH1D("mc_E_p_e_fin","FINAL  E/p for electrons(mc)", 100, 0, 1.1);
  his1[139] = new TH1D("mc_E_p_oth_fin","FINAL  E/p for others(mc)", 100, 0, 1.1);


  his1[151] = new TH1D("mc_pidec_sum", "#pi 0,+/-1 - all, +/-10 -decayed", 41,-20.5,20.5);
  his2[152] = new TH2D("mc_pidec_2", "xy vs z for #pi decay",100,-200.,200.,100.,0.,150.);
  his1[153] = new TH1D("mc_pidec_sf", "status flag of pi",32,-0.5,31.5);
  his1[156] = new TH1D("mc_pidec_Pini", "Momentum of #pi ini.", 150,0.,1.5);
  his1[157] = new TH1D("mc_pidec_betaini", "#beta of #pi ini.", 100,0.,1.0);
  his1[158] = new TH1D("mc_pidec_Pdec", "Momentum of #pi decayed", 150,0.,1.5);
  his1[159] = new TH1D("mc_pidec_betadec", "#beta of #pi decayed", 100,0.,1.0);

  his2[114] = new TH2D("mc_depth_pi", "Depth_{#pi^{-}} vs Depth_{#pi^{+}} when #pi^{-} && #pi^{+} ", 100, -10., 60., 100, -10., 60.);
  his2[115] = new TH2D("mc_depth_mu", "Depth_{#mu^{-}} vs Depth_{#mu^{+}} when #mu^{-} || #mu^{+}", 100, -10., 60., 100, -10., 60.);
  his2[116] = new TH2D("mc_depth_others", "Depth_{-} vs Depth_{+} without #pi & #mu", 100, -10., 60., 100, -10., 60.);

  his2[117] = new TH2D("ch2_Mgg","ch2: 30-100 vs  Mgg +- (1-5)sigma",8,-0.5,7.5,5,-0.5,4.5);
  his2[118] = new TH2D("ch2_Mgg_152","ch2: 30-100 vs  Mgg +- (1-5)sigma for decCode=152",8,-0.5,7.5,5,-0.5,4.5);
  his2[119] = new TH2D("mc_chi2_Mgg_ratio", "ratio 152/all", 8,-0.5,7.5,5,-0.5,4.5);
  his2[119]->Sumw2();

  his1[140] = new TH1D("sys_chi2", "#chi^{2} for 4C-Fit", 100, 0, 100 );
  his1[141] = new TH1D("Mgg_fromChi2", "M_{#gamma#gamma}", 100, 0.45, 0.65);
  his1[182] = new TH1D("MKK_3SigmaEta", "M_{KK} (3#sigma_{#eta})", 120, 0.96, 1.08);
  his1[187] = new TH1D("MF1_3SigmaEta", "M_{f1} (3#sigma_{#eta})", 600, 1.0, 1.6);
  his1[387] = new TH1D("MF1_4_3SigmaEta", "M_{#pi^{+}#pi^{-}#eta} (3#sigma_{#eta})", 600, 1.0, 1.6);
  his1[388] = new TH1D("MF1_2_3SigmaEta", "M_{#pi^{0}#pi^{0}#eta} (3#sigma_{#eta})", 600, 1.0, 1.6);

  his1[183] = new TH1D("MKK_3S", "M_{KK} 3 sigma", 120, 0.96, 1.08);
  his1[184] = new TH1D("MKK_SB", "M_{KK} SB", 120, 0.96, 1.8);
  his1[185] = new TH1D("Mgg_3S_eta", "M_{#gamma#gamma} 3 SIGMA", 120, 0.49, 0.61);
  his1[186] = new TH1D("Mgg_SB_eta", "M_{#gamma#gamma} SB", 120, 0.49, 0.61);
  his2[120] = new TH2D("eta958_phi_SBc", "M_{#pi^{+}#pi^{-}#gamma#gamma}-x vs M_{KK}-y SB-c",1200, 0.9, 1.02, 1200, 0.96, 1.08);
  his2[121] = new TH2D("eta958_phi_SBxy", "M_{#pi^{+}#pi^{-}#gamma#gamma}-x vs M_{KK}-y  SB-xy ",1200, 0.9, 1.02, 1200, 0.96, 1.08);
  his2[122] = new TH2D("eta958_phi_SBd", "M_{#pi^{+}#pi^{-}#gamma#gamma}-x vs M_{KK}-y  SB-d ",1200, 0.9, 1.02, 1200, 0.96, 1.08);

  his1[180] = new TH1D("Mgg_3S_pi0", "M_{#gamma#gamma} 3 SIGMA", 100, 0.090, 0.190);
  his1[181] = new TH1D("Mgg_SB_pi0", "M_{#gamma#gamma} SB", 100, 0.090, 0.190);

  his1[35] = new TH1D("final_Etot","E_{KK#pi0} final",50,2.5,3.5);
  his1[36] = new TH1D("final_ Ptot","P_{KK#pi0} final", 50,0.,1.);
  his1[142] = new TH1D("final_cos_theta_kp","cos(#theta_{K^{+}})", 40, -1., 1.);
  his1[143] = new TH1D("final_cos_theta_km","cos(#theta_{K^{-}})", 40, -1., 1.);
  his1[144] = new TH1D("final_cos_theta_pi0","cos(#theta_{#pi^{0}})", 40, -1., 1.);
  his1[145] = new TH1D("final_Pkp","Final P_{K^{+}} rec.", 40, 0., 2.);
  his1[146] = new TH1D("final_Pkm","Final P_{K^{-}} rec.", 40, 0., 2.);
  his1[147] = new TH1D("final_Ppi0","Final P_{#pi^{0}} rec.", 40, 0., 2.);

  his1[160] = new TH1D("mc_Pmu_all_ini", "Momentum of #mu all init", 80, 0., 4.);
  his1[161] = new TH1D("mc_Pmu_all_fin", "Momentum of #mu all final", 80, 0., 4.);
  his1[162] = new TH1D("mc_Pmu_nonJpsi_ini", "Momentum of #mu non J/#psi init", 80, 0., 4.);
  his1[163] = new TH1D("mc_Pmu_nonJpsi_fin", "Momentum of #mu non J/#psi final", 80, 0., 4.);
  his1[164] = new TH1D("mc__Pmu_Jpsi_ini", "Momentum of #mu from J/#psi init", 80, 0., 4.);
  his1[165] = new TH1D("mc__Pmu_Jpsi_fin", "Momentum of #mu from J/#psi final", 80, 0., 4.);

  his1[174] = new TH1D("E_p_mu_nonJpsi","#frac{Energy}{p} #mu non J/#psi", 100, 0, 1.1);
  his1[175] = new TH1D("E_p_mu_Jpsi","#frac{Energy}{p} #mu from J/#psi", 100, 0, 1.1);
  his1[176] = new TH1D("E_p_e_Jpsi","#frac{Energy}{p} e from J/#psi", 100, 0, 1.1);

  his1[500] = new TH1D("mc__Pee_Jpsi_ini", "Momentum of e from J/#psi init", 80, 0., 4.);
  his1[501] = new TH1D("mc__Pee_Jpsi_fin", "Momentum of e from J/#psi final", 80, 0., 4.);
  his1[502] = new TH1D("mc__PJpsi_ini", "Momentum of J/#psi init", 80, 0., 4.);
  his1[503] = new TH1D("mc__PJpsi_fin", "Momentum of J/#psi final", 80, 0., 4.);
  his1[504] = new TH1D("mc_Emu_nonJpsi_ini", "Energy of #mu non J/#psi init", 80, 0., 4.);
  his1[505] = new TH1D("mc_Emu_nonJpsi_fin", "Energy of #mu non J/#psi final", 80, 0., 4.);
  his1[506] = new TH1D("mc__Emu_Jpsi_ini", "Energy of #mu from J/#psi init", 80, 0., 4.);
  his1[507] = new TH1D("mc__Emu_Jpsi_fin", "Energy of #mu from J/#psi final", 80, 0., 4.);
  his1[508] = new TH1D("mc__Eee_Jpsi_ini", "Energy of e from J/#psi init", 80, 0., 4.);
  his1[509] = new TH1D("mc__Eee_Jpsi_fin", "Energy of e from J/#psi final", 80, 0., 4.);
  his1[510] = new TH1D("mc_Petot","MC: E_{tot}",180, 3.7, 5.5);
  his1[511] = new TH1D("mc_Pxtot","MC: P_{x}^{tot}", 200,-100.,100.);
  his1[512] = new TH1D("mc_Pytot","MC: P_{y}^{tot}", 200,-100.,100.);
  his1[513] = new TH1D("mc_Pztot","MC: P_{z}^{tot}", 200,-100.,100.);
  his1[514] = new TH1D("PJpsi_rec", "Momentum of J/#psi rec.", 80, 0., 4.);
  his1[515] = new TH1D("iMC_JpsiISR_Prec", "iMC(J/#psi+ISR): Momentum of J/#psi rec. ", 80, 0., 4.);

  his1[522] = new TH1D("mc__Ppsi2S_ini", "Momentum of #psi(2S) init", 80, 0., 4.);
  his1[523] = new TH1D("mc__Ppsi2S_fin", "Momentum of #psi(2S) final", 80, 0., 4.);
  his1[524] = new TH1D("Ppsi2S_rec", "Momentum of #psi(2S) rec.", 80, 0., 4.);
  his1[525] = new TH1D("M2mu_1C", "M_{#mu^{+}#mu^{-}} after 1C-Fit", 200, 3.0, 3.2);
  his1[526] = new TH1D("M_Jpsi2pi_1C", "1C: M_{J/#psi#pi^{+}#pi^{-}} for #psi(2S)", 2000, 3., 5.);
  his1[527] = new TH1D("MJpsiGamma_1C", "1C: M_{J/#psi#gamma}", 1700, 2.8, 4.5);
  his1[528] = new TH1D("Ppsi2S_rec_1C", "1C: Momentum of #psi(2S) rec.", 80, 0., 4.);
  his1[529] = new TH1D("M_Jpsi2k_1C", "1C: M_{J/#psiK^{+}K^{-}} if #]i<-> K from M_{J/#psi#pi^{+}#pi^{-}} 3.72-3.9 GeV", 2000, 3., 5.);
  his1[530] = new TH1D("mc_PJpsiAll_rec", "MC (all J/#psi): Momentum of J/#psi rec.", 80, 0., 4.);
  his1[531] = new TH1D("mc_Jpsi_phi_ini", "MC: phi for  J/#psi init.", 720, -360., 360.);
  his1[532] = new TH1D("mc_Jpsi_phi_fin", "MC: phi for  J/#psi fin.", 720, -360., 360.);
  his1[533] = new TH1D("mc_Jpsi_phi_rec", "MC: phi for  J/#psi rec.", 720, -360., 360.);
  his1[534] = new TH1D("mc_Jpsi_cosTheta_ini", "MC: cos(#theta) for  J/#psi init.", 200, -1., 1.);
  his1[535] = new TH1D("mc_Jpsi_cosTheta_fin", "MC: cos(#theta) for  J/#psi fin.", 200, -1., 1.);
  his1[536] = new TH1D("mc_Jpsi_cosTheta_rec", "MC: cos(#theta) for  J/#psi rec.", 200, -1., 1.);
  his1[537] = new TH1D("mc_psi2S_phi_rec", "MC: phi for  #psi(2S) rec.", 720, -360., 360.);
  his1[538] = new TH1D("mc_psi2S_cosTheta_rec", "MC: cos(#theta) for  #psi(2S) rec.", 200, -1., 1.);
  his1[539] = new TH1D("mc_chic1_phi_rec", "MC: phi for  #chi_{c1} rec.", 720, -360., 360.);
  his1[540] = new TH1D("mc_chic1_cosTheta_rec", "MC: cos(#theta) for  #chi_{c1} rec.", 200, -1., 1.);
  his1[541] = new TH1D("mc_chic2_phi_rec", "MC: phi for  #chi_{c2} rec.", 720, -360., 360.);
  his1[542] = new TH1D("mc_chic2_cosTheta_rec", "MC: cos(#theta) for  #chi_{c2} rec.", 200, -1., 1.);
  his1[545] = new TH1D("mc_chic1_phi_ini", "MC: phi for  #chi_{c1} init.", 720, -360., 360.);
  his1[546] = new TH1D("mc_chic1_cosTheta_ini", "MC: cos(#theta) for  #chi_{c1} init.", 200, -1., 1.);
  his1[547] = new TH1D("mc_chic2_phi_ini", "MC: phi for  #chi_{c2} init.", 720, -360., 360.);
  his1[548] = new TH1D("mc_chic2_cosTheta_ini", "MC: cos(#theta) for  #chi_{c2} init.", 200, -1., 1.);
  his1[549] = new TH1D("mc_E_Z0","iMC(J/#psi + ISR): E_{Z^{0}} ", 500, 0., 5.);
  his1[550] = new TH1D("PJpsi_rec_withoutISRcut", "Momentum of J/#psi rec. without ISR cut", 80, 0., 4.);
  his1[551] = new TH1D("cosTheta_mis","(J/#psi#gamma_{max})_{missing}: cos(#theta)", 200, -1., 1.);
  his1[552] = new TH1D("M2_mis_up95", "M^{2}-missing |cos(#theta)| > 0.95", 2000, -1., 1.);
  his1[553] = new TH1D("M2_mis_bel95", "M^{2}-missing |cos(#theta)| < 0.95", 2000, -1., 1.);
  his1[554] = new TH1D("Jpsi2g_mis_M2", "(J/#psi2#gamma)_{missing}: M^{2}", 2000, -1., 1.);
  his1[555] = new TH1D("GammaPair_M", "2#gamma M, GeV", 2000, 0., 2.);
  his1[556] = new TH1D("M_Jpsi2pi_1C_ISR", "1C: M_{J/#psi#pi^{+}#pi^{-}} for ISR #psi(2S)", 2000, 3., 5.);
  his1[557] = new TH1D("MllmuEP_ISR", "M_{#mu^{+}#mu^{-}}  for ISR J/#psi", 120, 2.8, 4.0);
  his1[558] = new TH1D("M_Jpsi2gamma_chic", "1C: M_{J/#psi2#gamma} for #psi'#rightarrow #gamma#chi_{c1,2}", 2000, 3., 5.);
  his1[559] = new TH1D("M_2gamma_psichi", "1C: M_{2#gamma} for #psi'#rightarrow #gamma#chi_{c1,2}", 1000, 0., 1.);
  his1[560] = new TH1D("M_Jpsi2gamma_pi0", "1C: M_{J/#psi2#gamma} with #pi^{0}", 2000, 3., 5.);
  his1[561] = new TH1D("M_Jpsi2gamma_nopi0", "1C: M_{J/#psi2#gamma} without #pi^{0}", 2000, 3., 5.);
  his1[562] = new TH1D("M_2gamma_psiONLY", "1C: M_{2#gamma} for #psi' only", 1000, 0., 1.);
  his1[563] = new TH1D("MJpsiGamma_psi", "1C: M_{J/#psi#gamma} for psi'", 1700, 2.8, 4.5);
  his1[564] = new TH1D("MlleEP_ISR", "M_{e^{+}e^{-}}  for ISR J/#psi", 120, 2.8, 4.0);






  his1[166] = new TH1D("mc_pdg_Ppip","Final P_{#pi^{+}} rec.", 40, 0., 2.);
  his1[167] = new TH1D("mc_pdg_Ppim","Final P_{#pi^{-}} rec.", 40, 0., 2.);
  his1[168] = new TH1D("mc_pdg_Pmup","Final P_{#mu^{+}} rec.", 40, 0., 2.);
  his1[169] = new TH1D("mc_pdg_Pmum","Final P_{#mu^{-}} rec.", 40, 0., 2.);
  his1[170] = new TH1D("mc_pdg_Pep","Final P_{#e^{+}} rec.", 40, 0., 2.);
  his1[171] = new TH1D("mc_pdg_Pem","Final P_{#e^{-}} rec.", 40, 0., 2.);
  his1[172] = new TH1D("mc_pdg_Potherp","Final P_other_{+} rec.", 40, 0., 2.);
  his1[173] = new TH1D("mc_pdg_Potherm","Final P_other_{-} rec.", 40, 0., 2.);

  his1[208] = new TH1D("NvarMW","N good variants for Mass-Windows  in event", 5, -0.5, 4.5);
  his1[209] = new TH1D("nEvents_chisq","N of events (#chi^{2} < '90' - 1, '80' - 2, '70' - 3, '60' - 4, '50'- 5)", 11, -0.5, 10.5);
  his1[210] = new TH1D("mc_deccode90", "decCode if #chi^{2} < 90",258,-1.5,256.5);
  his1[211] = new TH1D("mc_deccode80", "decCode if #chi^{2} < 80",258,-1.5,256.5);
  his1[212] = new TH1D("mc_deccode70", "decCode if #chi^{2} < 70",258,-1.5,256.5);
  his1[213] = new TH1D("mc_deccode60", "decCode if #chi^{2} < 60",258,-1.5,256.5);
  his1[214] = new TH1D("mc_deccode50", "decCode if #chi^{2} < 50",258,-1.5,256.5);

  //ISR
  his1[251] = new TH1D("a4C_chi2","After 4C chi^2", 200,0.,200.);
  his1[252] = new TH1D("a4C_Eisr","After 4C Eisr", 1000, 2.2, 3.2);
  his1[253] = new TH1D("a4C_PZisr","After 4C PZisr", 1000, -3.2, 3.2);
  his1[254] = new TH1D("a4C_PTisr","After 4C PTisr", 200, 0., 0.1);
  his1[255] = new TH1D("a4C_Etot","After 4C Etot", 1000, 3.0, 3.2);
  his1[256] = new TH1D("a4C_Ptot","After 4C Ptot", 1000, 0., 0.5);
  his1[257] = new TH1D("a4C_Esum","After 4C Etot+Eisr", 1000, 3.0, 3.2);
  his1[258] = new TH1D("a4C_Psum","After 4C Ptot+Pisr", 1000, 0., 0.1);
  his2[251] = new TH2D("a4C_Et_Ei","After 4C Etot+Eisr vs Eisr", 1000, 0., 1.,1000, 3.0,3.2);

  his1[260] = new TH1D("Pz_mis", "Pz-missing", 300,-1.5, 1.5);
  his1[261] = new TH1D("Pt_mis", "Pt-missing", 150, 0., 1.5);
  his1[262] = new TH1D("M_mis", "M-missing", 3000, -1.5, 1.5);
  his1[263] = new TH1D("mc_Pz_mis_70", "Pz-missing, iMC signal", 300,-1.5, 1.5);
  his1[264] = new TH1D("mc_Pt_mis_70", "Pt-missing, iMC signal", 150, 0., 1.5);
  his1[265] = new TH1D("mc_M_mis_70", "M-missing, iMC signal", 2000, -1., 1.);
  his1[266] = new TH1D("mc_Pz_mis_bkgr", "Pz-missing, iMC bkgr", 300,-1.5, 1.5);
  his1[267] = new TH1D("mc_Pt_mis_bkgr", "Pt-missing, iMC bkgr", 150, 0., 1.5);
  his1[268] = new TH1D("mc_M_mis_bkgr", "M-missing, iMC bkgr", 2000, -1., 1.);
  his1[269] = new TH1D("mc_deccode_Mmis008", "decCode, M-missing < 0.08",258,-1.5,256.5);
  his1[270] = new TH1D("mc_deccode_Mmis010", "decCode, M-missing < 0.1",258,-1.5,256.5);
  his1[271] = new TH1D("mc_deccode_Mmis012", "decCode, M-missing < 0.12",258,-1.5,256.5);
  his1[272] = new TH1D("mc_deccode_Mmis006", "decCode, M-missing < 0.06",258,-1.5,256.5);
  his1[273] = new TH1D("mc_deccode_Mmis005", "decCode, M-missing < 0.05",258,-1.5,256.5);
  his1[274] = new TH1D("mc_deccode_Mmis007", "decCode, M-missing < 0.07",258,-1.5,256.5);
  his1[275] = new TH1D("Nevents_Mmis","N_{events} for M-missing cuts 5-12 (*10 MeV)", 11, 2.5, 13.5);
  his1[276] = new TH1D("M2_mis", "M^{2}-missing", 3000, -1.5, 1.5);
  his1[277] = new TH1D("M2_mis_cut", "M^{2}-missing, Pz-cut", 3000, -1.5, 1.5);
  his1[278] = new TH1D("mc_M2_mis_bkgr", "M^{2}-missing, iMC bkgr", 2000, -1., 1.);

  his1[300] = new TH1D("MlleEP_Theta", "M_{e^{+}e^{-}} selected by E/p (cos(#Theta^{+,-}) > -0.98)", 120, 2.8, 4.0);


   his2[300] = new TH2D("Mpm0_Mpm1", "M_0^(trk{+}trk^{-}) vs M_1^{2}(trk^{+}trk^{-})", 1000, 2.5, 3.5, 1000, 2.5, 3.5);
   his2[301] = new TH2D("M2pm0_M2pm1", "M2_0^(trk{+}trk^{-}) vs M2_1^{2}(trk^{+}trk^{-})", 1000, 8.5, 9.5, 1000, 8.5, 9.5);
   his2[302]  = new TH2D("eMC_p_vs_cosTheta_Jpsi_ini","MC: p vs cos(#Theta) J/#psi init.; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
   his2[303]  = new TH2D("eMC_p_vs_cosTheta_Jpsi_fin","MC: p vs cos(#Theta) J/#psi fin. ; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
   his2[304]  = new TH2D("eMC_p_vs_cosTheta_Jpsi_aft","MC: p vs cos(#Theta) J/#psi after ISR cuts ; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);



   his1[400] = new TH1D("Mgg", "4C-fit: M_{#gamma#gamma}(for #pi_0)", 100, 0.080, 0.18);
   his1[401] = new TH1D("Mlep", "4C-fit: M_{lep+lep-}(for J/#psi)", 100, 3.050, 3.150);
   his1[402] = new TH1D("MJpsiWind", "M_{trk+trk-}(for J/#psi)", 200, 3.087, 3.107);
   his1[403] = new TH1D("Mmus", "M_{#mu^{+}#mu^{-}}", 500, 0., 5.);
   his1[404] = new TH1D("M4", "M_{tot}  (PID: one #mu^{+/-})", 150, 3., 4.5);
   his1[405] = new TH1D("M4all", "M_{tot}  (PID: at least one #mu^{+/-})", 150, 3., 4.5);
   his1[406] = new TH1D("Nmucomb","N_{mu comb.} in event", 5, -0.5, 4.5);
   his1[407] = new TH1D("E4Jpsi", "E_{tot} (PID: one #mu^{+/-})", 150, 3., 4.5);
   his1[408] = new TH1D("E4allJpsi", "E_{tot} (PID: at least one #mu^{+/-})", 150, 3., 4.5);
   his1[409] = new TH1D("M4ee", "M_{tot}  (PID: at least one e^{+/-})", 150, 3., 4.5);
   his1[410] = new TH1D("E4eeJpsi", "E_{tot} (PID: at least one e^{+/-})", 150, 3., 4.5);

   his1[411] = new TH1D("PmuJpsi","Momentum P for selected #mu pair", 250, 0., 2.5);
   his1[412] = new TH1D("EmuJpsi","Energy(EMC) for selected #mu pair", 450, 0., 4.5);
   his1[413] = new TH1D("Mes", "M_{e^{+}e^{-}}", 500, 0., 5.);
   his1[414] = new TH1D("PeJpsi","Momentum P for selected e^{+}e^{-}", 250, 0, 2.5);
   his1[415] = new TH1D("EeJpsi","Energy(EMC) for selected e^{+}e^{-}", 450, 0., 4.5);
   his1[416] = new TH1D("Necomb","N_{ee comb.} in event", 5, -0.5, 4.5);
   his1[417] = new TH1D("EoverPall","#frac{Energy_{EMC}}{p_{MDC}} for all trks", 100, 0, 1.1);
   his1[418] = new TH1D("EoverPmu","#frac{Energy_{EMC}}{p_{MDC}} for selected #mu pairs", 200, -1.1, 1.1);
   his1[419] = new TH1D("M4ee08", "M_{tot}  (PID: at least one e^{+/-}, E/p > 0.8)", 150, 3., 4.5);
   his1[420] = new TH1D("M4eenon08", "M_{tot}  (PID: at least one e^{+/-}, E/p <= 0.8)", 150, 3., 4.5);
   his1[421] = new TH1D("NeEP","N_{ee comb. E/p > 0.8} in event", 5, -0.5, 4.5);
   his1[422] = new TH1D("MlleEP", "M_{e^{+}e^{-}} selected by E/p", 120, 2.8, 4.0);
   his1[423] = new TH1D("MllmuEP", "M_{#mu^{+}#mu^{-}} selected by E/p", 120, 2.8, 4.0);
   his1[424] = new TH1D("NlleEPcomb","N_{ee comb.} selected by E/p in event", 5, -0.5, 4.5);
   his1[425] = new TH1D("NllmuEPcomb","N_{#mu pairs comb.} selected in 3-#sigma", 5, -0.5, 4.5);
   his1[426] = new TH1D("cos_e_Jpsi","cos(#theta) for e pairs from J/#psi", 200, -1., 1.);
   his1[427] = new TH1D("cos_mu_Jpsi","cos(#theta) for mu pairs from J/#psi", 200, -1., 1.);
   his1[428] = new TH1D("MlleEP097", "M_{e^{+}e^{-}} selected by E/p (P<0.97 Ebeam)", 60, 2.8, 3.4);
   his1[429] = new TH1D("MllmuEP097", "M_{#mu^{+}#mu^{-}} selected by E/p (P<0.97 Ebeam)", 60, 2.8, 3.4);
   his1[430] = new TH1D("PmuJpsi097","Momentum P for selected #mu pair (P<0.97 Ebeam)", 250, 0., 2.5);
   his1[431] = new TH1D("EmuJpsi097","Energy(EMC) for selected #mu pair (P<0.97 Ebeam)", 450, 0., 4.5);
   his1[432] = new TH1D("PeJpsi097","Momentum P for selected e^{+}e^{-} (P<0.97 Ebeam)", 250, 0, 2.5);
   his1[433] = new TH1D("EeJpsi097","Energy(EMC) for selected e^{+}e^{-} (P<0.97 Ebeam)", 450, 0., 4.5);
   his1[434] = new TH1D("cos_ee097","cos(#theta) for ee (P<0.97 Ebeam)", 100, 0., 1.);
   his1[435] = new TH1D("cos_mu097","cos(#theta) for mu pairs (P<0.97 Ebeam)", 100, 0., 1.);
   his1[436] = new TH1D("NGoodHits","N_{Good hits Dedx}", 40, -0.5, 39.5);
   his1[437] = new TH1D("Etrk","Energy(EMC) for all tracks", 450, 0., 4.5);
   his1[438] = new TH1D("P4C","Momentum P for 4 tracks (after 4C)", 250, 0., 2.5);
   his1[439] = new TH1D("E4C","Energy(EMC) for 4 tracks (after 4C)", 450, 0., 4.5);
   his1[440] = new TH1D("EoverP4C","#frac{Energy_{EMC}}{p_{MDC}} for 4  tracks (after 4C)", 100, 0, 1.1);
   his1[441] = new TH1D("chi2_4C_l08", "min #chi^{2} for 4C-Fit (E/p < 0.8)", 300, 0, 300 );
   his1[442] = new TH1D("chi2_4C_Epl08_pl1", "min #chi^{2} for 4C-Fit (E/p < 0.8, p < 1.)", 300, 0, 300 );
   his1[443] = new TH1D("NpartIDDedx","N_{partID Dedx}", 10, -0.5, 9.5);
   his1[444] = new TH1D("DedxP11","Dedx, P(0.95, 1.07)", 300, 0., 3000.);
   his1[445] = new TH1D("M_mis2", "M-missing (nGoog = 2)", 4000, -2., 2.);
   his1[446] = new TH1D("M_mis4", "M-missing (nGoog = 4)", 4000, -2., 2.);
   his1[447] = new TH1D("M_mis6", "M-missing (nGoog = 6)", 4000, -2., 2.);
   his1[448] = new TH1D("M2_mis2", "M^{2}-missing (nGood = 2)", 4000, -2., 2.);
   his1[449] = new TH1D("M2_mis4", "M^{2}-missing (nGood = 4)", 4000, -2., 2.);
   his1[450] = new TH1D("M2_mis6", "M^{2}-missing (nGood = 6)", 4000, -2., 2.);
   his1[451] = new TH1D("MJpsiGG2", "M_{J/#psi#gamma#gamma} (nGood = 2)", 170, 2.8, 4.5);
   his1[452] = new TH1D("MJpsiGG4", "M_{J/#psi#gamma#gamma} (nGood = 4)", 170, 2.8, 4.5);
   his1[453] = new TH1D("MJpsiGG6", "M_{J/#psi#gamma#gamma} (nGood = 6)", 170, 2.8, 4.5);
   his1[454] = new TH1D("M_mis_Jpsieta", "M-missing for J/#psi#eta", 4000, -2., 2.);
   his1[455] = new TH1D("M_gg_eta", "M_{#gamma#gamma} for #eta", 300, 0.4, 0.7);
   his1[456] = new TH1D("M_mis_Jpsi4g", "M-missing for J/#psi if 4#gamma", 4000, -2., 2.);
   his1[457] = new TH1D("M_mis_Jpsiplus4g", "M-missing for J/#psi4#gamma", 4000, -2., 2.);
   his1[458] = new TH1D("Etot_dd","sum Energyfor 2 selected tracks (dd)", 450, 0., 4.5);
   his1[459] = new TH1D("Etot_4p","sum Energy for 4 selected tracks (4p, Kalman)", 450, 0., 4.5);
   his1[460] = new TH1D("Etot_4p_mdc","sum Energy for 4 selected tracks (4p, MDC)", 450, 0., 4.5);
   his1[461] = new TH1D("M_mis_JpsiEta", "M-missing J/#psi #eta", 4000, -2., 2.);
   his1[462] = new TH1D("Ncomb_JpsiEta","N_{#gamma#gamma} for #eta", 11, -0.5, 10.5);
   his1[463] = new TH1D("M_KK", "M_{KK} if #pi<->K from M_{J/#psi#pi^{+}#pi^{-}}  3.72-3.9 GeV", 400, 0.8, 1.2);
   his1[464] = new TH1D("M_Ks", "M_{#pi^{+}#pi^{-}} for K_{s}", 600, 0.2, 0.8);
   his1[465] = new TH1D("M_Lambda", "M_{p#pi^{-}} for #Lambda", 500, 0.9, 1.4);
   his1[466] = new TH1D("M_aLambda", "M_{pbar#pi^{+}} for #Lambda bar", 2000, 0.9, 1.4);
   his1[467] = new TH1D("svtx_chi2", "second vertex fit #chi^2 #Lambda", 300, 0., 300.);
   his1[468] = new TH1D("DecayLoverErr", "Decay length/Err for #Lambda", 500, -10., 40.);
   his1[469] = new TH1D("DecayLoverErr_3sig", "Decay length/Err within 3-#sigma #Lambda", 500, -10., 40.);
   his1[470] = new TH1D("DecayLoverErr_a", "Decay length/Err for #Lambda bar", 500, -10., 40.);
   his1[471] = new TH1D("DecayLoverErr_3sig_a", "Decay length/Err within 3-#sigma #Lambda bar", 500, -10., 40.);
   his1[472] = new TH1D("M_Lambda_Ls", "M_{p#pi^{-}} for #Lambda (Decay length/Err > 2)", 500, 0.9, 1.4);
   his1[473] = new TH1D("M_aLambda_Ls", "M_{pbar#pi^{+}} for #Lambda bar (Decay length/Err > 2)", 2000, 0.9, 1.4);
   his1[474] = new TH1D("M_LLbar", "M_{#Lambda#Lambda bar} (Decay length/Err > 2)", 2500, 2., 4.5);
   his1[475] = new TH1D("mc_deccode_pbar", "decCode (pbar, >3 chtrk) ",258,-1.5,256.5);
   his1[476] = new TH1D("mc_deccode_ppbar", "decCode (ppbar, >3 chtrk) ",258,-1.5,256.5);
   his1[477] = new TH1D("mc_deccode_2p2pi", "decCode (ppbar+#pi^{+}#pi{-}, >3 chtrk) ",258,-1.5,256.5);
   his1[478] = new TH1D("Npart","N of  pbar(0), K^{+/-} (1/2), #pi^{+/-} (3/4), e^{+} barrel (5)", 11, -0.5, 10.5);
   his1[479] = new TH1D("M_aXiPlus", "M_{\\=p#pi^{+}#pi^{+}} for \\={#Xi^{+}}", 500, 1.1, 1.6);
   his1[480] = new TH1D("M_aOmegaPlus", "M_{\\=p#pi^{+}K^{+}} for \\={#Omega^{+}}", 500, 1.5, 2.);
   his1[481] = new TH1D("M_gg_pi0eta", "M_{#gamma#gamma} for #pi^{0}, #eta", 1000, 0., 1.);
   his1[482] = new TH1D("M_eta2pi", "M_{#eta#pi^{+}#pi^{-}} for #eta'", 500, 0.7, 1.2);
   his1[483] = new TH1D("M_aSigmaMinus", "M_{\\=p#pi^{0}} for \\={#Sigma^{-}}", 500, 1., 1.5);
   his1[484] = new TH1D("M_Jpsi2pi", "M_{J/#psi#pi^{+}#pi^{-}} for #psi(2S)", 2000, 3., 5.);
   his1[485] = new TH1D("NtrkIP","N_{charged tracks} in event (IP); N_{events} with Ntrk > 3 (-1)", 12, -1.5, 10.5);
   his1[486] = new TH1D("NtrkIP_phi","N_{charged tracks} in event (IP) with  #phi", 12, -1.5, 10.5);
   his1[487] = new TH1D("NtrkIP_eta","N_{charged tracks} in event (IP) with  #eta", 12, -1.5, 10.5);
   his1[488] = new TH1D("Ngamma_eta","N_{gamma} in event with  #eta", 11, -0.5, 10.5);
   his1[489] = new TH1D("Ngamma_phi","N_{gamma} in event with  #phi", 11, -0.5, 10.5);
   his1[490] = new TH1D("M_D0", "M_{#pi^{0}#pi^{+}K^{-}} for D_{0}", 600, 1.6, 2.2);
   his1[491] = new TH1D("M_Ds", "M_{#pi^{+/-}K^{+}K^{-}} for D_{s}^{+/-}", 600, 1.7, 2.3);
   his1[492] = new TH1D("M_Dpm", "M_{#pi^{+}#pi^{+}K^{-} + cc.} for D^{+/-}", 600, 1.6, 2.2);
   his1[493] = new TH1D("M_aLambda_all", "M_{pbar#pi^{+}} for #Lambda bar (without PID)", 2000, 0.9, 1.4);
   his1[494] = new TH1D("MJpsiGamma", "M_{J/#psi#gamma}", 1700, 2.8, 4.5);
   his1[495] = new TH1D("cosQ_mup_Jpsi","cos(#theta_{#mu^{+}}) (J/#psi frame, P_{#mu} - raw)", 200, -1., 1.);
   his1[496] = new TH1D("cosQ_mum_Jpsi","cos(#theta_{#mu^{-}}) (J/#psi frame, P_{#mu} - raw)", 200, -1., 1.);
   his1[497] = new TH1D("cosQ_mup_Jpsi_1C","cos(#theta_{#mu^{+}}) (J/#psi frame, P_{#mu} from 1C)", 200, -1., 1.);
   his1[498] = new TH1D("cosQ_mum_Jpsi_1C","cos(#theta_{#mu^{-}}) (J/#psi frame, P_{#mu} from 1C)", 200, -1., 1.);
   his1[499] = new TH1D("Q_mup_Jpsi","#theta_{#mu^{+}} (J/#psi frame, P_{#mu} - raw); degrees", 180, 0., 180.);
   his1[600] = new TH1D("Q_mum_Jpsi","#theta_{#mu^{-}} (J/#psi frame, P_{#mu} - raw); degrees", 180, 0., 180.);
   his1[601] = new TH1D("Q_mup_Jpsi_1C","#theta_{#mu^{+}} (J/#psi frame, P_{#mu} from 1C); degrees", 180, 0., 180.);
   his1[602] = new TH1D("Q_mum_Jpsi_1C","#theta_{#mu^{-}} (J/#psi frame, P_{#mu} from 1C); degrees", 180, 0., 180.);
   his1[603] = new TH1D("mc_cosQ_mup_Jpsi","cos(#theta_{#mu^{+}}) (J/#psi frame, P_{#mu} - raw)", 200, -1., 1.);
   his1[604] = new TH1D("mc_Q_mup_Jpsi","#theta_{#mu^{+}} (J/#psi frame, P_{#mu} - raw); degrees", 180, 0., 180.);
   his1[605] = new TH1D("mc__Pchic1_ini", "Momentum of #chi_{c1} init", 80, 0., 4.);
   his1[606] = new TH1D("mc__Pchic1_fin", "Momentum of #chi_{c1} final", 80, 0., 4.);
   his1[607] = new TH1D("mc__Pchic2_ini", "Momentum of #chi_{c2} init", 80, 0., 4.);
   his1[608] = new TH1D("mc__Pchic2_fin", "Momentum of #chi_{c2} final", 80, 0., 4.);
   his1[609] = new TH1D("Pchic1_rec_1C", "1C: Momentum of #chi_{c1} rec.", 80, 0., 4.);
   his1[610] = new TH1D("Pchic2_rec_1C", "1C: Momentum of #chi_{c2} rec.", 80, 0., 4.);
   his1[611] = new TH1D("mc_muQP_R1","MC(J/#psi + ISR): Charged momentum (Q*P) for #mu from Region 1", 1000, -5.0, 5.0);
   his1[612] = new TH1D("mc_muQP_R2","MC(J/#psi + ISR): Charged momentum (Q*P) for #mu from Region 2", 1000, -5.0, 5.0);
   his1[613] = new TH1D("mc_muQP_R3","MC(J/#psi + ISR): Charged momentum (Q*P) for #mu from Region 3", 1000, -5.0, 5.0);
   his1[614] = new TH1D("mc_muQP_R4","MC(J/#psi + ISR): Charged momentum (Q*P) for #mu from Region 4", 1000, -5.0, 5.0);
   his1[615] = new TH1D("mc_E_ecm_mu","MC: Energy(EMC) for #mu^{+/-}", 450, 0., 4.5);
   his1[616] = new TH1D("mc_E_ecm_mu_slct","Energy(EMC) for #mu^{+/-} (selected)", 450, 0., 4.5);
   his1[617] = new TH1D("mc_chiC", "MC: 1/2 - #chi_{c1/2}, -1/-2 - J/#psi#gamma", 5, -2.5, 2.5);
   his1[618] = new TH1D("imc_Ppsi2S_rec_1C", "1C: Momentum of #psi(2S) rec. (iMC)", 80, 0., 4.);
   his1[619] = new TH1D("imc_Pchic1_rec_1C", "1C: Momentum of #chi_{c1} rec. (iMC)", 80, 0., 4.);
   his1[620] = new TH1D("imc_Pchic2_rec_1C", "1C: Momentum of #chi_{c2} rec. (iMC)", 80, 0., 4.);
   his1[621] = new TH1D("Ntrks_Jpsi","J/#psi rec.: N_{charged tracks} in event", 11, -0.5, 10.5);
   his1[622] = new TH1D("Ntrks_psi2S","#psi(2S) rec.: N_{charged tracks} in event", 11, -0.5, 10.5);
   his1[623] = new TH1D("mc_PJpsi_bins_ini", "Momentum of J/#psi 0-0.1, 0.1-0.2 ... GeV (init)", 15, -0.5, 14.5);
   his1[624] = new TH1D("mc_PJpsi_bins_fin", "Momentum of J/#psi 0-0.1, 0.1-0.2 ... GeV (fin)", 15, -0.5, 14.5);
   his1[625] = new TH1D("mc_PJpsi_0", "MC: Momentum of J/#psi without ISR cut", 80, 0., 4.);
   his1[626] = new TH1D("mc_Ppsi2S_0", "MC: Momentum of #psi(2S) without ISR cut", 80, 0., 4.);
   his1[627] = new TH1D("mc_Pchic1_0", "Momentum of #chi_{c1} without ISR cut", 80, 0., 4.);
   his1[628] = new TH1D("mc_Pchic2_0", "Momentum of #chi_{c2}  without ISR cut", 80, 0., 4.);
   his1[629] = new TH1D("Ntrks_Jpsi_SB","J/#psi rec.: N_{charged tracks} in event (SideBand)", 11, -0.5, 10.5);
   his1[630] = new TH1D("Ntrks_psi2S_SB","#psi(2S) rec.: N_{charged tracks} in event (SideBand)", 11, -0.5, 10.5);
   his1[631] = new TH1D("mc_pionsFromPsi2S", "1 - match #pi from #psi(2S), -1 - npt match", 3, -1.5, 1.5);
   his1[632] = new TH1D("PiEN","EN: #pi^{+/-} <-> N trks", 51, -50.5, 50.5);
   his1[633] = new TH1D("PiTrkEN","EN: #pi^{+/-} Track <-> N other particles", 51, -50.5, 50.5);
   his1[634] = new TH1D("mc_Eg_ISRregion", "ISR cut region: Energy of one gamma ISR (MeV)", 5000,0.,5000.);
   his1[635] = new TH1D("EoverPslct","#frac{Energy_{EMC}}{p_{MDC}} for selected #mu^{+/-}", 100, 0, 1.1);
   his1[636] = new TH1D("MllmuEoP5", "M_{#mu^{+}#mu^{-}} selected by E/p < 0.05", 120, 2.8, 4.0);
   his1[637] = new TH1D("mc_E_ecm_mu_nu","MC: Energy(EMC) for #mu^{+/-}", 450, 0., 4.5);
   his1[638] = new TH1D("p_Psi2S_ISR_slct","p_{#psi'} ISR selected", 80, 0., 4.);
   his1[639] = new TH1D("cosTheta_PsiISR_slct","cos(#theta) #psi' ISR selected", 200, -1., 1.);
   his1[640] = new TH1D("phi_PsiISR_slct","#phi  #psi' ISR selected",720, -360., 360.);
   his1[641] = new TH1D("mc_Ppsi_bins_ini", "Momentum of #psi' 0-0.1, 0.1-0.2 ... GeV (init)", 15, -0.5, 14.5);
   his1[642] = new TH1D("mc_Ppsi_bins_fin", "Momentum of #psi' 0-0.1, 0.1-0.2 ... GeV (fin)", 15, -0.5, 14.5);
   his1[643] = new TH1D("mc_psi_phi_ini", "MC: phi for  #psi' init.", 720, -360., 360.);
   his1[644] = new TH1D("mc_psi_phi_fin", "MC: phi for  #psi' fin", 720, -360., 360.);
   his1[645] = new TH1D("mc_psi_cosTheta_ini", "MC: cos(#theta) for  #psi' init.", 200, -1., 1.);
   his1[646] = new TH1D("mc_psi_cosTheta_fin", "MC: cos(#theta) for  #psi' fin", 200, -1., 1.);
   his1[647] = new TH1D("M2mutrk", "M_{#mu_{trk}^{+}#mu_{trk}^{-}}  for J/#psi", 500, 0. , 5.0);
   his1[648] = new TH1D("NpsiInt","Data: 0 - #psi'_{ISR}, 1 - #psi' X", 2,-0.5,1.5);
   his1[649] = new TH1D("NJpsiInt","Data: 0 - J/#psi'_{ISR}, 1 - J/#psi' X", 2,-0.5,1.5);
   his1[650] = new TH1D("NchicInt","Data: 0 - #chi_{c1}, 1 - #chi_{c2}", 2,-0.5,1.5);
   his1[651] = new TH1D("AngTrk2","The angle squared between the photon & the nearest charged track, degrees^{2}", 288, 0., 14400.);
   his1[652] = new TH1D("mc_gchic12", "1 - match #gamma from #chi_{c12}, -1 - npt match", 3, -1.5, 1.5);
   his1[653] = new TH1D("gTrkEN","EN: #gamma Track <-> N other particles", 51, -50.5, 50.5);
   his1[654] = new TH1D("mc_gam_dp", "delta momentum #gamma rec-MC shower",100,-0.2,0.2);
   his1[655] = new TH1D("qq_MllmuEP", "iMC: M_{#mu^{+}#mu^{-}}", 120, 2.8, 4.0);
   his1[656] = new TH1D("qq_MllmuEP_ISR", "iMC: M_{#mu^{+}#mu^{-}}  for ISR J/#psi", 120, 2.8, 4.0);
   his1[657] = new TH1D("qq_M_Jpsi2pi_1C", "iMC: M_{J/#psi#pi^{+}#pi^{-}} for #psi(2S)", 2000, 3., 5.);
   his1[658] = new TH1D("qq_M_Jpsi2pi_1C_ISR", "iMC: M_{J/#psi#pi^{+}#pi^{-}} for ISR #psi(2S)", 2000, 3., 5.);
   his1[659] = new TH1D("qq_MJpsiGamma_1C", "iMC: M_{J/#psi#gamma}", 1700, 2.8, 4.5);
   his1[660] = new TH1D("data_PJpsi_sel", "Data: Momentum of J/#psi selected", 80, 0., 4.);
   his1[661] = new TH1D("data_PtJpsi_sel", "Data: Transverse momentum of J/#psi selected", 80, 0., 4.);
   his1[662] = new TH1D("data_Jpsi_phi_sel", "Data: phi for  J/#psi selected", 720, -360., 360.);
   his1[663] = new TH1D("data_Jpsi_cosTheta_sel", "Data: cos(#theta) for  J/#psi selected", 200, -1., 1.);
   his1[664] = new TH1D("mc_PJpsi_sel", "MC: Momentum of J/#psi selected", 80, 0., 4.);
   his1[665] = new TH1D("mc_PtJpsi_sel", "MC: Transverse momentum of J/#psi selected", 80, 0., 4.);
   his1[666] = new TH1D("mc_Jpsi_phi_sel", "MC: phi for  J/#psi selected", 720, -360., 360.);
   his1[667] = new TH1D("mc_Jpsi_cosTheta_sel", "MC: cos(#theta) for  J/#psi selected", 200, -1., 1.);
   his1[668] = new TH1D("JpsiX_Ntrks","JpsiX: N_{charged tracks}", 11, -0.5, 10.5);
   his1[669] = new TH1D("JpsiX_mom_trk", "JpsiX: Momentum of charged tracks", 80, 0., 4.);
   his1[670] = new TH1D("JpsiX_cos_trk","JpsiX: cos(#theta) of charged tracks", 200, -1., 1.);
   his1[671] = new TH1D("JpsiX_Ngam","JpsiX: N_{#gamma}", 11, -0.5, 10.5);
   his1[672] = new TH1D("JpsiX_mom_gam", "JpsiX: Momentum of photons", 80, 0., 4.);
   his1[673] = new TH1D("JpsiX_cos_gam","JpsiX: cos(#theta) of photons", 200, -1., 1.);



   his2[400] = new TH2D("Lep_pVSm","N_{#mu^{+}} vs N_{#mu^{-}}", 5, 0, 4, 5, 0, 4 );
   his2[401] = new TH2D("DedxVSpm10","Dedx vs p, N_{Good hits Dedx} > 10", 6000, -3., 3. , 1000, 0., 10000);
   his2[402] = new TH2D("DedxVSpl10","Dedx vs p, N_{Good hits Dedx} <= 10",  600, -3., 3. , 1000, 0., 10000 );
   his2[403] = new TH2D("EVSp","Energy vs p for all tracks",  600, -3., 3. , 450, 0., 4.5 );
   his2[404] = new TH2D("DedxVSp0","Dedx vs p, N_{Good hits Dedx} > 10, PID = 0", 6000, -3., 3. , 1000, 0., 10000);
   his2[405] = new TH2D("DedxVSp1","Dedx vs p, N_{Good hits Dedx} > 10, PID = 1", 6000, -3., 3. , 1000, 0., 10000);
   his2[406] = new TH2D("DedxVSp2","Dedx vs p, N_{Good hits Dedx} > 10, PID = 2", 6000, -3., 3. , 1000, 0., 10000);
   his2[407] = new TH2D("DedxVSp3","Dedx vs p, N_{Good hits Dedx} > 10, PID = 3", 6000, -3., 3. , 1000, 0., 10000);
   his2[408] = new TH2D("DedxVSp4","Dedx vs p, N_{Good hits Dedx} > 10, PID = 4", 6000, -3., 3. , 1000, 0., 10000);
   his2[409] = new TH2D("MmisVSEmis_2","M_{mis} vs E_{mis} (nGood = 2, mis: P{cm} - P_{J/#psi})", 400, -2., 2., 2250, 0, 4.5);
   his2[410] = new TH2D("MmisVSEmis_4","M_{mis} vs E_{mis} (nGood = 4, mis: P{cm} - P_{J/#psi})", 400, -2., 2., 2250, 0, 4.5);
   his2[411] = new TH2D("MmisVSEmis_6","M_{mis} vs E_{mis} (nGood = 6, mis: P{cm} - P_{J/#psi})", 400, -2., 2., 2250, 0, 4.5);
   his2[412] = new TH2D("MJpsi1gVSEmisJpsi1g_2","M_{J/#psi1#gamma} vs Emis_{J/#psi1#gamma} (nGood = 2)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[413] = new TH2D("MJpsi1gVSEmisJpsi1g_4","M_{J/#psi1#gamma} vs Emis_{J/#psi1#gamma} (nGood = 4)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[414] = new TH2D("MJpsi1gVSEmisJpsi1g_6","M_{J/#psi1#gamma} vs Emis_{J/#psi1#gamma} (nGood = 6)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[415] = new TH2D("MJpsi2gVSEmisJpsi2g_2","M_{J/#psi2#gamma} vs Emis_{J/#psi2#gamma} (nGood = 2)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[416] = new TH2D("MJpsi2gVSEmisJpsi2g_4","M_{J/#psi2#gamma} vs Emis_{J/#psi2#gamma} (nGood = 4)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[417] = new TH2D("MJpsi2gVSEmisJpsi2g_6","M_{J/#psi2#gamma} vs Emis_{J/#psi2#gamma} (nGood = 6)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[418] = new TH2D("MJpsi3gVSEmisJpsi3g_2","M_{J/#psi3#gamma} vs Emis_{J/#psi3#gamma} (nGood = 2)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[419] = new TH2D("MJpsi3gVSEmisJpsi3g_4","M_{J/#psi3#gamma} vs Emis_{J/#psi3#gamma} (nGood = 4)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[420] = new TH2D("MJpsi3gVSEmisJpsi3g_6","M_{J/#psi3#gamma} vs Emis_{J/#psi3#gamma} (nGood = 6)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[421] = new TH2D("MJpsi4gVSEmisJpsi4g_2","M_{J/#psi4#gamma} vs Emis_{J/#psi4#gamma} (nGood = 2)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[422] = new TH2D("MJpsi4gVSEmisJpsi4g_4","M_{J/#psi4#gamma} vs Emis_{J/#psi4#gamma} (nGood = 4)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[423] = new TH2D("MJpsi4gVSEmisJpsi4g_6","M_{J/#psi4#gamma} vs Emis_{J/#psi4#gamma} (nGood = 6)", 200, 2.8, 4.8, 2250, 0, 4.5);
   his2[424] = new TH2D("MJpsi2gVSM2g_2","M_{J/#psi2#gamma} vs M_{2#gamma} (nGood = 2)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[425] = new TH2D("MJpsi2gVSM2g_4","M_{J/#psi2#gamma} vs M_{2#gamma} (nGood = 4)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[426] = new TH2D("MJpsi2gVSM2g_6","M_{J/#psi2#gamma} vs M_{2#gamma} (nGood = 6)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[427] = new TH2D("MJpsi3gVSM3g_2","M_{J/#psi3#gamma} vs M_{3#gamma} (nGood = 2)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[428] = new TH2D("MJpsi3gVSM3g_4","M_{J/#psi3#gamma} vs M_{3#gamma} (nGood = 4)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[429] = new TH2D("MJpsi3gVSM3g_6","M_{J/#psi3#gamma} vs M_{3#gamma} (nGood = 6)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[430] = new TH2D("MJpsi4gVSM4g_2","M_{J/#psi4#gamma} vs M_{4#gamma} (nGood = 2)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[431] = new TH2D("MJpsi4gVSM4g_4","M_{J/#psi4#gamma} vs M_{4#gamma} (nGood = 4)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[432] = new TH2D("MJpsi4gVSM4g_6","M_{J/#psi4#gamma} vs M_{4#gamma} (nGood = 6)", 200, 2.8, 4.8, 750, 0, 1.5);
   his2[433] = new TH2D("EVSp_dd","Energy vs p for dd",  600, -3., 3. , 450, 0., 4.5 );
   his2[434] = new TH2D("EVSp_4p","Energy vs p for 4p",  600, -3., 3. , 450, 0., 4.5 );
   his2[435] = new TH2D("DedxVSp_dd","Dedx vs p, N_{Good hits Dedx} > 10, PID = 0 (dd)", 6000, -3., 3. , 1000, 0., 10000);
   his2[436] = new TH2D("DedxVSp_4p","Dedx vs p, N_{Good hits Dedx} > 10, PID = 4 (4p)", 6000, -3., 3. , 1000, 0., 10000);
   his2[437] = new TH2D("M2gVSM2g_2","M_{2#gamma} vs M_{2#gamma} (nGood = 2)", 750, 0, 1.5, 750, 0, 1.5);
   his2[438] = new TH2D("M2gVSM2g_4","M_{2#gamma} vs M_{2#gamma} (nGood = 4)", 750, 0, 1.5, 750, 0, 1.5);
   his2[439] = new TH2D("M2gVSM2g_6","M_{2#gamma} vs M_{2#gamma} (nGood = 6)", 750, 0, 1.5, 750, 0, 1.5);
   his2[440] = new TH2D("MmisKK","M_{mis} vs M_{K^{+}K^{-}} (mis: P{cm} - P_{K^{+}K^{-}})", 800, -4., 4., 400, 0.8, 1.2);
   his2[441] = new TH2D("Mmis2pi","M_{mis} vs M_{#pi^{+}#pi^{-}} (mis: P{cm} - P_{#pi^{+}#pi^{-}})", 860, -4.3, 4.3, 600, 0.2, 0.8);
   his2[442] = new TH2D("MmisLambda","M_{mis} vs M_{pbar#pi^{+}} (mis: P{cm} - P_{pbar#pi^{+}})", 720, -3.6, 3.6, 2000, 0.9, 1.4);
   his2[443] = new TH2D("Mmis2gamma","M_{mis} vs M_{#gamma#gamma} (mis: P{cm} - P_{#gamma#gamma})", 920, -4.6, 4.6, 1000, 0., 1.);
   his2[444] = new TH2D("P_KK","p_{K^{+}K^{-}} vs M_{K^{+}K^{-}}", 3600, 0., 3.6, 400, 0.8, 1.2);
   his2[445] = new TH2D("P_2pi","p_{#pi^{+}#pi^{-}} vs M_{#pi^{+}#pi^{-}}", 4100, 0., 4.1, 600, 0.2, 0.8);
   his2[446] = new TH2D("P_aLambda","p_{pbar#pi^{+}} vs M_{pbar#pi^{+}}", 3600, 0., 3.6, 500, 0.9, 1.4);
   his2[447] = new TH2D("P_2gamma","p_{#gamma#gamma} vs M_{#gamma#gamma}", 4600, 0., 4.6, 1000, 0., 1.);
   his2[448] = new TH2D("Npart_phi","N of p/pbar(1), #pi^{+/-} (2), e^{+/-} barrel (3), #mu^{+/-} (4), K^{+/-} (5) in event with #phi", 11, -5.5, 5.5, 20, 0.5, 19.5);
   his2[449] = new TH2D("M2pivsMkk","M_{#pi^{+}#pi^{-}} vs M_{K^{+}K^{-}} (2#pi in event)", 2000, 0.2, 2.2, 400, 0.8, 1.2);
   his2[450] = new TH2D("M22pivsMkk","M_{#pi^{+}#pi^{-}} vs M_{K^{+}K^{-}} (4#pi in event)", 2000, 0.2, 2.2, 400, 0.8, 1.2);
   his2[451] = new TH2D("MkkvsMkk","M_{K^{+}K^{-}} vs M_{K^{+}K^{-}}", 400, 0.8, 1.2, 400, 0.8, 1.2);
   his2[452] = new TH2D("M2pivsM22k","M_{#pi^{+}#pi^{-}} vs M_{K^{+}K^{-}} (4K in event)", 2000, 0.2, 2.2, 400, 0.8, 1.2);
   his2[453] = new TH2D("PvsE_pbar","p vs Energy for pbar",  370, 0, 3.7 , 460, 0., 4.6 );
   his2[454] = new TH2D("PvsE_kaons","p vs Energy for kaons",  840, -4.2, 4.2 , 460, 0., 4.6 );
   his2[455] = new TH2D("E_vs_CosTheta_JpsiGammaMIS","(J/#psi#gamma_{max})_{missing}: E  vs cos(#Theta); E, GeV; cos(#Theta)", 2000., 0., 2., 200, -1., 1.);
   his2[456] = new TH2D("M2_vs_CosTheta_JpsiGammaMIS","(J/#psi#gamma_{max})_{missing}: M^{2}  vs cos(#Theta); E, GeV; cos(#Theta)", 2000., -1., 1., 200, -1., 1.);
   his2[457]  = new TH2D("P_vs_CosTheta_Jpsi_ISR_slct","p vs cos(#Theta) J/#psi ISR selected; p J/#psi; cos(#Theta) J/#psi", 3000, 0., 3., 200, -1., 1.);
   his2[458]  = new TH2D("P_vs_CosTheta_psi_ISR_slct","p vs cos(#Theta) #psi' ISR selected; p #psi'; cos(#Theta) #psi'", 3000, 0., 3., 200, -1., 1.);
   his2[459] = new TH2D("M2_psi_chic", "M_{J/#psi2#gamma}^{2} vs M_{J/#psi#gamma}^{2} (#gamma from #chi_{c}); M_{J/#psi2#gamma}^{2}, GeV^{2}; M_{J/#psi#gamma}^{2}, GeV^{2}", 700, 9., 16., 700, 9., 16.);
   his2[460] = new TH2D("M2_psi_NONchic", "M_{J/#psi2#gamma}^{2} vs M_{J/#psi#gamma}^{2} (other #gamma); M_{J/#psi2#gamma}^{2}, GeV^{2};  M_{J/#psi#gamma}^{2}, GeV^{2}", 700, 9., 16., 700, 9., 16.);
   his2[461] = new TH2D("M2_psi_2g", "M_{J/#psi2#gamma}^{2} vs M_{2#gamma}^{2}; M_{J/#psi2#gamma}^{2}, GeV^{2};  M_{2#gamma}^{2}, GeV^{2}", 700, 9., 16., 1000, 0., 1.);
   his2[462] = new TH2D("mc_E_ecm_mu_cosTheta","MC: Energy(EMC) vs. cos(#Theta) for #mu^{+/-}", 450, 0., 4.5, 200, -1., 1.);
   his2[463] = new TH2D("M_Jpsi2pi_2pi_1C", "1C: M_{J/#psi#pi^{+}#pi^{-}} vs. M_{#pi^{+}#pi^{-}} for #psi(2S); M_{J/#psi#pi^{+}#pi^{-}}, GeV^{2}; M_{#pi^{+}#pi^{-}}, GeV^{2};", 2000, 3., 5., 5000., 0., 5.);
   his2[464] = new TH2D("M_Jpsi2pi_2pi_1C_ISR", "1C: M_{J/#psi#pi^{+}#pi^{-}} vs. M_{#pi^{+}#pi^{-}} for ISR #psi(2S); M_{J/#psi#pi^{+}#pi^{-}}, GeV^{2}; M_{#pi^{+}#pi^{-}}, GeV^{2};", 2000, 3., 5., 5000., 0., 5.);
   his2[465] = new TH2D("Pz_M2_mis", "Pz vs. M^{2}  missing; Pz ; M^{2}", 4000,-2., 2., 4000, -1., 3.);
   his2[466] = new TH2D("Pz_Pt_mis", "Pz vs. Pt  missing; Pz ; Pt", 4000,-2., 2., 1500, 0., 1.5);
   his2[467] = new TH2D("Pz_M2_mis_2g", "Pz vs. M^{2}  missing; Pz ; M^{2}", 4000,-2., 2., 4000, -1., 3.);
   his2[468] = new TH2D("Pz_Pt_mis_2g", "Pz vs. Pt  missing; Pz ; Pt", 4000,-2., 2., 1500, 0., 1.5);
   his2[469]  = new TH2D("MCpsi2pi_psiISR_slct","p vs cos(#Theta) #psi' ISR selected; p #psi'; cos(#Theta) #psi'", 3000, 0., 3., 200, -1., 1.);
   his2[470] = new TH2D("Eg_chic1", "#chi_{c1}: E_{#gamma}^{rec.} vs. E_{#gamma}^{true}; E_{#gamma}^{rec.}, MeV; E_{#gamma}^{true}, MeV", 1000,0.,1000., 1000,0.,1000.);
   his2[471] = new TH2D("Eg_chic2", "#chi_{c2}: E_{#gamma}^{rec.} vs. E_{#gamma}^{true}; E_{#gamma}^{rec.}, MeV; E_{#gamma}^{true}, MeV", 1000,0.,1000., 1000,0.,1000.);
   his2[472] = new TH2D("mom_vs_Ntrk","JpsiX: Momentum vs N_{charged tracks}", 80, 0., 4., 11, -0.5, 10.5);
   his2[473] = new TH2D("cos_vs_Ntrk", "JpsiX: cos(#theta) vs N_{charged tracks}", 200, -1., 1., 11, -0.5, 10.5);
   his2[474] = new TH2D("mom_vs_Ng","JpsiX: Momentum vs N_{#gamma}", 80, 0., 4., 11, -0.5, 10.5);
   his2[475] = new TH2D("cos_vs_Ng", "JpsiX: cos(#theta) vs N_{#gamma}", 200, -1., 1., 11, -0.5, 10.5);


  // register in selector to save in given directory
  VecObj his1o(his1.begin(),his1.end());
  selector->RegInDir(his1o,"SelectJpsiX");

  VecObj his2o(his2.begin(),his2.end());
  selector->RegInDir(his2o,"SelectJpsiX");
}




//-----------------------------------------------------------------------------
  static Hep3Vector getVertexOrigin(int runNo, HepSymMatrix& ErrMatrixOrigin)
//-----------------------------------------------------------------------------
{
  static int save_runNo = 0;
  static Hep3Vector xorigin;

  if ( runNo == save_runNo ) return xorigin;

  // update vertex for new run
  xorigin.set(0.,0.,0.);
  VertexDbSvc* vtxsvc = VertexDbSvc::instance();


  if(abs(runNo) > 65200){

    //3-09-2021
    vtxsvc->SetBossVer("7.0.7");
  }else if ( abs(runNo) > 63000 && abs(runNo) < 65100 ){

    //17-08-2021
    vtxsvc->SetBossVer("7.0.6");
  }else if ( abs(runNo) > 59000 && abs(runNo) < 63000 ){

    //6-06-20
    vtxsvc->SetBossVer("7.0.5");
  }else{

    //1-02-2019
    vtxsvc->SetBossVer("7.0.3");
  }



  //================================
  /*
  if ( abs(runNo) < 29000 ) {

    //4009peak
    vtxsvc->SetBossVer("6.6.4");

  }else if( abs(runNo)  > 40000){


    //new data
    vtxsvc->SetBossVer("7.0.3");


  }else{

    vtxsvc->SetBossVer("6.6.4.p01");
  }
  */
  //=========================================







  vtxsvc->handle(runNo);

  if ( vtxsvc->isVertexValid() ) {
    double* dbv = vtxsvc->PrimaryVertex();
    xorigin.set(dbv[0],dbv[1],dbv[2]);
    double* dbvs = vtxsvc->SigmaPrimaryVertex();
    ErrMatrixOrigin[0][0] = dbvs[0];
    ErrMatrixOrigin[1][1] = dbvs[1];
    ErrMatrixOrigin[2][2] = dbvs[2];
    // cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
  } else {
    cout << "Cannot obtain vertex information for run#" << runNo << endl;
    ERROR_WARNING++;
  }

  save_runNo = runNo;
  return xorigin;
}

//-----------------------------------------------------------------------------
static void SavePionsMom(string file, const Select& S3pi)
//-----------------------------------------------------------------------------
{
  static ofstream ofs;
  static bool open_file = false;
  if ( !open_file ) {
    ifstream tmp(file.c_str());
    if ( tmp ) {
       cerr << " SavePionsMom::ERROR: File already exists:" << endl
            << file << endl;
       exit(1);
    }
    ofs.open(file.c_str());
    if ( !ofs.is_open() ) {
       cerr << " SavePionsMom::ERROR: Can not open file:" << endl
            << file << endl;
       exit(1);
    }
    open_file = true;
  }

  ofs << "0 1 0 0" << endl
//      << S3pi.Pgg << endl
//      << S3pi.Ppip << endl
//      << S3pi.Ppim << endl;
      << S3pi.Pgg.px() << " " << S3pi.Pgg.py() << " "
      << S3pi.Pgg.pz() << " " << S3pi.Pgg.e() << endl
      << S3pi.Ppip.px() << " " << S3pi.Ppip.py() << " "
      << S3pi.Ppip.pz() << " " << S3pi.Ppip.e() << endl
      << S3pi.Ppim.px() << " " << S3pi.Ppim.py() << " "
      << S3pi.Ppim.pz() << " " << S3pi.Ppim.e() << endl;
}

//-----------------------------------------------------------------------------
  static void SaveInfoMC(string file, const ReadDst* selector, int region)
//-----------------------------------------------------------------------------
{
  //for file
  static ofstream ofs;
  static bool open_file = false;
  if ( !open_file ) {
    ifstream tmp(file.c_str());
    if ( tmp ) {
      cerr << " SavePionsMom::ERROR: File already exists:" << endl
           << file << endl;
      exit(1);
    }
    ofs.open(file.c_str());
    if ( !ofs.is_open() ) {
      cerr << " SavePionsMom::ERROR: Can not open file:" << endl
           << file << endl;
      exit(1);
    }
    open_file = true;
  }
  //---------------------------------------------------------------

  vector<HepLorentzVector> P_ISRgamma;
  vector<HepLorentzVector> P_Jpsi;
  vector<HepLorentzVector> P_muon;
  vector<double> cos_ISRgamma;
  vector<double> cos_Jpsi;
  vector<double> cos_muon;
  long int Jpsi_mother = -99;

  /*
  if( region == 1) { ofs << "\n\n=== New Event from Region:  p_Jpsi < 1 GeV, |cosTheta| < 0.8  ===\n" << endl;}
  if( region == 2) { ofs << "\n\n=== New Event from Region:  p_Jpsi < 1 GeV, |cosTheta| > 0.8  ===\n" << endl;}
  if( region == 3) { ofs << "\n\n=== New Event from Region:  p_Jpsi > 1 GeV, |cosTheta| < 0.8  ===\n" << endl;}
  if( region == 4) { ofs << "\n\n=== New Event from Region:  p_Jpsi > 1 GeV, |cosTheta| > 0.8  ===\n" << endl;}
  */
  ofs << "\n\n=== New Event ===\n\n" << endl;
  ofs << "Particles: ";



  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  int ISR_count = 0;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();
    ofs <<" "<< part_pdg;

    /*
    if ( part_pdg==22 ) {
      if ( part->getMother() == -99 ) { // from primary vertex
        ISR_count++;

        Hep3Vector P4_3p( part->getInitialMomentumX(),
                          part->getInitialMomentumY(),
                          part->getInitialMomentumZ() );
        double P4_E = part->getInitialMomentumE();

        HepLorentzVector P4_Vector = HepLorentzVector( P4_3p, P4_E );

        P_ISRgamma.push_back(P4_Vector);
        cos_ISRgamma.push_back(P4_Vector.cosTheta());

      }
    }

    if( part_pdg==13 || part_pdg==-13 ){

        Hep3Vector P4_3p( part->getInitialMomentumX(),
                          part->getInitialMomentumY(),
                          part->getInitialMomentumZ() );
        double P4_E = part->getInitialMomentumE();

        HepLorentzVector P4_Vector = HepLorentzVector( P4_3p, P4_E );

        P_muon.push_back(P4_Vector);
        cos_muon.push_back(P4_Vector.cosTheta());

    }//end if(mu)



    if( part_pdg==443 ){

      if ( part->getMother() != -99 ) {//not from IP

         const TMcParticle *Jpsi_mother_pointer = m_TMcEvent->getMcParticle (part->getMother());

         Jpsi_mother = Jpsi_mother_pointer->getParticleID ();
       }

        ofs <<"\n J/psi -> : ";

        Hep3Vector P4_3p( part->getInitialMomentumX(),
                          part->getInitialMomentumY(),
                          part->getInitialMomentumZ() );
        double P4_E = part->getInitialMomentumE();

        HepLorentzVector P4_Vector = HepLorentzVector( P4_3p, P4_E );


        P_Jpsi.push_back(P4_Vector);
        cos_Jpsi.push_back(P4_Vector.cosTheta());

        vector<Int_t> dref = part->getDaughters ();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          ofs <<" "<< daughter_pdg;



      }//end for(Daughters)


      ofs <<"\n Particles continue: ";
    }//end if(J/psi)
    */


     if( part_pdg==100443 ){


        ofs <<"\n psi(2S) -> : ";


        vector<Int_t> dref = part->getDaughters ();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          ofs <<" "<< daughter_pdg;

      }//end for(Daughters)


      ofs <<"\n Particles continue: ";
    }//end if(psi2S)


     // chi_c1 = 20443, chi_c2 = 445

  if( part_pdg==20443 ){


        ofs <<"\n chi_c1 -> : ";


        vector<Int_t> dref = part->getDaughters ();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          ofs <<" "<< daughter_pdg;

      }//end for(Daughters)


      ofs <<"\n Particles continue: ";
    }//end if(chi_c1)

  if( part_pdg==445 ){


        ofs <<"\n chi_c2 -> : ";


        vector<Int_t> dref = part->getDaughters ();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          ofs <<" "<< daughter_pdg;

      }//end for(Daughters)


      ofs <<"\n Particles continue: ";
    }//end if(chi_c2)


  }//end while()


  /*
  ofs <<"\n J/psi mother:  "<< Jpsi_mother;


  ofs << "\n ISR gamma: P4, cosTheta " << endl;

  for(int i = 0; i < P_ISRgamma.size(); i++){

    ofs <<"( "<<P_ISRgamma[i].px() << ", " << P_ISRgamma[i].py() << ", "
        << P_ISRgamma[i].pz() << ", " << P_ISRgamma[i].e() <<" ) "<<  cos_ISRgamma[i]<< endl;

  }//end for(ISR)

  ofs << "\n J/psi: P4, cosTheta " << endl;
  ofs <<"( "<<P_Jpsi[0].px() << ", " << P_Jpsi[0].py() << ", "
      << P_Jpsi[0].pz() << ", " << P_Jpsi[0].e() <<" ) "<<  cos_Jpsi[0]<< endl;


  ofs << "\n Muon: P4, cosTheta " << endl;

  for(int i = 0; i < P_muon.size(); i++){

    ofs <<"( "<<P_muon[i].px() << ", " << P_muon[i].py() << ", "
        << P_muon[i].pz() << ", " << P_muon[i].e() <<" ) "<<  cos_muon[i]<< endl;

  }//end for(ISR)
  */




  return;
}


//-------------------------------------------------------------------------------- Loop events ----------------------------------------------------------------


  //----------------------------------------------------------------------------
static void chic12MomentumHistMC(const ReadDst* selector, bool ini_fin_flag)                  // flag  0 - ini, 1 - fin
//----------------------------------------------------------------------------
{
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;


  //J/psi -> mu+ mu-
  if(!MCJpsiTo2mu) return;


  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    //PDG:  e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443, chi_c1 = 20443, chi_c2 = 445




    // chi_c1
    if(part_pdg==20443){

      Hep3Vector chi_3p( part->getInitialMomentumX(),
                         part->getInitialMomentumY(),
                         part->getInitialMomentumZ() );

      double chi_E = part->getInitialMomentumE();
      double chi_P = chi_3p.mag();
      HepLorentzVector chi_4Vector = HepLorentzVector( chi_3p, chi_E );

      if ( !ini_fin_flag ) {

        //ISR cut ------------------------------------------------------
        //partDaughters
        vector<Int_t> dref = part->getDaughters();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          if(daughter_pdg==443){

            Hep3Vector Jpsi_3p( daughter->getInitialMomentumX(),
                                daughter->getInitialMomentumY(),
                                daughter->getInitialMomentumZ() );

            double Jpsi_E = daughter->getInitialMomentumE();
            double Jpsi_P = Jpsi_3p.mag();
            HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );

            double pISR = (Ecms*Ecms - mJpsi*mJpsi )/( 2.*Ecms );
            if(  (Jpsi_P < (pFracISR*pISR))  &&  (fabs(Jpsi_4Vector.cosTheta()) < cosThetaISR) ){

              his1[605]->Fill( chi_P );

            }//end if(ISR)
          }//end if(Jpsi)
        }//end ISR cut -------------------------------------------------


        his1[627]->Fill( chi_P );
        his1[545]->Fill( RtoD( chi_4Vector.phi() )  );
        his1[546]->Fill( chi_4Vector.cosTheta() );
        his2[88]->Fill(  RtoD( chi_4Vector.phi() ) , chi_4Vector.cosTheta() );

      }else {
        his1[606]->Fill( chi_P );
        his2[89]->Fill(  RtoD( chi_4Vector.phi() ) , chi_4Vector.cosTheta() );
      }
    }// chi_c1

      // chi_c2
      if(part_pdg==445){

        Hep3Vector chi_3p( part->getInitialMomentumX(),
                           part->getInitialMomentumY(),
                           part->getInitialMomentumZ() );

        double chi_E = part->getInitialMomentumE();
        double chi_P = chi_3p.mag();
        HepLorentzVector chi_4Vector = HepLorentzVector( chi_3p, chi_E );

        if ( !ini_fin_flag ) {

          //ISR cut ------------------------------------------------------
          //partDaughters
          vector<Int_t> dref = part->getDaughters();
          for (unsigned int i = 0; i < dref.size(); i++) {
            const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
            int daughter_pdg = daughter->getParticleID ();

            if(daughter_pdg==443){

              Hep3Vector Jpsi_3p( daughter->getInitialMomentumX(),
                                  daughter->getInitialMomentumY(),
                                  daughter->getInitialMomentumZ() );

              double Jpsi_E = daughter->getInitialMomentumE();
              double Jpsi_P = Jpsi_3p.mag();
              HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );

              double pISR = (Ecms*Ecms - mJpsi*mJpsi )/( 2.*Ecms );
              if(  (Jpsi_P < (pFracISR*pISR))  &&  (fabs(Jpsi_4Vector.cosTheta()) < cosThetaISR) ){

                his1[607]->Fill( chi_P );

              }//end if(ISR)
            }//end if(Jpsi)
          }//end ISR cut -------------------------------------------------


          his1[628]->Fill( chi_P );
          his1[547]->Fill( RtoD( chi_4Vector.phi() )  );
          his1[548]->Fill( chi_4Vector.cosTheta() );
          his2[90]->Fill(  RtoD( chi_4Vector.phi() ) , chi_4Vector.cosTheta() );

        }else {
          his1[608]->Fill( chi_P );
          his2[91]->Fill(  RtoD( chi_4Vector.phi() ) , chi_4Vector.cosTheta() );
        }
      }// chi_c2

    //=========================================


  }//end while()

  return;
}



//----------------------------------------------------------------------------
static void Psi2SMomentumHistMC(const ReadDst* selector, bool ini_fin_flag)                  // flag  0 - ini, 1 - fin
//----------------------------------------------------------------------------
{
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    //PDG:  e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443




    // psi(2S)
    if(part_pdg==100443){

      Hep3Vector psi2S_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

      double psi2S_E = part->getInitialMomentumE();
      double psi2S_P = psi2S_3p.mag();
      double psi2S_EoverP = psi2S_E/psi2S_P;
      HepLorentzVector psi2S_4Vector = HepLorentzVector( psi2S_3p, psi2S_E );

      if ( !ini_fin_flag ) {


        //ISR cut ------------------------------------------------------
        //partDaughters
        vector<Int_t> dref = part->getDaughters();
        for (unsigned int i = 0; i < dref.size(); i++) {
          const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
          int daughter_pdg = daughter->getParticleID ();

          if(daughter_pdg==443){

            Hep3Vector Jpsi_3p( daughter->getInitialMomentumX(),
                                daughter->getInitialMomentumY(),
                                daughter->getInitialMomentumZ() );

            double Jpsi_E = daughter->getInitialMomentumE();
            double Jpsi_P = Jpsi_3p.mag();
            HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );

            double pISR = (Ecms*Ecms - mJpsi*mJpsi )/( 2.*Ecms );
            if(  (Jpsi_P < (pFracISR*pISR))  &&  (fabs(Jpsi_4Vector.cosTheta()) < cosThetaISR) ){

              his1[522]->Fill( psi2S_P );

            }//end if(ISR)
          }//end if(Jpsi)
        }//end ISR cut -------------------------------------------------


         //pPsi bins init.
        for(int i = 0; i < 15; i++){
          if( (psi2S_P >= i*0.1)  &&  (psi2S_P < (i+1)*0.1) ) his1[641]->Fill(i);
        }//end for(i)


        his1[626]->Fill( psi2S_P );
        his1[643]->Fill( RtoD( psi2S_4Vector.phi() )  );
        his1[645]->Fill( psi2S_4Vector.cosTheta() );
        his2[86]->Fill(  RtoD( psi2S_4Vector.phi() ) , psi2S_4Vector.cosTheta() );

      }else {
        his1[523]->Fill( psi2S_P );
        his2[87]->Fill(  RtoD( psi2S_4Vector.phi() ) , psi2S_4Vector.cosTheta() );
        his1[644]->Fill( RtoD( psi2S_4Vector.phi() )  );
        his1[646]->Fill( psi2S_4Vector.cosTheta() );

        //pPsi bins fin.
        for(int i = 0; i < 15; i++){
          if( (psi2S_P >= i*0.1)  &&  (psi2S_P < (i+1)*0.1) ) his1[642]->Fill(i);
        }//end for(i)


      }

      /*
      //partDaughters
      vector<Int_t> dref = part->getDaughters();
      for (unsigned int i = 0; i < dref.size(); i++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
        int daughter_pdg = daughter->getParticleID ();

        //J/psi -> mu+mu-
        if(daughter_pdg==13 || daughter_pdg==-13){

          Hep3Vector mu_3p( daughter->getInitialMomentumX(),
                            daughter->getInitialMomentumY(),
                            daughter->getInitialMomentumZ() );

          double mu_E = daughter->getInitialMomentumE();
          double mu_P = mu_3p.mag();
          double mu_EoverP = mu_E/mu_P;
          HepLorentzVector mu_4Vector = HepLorentzVector( mu_3p, mu_E );

          if ( !ini_fin_flag ) {
            his1[164]->Fill( mu_P );
            his2[64]->Fill(  RtoD( mu_4Vector.phi() ) , mu_4Vector.cosTheta() );

            //other
            his1[506]->Fill( mu_E );
            //eff
            his2[73]->Fill(  mu_P , mu_4Vector.cosTheta() );

          }else {
            his1[165]->Fill( mu_P );
            his2[65]->Fill(  RtoD( mu_4Vector.phi() ) , mu_4Vector.cosTheta() );

            //other
            his1[507]->Fill( mu_E );
            //eff
            his2[74]->Fill(  mu_P , mu_4Vector.cosTheta() );
          }


        }//end if(mu+mu-)

        //J/psi -> e+e-
        if(daughter_pdg==11 || daughter_pdg==-11){

          Hep3Vector e_3p( daughter->getInitialMomentumX(),
                           daughter->getInitialMomentumY(),
                           daughter->getInitialMomentumZ() );

          double e_E = daughter->getInitialMomentumE();
          double e_P = e_3p.mag();
          double e_EoverP = e_E/e_P;
          HepLorentzVector e_4Vector = HepLorentzVector( e_3p, e_E );

          if ( !ini_fin_flag ) {
            his1[500]->Fill( e_P );
            his2[66]->Fill(  RtoD( e_4Vector.phi() ) , e_4Vector.cosTheta() );

            //other
            his1[508]->Fill( e_E );

          }else {
            his1[501]->Fill( e_P );
            his2[67]->Fill(  RtoD( e_4Vector.phi() ) , e_4Vector.cosTheta() );
            his1[509]->Fill( e_E );
          }
        }//end if(e+e-)

      }//end for(daughters)
      */

    }//end psi(2S)



    //=========================================


  }//end while()

  return;
}



//----------------------------------------------------------------------------
static void MuonsMomentumHistMC(const ReadDst* selector, bool ini_fin_flag)                  // flag  0 - ini, 1 - fin
//----------------------------------------------------------------------------
{
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    //PDG:  e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443


    //--------------------------------------
    if( !ini_fin_flag && ISRcutRegion ){
      if ( part_pdg==22 ) {

        double Eg = part->getInitialMomentumE() * 1e3; // MeV
        his1[634]->Fill(Eg);

      }
    }//--------------------------------------





    //===========================================

    //J/psi
    if(part_pdg==443){

      Hep3Vector Jpsi_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

      double Jpsi_E = part->getInitialMomentumE();
      double Jpsi_P = Jpsi_3p.mag();
      double Jpsi_EoverP = Jpsi_E/Jpsi_P;
      HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );

      if ( !ini_fin_flag ) {


        //ISR cut -----------------------------------------------------------------------
        double pISR = (Ecms*Ecms - mJpsi*mJpsi )/( 2.*Ecms );
        if(  (Jpsi_P < (pFracISR*pISR))  &&  (fabs(Jpsi_4Vector.cosTheta()) < cosThetaISR) ){

          his1[502]->Fill( Jpsi_P );
          his2[304]->Fill( Jpsi_P , Jpsi_4Vector.cosTheta() );

          //pJpsi bins init.
          for(int i = 0; i < 15; i++){
            if( (Jpsi_P >= i*0.1)  &&  (Jpsi_P < (i+1)*0.1) ) his1[623]->Fill(i);
          }//end for(i)

        }//end ISR cut -----------------------------------------------------------------


        his1[625]->Fill( Jpsi_P );
        his1[531]->Fill( RtoD( Jpsi_4Vector.phi() )  );
        his1[534]->Fill( Jpsi_4Vector.cosTheta() );
        his2[68]->Fill(  RtoD( Jpsi_4Vector.phi() ) , Jpsi_4Vector.cosTheta() );
        his2[302]->Fill( Jpsi_P , Jpsi_4Vector.cosTheta() );



      }else {
        his1[503]->Fill( Jpsi_P );
        his1[532]->Fill( RtoD( Jpsi_4Vector.phi() )  );
        his1[535]->Fill( Jpsi_4Vector.cosTheta() );
        his2[69]->Fill(  RtoD( Jpsi_4Vector.phi() ) , Jpsi_4Vector.cosTheta() );
        his2[303]->Fill( Jpsi_P , Jpsi_4Vector.cosTheta() );

        //pJpsi bins fin.
        for(int i = 0; i < 15; i++){
          if( (Jpsi_P >= i*0.1)  &&  (Jpsi_P < (i+1)*0.1) ) his1[624]->Fill(i);
        }//end for(i)


        if(JpsiSelectedForDataMCcompar){

          his1[664]->Fill( Jpsi_P );
          his1[665]->Fill( Jpsi_3p.perp() );
          his1[666]->Fill( RtoD( Jpsi_4Vector.phi() ) );
          his1[667]->Fill( Jpsi_4Vector.cosTheta() );

        }//end if(JpsiSelectedForDataMCcompar)

      }//end if(flags)


      //partDaughters
      vector<Int_t> dref = part->getDaughters();
      for (unsigned int i = 0; i < dref.size(); i++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
        int daughter_pdg = daughter->getParticleID ();

        //J/psi -> mu+mu-
        if(daughter_pdg==13 || daughter_pdg==-13){

          Hep3Vector mu_3p( daughter->getInitialMomentumX(),
                            daughter->getInitialMomentumY(),
                            daughter->getInitialMomentumZ() );

          double mu_E = daughter->getInitialMomentumE();
          double mu_P = mu_3p.mag();
          double mu_EoverP = mu_E/mu_P;
          HepLorentzVector mu_4Vector = HepLorentzVector( mu_3p, mu_E );

          if ( !ini_fin_flag ) {
            his1[164]->Fill( mu_P );
            his2[64]->Fill(  RtoD( mu_4Vector.phi() ) , mu_4Vector.cosTheta() );

            //other
            his1[175]->Fill( mu_EoverP );
            his1[506]->Fill( mu_E );
            //eff
            his2[73]->Fill(  mu_P , mu_4Vector.cosTheta() );

          }else {
            his1[165]->Fill( mu_P );
            his2[65]->Fill(  RtoD( mu_4Vector.phi() ) , mu_4Vector.cosTheta() );

            //other
            his1[507]->Fill( mu_E );
            //eff
            his2[74]->Fill(  mu_P , mu_4Vector.cosTheta() );
          }


        }//end if(mu+mu-)

        //J/psi -> e+e-
        if(daughter_pdg==11 || daughter_pdg==-11){

          Hep3Vector e_3p( daughter->getInitialMomentumX(),
                           daughter->getInitialMomentumY(),
                           daughter->getInitialMomentumZ() );

          double e_E = daughter->getInitialMomentumE();
          double e_P = e_3p.mag();
          double e_EoverP = e_E/e_P;
          HepLorentzVector e_4Vector = HepLorentzVector( e_3p, e_E );

          if ( !ini_fin_flag ) {
            his1[500]->Fill( e_P );
            his2[66]->Fill(  RtoD( e_4Vector.phi() ) , e_4Vector.cosTheta() );

            //other
            his1[176]->Fill( e_EoverP );
            his1[508]->Fill( e_E );

          }else {
            his1[501]->Fill( e_P );
            his2[67]->Fill(  RtoD( e_4Vector.phi() ) , e_4Vector.cosTheta() );
            his1[509]->Fill( e_E );
          }
        }//end if(e+e-)

      }//end for(daughters)
    }//end J/psi



    //=========================================


  }//end while()

  return;
}




//----------------------------------------------------------------------------
static void MuonsMomentumHist_iMC(const ReadDst* selector, bool ini_fin_flag)                  // flag  0 - ini, 1 - fin
//----------------------------------------------------------------------------
{
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  int mu_count=0;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    //PDG:  e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443



    //===========================================

    //Z0
    if(part_pdg==23 && iMCJpsiISR ){

      double Z_E = part->getInitialMomentumE();

      if ( !ini_fin_flag ) {

        his1[549]->Fill(Z_E);
      }
    }
    //--------------------------------



    //J/psi
    if(part_pdg==443 && iMCJpsiISR){

      Hep3Vector Jpsi_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

      double Jpsi_E = part->getInitialMomentumE();
      double Jpsi_P = Jpsi_3p.mag();
      double Jpsi_EoverP = Jpsi_E/Jpsi_P;
      HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );

      pJpsi_iMC = Jpsi_P;

      if ( !ini_fin_flag ) {
        his2[106]->Fill( Jpsi_P,  Jpsi_4Vector.cosTheta() );

      }else {
        his2[107]->Fill( Jpsi_P,  Jpsi_4Vector.cosTheta() );
      }


      //partDaughters
      vector<Int_t> dref = part->getDaughters();
      for (unsigned int i = 0; i < dref.size(); i++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
        int daughter_pdg = daughter->getParticleID ();

        //J/psi -> mu+mu-
        if(daughter_pdg==13 || daughter_pdg==-13){
          mu_count++;
        }//end if(mu+mu-)

      }//end for(daughters)


      //--------------------------------------------------------
      if(mu_count == 2){

        if ( !ini_fin_flag ) {
          his2[108]->Fill( Jpsi_P,  Jpsi_4Vector.cosTheta() );

        }else {
          his2[109]->Fill( Jpsi_P,  Jpsi_4Vector.cosTheta() );
        }
      }
      //-------------------------------------------------------

    }//end J/psi



    //=========================================


  }//end while()


  return;
}



//-----------------------------------------------------------------------------
static int ParticleTypesMC(const ReadDst* selector, Select& S2K2piG)
//-----------------------------------------------------------------------------
{
  if ( !isMC ) return -1;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return -1;

   //------------------------------------- For Event Navigator check
  // cout << "\n\n MC Particles: " <<  mcParticles->GetEntries()  << endl;
  //--------------------------------------

  int pi_count=0, rho_count=0, PiFromRho_count=0, other_count=0;

  int mupJpsi_count = 0, mumJpsi_count = 0, epJpsi_count = 0, emJpsi_count = 0, mupPsi2S_count = 0, mumPsi2S_count = 0, epPsi2S_count = 0, emPsi2S_count = 0, Jpsi_count = 0, Psi2S_count = 0, Psi3770_count = 0, Psi4040_count = 0, ISR_count = 0, JpsiPsi2S_count = 0, piPsi2S_count = 0, chiC1_count = 0, chiC2_count = 0, JpsichiC1_count = 0, GammachiC1_count = 0, JpsichiC2_count = 0, GammachiC2_count = 0;

  HepLorentzVector mc_pi[3]; // pi+, pi-, pi0
  int n_mcpi[3] = { 0, 0, 0};
  double Eg_sum = 0.;

  HepLorentzVector mup_4Vector, mum_4Vector;
  bool JpsiToMu = false;
  iMCJpsiISR = false;
  iMCJpsi = false;
  MCJpsiTo2mu = false;
  iMCpsi2SToJpsi2pi = false;
  iMCchiC1ToJpsiGamma = false;
  iMCchiC2ToJpsiGamma = false;
  iMCpsi2SInEvent = false;
  iMCchiC1InEvent = false;
  iMCchiC2InEvent = false;
  ISRcutRegion=false;

  isJpsiX=false;
  isJpsiISR=false;
  ispsiX=false;
  ispsiISR=false;
  ischic1X=false;
  ischic2X=false;


  double pJpsi = -100., phimup = -500., phimum = -500.;

  HepLorentzVector P4tot = {0., 0., 0., 0.};

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part =
                static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();
    if ( part_pdg==22 ) {
      if ( part->getMother() == -99 ) { // from primary vertex
        double Eg = part->getInitialMomentumE() * 1e3; // MeV
        his1[111]->Fill(Eg);
        Eg_sum += Eg;
        ISR_count++;
      }
      continue; // ignore gammas (ISR)
    }


    //PDG:  gamma = 22, e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443, psi(3770) = 30443, psi(4040) = 9000443,   chi_c1 = 20443, chi_c2 = 445

    if ( part_pdg== 30443) {  Psi3770_count++; }

     if ( part_pdg== 900044) {  Psi4040_count++;  }


    //P4tot
    if(part_pdg!=443){



      Hep3Vector part_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

      double part_E = part->getInitialMomentumE();
      HepLorentzVector part_4Vector = HepLorentzVector( part_3p, part_E );

      P4tot = P4tot + part_4Vector;

    }//end if(P4tot)




    //J/psi --------------------------------------------------------------------------------------------
    if(part_pdg == 443){

      Hep3Vector Jpsi_3p( part->getInitialMomentumX(),
                          part->getInitialMomentumY(),
                          part->getInitialMomentumZ() );

      double Jpsi_E = part->getInitialMomentumE();
      double Jpsi_P = Jpsi_3p.mag();
      HepLorentzVector Jpsi_4Vector = HepLorentzVector( Jpsi_3p, Jpsi_E );


      double pISR = (Ecms*Ecms - mJpsi*mJpsi )/( 2.*Ecms );
      if(  (Jpsi_P < (pFracISR*pISR))  &&  (fabs(Jpsi_4Vector.cosTheta()) < cosThetaISR) ){ ISRcutRegion=true; }

    }//end if(J/psi)  -------------------------------------------------------------------------------------




    /*
    if ( part_pdg==11 || part_pdg==-11 || part_pdg==13 || part_pdg==-13 ) {


      Hep3Vector pi_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );
      double pi_E = part->getInitialMomentumE();
      double pi_P = pi_3p.mag();
      double pi_beta = pi_P/pi_E;

      if ( part_pdg==111 ) {
        his1[151]->Fill(0);
        mc_pi[2] = HepLorentzVector( pi_3p, pi_E );
        n_mcpi[2]++;
      } else {
        his1[151]->Fill(part_pdg/211.);
        his1[156]->Fill(pi_P);
        his1[157]->Fill(pi_beta);
        if ( part_pdg==211 ) {
          mc_pi[0] = HepLorentzVector( pi_3p, pi_E );
          n_mcpi[0]++;
        }
        if ( part_pdg==-211 ) {
          mc_pi[1] = HepLorentzVector( pi_3p, pi_E );
          n_mcpi[1]++;
        }
      }

      vector<Int_t> dref = part->getDaughters ();
      if ( dref.empty() ) continue;


      if ( abs(part_pdg) == 211 ) {
        his1[151]->Fill(part_pdg/21.);
        his1[158]->Fill(pi_P);
        his1[159]->Fill(pi_beta);
        TVector3 xyz( part->getFinalPositionX(),
                      part->getFinalPositionY(),
                      part->getFinalPositionZ() );

        cout << " pi-decay: flag= " << part->getStatusFlags();
        cout << " ini (t)= "
             << part->getInitialPositionT() << endl;
        cout << " fin (x,y,z,t)= "
             << xyz.X() <<", "<<xyz.Y()<<", "<<xyz.Z()
             <<", "<< part->getFinalPositionT() << endl;
        selector->PrintMcDecayTree(part->getTrackIndex(),1);

        his2[152]->Fill(xyz.Z(),xyz.Perp());
        his1[153]->Fill(part->getStatusFlags());
      }


      continue;
    } // end if(mu/e)
    */



    //J/psi= 443, psi(2S) = 100443, chi_c1 = 20443, chi_c2 = 445
    if ( part_pdg==443 || part_pdg==100443 || part_pdg==20443 ||  part_pdg==445 ) {

      Hep3Vector part_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );


      if( part_pdg==443 ){



        pJpsi = part_3p.mag();

        his1[104]->Fill(0);
        Jpsi_count++;
      }else if( part_pdg==100443 ){
        iMCpsi2SInEvent = true;
        his1[105]->Fill(0);
        Psi2S_count++;
      }else if( part_pdg==20443 ){
        iMCchiC1InEvent = true;
        his1[617]->Fill(1);
        chiC1_count++;
      }else if( part_pdg==445 ){
        iMCchiC2InEvent = true;
        his1[617]->Fill(2);
        chiC2_count++;
      }

      vector<Int_t> dref = part->getDaughters ();
      for (unsigned int i = 0; i < dref.size(); i++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
        int daughter_pdg = daughter->getParticleID ();

        //mu 4-vector
        Hep3Vector mu_3p( daughter->getInitialMomentumX(),
                          daughter->getInitialMomentumY(),
                          daughter->getInitialMomentumZ() );
        double mu_E = daughter->getInitialMomentumE();
        double mu_P = mu_3p.mag();
        HepLorentzVector mu_4Vector = HepLorentzVector( mu_3p, mu_E );

        //---------------------------------------------------------------------------
        /*
          HepDouble m     = v4.rho();      // get length of the vector component
          HepDouble t     = v4.theta();    // get polar angle
          HepDouble cost  = v4.cosTheta(); // get cosine of polar angle
          HepDouble p     = v4.phi();      // get azimuth angle
         */
        //--------------------------------------------------------------------------

        if( part_pdg==443){
          if (daughter_pdg==11){
            JpsiToMu = false;
            his1[104]->Fill(-2);
            emJpsi_count++;
          }else if(daughter_pdg== -11){
            JpsiToMu = false;
            his1[104]->Fill(2);
            epJpsi_count++;
          }else if(daughter_pdg==13){

            MCJpsiTo2mu = true;
            phimum = RtoD( mu_4Vector.phi() );

            mum_4Vector  =  mu_4Vector;
            JpsiToMu = true;

            his1[104]->Fill(-1);
            his1[109]->Fill(-mu_P);
            his1[60]->Fill( mu_4Vector.cosTheta() );
            his1[62]->Fill( -RtoD( mu_4Vector.phi() ) );
            mumJpsi_count++;
          }else if(daughter_pdg== -13){

            phimup = RtoD( mu_4Vector.phi() );

            mup_4Vector  =  mu_4Vector;
            JpsiToMu = true;

            his1[104]->Fill(1);
            his1[109]->Fill(mu_P);
            his1[60]->Fill( mu_4Vector.cosTheta() );
            his1[62]->Fill( RtoD( mu_4Vector.phi() ) );
            mupJpsi_count++;


          }
        }//end if(J/psi -> ...)

        if( part_pdg==100443){
          if (daughter_pdg==11){
            his1[105]->Fill(-2);
            emPsi2S_count++;
          }else if(daughter_pdg== -11){
            his1[105]->Fill(2);
            epPsi2S_count++;
          }else if(daughter_pdg==13){
            his1[105]->Fill(-1);
            his1[110]->Fill(-mu_P);
            his1[61]->Fill( mu_4Vector.cosTheta() );
            his1[63]->Fill( -RtoD( mu_4Vector.phi() ) );
            mumPsi2S_count++;
          }else if(daughter_pdg== -13){
            his1[105]->Fill(1);
            his1[110]->Fill(mu_P);
            his1[61]->Fill( mu_4Vector.cosTheta() );
            his1[63]->Fill( RtoD( mu_4Vector.phi() ) );
            mupPsi2S_count++;
          }else if(daughter_pdg==443){
            JpsiPsi2S_count++;
          }else if(abs(daughter_pdg)==211){
            piPsi2S_count++;
          }
        }//end if(psi(2S) -> ...)


        if( part_pdg==20443 ){
          if (daughter_pdg==443){
            JpsichiC1_count++;
          }else if(daughter_pdg==22){
            GammachiC1_count++;
            S2K2piG.Eg_chic1 = daughter->getInitialMomentumE() * 1e3; // MeV
          }
        }


        if( part_pdg==445 ){
          if (daughter_pdg==443){
            JpsichiC2_count++;
          }else if(daughter_pdg==22){
            GammachiC2_count++;
            S2K2piG.Eg_chic2 = daughter->getInitialMomentumE() * 1e3; // MeV
          }
        }


      }//end for(daughters)

      //-----------------------------------{  J/psi polarization
      if( (part_pdg==443)  && JpsiToMu  ){

        HepLorentzVector PJpsi = mup_4Vector + mum_4Vector;
        Hep3Vector beta = PJpsi.boostVector();
        HepLorentzRotation L(-beta); // Lorentz boost
        HepLorentzVector P_mup_Jpsi, P_mum_Jpsi;
        P_mup_Jpsi = L * mup_4Vector;
        P_mum_Jpsi = L * mum_4Vector;
        his1[603]->Fill( P_mup_Jpsi.cosTheta() );
        his1[604]->Fill(  RtoD(P_mup_Jpsi.theta()) );

    }
    //----------------------------------  J/psi polarization }


      continue;
    }//end if(J/psi  ||  psi )






    other_count++;

  } //-----------------------------------------------End  while()


  his1[510]->Fill( P4tot.e() );
  his1[511]->Fill( P4tot.px() );
  his1[512]->Fill( P4tot.py() );
  his1[513]->Fill( P4tot.pz() );

  his2[71]->Fill(  phimup,  pJpsi );
  his2[71]->Fill(  phimum,  pJpsi );

  if ( selector->Verbose() ) {
    selector->PrintMCParticle(mcParticles);
    cout << " pi_count= " << pi_count << " rho_count= " << rho_count
         << " PiFromRho= " << PiFromRho_count << " other_count= "
         << other_count << endl;
  }


  /*
  // fill histograms for  pi+, pi-, pi0
  if ( n_mcpi[0] == 1 && n_mcpi[1] == 1 && n_mcpi[2] == 1 ) {
    his1[112]->Fill(Eg_sum);
    HepLorentzVector P3pi = mc_pi[0] + mc_pi[1] + mc_pi[2];
    double E3pi = P3pi.e() * 1e3;
    his1[113]->Fill( E3pi );
    his1[114]->Fill( E3pi+Eg_sum ); // Etot?
    if ( selector->Verbose() ) {
      cout << " pi+, pi-, pi0: P3pi= " << P3pi << endl;
      cout << " Eg_sum= " << Eg_sum << " Etot= " << E3pi+Eg_sum << endl;
    }
    HepLorentzVector Ppm = mc_pi[0] + mc_pi[1];
    HepLorentzVector Pp0 = mc_pi[0] + mc_pi[2];
    HepLorentzVector Pm0 = mc_pi[1] + mc_pi[2];
    his2[111]->Fill(Ppm.m2(),Pp0.m2());
    his2[112]->Fill(Ppm.m2(),Pm0.m2());
    his2[113]->Fill(Pp0.m2(),Pm0.m2());
  }

  */


  if( Jpsi_count > 0 ){ iMCJpsi = true; }

  if( ( ISR_count > 0) && (Jpsi_count > 0) ) { iMCJpsiISR = true; }

  if( (JpsiPsi2S_count==1) && (piPsi2S_count==2) ){ iMCpsi2SToJpsi2pi = true;  his1[105]->Fill(3); }

  if( (JpsichiC1_count==1) && (GammachiC1_count==1) ){  iMCchiC1ToJpsiGamma = true; his1[617]->Fill(-1);  }

   if( (JpsichiC2_count==1) && (GammachiC2_count==1) ){  iMCchiC2ToJpsiGamma = true; his1[617]->Fill(-2);  }


   //bkgr study
   if( (mupJpsi_count==1) && (mumJpsi_count==1) ){

     //JpsiX
     isJpsiX=true;

     //Jpsi ISR
     if( (ISR_count>0) && (Psi2S_count==0) && (chiC1_count==0) && (chiC2_count==0) ){ isJpsiISR=true; }

     //psiX
     if( (JpsiPsi2S_count==1) && (piPsi2S_count==2) ){ ispsiX=true; }

     //psi ISR
     if( ( ISR_count > 0) &&  (JpsiPsi2S_count==1) && (piPsi2S_count==2) ){ ispsiISR=true; }

     //chiC1 X
     if( (JpsichiC1_count==1) && (GammachiC1_count==1) ){ ischic1X=true; }

     //chiC2 X
     if( (JpsichiC2_count==1) && (GammachiC2_count==1) ){ ischic2X=true; }

   }//end if(J/psi -> mu+mu-)

  return 0;
}














//-----------------------------------------------------------------------------
static bool MatchRecMcTrks_2pi(ReadDst* selector, Select& S2K2piG)
//-----------------------------------------------------------------------------
{
  //const static double maxdp = 0.1;
  //  const static double maxangle = 0.035; // ~ 2deg.

  const static double maxdp = 0.09;
  const static double maxangle = 0.08;

  int MCpipTrkIndex = -1, MCpimTrkIndex = -1;
  int ENpipTrkID = -1, ENpimTrkID = -1;

  double p3pipmc = 0, p3pimmc = 0;



  // Match charged reconstructed tracks with MC particles
  if ( !isMC ) return false;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();


  //obtaining pion TrackIndex from psi(2S) decay
  for(int i = 0; i < mcParticles->GetEntries(); i++){
    TMcParticle * part =(TMcParticle *) mcParticles->At(i);

    long int part_pdg = part->getParticleID();

    //psi(2S)
    if( part_pdg == 100443 ){

      vector<Int_t> dref = part->getDaughters();
      for (unsigned int j = 0; j < dref.size(); j++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[j]);
        int daughter_pdg = daughter->getParticleID();

        //pi+
        if( daughter_pdg == 211){

          MCpipTrkIndex = daughter->getTrackIndex();

          Hep3Vector p3mc( daughter->getInitialMomentumX(),
                           daughter->getInitialMomentumY(),
                           daughter->getInitialMomentumZ() );
          p3pipmc = p3mc.mag();

        }//end if(pi+)

        //pi-
        if( daughter_pdg ==-211){

          MCpimTrkIndex = daughter->getTrackIndex();

          Hep3Vector p3mc( daughter->getInitialMomentumX(),
                           daughter->getInitialMomentumY(),
                           daughter->getInitialMomentumZ() );
          p3pimmc = p3mc.mag();

        }//end if(pi-)

      }//end for(j)
    }//end if( psi(2S) )
  }//end for(i)




  //Event Navigator map
  const TEvtNavigator* m_TEvtNavigator = selector->GetEventNavigator();
  std::multimap <int, int>::const_iterator it;
  const std::multimap <int, int>&  index = m_TEvtNavigator->m_mcMdcTracks;



  for(it = index.begin(); it != index.end(); it++)
    {
      // (*it).first = McParticle id; (*it).second = MdcTrack id


      //Selection
      if( (*it).first == MCpipTrkIndex ) { S2K2piG.qpipENsort.push_back((*it).second); }

      if( (*it).first == MCpimTrkIndex ) { S2K2piG.qpimENsort.push_back((*it).second); }

    }//end for(it)


  if(S2K2piG.qpipENsort.size() != 0){
    //The best pi+ track
    int tmpENpipTrkID = S2K2piG.qpipENsort[0];
    int tmppipTrkCount = 1;
    int pipTrkCount = 0;

    if(S2K2piG.qpipENsort.size() > 1){
      for(int i = 1; i < S2K2piG.qpipENsort.size(); i++){

        if(tmpENpipTrkID != S2K2piG.qpipENsort[i]){

          if( tmppipTrkCount > pipTrkCount ){ ENpipTrkID = tmpENpipTrkID; pipTrkCount = tmppipTrkCount; }

          tmpENpipTrkID = S2K2piG.qpipENsort[i];
          tmppipTrkCount = 1;
        }else{ tmppipTrkCount++; }


      }//end for(i)
    }//end if(S2K2piG.qpipENsort.size() > 1)

    //after for
    if( tmppipTrkCount > pipTrkCount ){ ENpipTrkID = tmpENpipTrkID;  pipTrkCount = tmppipTrkCount; }

     //Choose NhitCut
    if( pipTrkCount  < 5 ){  his1[631]->Fill(-1);  return false;  }

  }//end if(S2K2piG.qpipENsort.size())


  if(S2K2piG.qpimENsort.size() != 0){

    //The best pi- track
    int tmpENpimTrkID = S2K2piG.qpimENsort[0];
    int tmppimTrkCount = 1;
    int pimTrkCount = 0;

    if(S2K2piG.qpimENsort.size() > 1){
      for(int i = 1; i < S2K2piG.qpimENsort.size(); i++){

        if(tmpENpimTrkID != S2K2piG.qpimENsort[i]){

          if( tmppimTrkCount > pimTrkCount ){ ENpimTrkID = tmpENpimTrkID; pimTrkCount = tmppimTrkCount; }

          tmpENpimTrkID = S2K2piG.qpimENsort[i];
          tmppimTrkCount = 1;
        }else{ tmppimTrkCount++; }

      }//end for(i)
    }//end if(S2K2piG.qpimENsort.size() > 1)


    //after for
    if( tmppimTrkCount > pimTrkCount ){ ENpimTrkID = tmpENpimTrkID;  pimTrkCount = tmppimTrkCount; }

    //Choose NhitCut
    if( pimTrkCount < 5){  his1[631]->Fill(-1);  return false;  }

  }//end if(S2K2piG.qpimENsort.size())






  /*
  //----------------------------------------------- To check Event Navigator
  his1[632]->Fill(pipCount);
  his1[632]->Fill(-pimCount);
  his1[633]->Fill(pipTrkCount);
  his1[633]->Fill(-pimTrkCount);


  if( S2K2piG.qpipENsort.size() > 1 ){

    cout << "\n qpipENsort.size() & ID:    " <<  S2K2piG.qpipENsort.size() << "\n " << endl;
    for(int k = 0; k < S2K2piG.qpipENsort.size(); k++){cout <<S2K2piG.qpipENsort[k] << ", " << endl; }
    cout << "\n qtrksPos.size():    " <<  S2K2piG.qtrksPos.size()  << endl;


    cout << "\nEN: pi+ (p_MC - p_rec):    " << endl;

    const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
    for(int k = 0; k < S2K2piG.qpipENsort.size(); k++){
      for(int i = 0; i < S2K2piG.qtrksPos.size(); i++){
        DstEvtRecTracks* itTrkp=(DstEvtRecTracks*) evtRecTrkCol->At( S2K2piG.qtrksPos[i] );
        RecMdcTrack *mdcTrkp = itTrkp->mdcTrack();

        if( (mdcTrkp->trackId()) == S2K2piG.qpipENsort[k]){

          RecMdcKalTrack* mdcKalTrkp = itTrkp->mdcKalTrack();
          if ( !mdcKalTrkp ) {
            cout << " No Kalman track !"<< endl;
            continue; }

          mdcKalTrkp->setPidType(RecMdcKalTrack::pion);
          Hep3Vector p3pip(mdcKalTrkp->px(), mdcKalTrkp->py(), mdcKalTrkp->pz());

          cout << (p3pipmc - p3pip.mag()) << ",  ";

        }//end if(trks)

      }//end for(i)
    }//end for(k)

  }//end if(sort > 1)

  //-------------------------------------------------------------------

  */



  //MC & Rec. comparison
  for(int k = 0; k < S2K2piG.qpipPsi2StrkID.size(); k++){

    cout << "\n Rec. TrackID: pipID,  pimID  " <<  S2K2piG.qpipPsi2StrkID[k]  << ",  " <<  S2K2piG.qpimPsi2StrkID[k]  << endl;

    if( S2K2piG.qpipPsi2StrkID[k] == ENpipTrkID  && S2K2piG.qpimPsi2StrkID[k] == ENpimTrkID ){

      Hep3Vector p3pip = S2K2piG.PpipsFromPsi2S[k].vect();
      Hep3Vector p3pim = S2K2piG.PpimsFromPsi2S[k].vect();
      his1[122]->Fill(  p3pipmc - p3pip.mag()  );
      his1[122]->Fill(  p3pimmc - p3pim.mag()  );



      his1[631]->Fill(1);
      return true;
    }//end if(trks)

  }//end for(k)


  /*
? в dst есть объект TEvtNavigator

внутри лежит std::multimap <int, int> m_mcMdcTracks; // McParticle id <-> RecMdcTrack (RecMdcKalTrack) id
т.е. таблица соответствий id трека и id McParticle

так что можно сначала сделать всмопогательные массивы id<->указатель для обоих типов, чтобы сразу иметь указатель

но это для ускорения и упрощения

проблемы

1) надо понять как взять указатель из dst

2) надо получить доступ к multimap. Я посмотрел, в классе он почему-то private

3) надо проверить что он не пустой, а если пустой - перегенерить

в boss пакет Event/EventNavigator/

cd $EVENTNAVIGATORROOT на китайской ферме

реально там только EventNavigator.cxx и NavigationTestAlg.cxx представляют интерес

-----------------------------------------------------
IndexMap::const_iterator i;
  IndexMap& index = m_navigator->getMcMdcTracksIdx();
  for(i = index.begin(); i != index.end(); i++)
    {
      // i.first = McParticle id; i.second = MdcTrack id

где IndexMap это typedef  std::multimap <int, int>, чтобы пальчики не уставали печатать


Alexey 14:48
ну, естественно index=m_TEvtNavigator->m_mcMdcMcHits; в твоем случае


//  const Int_t     trackId()      const { return  m_trackId;  }  //TMdcTrack.h from  RootEventData


// Get the track id
//Int_t getTrackIndex() const {return m_trackIndex; } //TMcParticle.h  from  RootEventData

// const TMcEvent*      GetMcEvent()      const {return m_TMcEvent;} //DstFormat.h
// const TEvtNavigator* GetEventNavigator()     const {return m_TEvtNavigator;} //DstFormat.h


//from TMcEvent.cxx
const TMcParticle*  TMcEvent::getMcParticle(Int_t i) const {
        if(Int_t(i) >=m_mcParticleCol->GetEntries())
                return 0;
        return (TMcParticle*) m_mcParticleCol->At(i);
}


//typedef         BMdcTrack       RecMdcTrack;
//class BMdcTrack : public TMdcTrack



т.е. количество записей в мапе получается равно числу Rec хитов
надо для начала сосчитать идентичные пары, и выкинуть из рассмотрения те которые встречаются скажем меньше 5 раз


---------------------------------------------------


  */


  /*

  //2-pions combinations from 3-sigma psi(2S)
  for(int j = 0; j < S2K2piG.qpipsFromPsi2S.size(); j++){

    double maxdp_p = maxdp,maxdp_m = maxdp;
    double maxa_p = maxangle, maxa_m = maxangle;

    int ipip = -1, ipim = -1;
    int pip_ID = 0, pim_ID = 0;

    Hep3Vector p3pip = S2K2piG.PpipsFromPsi2S[j].vect();
    Hep3Vector p3pim = S2K2piG.PpimsFromPsi2S[j].vect();

    for(int i = 0; i < mcParticles->GetEntries(); i++){
      TMcParticle * part =(TMcParticle *) mcParticles->At(i);

      int pId = part->getParticleID();
      TParticlePDG* pdgParticle = TDatabasePDG::Instance()->GetParticle(pId);
      if( pdgParticle==0  || pdgParticle->Charge() == 0 ) {
        continue;
      }

      Hep3Vector p3mc( part->getInitialMomentumX(),
                       part->getInitialMomentumY(),
                       part->getInitialMomentumZ() );

      double angle_p = p3pip.angle(p3mc);
      double angle_m = p3pim.angle(p3mc);
      double dP_p = abs( p3mc.mag() - p3pip.mag() );
      double dP_m = abs( p3mc.mag() - p3pim.mag() );

      if ( dP_p < maxdp && angle_p < maxa_p && pId == 211) {
      //   if ( (dP_p < maxdp && angle_p < 0.035 && pId == 211) ) {

        maxdp_p = dP_p;
        maxa_p  = angle_p;
        ipip = i;
        pip_ID = pId;
      }

      if ( dP_m < maxdp && angle_m < maxa_m && pId == -211) {
    //  if ( (dP_m < maxdp && angle_m < 0.035 && pId == -211) ) {
        maxdp_m = dP_m;
        maxa_m  = angle_m;
        ipim = i;
        pim_ID = pId;
      }


    }//----------------------------------------------------End for(i)


    // check matching:
    his1[121]->Fill( maxa_p );
    his1[121]->Fill( -maxa_m );
    his1[122]->Fill( maxdp_p );
    his1[122]->Fill( -maxdp_m );


    //check matching pions
    if( (pip_ID == 211) && (pim_ID == -211) ){

      TMcParticle * part_pip =(TMcParticle *) mcParticles->At(ipip);
      TMcParticle * part_pim =(TMcParticle *) mcParticles->At(ipim);

      //pions from psi(2S)
      if( ( part_pip->getMother() != -99) && (part_pim->getMother() != -99)) { his1[631]->Fill(1); return true; }

    }//end if(pions)

  }//end for(j)

  */


  cout << "\n Event Navigator: TrkpipID,  TrkpimID  " <<  ENpipTrkID << ",  " <<  ENpimTrkID  << endl;

  his1[631]->Fill(-1);
  return false;
}







//-----------------------------------------------------------------------------
static bool MatchRecMcEMC_gamma(ReadDst* selector, Select& S2K2piG)
//-----------------------------------------------------------------------------
{

  const static double maxdp = 0.09;
  const static double maxangle = 0.08;

  int MCTrkIndex = -1;
  int ENTrkID = -1;

  double p3gmc = 0;



  // Match reconstructed tracks with MC particles
  if ( !isMC ) return false;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();


  //PDG:  gamma = 22, e- = 11, mu- = 13, J/psi = 443, psi(2S) = 100443, psi(3770) = 30443, psi(4040) = 9000443,   chi_c1 = 20443, chi_c2 = 445


  //obtaining gamma TrackIndex from chic1,2 decay
  for(int i = 0; i < mcParticles->GetEntries(); i++){
    TMcParticle * part =(TMcParticle *) mcParticles->At(i);

    long int part_pdg = part->getParticleID();

    //chic1,2
    if( (part_pdg == 20443) || (part_pdg == 445)  ){

      vector<Int_t> dref = part->getDaughters();
      for (unsigned int j = 0; j < dref.size(); j++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[j]);
        int daughter_pdg = daughter->getParticleID();

        //gamma
        if( daughter_pdg == 22){

          MCTrkIndex = daughter->getTrackIndex();

          Hep3Vector p3mc( daughter->getInitialMomentumX(),
                           daughter->getInitialMomentumY(),
                           daughter->getInitialMomentumZ() );
          p3gmc = p3mc.mag();

        }//end if(gamma)
      }//end for(j)
    }//end if( chic1,2 )
  }//end for(i)




  //Event Navigator map
  const TEvtNavigator* m_TEvtNavigator = selector->GetEventNavigator();
  std::multimap <int, int>::const_iterator it;
  const std::multimap <int, int>&  index = m_TEvtNavigator->m_mcEmcRecShowers;



  for(it = index.begin(); it != index.end(); it++)
    {
      // (*it).first = McParticle id; (*it).second = EmcRecShower id

      //Selection
      if( (*it).first == MCTrkIndex ) {

        S2K2piG.qgENsort.push_back((*it).second);

        //      cout << "\n McParticle id: EmcRecShower id   " << (*it).first << ": " << (*it).second << endl;
      } ///end if(index)

    }//end for(it)




  if(S2K2piG.qgENsort.size() != 0){
    //The best gamma track
    int tmpENTrkID = S2K2piG.qgENsort[0];
    int tmpTrkCount = 1;
    int TrkCount = 0;

    if(S2K2piG.qgENsort.size() > 1){
      for(int i = 1; i < S2K2piG.qgENsort.size(); i++){

        if(tmpENTrkID != S2K2piG.qgENsort[i]){


          if( tmpTrkCount > TrkCount ){ ENTrkID = tmpENTrkID; TrkCount = tmpTrkCount; }

          tmpENTrkID = S2K2piG.qgENsort[i];
          tmpTrkCount = 1;
        }else{ tmpTrkCount++; }


      }//end for(i)
    }//end if(S2K2piG.qgENsort.size() > 1)

    //after for
    if( tmpTrkCount > TrkCount ){ ENTrkID = tmpENTrkID;  TrkCount = tmpTrkCount; }


    //Choose NhitCut
    // if( TrkCount  < 5 ){  his1[652]->Fill(-1);  return false;  }


    his1[653]->Fill(TrkCount); //fot EN test

  }//end if(S2K2piG.qgENsort.size())








  //MC & Rec. comparison
  for(int k = 0; k < S2K2piG.qgChiCshwrID.size(); k++){

    // cout << "\n Rec. Shower TrackID : Event Navigator ID  " <<  S2K2piG.qgChiCshwrID[k] <<" : " << ENTrkID  << endl;

    if( S2K2piG.qgChiCshwrID[k] == ENTrkID ){


      Hep3Vector p3gam = S2K2piG.PgFromChiC[k].vect();
      his1[654]->Fill(  p3gmc - p3gam.mag()  );


      his1[652]->Fill(1);
      return true;
    }//end if(trks)

  }//end for(k)



  return false;
}






//-----------------------------------------------------------------------------
static double GetEcms(int run, double& Lumi)
//-----------------------------------------------------------------------------
{
  struct runInfo {
    int runS, runE; // first and last runs in period
    double Ebeam;   // energy of beam (MeV)
    double Spread;  // beam spread (KeV)
    double Lumi;    // luminosity (pb^-1)
  };

  //Ecm via di-muon (BESIII)
  static runInfo ListRuns[] =  {
        { 33490, 33556,     3807.65,   1088,     50.54 },
        { 33572, 33657,     3896.24,   1025,     52.61 },
        { 33659, 33719,     4085.45,    722,     52.63 },
        { 30372, 30437,     4188.59,    644,     43.09 },
        { 32046, 32140,     4217.13,   1497,     54.13 },
        { 30438, 30491,     4226.26,   1025,     44.40 }, //4230 scan
        { 32141, 32226,     4241.66,   1583,     55.59 },
        { 30492, 30557,     4307.89,   1071,     44.90 },
        { 31281, 31325,     4387.40,   1577,     55.18 },
                        };

  int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

  int absrun = abs(run);

  // we need some values any way
  double Ecms   = 4.260;
  double Spread = 1000.;
  Lumi          = 0;

  int i = 0;
  for(; i < Np; i++) {
     if ( absrun >= ListRuns[i].runS && absrun <= ListRuns[i].runE ) {
       Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV
       Spread = ListRuns[i].Spread * 1.e-6;  // KeV -> GeV
       Lumi   = ListRuns[i].Lumi;
       break;
     }
  }


    // run not in the list
    if ( absrun >= 31327 && absrun <= 31390 ) {
      Ecms =  4.41558; // 4420 scan
      Lumi =  44.67;
    }else if ( absrun >= 36773 && absrun <= 38140 ) {
      Ecms =  4.41558; // 4420 peak
      Lumi =  1028.89;
    } else if ( absrun >= 31983 && absrun <= 32045 ) {
      Ecms =  4.20773; // 4210
      Lumi =  54.55;
    } else if ( absrun >= 32239 && absrun <= 33484 ) {
      Ecms =  4.22626; // 4230
      Lumi =  1047.34;
    } else if ( absrun >= 29677 && absrun <= 29805 ) {
      Ecms =  4.25797; // 4260-1
      Lumi =  523.74;  // (1) + (2)
    } else if ( absrun >= 29822 && absrun <= 30367 ) {
      Ecms =  4.25797; // 4260-2
      Lumi =  523.74;  // (1) + (2)
    }  else if ( absrun >= 31561 && absrun <= 31981 ) {
      Ecms =  4.25797; // 4260-3
      Lumi =  301.93;  // (3)
    } else if ( absrun >= 30616 && absrun <= 31279 ) {
      Ecms =  4.35826; // 4360
      Lumi =  539.84;
    } else if ( absrun >= 36245 && absrun <= 36393 ) {
      Ecms =  4.46706; // 4470
      Lumi =  109.94;
    } else if ( absrun >= 36398 && absrun <= 36588 ) {
      Ecms =  4.52714; // 4530
      Lumi =  109.98;
    } else if ( absrun >= 36603 && absrun <= 36699 ) {
      Ecms =  4.57450; // 4575
      Lumi =  47.67;
    }else if ( absrun >= 35227 && absrun <= 36213 ) {
      Ecms =  4.59953; // 4600
      Lumi =  566.93;
    }else if ( absrun >= 23463 && absrun <= 24141 ) {
      Ecms =  4.00762; // 4009 or 4040
      Lumi =  481.96;
    } else {
//      cout << " GetEcms::WARNING unknown run# " << run << " use Ecms= " << Ecms << endl;
      ERROR_WARNING++;
    }


  return Ecms;
}






//-------------------------------------------------------------------------------- Select3piEvent()------------------------------------------------------------
bool SelectJpsiXEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                      THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << " SelectJpsiXEvent() " << endl;

  m_abscor->AbsorptionCorrection(selector);


  Select S2K2piG; // information about selected kaons, pions and gammas

  leptonISmuon = false;

  //flags  initialization  BhabhaInEvent=false;
  JpsiVertexFit = false;
  Jpsi1CFit = false;
  Jpsi2piFromPsi2S=false;
  JpsiGammaFromChiC1=false;
  JpsiGammaFromChiC2=false;
  MuonsFromJpsi=false;
  MuonsFromJpsi3sigma=false;
  MuonsFromJpsi3032=false;
  ElectronsFromJpsi=false;
  MuonsFromPsi2S=false;
  ElectronsFromPsi2S=false;
  isAntiproton = false;
  isAntideuteron = false;
  isDeuteron = false;
  isLambda = false;
  isAntiLambda = false;
  iMCJpsiISR = false;
  iMCJpsi = false;
  MCJpsiTo2mu = false;
  iMCpsi2SToJpsi2pi = false;
  iMCchiC1ToJpsiGamma = false;
  iMCchiC2ToJpsiGamma = false;
  iMCpsi2SInEvent = false;
  iMCchiC1InEvent = false;
  iMCchiC2InEvent = false;
  JpsiWithoutPcut = false;
  JpsiWithPcut = false;
  SBJpsi=false;
  SBpsi2S=false;
  ISRcutRegion=false;
  isISRevent=false;
  isTrk4psi=false;
  isGoodMCJpsi24pi = true;
  JpsiSelectedForDataMCcompar=false;

  //bkgr study
  isJpsiX=false;
  isJpsiISR=false;
  ispsiX=false;
  ispsiISR=false;
  ischic1X=false;
  ischic2X=false;

  HepLorentzVector P_missing;

  his1[107]->Fill(0);
//-------------------------------------------------------------------------------- Get event information-------------------------------------------------------
  int runNo = m_TEvtHeader->getRunId();
  his1[0]->Fill(abs(runNo));

  int eventNo = m_TEvtHeader->getEventId();


  //--------------------------------------------------Ecms
  double lum = 0;
  Ecms =  GetEcms(abs(runNo), lum);

  his1[30]->Fill(Ecms);
  //for kinematic fits & missing mass
  HepLorentzVector ecms(Ecms*sin(beam_angle), 0, 0, Ecms);
  //------------------------------------------------------



  isMC = (runNo < 0);
  if ( (isMC && m_TMcEvent->getMcParticleCol()->IsEmpty() ) ||
     ( !isMC && m_TMcEvent != 0 && !m_TMcEvent->getMcParticleCol()->IsEmpty() ) ) {
    cout << " ERROR: something wrong: isMC= " << isMC
         << " McParticles= " << m_TMcEvent->getMcParticleCol()->GetEntries()
         << endl;
    return false;
  }

  int mc_type = -1;
  int decCode = -1;
  if ( isMC ) {
    m_EventTagSvc->setMcEvent(m_TMcEvent);
    mc_type = ParticleTypesMC(selector, S2K2piG);
    decCode = ((m_EventTagSvc->getEventTag())&0xFF00) >> 8;
    his1[101]->Fill(mc_type);
    his1[102]->Fill(decCode);

    if( decCode==70 ) his1[108]->Fill(0);
    //   if( decCode==74 ||decCode==75 ||decCode==156 ) his1[109]->Fill(0);
  }



  //=================================================================
  /*  const TEvtNavigator* m_TEvtNavigator = selector->GetEventNavigator();

  if(!m_TEvtNavigator){cout << "  m_TEvtNavigator = NULL" << endl;}
  else{
    m_TEvtNavigator->Print();
    cout << "Hand read" << endl;
    cout << " McParticle id <-> RecMdcTrack (RecMdcKalTrack) id  size   "<<  m_TEvtNavigator-> m_mcMdcTracks.size() << endl;

    }*/
  //=================================================================


  MuonsMomentumHistMC( selector, 0);  // init MC  momentum  hists
  if(iMCpsi2SToJpsi2pi) { Psi2SMomentumHistMC( selector, 0); }   // init MC  momentum  hists
  if( iMCchiC1ToJpsiGamma || iMCchiC2ToJpsiGamma ) { chic12MomentumHistMC( selector, 0); }   // init MC  momentum  hists
  MuonsMomentumHist_iMC( selector, 0);  // init MC  momentum  hists

  //  if( (!iMCpsi2SToJpsi2pi && iMCpsi2SInEvent) || (!iMCchiC1ToJpsiGamma && iMCchiC1InEvent) || (!iMCchiC2ToJpsiGamma && iMCchiC2InEvent) ){SaveInfoMC(string("infoMC3mix_4600_irregPSICHI.dat"),selector, 1);}

  HepSymMatrix ErrMatrixOrigin(3, 0);
  Hep3Vector xorigin = getVertexOrigin(runNo, ErrMatrixOrigin);

  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();



  //------------------------------------- For Event Navigator check
  //cout << "\n\n Charged tracks: " << evtRecEvent->totalCharged() << endl;
  //   bool tmpflag = MatchRecMcTrks_2pi(selector, S2K2piG);
  //--------------------------------------




  //------------------------------------------------------------------------------ All charged tracks----------------------------------------------------------
  // 1) select 2 good tracks with opposite charges from primary vertex
  int Ncharge = 0;
  int nGoodTrkIP = 0;
  int Np= 0, Nap = 0;
  double pPQ[2] = {0., 0.}, apPQ[2] = {0., 0.}, dPQ = 0, adPQ = 0, pE[2] = {0., 0.}, apE[2] = {0., 0.}, dE = 0, adE = 0, pDedx[2] = {0., 0.}, apDedx[2] = {0., 0.}, dDedx = 0, adDedx = 0;
  // bool Pless2 = true;
  for(int i = 0; i < evtRecEvent->totalCharged(); i++){

    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if( !itTrk->isMdcTrackValid() ) continue;



    RecMdcTrack *mdcTrk = itTrk->mdcTrack();



    double theta = mdcTrk->theta();
    double cosTheta = cos(theta);




    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
    VFHelix helixip(point0,a,Ea);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction

    his1[2]->Fill(Rvxy0);
    his1[3]->Fill(Rvz0);
    his1[4]->Fill(fabs(cosTheta));


    if(fabs(cosTheta) >= 0.93) continue;
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  J/psi  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(fabs(Rvz0) >= 10.0) continue;
    if(fabs(Rvxy0) >= 1.0) continue;



    //MDC save
    S2K2piG.qtrks.push_back(i);
    S2K2piG.qRxy0.push_back(Rvxy0);
    S2K2piG.qRz0.push_back(Rvz0);
    Ncharge += mdcTrk->charge();
    his1[9]->Fill( mdcTrk->charge()*mdcTrk->p() );

    //MC vs Data study
    S2K2piG.qMom.push_back( mdcTrk->p() );
    S2K2piG.qCosTheta.push_back(cosTheta);


    //if Kalman -> save
    RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
    if ( !mdcKalTrk ) {
      // cout << " No Kalman track ! EventNo= "<<eventNo << endl;
      continue;
    }


    if(mdcTrk->charge() > 0) { S2K2piG.qtrksPos.push_back(i);   S2K2piG.qRxy0Pos.push_back(Rvxy0);  S2K2piG.qRz0Pos.push_back(Rvz0);}
    else {S2K2piG.qtrksNeg.push_back(i);   S2K2piG.qRxy0Neg.push_back(Rvxy0);  S2K2piG.qRz0Neg.push_back(Rvz0);}

  }
  //-------------------------------------------------End for(i)

  int nGood = S2K2piG.nGood();  //MDC qtrks.size()
  int nPos = S2K2piG.nPos(); //Kalman qtrksPos.size()
  int nNeg = S2K2piG.nNeg(); //Kalman qtrksNeg.size()
  his1[5]->Fill(nGood);
  his1[6]->Fill(Ncharge);





  //for preselection J/psi + anything
  if (  nPos == 0 || nNeg == 0 )  return false; // skip event






  his1[107]->Fill(1);


//------------------------------- Finish Good Charged Track Selection ------------------------












  //for preselection J/psi + anything
  int  nlleEPcomb = 0;
  int  nllmuEPcomb = 0;
  for(int i = 0; i < nPos; i++){
    for(int j = 0; j < nNeg; j++){

      DstEvtRecTracks* itTrkp=(DstEvtRecTracks*) evtRecTrkCol->At( S2K2piG.qtrksPos[i] );
      DstEvtRecTracks* itTrkm=(DstEvtRecTracks*) evtRecTrkCol->At( S2K2piG.qtrksNeg[j] );

      //-------------------------------------------------------
      RecMdcKalTrack* mdcKalTrkp = itTrkp->mdcKalTrack();
      if ( !mdcKalTrkp ) {
        cout << " No Kalman track ! EventNo= "<<eventNo << endl;
        return false;
      }

      RecMdcKalTrack* mdcKalTrkm = itTrkm->mdcKalTrack();
      if ( !mdcKalTrkm ) {
        cout << " No Kalman track ! EventNo= "<<eventNo << endl;
        return false;
      }
      //------------------------------------------------------



      RecMdcTrack *mdcTrkp = itTrkp->mdcTrack();
      double Ptrkp = mdcTrkp->p();
      double cosThetap = cos(mdcTrkp->theta());
      RecMdcTrack *mdcTrkm = itTrkm->mdcTrack();
      double Ptrkm = mdcTrkm->p();
      double cosThetam = cos(mdcTrkm->theta());

      if ( itTrkp->isEmcShowerValid() && itTrkm->isEmcShowerValid() ){

        RecEmcShower *emcTrkp = itTrkp->emcShower();
        double EnergyEMCtrkp = emcTrkp->energy();
        double EoverPp =  EnergyEMCtrkp / Ptrkp;
        RecEmcShower *emcTrkm = itTrkm->emcShower();
        double EnergyEMCtrkm = emcTrkm->energy();
        double EoverPm =  EnergyEMCtrkm / Ptrkm;

        his1[417]->Fill(EoverPp);
        his1[417]->Fill(EoverPm);


        his1[615]->Fill(EnergyEMCtrkp);
        his1[615]->Fill(EnergyEMCtrkm);

        his2[462]->Fill(EnergyEMCtrkp, cosThetap);
        his2[462]->Fill(EnergyEMCtrkm, cosThetam);



        if( isMC &&  MCJpsiTo2mu ){
          his1[637]->Fill(EnergyEMCtrkp);
          his1[637]->Fill(EnergyEMCtrkm);
        }








        //for mu+mu- candidates
        if( (EnergyEMCtrkp < 0.6) && (EnergyEMCtrkm < 0.6) ){



          //syst. study --------------------------------------------
//      if( (EnergyEMCtrkp < 0.5) && (EnergyEMCtrkm < 0.5) ){
        //      if( (EnergyEMCtrkp < 0.7) && (EnergyEMCtrkm < 0.7) ){
          //------------------------------------------------------

          mdcKalTrkp->setPidType(RecMdcKalTrack::muon);
          Hep3Vector p3mup(mdcKalTrkp->px(), mdcKalTrkp->py(), mdcKalTrkp->pz());
          double emup = sqrt( p3mup.mag2() + mmu*mmu );
          int tmp_imup = S2K2piG.qtrksPos[i];
          HepLorentzVector tmp_pmup =  HepLorentzVector( p3mup, emup );

          mdcKalTrkm->setPidType(RecMdcKalTrack::muon);
          Hep3Vector p3mum(mdcKalTrkm->px(), mdcKalTrkm->py(), mdcKalTrkm->pz());
          double emum = sqrt( p3mum.mag2() + mmu*mmu );
          int tmp_imum = S2K2piG.qtrksNeg[j];
          HepLorentzVector tmp_pmum =  HepLorentzVector( p3mum, emum );



          //for Bhabha check
          double cos2trks = ( p3mup.x()*p3mum.x() +p3mup.y()*p3mum.y() + p3mup.z()*p3mum.z() ) /( p3mup.mag()*p3mum.mag() );

          HepLorentzVector Pllmu = tmp_pmup + tmp_pmum;



        //total check
        //  his1[647]->Fill(Pllmu.m());


          //Jpsi
          if(Pllmu.m() > 2.8  &&  Pllmu.m() < 3.4){



            if(nGood!=2){his1[423]->Fill(Pllmu.m());}



            //selection flag
            MuonsFromJpsi=true;


            //bkgr study
            if( nGood!=2 && isJpsiX ){his1[655]->Fill(Pllmu.m());}

            //--------------------------------------- J/psi candidate -----------------------
            S2K2piG.qmups.push_back(tmp_imup);
            S2K2piG.qmums.push_back(tmp_imum);
            S2K2piG.Pmups.push_back(tmp_pmup);
            S2K2piG.Pmums.push_back(tmp_pmum);



            if(isMC){ his1[103]->Fill(decCode); }

            //-------------------------------------------------------------------------------




            //--------------------------------------------------------------------Data vs MC
            double sigma_Jpsi = 0.017;

            //Jpsi Side Band  NOT BE USED! 3-09-2021
            if( ( Pllmu.m() > (mJpsi - 7.*sigma_Jpsi) && Pllmu.m()  < (mJpsi - 4.*sigma_Jpsi) ) || ( Pllmu.m() < (mJpsi + 7.*sigma_Jpsi) && Pllmu.m()  > (mJpsi + 4.*sigma_Jpsi)  ) ){ SBJpsi=true; }

            if( Pllmu.m() > (mJpsi - 3.*sigma_Jpsi)  &&  Pllmu.m() < (mJpsi + 3.*sigma_Jpsi) ){

              JpsiSelectedForDataMCcompar=true;

              if(nGood!=2){

                his1[660]->Fill( Pllmu.vect().mag() );
                his1[661]->Fill( Pllmu.vect().perp() );
                his1[662]->Fill( RtoD( Pllmu.phi() ) );
                his1[663]->Fill( Pllmu.cosTheta() );

                his1[668]->Fill(nGood);
                for(int k=0; k<nGood; k++){

                  his1[669]->Fill(S2K2piG.qMom[k]);
                  his1[670]->Fill(S2K2piG.qCosTheta[k]);
                  his2[472]->Fill(S2K2piG.qMom[k], nGood);
                  his2[473]->Fill(S2K2piG.qCosTheta[k],nGood);
                }//end for(k)


              }//end if(nGood!=2)
            }//end if 3-sigma
            //==============================================================================









              //other info =========================================================================================================


            //rec
            his2[75]->Fill(  S2K2piG.pmup.vect().mag() , S2K2piG.pmup.cosTheta() );
            his2[75]->Fill(  S2K2piG.pmum.vect().mag() , S2K2piG.pmum.cosTheta() );



            if( (EoverPp < 0.05) || (EoverPm < 0.05) ){ his1[636]->Fill(Pllmu.m()); }











            //---------J/psi: 3.0 - 3.2
            if( Pllmu.m() > 3.0  &&  Pllmu.m() < 3.2 ){



              //flag
              MuonsFromJpsi3sigma=true;



              //N integral
              if(nGood!=2){ his1[649]->Fill(1); }



              nllmuEPcomb++;









              //------------------------------------------------ iMC
              if(iMCJpsiISR && isMC){
                his2[101]->Fill( Pllmu.vect().mag(), Pllmu.m()  );
                his2[102]->Fill( Pllmu.vect().mag(),  nGood);
                his2[103]->Fill( Pllmu.vect().mag(),  Pllmu.cosTheta() );
                his2[110]->Fill( Pllmu.vect().mag(), pJpsi_iMC);

                //if( (Pllmu.vect().mag() < 1.) && (fabs(Pllmu.cosTheta()) < 0.8) ){ SaveInfoMC(string("infoMC_JpsiISR_4600_CentreNonISR.dat"),selector, 1);}
                //if( (Pllmu.vect().mag() < 1.) && (fabs(Pllmu.cosTheta()) > 0.8) ){ SaveInfoMC(string("infoMC_JpsiISR_4600_NonCentreNonISR.dat"),selector, 2);}
                //if( (Pllmu.vect().mag() > 1.) && (fabs(Pllmu.cosTheta()) < 0.8) ){ SaveInfoMC(string("infoMC_JpsiISR_4600_CentreISR.dat"),selector, 3);}
                //if( (Pllmu.vect().mag() > 1.) && (fabs(Pllmu.cosTheta()) > 0.8) ){ SaveInfoMC(string("infoMC_JpsiISR_4600_NonCentreISR.dat"),selector, 4);}

                //test
                /*      if( (Pllmu.vect().mag() < 1.) && (fabs(Pllmu.cosTheta()) < 0.8) ){ SaveInfoMC(string("infoMC_4600_test.dat"),selector, 1);}
                if( (Pllmu.vect().mag() < 1.) && (fabs(Pllmu.cosTheta()) > 0.8) ){ SaveInfoMC(string("infoMC_4600_test.dat"),selector, 2);}
                if( (Pllmu.vect().mag() > 1.) && (fabs(Pllmu.cosTheta()) < 0.8) ){ SaveInfoMC(string("infoMC_4600_test.dat"),selector, 3);}
                if( (Pllmu.vect().mag() > 1.) && (fabs(Pllmu.cosTheta()) > 0.8) ){ SaveInfoMC(string("infoMC_4600_test.dat"),selector, 4);}
                */
              }
              //-------------------------------------------------




              his2[79]->Fill( Pllmu.vect().mag(), Pllmu.m()  );



              //nGood
              if(nGood!=2){
                his2[76]->Fill( Pllmu.vect().mag(), Pllmu.m()  );
                his2[77]->Fill( Pllmu.vect().mag(),  nGood);
                his2[78]->Fill( Pllmu.m() ,  nGood);
                his2[80]->Fill( Pllmu.vect().mag(),  Pllmu.cosTheta() );
                his1[550]->Fill(Pllmu.vect().mag());
                his1[616]->Fill(EnergyEMCtrkp);
                his1[616]->Fill(EnergyEMCtrkm);
                his1[635]->Fill(EoverPp);
                his1[635]->Fill(EoverPm);
              }



              //PJpsi_rec
              his1[514]->Fill( Pllmu.vect().mag() );
              his1[533]->Fill( RtoD( Pllmu.phi() ) );
              his1[536]->Fill( Pllmu.cosTheta() );
              his2[105]->Fill( RtoD( Pllmu.phi() ),  Pllmu.cosTheta() );

              //---------------------------------------!!!
              if(iMCJpsiISR && isMC){ his1[515]->Fill( Pllmu.vect().mag() ); }

              if( iMCJpsi && isMC ) {
                his1[530]->Fill( Pllmu.vect().mag() );

                // SaveInfoMC(string("info_iMC_3mix_4600_JpsiONLY_rec.dat"),selector, 1);
              }
              //----------------------------------------


            }//------------------------------------------------------end if(3-sigma Jpsi)
            //===============================================




            //for Bhabha check
            his1[411]->Fill(Ptrkp);
            his1[411]->Fill(Ptrkm);
            his1[412]->Fill(EnergyEMCtrkp);
            his1[412]->Fill(EnergyEMCtrkm);
            his1[427]->Fill( cos2trks  );
          }//end if(ll inv mass Jpsi)



        }//end if(mumu)
      }//end if(isEmcShowerValid)

    }//end for(nNeg)
  }//end for(nPos)

  his1[425]->Fill( nllmuEPcomb);



   //reject events without J/psi -> mu+mu-/e+e-
  //  if( !ElectronsFromJpsi && !MuonsFromJpsi && !ElectronsFromPsi2S && !MuonsFromPsi2S) return false; //skip event





   //3-sigma J/psi: 3.0 - 3.2
  if(  MuonsFromJpsi3sigma ) {

    his1[107]->Fill(2);


    if(nGood!=2){
      his1[621]->Fill(nGood);
      MuonsMomentumHistMC( selector, 1);
      MuonsMomentumHist_iMC( selector, 1);

    }//end if(nGood)
  }//end if( MuonsFromJpsi3sigma)



  // J/psi Side Band
  if( !MuonsFromJpsi3sigma && SBJpsi && (nGood!=2) ){ his1[629]->Fill(nGood); }









    //extract  J/psi -> mu+mu-   ONLY!
  if(  !MuonsFromJpsi ) return false; //skip event
  //  if(  !MuonsFromJpsi && !ElectronsFromJpsi ) return false; //skip event






  //------------------------------------------------
  if (init_selection) return true;
  //------------------------------------------------






//------------------------------------- All neutral tracks----------------------------------------
  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if(!itTrk->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = itTrk->emcShower();

    // 1) good EMC time:
    if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) continue;

//-------------------------------The energy deposited in the barrel (endcap) part of EMC---------
    double eraw = emcTrk->energy();
    double absCosTheta = fabs(  cos(emcTrk->theta()) );

    // 2) good EMC energy
    bool GoodCluster=false;
    if ( absCosTheta < 0.8 && eraw > 25E-3 ) GoodCluster=true;                          //barrel
    if ( absCosTheta > 0.85 && absCosTheta < 0.92 && eraw >  50E-3 ) GoodCluster=true;  //endcap
    if ( !GoodCluster ) continue;

    // 3) the nearest charged track is far from cluster  --------------------------
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    double dang = 200.; // min angle between cluster and track

    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      DstEvtRecTracks* jtTrk = (DstEvtRecTracks*) evtRecTrkCol->At(j);
      if ( !jtTrk->isExtTrackValid() ) continue;
      RecExtTrack *extTrk = jtTrk->extTrack();
      if ( extTrk->emcVolumeNumber() == -1 ) continue; // track does not hit EMC
      Hep3Vector extpos = extTrk->emcPosition();
      double angd = extpos.angle(emcpos);

      if ( angd < dang ) dang = angd;
    }
    //--------------------------------------------End for(j)

    his1[7]->Fill(fabs(RtoD(dang)));
    his1[651]->Fill( RtoD(dang)*RtoD(dang) );

    if ( fabs(RtoD(dang)) <= 20. )  continue;

    S2K2piG.gammas.push_back(i);

    //---------------------------------------------------------------------------- (E,Px,Py,Pz) for good gammas -----------------------------------------------
    Hep3Vector p3 = emcpos - xorigin; // position emc cluster wrt of vertex
    p3 *= eraw / p3.mag();            // normalization on energy
    S2K2piG.Pg.push_back( HepLorentzVector(p3,eraw) );
  }
//-----------------------------------------------End for(i)

  int nGam = S2K2piG.nGam();
  his1[10]->Fill(nGam);
  //  if ( nGam < 1) return false; // skip event

  his1[107]->Fill(3);
  if ( isMC && (decCode==70) )  his1[108]->Fill(3);

//------------------------------------- Finish Good Photon Selection





  //TEST
  for(int i = 0; i < nGam; i++){

#ifdef  MC_chi_c1
    his2[470]->Fill(S2K2piG.Pg[i].e() *1.e3, S2K2piG.Eg_chic1);
#endif

#ifdef  MC_chi_c2
    his2[471]->Fill(S2K2piG.Pg[i].e() *1.e3, S2K2piG.Eg_chic2);
#endif
}
  //====================================================





//=================================================================== J/psi -> mu+mu-
  if(MuonsFromJpsi){
    //===============================================================  2-track events:
    if( nGood==2 && nGam > 1 ){


      HepLorentzVector Pllmu = S2K2piG.Pmups[0] +  S2K2piG.Pmums[0];
      his1[423]->Fill(Pllmu.m());


      //3-sigma
      if(  MuonsFromJpsi3sigma ) {

        //N integral
        his1[649]->Fill(1);

        //J/psi signal
        his1[621]->Fill(nGood);
        MuonsMomentumHistMC( selector, 1);
        MuonsMomentumHist_iMC( selector, 1);



        //other info
        his2[76]->Fill( Pllmu.vect().mag(), Pllmu.m()  );
        his2[77]->Fill( Pllmu.vect().mag(),  nGood);
        his2[78]->Fill( Pllmu.m() ,  nGood);
        his2[80]->Fill( Pllmu.vect().mag(),  Pllmu.cosTheta() );
        his1[550]->Fill(Pllmu.vect().mag());


      }//end if(MuonsFromJpsi3sigma)


      //bkgr study
      if(isJpsiX){ his1[655]->Fill(Pllmu.m()); }



      if(JpsiSelectedForDataMCcompar){

        his1[660]->Fill( Pllmu.vect().mag() );
        his1[661]->Fill( Pllmu.vect().perp() );
        his1[662]->Fill( RtoD( Pllmu.phi() ) );
        his1[663]->Fill( Pllmu.cosTheta() );


        his1[668]->Fill(nGood);
        for(int k=0; k<nGood; k++){

          his1[669]->Fill(S2K2piG.qMom[k]);
          his1[670]->Fill(S2K2piG.qCosTheta[k]);
          his2[472]->Fill(S2K2piG.qMom[k], nGood);
          his2[473]->Fill(S2K2piG.qCosTheta[k],nGood);

        }//end for(k)

        his1[671]->Fill(nGam);
        for(int k=0; k<nGam; k++){

          his1[672]->Fill( S2K2piG.Pg[k].vect().mag() );
          his1[673]->Fill( S2K2piG.Pg[k].cosTheta() );
          his2[474]->Fill( S2K2piG.Pg[k].vect().mag(), nGam);
          his2[475]->Fill( S2K2piG.Pg[k].cosTheta(), nGam);
        }//end for(k)


      }//end if(JpsiSelectedForDataMCcompar)



      // J/psi Side Band
      if( !MuonsFromJpsi3sigma && SBJpsi ){ his1[629]->Fill(nGood); }

    }//end if( 2-track events)





    //ISR Jpsi measurement
    if(nGood==2 && nGam < 2){

      HepLorentzVector Pllmu = S2K2piG.Pmups[0] + S2K2piG.Pmums[0];
      his1[557]->Fill(Pllmu.m());

      //Jpsi ISR efficiency
      if(MuonsFromJpsi3sigma){
        his2[457]->Fill( Pllmu.vect().mag(),  Pllmu.cosTheta() );

        //N integral
        his1[649]->Fill(0);
      }//end if(3-sigma)

      //bkgr study
      if(isJpsiISR){his1[656]->Fill(Pllmu.m());}

    }//end if(ISR Jpsi)
    //====================================================================


    //MC vs Data JpsiX
    if( (nGood!=2) && JpsiSelectedForDataMCcompar){

      his1[671]->Fill(nGam);
      for(int k=0; k<nGam; k++){

        his1[672]->Fill( S2K2piG.Pg[k].vect().mag() );
        his1[673]->Fill( S2K2piG.Pg[k].cosTheta() );
        his2[474]->Fill( S2K2piG.Pg[k].vect().mag(), nGam);
        his2[475]->Fill( S2K2piG.Pg[k].cosTheta(), nGam);
      }//end for(k)
    }//end if(nGood!=2)

  }//end if(MuonsFromJpsi)
















  //--------------------------------------------------------------------------------------------- { J/psi  Vertex + 1C Fit
  VertexFit* vtxfit = VertexFit::instance();
  KinematicFit * kmfit = KinematicFit::instance();
  WTrackParameter wmup, wmum;
  double chisq = 9999.;

  if(MuonsFromJpsi){

    for(int i = 0; i < S2K2piG.qmups.size(); i++){

      JpsiVertexFit = false;

      HepLorentzVector PJpsi = S2K2piG.Pmups[i] + S2K2piG.Pmums[i];




      if( PJpsi.m() > 3.0  &&  PJpsi.m() < 3.2 ){


        RecMdcKalTrack *mupTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(S2K2piG.qmups[i]))->mdcKalTrack();
        RecMdcKalTrack *mumTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(S2K2piG.qmums[i]))->mdcKalTrack();

        WTrackParameter wvmupTrk, wvmumTrk;
        wvmupTrk = WTrackParameter(mmu, mupTrk->getZHelix(), mupTrk->getZError());
        wvmumTrk = WTrackParameter(mmu, mumTrk->getZHelix(), mumTrk->getZError());

        HepPoint3D vx(0., 0., 0.);
        HepSymMatrix Evx(3, 0);
        double bx = 1E+6;
        double by = 1E+6;
        double bz = 1E+6;
        Evx[0][0] = bx*bx;
        Evx[1][1] = by*by;
        Evx[2][2] = bz*bz;

        VertexParameter vxpar;
        vxpar.setVx(vx);
        vxpar.setEvx(Evx);


        vtxfit->init();
        vtxfit->AddTrack(0,  wvmupTrk);
        vtxfit->AddTrack(1,  wvmumTrk);
        vtxfit->AddVertex(0, vxpar,0, 1);
        bool ok = vtxfit->Fit(0);
        his1[14]->Fill( double(ok) );

        if ( ok ){

          JpsiVertexFit = true;

          vtxfit->BuildVirtualParticle(0);
          his1[15]->Fill(vtxfit->chisq(0));
          vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0

          wmup = vtxfit->wtrk(0);
          wmum = vtxfit->wtrk(1);
        }


        if(JpsiVertexFit){

          kmfit->init();
          kmfit->AddTrack(0, wmup);
          kmfit->AddTrack(1, wmum);
          kmfit->AddResonance(0, mJpsi, 0, 1);
          bool oksq = kmfit->Fit();

          if ( oksq ) {
            double chi2 = kmfit->chisq();
            if ( chi2 < chisq ) {
              chisq = chi2;
              S2K2piG.imup = S2K2piG.qmups[i];
              S2K2piG.imum = S2K2piG.qmums[i];
              S2K2piG.pmup = S2K2piG.Pmups[i];
              S2K2piG.pmum = S2K2piG.Pmums[i];
            }//end if(chi2)
          }//end if(oksq)
        }//if(JpsiVertexFit)
      }//end if(Jpsi 3 Sigma)
    }//end for(i)



    //min chisq
    his1[31]->Fill(chisq);


    if (chisq < 50.){

      Jpsi1CFit = true;

      S2K2piG.pmup_1C = kmfit->pfit(0);
      S2K2piG.pmum_1C = kmfit->pfit(1);

      HepLorentzVector P_Jpsi_1C = S2K2piG.pmup_1C + S2K2piG.pmum_1C;
      his1[525]->Fill( P_Jpsi_1C.m() );

      //------------------------------------------------------------{  J/psi polarization (1C)
      Hep3Vector beta = P_Jpsi_1C.boostVector();
      HepLorentzRotation L(-beta); // Lorentz boost
      HepLorentzVector P_mup_Jpsi, P_mum_Jpsi;
      P_mup_Jpsi = L * S2K2piG.pmup_1C;
      P_mum_Jpsi = L * S2K2piG.pmum_1C;
      his1[497]->Fill( P_mup_Jpsi.cosTheta() );
      his1[498]->Fill( P_mum_Jpsi.cosTheta() );
      his1[601]->Fill(  RtoD(P_mup_Jpsi.theta()) );
      his1[602]->Fill(  RtoD(P_mum_Jpsi.theta()) );
      //============================================================= J/psi polarization (1C)}

    }//if(chisq<50)







    //======================================================================= 1C-Fit: J/psi pi+ pi-
    HepLorentzVector invJpsi2pi;
    if(Jpsi1CFit){
      if( (nPos > 1) && (nNeg > 1) ){

        for(int i = 0; i < nPos; i++){
          for(int j = 0; j < nNeg; j++){

            if( (S2K2piG.qtrksPos[i] != S2K2piG.imup) && ( S2K2piG.qtrksNeg[j] != S2K2piG.imum) ){

              DstEvtRecTracks* itTrkp=(DstEvtRecTracks*) evtRecTrkCol->At( S2K2piG.qtrksPos[i] );
              DstEvtRecTracks* itTrkm=(DstEvtRecTracks*) evtRecTrkCol->At( S2K2piG.qtrksNeg[j] );

              //-------------------------------------------------------
              RecMdcKalTrack* mdcKalTrkp = itTrkp->mdcKalTrack();
              if ( !mdcKalTrkp ) {
                cout << " No Kalman track ! EventNo= "<<eventNo << endl;
                return false;
              }

              RecMdcKalTrack* mdcKalTrkm = itTrkm->mdcKalTrack();
              if ( !mdcKalTrkm ) {
                cout << " No Kalman track ! EventNo= "<<eventNo << endl;
                return false;
              }
              //------------------------------------------------------



                mdcKalTrkp->setPidType(RecMdcKalTrack::pion);
                Hep3Vector p3pip(mdcKalTrkp->px(), mdcKalTrkp->py(), mdcKalTrkp->pz());
                double epip = sqrt( p3pip.mag2() + mpi*mpi );
                HepLorentzVector p4pip = HepLorentzVector(p3pip,epip);

                mdcKalTrkm->setPidType(RecMdcKalTrack::pion);
                Hep3Vector p3pim(mdcKalTrkm->px(), mdcKalTrkm->py(), mdcKalTrkm->pz());
                double epim = sqrt( p3pim.mag2() + mpi*mpi );
                HepLorentzVector p4pim = HepLorentzVector(p3pim,epim);


                HepLorentzVector inv2pi = p4pip + p4pim;


                invJpsi2pi = S2K2piG.pmup_1C + S2K2piG.pmum_1C  +  p4pip + p4pim;

                //non-ISR  psi'
                if(( nGood==4 && nGam > 1 ) || (nGood!=4)){

                  his1[526]->Fill( invJpsi2pi.m() );
                  his2[463]->Fill( invJpsi2pi.m(), inv2pi.m() );

                  //N integral
                  if((invJpsi2pi.m()  > 3.65) && (invJpsi2pi.m() < 3.72)){ his1[648]->Fill(1); }

                  //bkgr study
                  if(ispsiX){ his1[657]->Fill( invJpsi2pi.m() ); }

                }//end if(non-ISR  psi')


                //ISR psi'
                if( nGood==4 && nGam < 2 ){

                  his1[556]->Fill( invJpsi2pi.m() );
                  his2[464]->Fill( invJpsi2pi.m(), inv2pi.m() );

                  //N integral
                  if((invJpsi2pi.m()  > 3.59) && (invJpsi2pi.m() < 3.77)){ his1[648]->Fill(0); }

                  //bkgr study
                  if(ispsiISR){ his1[658]->Fill( invJpsi2pi.m() ); }

                }//end if(ISR  psi')




                //psi'X
                if( (invJpsi2pi.m()  > 3.65)  &&  (invJpsi2pi.m() < 3.72) ){



                  Jpsi2piFromPsi2S=true;



                  if( nGood==4 && nGam < 2){ isTrk4psi=true; }

                    //psi(2S) rec.
                    his1[528]->Fill( invJpsi2pi.vect().mag() );
                    his1[537]->Fill( RtoD( invJpsi2pi.phi() ) );
                    his1[538]->Fill( invJpsi2pi.cosTheta() );


                    //nGood
                    his2[81]->Fill( invJpsi2pi.vect().mag(), invJpsi2pi.m()  );
                    his2[82]->Fill( invJpsi2pi.vect().mag(),  nGood);
                    his2[83]->Fill( invJpsi2pi.m() ,  nGood);
                    his2[85]->Fill( invJpsi2pi.vect().mag(),  invJpsi2pi.cosTheta() );

                    if( isMC && iMCpsi2SToJpsi2pi ){ his1[618]->Fill( invJpsi2pi.vect().mag() ); }

                    //save pions info for eMC psi(2S) pi+ pi- near 4260
                    S2K2piG.qpipsFromPsi2S.push_back(S2K2piG.qtrksPos[i]);
                    S2K2piG.qpimsFromPsi2S.push_back(S2K2piG.qtrksNeg[j]);
                    S2K2piG.PpipsFromPsi2S.push_back(p4pip);
                    S2K2piG.PpimsFromPsi2S.push_back(p4pim);

                    //for Event Navigator
                    RecMdcTrack *mdcTrkp = itTrkp->mdcTrack();
                    RecMdcTrack *mdcTrkm = itTrkm->mdcTrack();
                    S2K2piG.qpipPsi2StrkID.push_back( mdcTrkp->trackId() );
                    S2K2piG.qpimPsi2StrkID.push_back( mdcTrkm->trackId() );



                }//end if(psi' X)


            }//end if(different tracks)
          }//end for(j)
        }//end for(i)
      }//if(nGood > 3)
    }//if(Jpsi1CFit)



    //result psi' X
    if(Jpsi2piFromPsi2S){

      bool isGoodMCevent = true;

#ifdef  MC_psi2S_2pi
      isGoodMCevent = MatchRecMcTrks_2pi(selector, S2K2piG); // !!! ONLY for eMC psi(2S) pi+ pi-
#endif


      if(isGoodMCevent) {
        if((nGood==4 && nGam > 1) || (nGood!=4)){Psi2SMomentumHistMC( selector, 1); }
      }//end if(isGoodMCevent)

      his1[107]->Fill(5);

      if((nGood==4 && nGam > 1) || (nGood!=4)){ his1[622]->Fill(nGood);}




      if( nGood==4 && nGam < 2 ){

        //MC psi' 2pi ONLY
        if(isGoodMCevent) {his2[469]->Fill( invJpsi2pi.vect().mag(),  invJpsi2pi.cosTheta() );}

      }//end if(2-2<2)

    }//end if(Jpsi2piFromPsi2S)







     //ISR psi'
    if( nGood==4 && nGam < 2 ){

      //mass cut
      if( (invJpsi2pi.m()  > 3.59)  &&  (invJpsi2pi.m() < 3.77) ){
        his2[458]->Fill( invJpsi2pi.vect().mag(),  invJpsi2pi.cosTheta() );
      }



      //to check
      his1[638]->Fill( invJpsi2pi.vect().mag() );
      his1[639]->Fill( invJpsi2pi.cosTheta() );
      his1[640]->Fill( invJpsi2pi.phi() );
    }


//=================================================================  J/psi pi+ pi- }



    //------------------------------------------------------------{  J/psi polarization
    if(Jpsi1CFit){

      /*
        TLorentzVector LV_DIMUON=LV_MU_1+LV_MU_2; - тут  LV_MU_1 и LV_MU_2 -  в лаб. системе
        TVector3 boost1=LV_DIMUON.BoostVector();
        TLorentzRotation l5;
        l5.Boost(-boost1);
        TLorentzVector LV_mu1,LV_mu2;
        LV_mu1=l5*LV_MU_1;
        LV_mu2=l5*LV_MU_2;  тут  LV_mu1 и LV_mu2 - в системе покоя димюона.

*/

      HepLorentzVector PJpsi = S2K2piG.pmup + S2K2piG.pmum;

      Hep3Vector beta = PJpsi.boostVector();
      HepLorentzRotation L(-beta); // Lorentz boost
      HepLorentzVector P_mup_Jpsi, P_mum_Jpsi;
      P_mup_Jpsi = L * S2K2piG.pmup;
      P_mum_Jpsi = L * S2K2piG.pmum;
      his1[495]->Fill( P_mup_Jpsi.cosTheta() );
      his1[496]->Fill( P_mum_Jpsi.cosTheta() );
      his1[499]->Fill(  RtoD(P_mup_Jpsi.theta()) );
      his1[600]->Fill(  RtoD(P_mum_Jpsi.theta()) );


    }//end if(Jpsi1CFit)
    //============================================================= J/psi polarization }


  }//if(MuonsFromJpsi)
  //===============================================================  J/psi  Vertex + 1C Fit }












//=================================================================== J/psi -> mu+mu-
  if(MuonsFromJpsi){

    //chi_c1,2
    if( (nGood==2 && nGam > 1)  || (nGood!=2  && nGam > 0)   ){


      //1C
      if(Jpsi1CFit){
        for(int i = 0; i < nGam; i++){
          HepLorentzVector invJpsiGamma = S2K2piG.Pg[i] + S2K2piG.pmup_1C  +  S2K2piG.pmum_1C;
          his1[527]->Fill( invJpsiGamma.m() );

          //bkgr study
          if( (ischic1X) || (ischic2X) ){ his1[659]->Fill( invJpsiGamma.m() ); }


          //chi_c1  -------------------------------------
          if( invJpsiGamma.m()  > 3.47  &&  invJpsiGamma.m() < 3.53343 ){

            JpsiGammaFromChiC1=true;

            his1[650]->Fill(0);




#ifdef  MC_chi_c1
DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(S2K2piG.gammas[i]);
RecEmcShower *emcTrk = itTrk->emcShower();
S2K2piG.qgChiCshwrID.push_back( emcTrk->cellId() );
S2K2piG.PgFromChiC.push_back(S2K2piG.Pg[i]);
#endif







            //chi_c1 rec.
            his1[609]->Fill( invJpsiGamma.vect().mag() );
            his1[539]->Fill( RtoD( invJpsiGamma.phi() ) );
            his1[540]->Fill( invJpsiGamma.cosTheta() );

            //nGood
            his2[92]->Fill( invJpsiGamma.vect().mag(), invJpsiGamma.m()  );
            his2[93]->Fill( invJpsiGamma.vect().mag(),  nGood);
            his2[94]->Fill( invJpsiGamma.m() ,  nGood);
            his2[95]->Fill( invJpsiGamma.vect().mag(),  invJpsiGamma.cosTheta() );

            if( isMC && iMCchiC1ToJpsiGamma ){ his1[619]->Fill( invJpsiGamma.vect().mag() ); }

          }
          //===============================================

          //chi_c2  -------------------------------------
          if( invJpsiGamma.m()  > 3.53343  &&  invJpsiGamma.m() < 3.6  ){

            JpsiGammaFromChiC2=true;

            his1[650]->Fill(1);



#ifdef  MC_chi_c2
DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(S2K2piG.gammas[i]);
RecEmcShower *emcTrk = itTrk->emcShower();
S2K2piG.qgChiCshwrID.push_back( emcTrk->cellId() );
S2K2piG.PgFromChiC.push_back(S2K2piG.Pg[i]);
#endif




            //chi_c2 rec.
            his1[610]->Fill( invJpsiGamma.vect().mag() );
            his1[541]->Fill( RtoD( invJpsiGamma.phi() ) );
            his1[542]->Fill( invJpsiGamma.cosTheta() );

            //nGood
            his2[96]->Fill( invJpsiGamma.vect().mag(), invJpsiGamma.m()  );
            his2[97]->Fill( invJpsiGamma.vect().mag(),  nGood);
            his2[98]->Fill( invJpsiGamma.m() ,  nGood);
            his2[99]->Fill( invJpsiGamma.vect().mag(),  invJpsiGamma.cosTheta() );


            if( isMC && iMCchiC2ToJpsiGamma ){ his1[620]->Fill( invJpsiGamma.vect().mag() ); }

          }
          //===============================================



        }//end for(i)
      }// if(Jpsi1CFit)
    }//end if(2-track)


  }//end if(MuonsFromJpsi)


#ifdef  MC_chi_c1
  if( JpsiGammaFromChiC1 ){

    bool isGoodMCevent = true;
    isGoodMCevent = MatchRecMcEMC_gamma(selector, S2K2piG);

    if(isGoodMCevent){ chic12MomentumHistMC( selector, 1); }
  }
#endif



#ifdef  MC_chi_c2
  if( JpsiGammaFromChiC2 ){

    bool isGoodMCevent = true;
    isGoodMCevent = MatchRecMcEMC_gamma(selector, S2K2piG);

    if(isGoodMCevent){ chic12MomentumHistMC( selector, 1); }
  }
#endif







  return true;
}


//-------------------------------------------------------------------------------- EndJob() -------------------------------------------------------------------
void SelectJpsiXEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << " SelectJpsiXEndJob() " << endl;

  his2[119]->Divide(his2[118],his2[117],1.,1.,"B");

  if ( ERROR_WARNING != 0 ) {
    cout << " ERROR_WARNING= " << ERROR_WARNING << " Check output!" << endl;
  }

  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }

/*
  // print efficiency
  if ( his1[107] ) {
     double Nini = his1[107]->GetBinContent(1);
     double Nfin = his1[107]->GetBinContent( his1[107]->GetNbinsX() );
     double eff = Nfin/Nini * 100;
     double deff = sqrt(Nfin)/Nini * 100;
     printf(" efficiency= %6.2f +/- %5.2f\n",eff,deff);
  }
*/
}
#ifdef __cplusplus
}
#endif
