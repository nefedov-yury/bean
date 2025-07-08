//==================================================================//
//                                                                  //
// SelectKKgg - search for e+ e- -> K+ K- 2gammas                   //
//                                                                  //
//==================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TNtupleD.h>
#include <TDatabasePDG.h>

// #include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "AbsCor/AbsCor.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "EventTag/EventTagSvc.h"
#include "EventTag/DecayTable.h"
#include "RscanDQ/RscanDQ.h"

#include "DstEvtRecTracks.h"
// #include "TofHitStatus.h"
#include "ReadDst.h"

using namespace std;

//--------------------------------------------------------------------
// {{{1 Structure to save variables for a single event
//--------------------------------------------------------------------
struct Select_KKgg {
   // Run-info
   int runNo;                      // run-number
   int event;                      // event-number
   HepLorentzVector LVcms;         // Momentum in center of mass sys.
   Hep3Vector xorig;               // interaction point from DB

   // MC information:
   int decJpsi;                    // decay codes for MC events
   int dec_eta;                    // 1 if eta->2gamma
   double Xisr;                    // Xisr = s'/s
   double mc_mkk;                  // true inv. mass of K+K-
   double mc_mkpeta;               // true inv. mass of K+eta

   // Kaons candidate
   vector<RecMdcKalTrack*> trk_Kp; // positive tracks
   vector<RecMdcKalTrack*> trk_Km; // negative tracks

   // gamma tracks
   vector<RecEmcShower*> gtrk;
   vector<double> angt;            // angle with closest charged track
   vector<HepLorentzVector> Pg;    // 4-momentum of gammas
   vector<RecEmcShower*> g4f;      // two best photons after 4C-fit

   Select_KKgg() {
      decJpsi = -1;
      dec_eta = 0;
      Xisr = 0.;
      mc_mkk = 0;
      mc_mkpeta = 0;
      g4f.resize(2,nullptr);
   }
};

typedef Select_KKgg Select;

//--------------------------------------------------------------------
// {{{1 Global variables
//--------------------------------------------------------------------
const static double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
static const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
const static double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
const static double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
const static double meta   = 0.547862; // 547.862   +/- 0.017   MeV
const static double momega = 0.78265;  // 782.65    +/- 0.12    MeV
const static double mk     = 0.493677; // 493.677   +/- 0.016   MeV
static const double mk0    = 0.497611; // 497.611   +/- 0.013   MeV
const static double mphi   = 1.019461; //1019.461   +/- 0.019   MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1*> hst;
static std::vector<TNtupleD*> m_tuple;

// decays classification
static DecayTable JpsiTbl;

// container for warnings
static map<string,int> warning_msg;

static bool isMC = false;

// {{{1 Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif

//--------------------------------------------------------------------
// static inline void Warning(const char* msg) {
//    warning_msg[string(msg)] += 1;
static inline void Warning(string msg) {
   warning_msg[msg] += 1;
}
//--------------------------------------------------------------------

//--------------------------------------------------------------------
static inline Double_t RtoD(Double_t ang) {
   return ang*180/M_PI;
}
//--------------------------------------------------------------------

//--------------------------------------------------------------------
static inline double SQ(double x) {
   return x*x;
}
//--------------------------------------------------------------------

// {{{1 StartJob, book histograms
//--------------------------------------------------------------------
BeanUserShared_EXPORT
void SelectKKggStartJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
      printf("Masses of particles from PDG\n");
      printf("M_J/Psi   = %f MeV\n", mjpsi*1e3);
      printf("M_pi^\\pm  = %f MeV\n", mpi*1e3);
      printf("M_pi^0    = %f MeV\n", mpi0*1e3);
      printf("M_eta     = %f MeV\n", meta*1e3);
      printf("Momega    = %f MeV\n", momega*1e3);
      printf("M_K^\\pm   = %f MeV\n", mk*1e3);
      printf("M_K^0     = %f MeV\n", mk0*1e3);
      printf("M_phi     = %f MeV\n", mphi*1e3);
   }

   hst.resize(500,nullptr);
   m_tuple.resize(10,nullptr);

   // initialize Absorption Correction -------------------------------
   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

   // initialize EventTag --------------------------------------------
   m_EventTagSvc = EventTagSvc::instance();
   if ( !m_EventTagSvc->IsInitialized() ) {
      // set paths to pdg & decayCodes files:
      string pdtFile("Analysis/EventTag/share/pdt_bean.table");
      m_EventTagSvc->setPdtFile( selector->AbsPath(pdtFile) );
      string dcFile("Analysis/EventTag/share/DecayCodes/"
            "dcode_charmonium.txt");
      m_EventTagSvc->setDecayTabsFile( selector->AbsPath(dcFile) );
      m_EventTagSvc->setIgnorePhotons(false); // ignore ISR & FSR
      if ( selector->Verbose() ) {
         m_EventTagSvc->setVerbose(1);
      }
      m_EventTagSvc->initialize();
   } else {
      cout << " WARNING in " << __func__ << ": "
           << "EventTagSvc has already been initialized" << endl;
      Warning("EventTagSvc has already been initialized");
   }

   // initialize DatabaseSvc -----------------------------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      string dirDB("Analysis/DatabaseSvc/dat");
      dbs->SetDBFilePath( selector->AbsPath(dirDB) );
   }

   // initialize Magnetic field --------------------------------------
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if ( (mf->GetPath()).empty() ) {
      // set path to directory with magnetic fields tables
      mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
      mf->UseDBFlag(false); // like in the boss program
      mf->RunMode(3); // like in the boss program
   } else {
      cout << " WARNING:"
           << " MagneticFieldSvc has already been initialized" << endl
           << "          path = " << mf->GetPath() << endl;
      Warning("MagneticFieldSvc has already been initialized");
   }

   // set path for ParticleID algorithm ------------------------------
   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   // Book histograms --------- --------------------------------------

   hst[0] = new TH1D("Runs", "run numbers", 1000,0.,50000);
   hst[1] = new TH1D("cuts_0","preselections cuts", 20,-0.5,19.5);
   hst[9] = new TH1D("Ecms","E_{cms}", 600, 2.85, 3.15);

   //  ChargedTracksKK:
   // angle between tracks of one charge
   hst[6] = new TH1D("ang_pip","cos(ang) between all pairs of K+",
         200,-1.,1.);
   hst[7] = new TH1D("ang_pim","cos(ang) between all pairs of K-",
         200,-1.,1.);

   hst[11] = new TH1D("Rxy","R_{xy}", 200,-2.,2.);
   hst[12] = new TH1D("Rz","R_{z}", 200,-20.,20.);
   hst[13] = new TH1D("cos_theta","cos(#theta)", 200,-1.,1.);
   hst[14] = new TH1D("theta","#theta", 180,0.,180.);
   hst[15] = new TH1D("PQ","Charged momentum (Q*P)", 600,-3.,3.);

   hst[21] = new TH1D("Pid_clK","lg(CL_{K})", 100,-4.,0.);
   hst[22] = new TH1D("Pid_isK","1(K),0(other) particle",2,-0.5,1.5);
   hst[23] = new TH1D("K_Q","(Q*P) for K", 300,-1.5,1.5);
   hst[24] = new TH2D("K_pm","N(K+) vs N(K-)",
         10,-0.5,9.5, 10,-0.5,9.5);
   hst[25] = new TH1D("K_N","N(K-)+N(K+)", 10,-0.5,9.5);
   hst[26] = new TH1D("N_o","N(other)", 10,-0.5,9.5);
   hst[27] = new TH1D("N_ok2","N(other) N(K)=2", 10,-0.5,9.5);

   // NeutralTracks:
   hst[31] = new TH1D("G_ang","angle with closest chg.trk",
         180,0.,180.);
   hst[32] = new TH1D("G_n","N_{#gamma} in event", 11,-0.5,10.5);
   hst[33] = new TH1D("K_PM","momentum K+ and K-", 240,-1.2,1.2);
   hst[34] = new TH1D("G_P","Momentum of gammas", 200,0.,2.);

   // VertKinFit:
   hst[41] = new TH1D("vtx_fit", "vertex fit: bad(0)/good(1)",
         2,-0.5,1.5);
   hst[42] = new TH1D("vtx_chi2", "vertex fit #chi^2", 100,0.,100.);
   hst[51] = new TH1D("f4c_chi2","4C-fit: min #chi^{2}", 300,0.,300.);
   hst[52] = new TH1D("fit_ch2_3delta",
                      "fit: #chi^{2}(2#gamma) - #chi^{2}(3#gamma)",
                      200,-100.,100.);
   hst[53] = new TH1D("fit_ch2_3","fit: #chi^{2}(3#gamma)",
         250,0.,250.);
   hst[54] = new TH1D("fit_ch2_2","fit: #chi^{2}(2#gamma)",
         250,0.,250.);

   hst[61] = new TH1D("f4c_ang","4C-fit: min angle gamma trk",
         180,0.,180.);
   hst[62] = new TH1D("fit_eg","fit: E(#gamma) rejected",
         500,0.,0.5);
   hst[63] = new TH1D("f4c_egm","4C-fit: E_{max}(#gamma) dropped",
         100,0.,2.);

   hst[71] = new TH2D("Mgg_Mkk", "M^{2}(K^{+}K^{-}) "
         "vs M^{2}(#gamma#gamma)", 220,0.,1.1, 200,0.,10.);
   hst[72] = new TH1D("Cphi_rest","cos #Theta(#phi) in J/Psi",
         200,-1.0,1.0);

   // Monte Carlo histograms:
   hst[100] = new TH1D("mc_dcj0", "dec J/psi nocut",300,-0.5,299.5);
   hst[101] = new TH1D("mc_dcj1", "dec J/psi cut#1",300,-0.5,299.5);
   hst[102] = new TH1D("mc_dcj2", "dec J/psi cut#2",300,-0.5,299.5);
   hst[104] = new TH1D("mc_dcj4", "dec J/psi cut#4",300,-0.5,299.5);

   // ISRhistoMC:
   hst[111] = new TH1D("mcisr_pdg", "PDG codes of all particles",
         2001,-1000.5,1000.5);
   hst[112] = new TH1D("mcisr_pdg0", "PDG of particles "
         "from primary vertex", 2001,-1000.5,1000.5);
   hst[113] = new TH1D("mcisr_Eg", "ISR: energy of one gamma (MeV)",
         1000,0.,1000.);
   hst[114] = new TH1D("mcisr_Ng", "ISR: N_{#gamma}", 10,-0.5,9.5);
   hst[115] = new TH1D("mcisr_Etot", "ISR: Eisr (MeV)",1000,0.,3000.);
   hst[116] = new TH1D("mcisr_xisr", "ISR: s'/s",100,0.,1.);
   hst[117] = new TH1D("mcisr_1x", "ISR: 1 - s'/s",100,0.,1.);
   hst[118] = new TH1D("mcisr_bspr", "MC: energy spread (MeV)",
         1000,-5.,5.);
   hst[119] = new TH2D("mcisr_Eg1_Eg2","Eg1 vs Eg2",
         1500,0.,1500., 1500,0.,1500.);

   // initial number of events for efficiency: Mkk(MC) for Xisr>0.9
   hst[121] = new TH1D("mc_Mkk_ini", "Minv(K^{+}K^{-})", 100,0.98,1.08);

   // FillHistoMC:
   hst[131] = new TH1D("mc_EtaP", "Momentum of #eta", 1000,0.,2.);
   hst[132] = new TH1D("mc_EtaPt","Pt of #eta", 1000,0.,2.);
   hst[133] = new TH1D("mc_EtaC", "cos(#Theta) of #eta", 100,-1.,1.);
   hst[135] = new TH1D("mc_PhiP", "Momentum of #phi", 1000,0.,2.);
   hst[136] = new TH1D("mc_PhiPt","Pt of #phi", 1000,0.,2.);
   hst[137] = new TH1D("mc_PhiC", "cos(#Theta) of #phi", 100,-1.,1.);

   hst[141] = new TH1D("mc_Kp", "Momentum of K^{+}", 1000,0.,2.);
   hst[142] = new TH1D("mc_KCp","cos(#Theta) of K^{+}", 100,-1.,1.);
   hst[143] = new TH1D("mc_Km", "Momentum of K^{-}", 1000,0.,2.);
   hst[144] = new TH1D("mc_KCm","cos(#Theta) of K^{-}", 100,-1.,1.);

   hst[146] = new TH1D("mc_Mkk", "Minv(K^{+}K^{-})", 140,0.98,1.12);
   hst[147] = new TH1D("mc_Mkpeta", "Minv(K^{+}eta)", 140,1.8,2.5);
   hst[148] = new TH2D("mc_M1M2", "Minv2(K+K-) vs Minv2(K+ eta)",
         140,0.,7.,140,0.,7.);
   hst[149] = new TH1D("mc_deceta", "1 if eta->2gamma", 2,-0.5,1.5);

   // ntuple for e+ e- -> K+ K- 2gammas
   m_tuple[0] = new TNtupleD("a4c","after 4C kinematic fit",
         "ch2:chsq3g:"     // chi^2 of 4C fit, chi^2 of 3gammas
         "Pkp:Ckp:phikp:"  // P,cos(Theta),phi of K+
         "Pkm:Ckm:phikm:"  // P,cos(Theta),phi of K-
         "Eg1:Cg1:phig1:"  // E,cos(Theta),phi of gamma-1
         "Eg2:Cg2:phig2:"  // E,cos(Theta),phi of gamma-2
         "Pgg:Cgg:phigg:"  // P,cos(Theta),phi of gamma,gamma
         "Mkk:Mgg:"        // invariant masses of K+K- and 2gammas
         "M2kpeta:M2kmeta" // M_inv^2( Keta )
         ":dec"            // MCtrue: decay codes of J/Psi
         ":xisr"           // MCtrue: Xisr = s'/s
         ":mcmkk"          // MCtrue: inv. mass of K+K-
         ":mcmkpet"        // MCtrue: inv. mass of K+eta
         );
//          "M2kpg1:M2kpg2:M2kmg1:M2kmg2:" // M_inv^2( Kg )

   // register in selector to save in given directory
   const char* SaveDir = "SelectKKgg";
   VecObj his1o(hst.begin(),hst.end());
   selector->RegInDir(his1o,SaveDir);
   VecObj ntuples(m_tuple.begin(),m_tuple.end());
   selector->RegInDir(ntuples,SaveDir);
}

// {{{1 getVertexOrigin() && GetEcms()
//--------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//--------------------------------------------------------------------
   // cache value for one run
   static int save_runNo = 0;
   static Hep3Vector xorigin;
   if ( runNo == save_runNo ) {
      return xorigin;
   }

   // update vertex for new run
   xorigin.set(0.,0.,0.);
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();

   int run = abs(runNo);
   if (
         (run >= 28241 && run <= 28266)  // 3080 no 6.6.4
         || (run >=  9947 && run <= 10878)  // J/Psi 2009 no 6.6.4
//          || (run >= 27147 && run <= 27233)  // 3082-old don't use it!
      ) {
      vtxsvc->SetBossVer("6.6.3");
   } else if ( run < 30000 ) {   // here: J/Psi scan 2012
      vtxsvc->SetBossVer("6.6.4");
   } else if ( run >=39355 && run <= 40069 ) { // R-scan (2015)
#if (BOSS_VER < 700)
      vtxsvc->SetBossVer("6.6.5.p01");
#else
      vtxsvc->SetBossVer("7.0.3");
#endif
   } else if ( run >=55060 && run <= 55109 ) { // J/Psi scan 2018
      vtxsvc->SetBossVer("7.0.4");
   } else if ( run >=59016 && run <= 59141 ) { // 3080 data 2019
      vtxsvc->SetBossVer("7.0.4");
   } else {
      cout << " WARNING:" << __func__ << " :" << " run=" << runNo
         << " We will use default Boss Version: "
         << vtxsvc->GetBossVer()
         << endl;
      Warning( string("getVertexOrigin: use default Boss Version:")
            + vtxsvc->GetBossVer() );
   }
   vtxsvc->handle(runNo);

   if ( vtxsvc->isVertexValid() ) {
      double* dbv = vtxsvc->PrimaryVertex();
      xorigin.set(dbv[0],dbv[1],dbv[2]);
      if( verbose ) {
         cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
      }
   } else {
      cout << " FATAL ERROR:"
           " Cannot obtain vertex information for run#"
           << runNo << endl;
      exit(1);
   }

   save_runNo = runNo;
   return xorigin;
}

// GetEcms: return negative value for bad runs
//--------------------------------------------------------------------
static double GetEcms(int runNo, bool verbose = false) {
//--------------------------------------------------------------------
   struct runInfo {
      int runS, runE; // first and last runs in period
      double Ebeam;   // energy of beam (MeV)
   };

   // List Runs for J/Psi-scan 2012 see arXiv:2206.13674v1
   // for 3080 take old -0.55MeV
   static const runInfo ListRuns[] =  {
      { 28312, 28346,     3049.642},
      { 28347, 28381,     3058.693},
      { 28241, 28266,     3079.645},
      { 28382, 28387,     3082.496},
      { 28466, 28469,     3082.496},
      { 28388, 28416,     3088.854},
      { 28472, 28475,     3088.854},
      { 28417, 28453,     3091.760},
      { 28476, 28478,     3091.760},
      { 28479, 28482,     3094.697},
      { 28487, 28489,     3095.430},
      { 28490, 28492,     3095.826},
      { 28493, 28495,     3097.213},
      { 28496, 28498,     3098.340},
      { 28499, 28501,     3099.042},
      { 28504, 28505,     3101.359},
      { 28506, 28509,     3105.580},
      { 28510, 28511,     3112.051},
      { 28512, 28513,     3119.878},
   };
   static const int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

   // List Runs for J/Psi scan 2018: 19.04-87006-bes3memo-v1.1.pdf
   static const runInfo ListJ2018[] =  {
      { 55060, 55065,     3087.593},
      { 55066, 55073,     3095.726},
      { 55074, 55083,     3096.203},
      { 55084, 55088,     3096.986},
      { 55089, 55091,     3097.226},
      { 55092, 55097,     3097.654},
      { 55098, 55103,     3098.728},
      { 55104, 55109,     3104.000},
      { 59016, 59141,     3080.   }, // 3080 data 2019
   };
   static const int Np3 = sizeof(ListJ2018)/sizeof(ListJ2018[0]);

   // List Runs for R-scan 2015 + RscanDQ
   static const runInfo ListRscan[] =  {
      { 39355, 39618,     3080.},
      { 39619, 39650,     2950.},
      { 39651, 39679,     2981.},
      { 39680, 39710,     3000.},
      { 39711, 39738,     3020.},
      { 39775, 40069,     2900.},
   };
   static const int Np2 = sizeof(ListRscan)/sizeof(ListRscan[0]);

   // cache value for one run
   static int save_runNo = 0;
   static double Ecms = 3.097;
   if ( runNo == save_runNo ) {
      return Ecms;
   }

   int absrun = abs(runNo);
   bool found = false;

   // search in J/Psi-scan 2012
   if ( !found ) {
      for(int i = 0; i < Np; i++) {
         if ( absrun >= ListRuns[i].runS &&
               absrun <= ListRuns[i].runE ) {
            Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV
            found  = true;
            break;
         }
      }
   }

   // search in J/Psi scan 2018
   if ( !found ) {
      for(int i = 0; i < Np3; i++) {
         if ( absrun >= ListJ2018[i].runS &&
               absrun <= ListJ2018[i].runE ) {
            Ecms   = ListJ2018[i].Ebeam  * 1.e-3;  // MeV -> GeV
            found  = true;
            break;
         }
      }
   }

   // search in R-scan
   if ( !found ) {
      for(int i = 0; i < Np2; i++) {
         if ( absrun >= ListRscan[i].runS &&
               absrun <= ListRscan[i].runE ) {
            Ecms   = ListRscan[i].Ebeam  * 1.e-3;  // MeV -> GeV
            found  = true;
            if ( runNo > 0 ) { // data
               RscanDQ rdq = RscanDQ(runNo);
               double Ebeam = rdq.getEbeam(); // beam energy in GeV
               int status = rdq.getStatus();

               if ( fabs(2*Ebeam-Ecms) > 1e-4 ) {
                  cout << " GetEcms::WARNING "
                     << "RscanDQ problem for run= " << runNo
                     << " Ebeam= " << Ebeam
                     << " Ecms= " << Ecms << endl;
                  Warning("RscanDQ problem");
               }

               if ( status != 1 ) {
                  cout << " GetEcms::RscanDQ bad run="
                     << runNo << endl;
                  Ecms = -Ecms;
               }
            }
            break;
         }
      }
   }

   if ( !found ) {
      // run not in the list
      if ( absrun >= 9947 && absrun <= 10878 ) { // J/Psi 2009
         Ecms =  3.097;
      } else {
         cout << " GetEcms::WARNING unknown run# " << runNo
              << " use Ecms= " << Ecms << endl;
         Warning("Unknown run");
      }
   }

   if( verbose ) {
      cout << " GetEcms: run= " << runNo << " -> "
           << " Ecms= " << Ecms << endl;
   }

   save_runNo = runNo;
   return Ecms;
}

/*
// GetEcms: return negative value for bad runs
//--------------------------------------------------------------------
static double GetEcmsOLD(int runNo, bool verbose = false) {
//--------------------------------------------------------------------
   struct runInfo {
      int runS, runE; // first and last runs in period
      double Ebeam;   // energy of beam (MeV)
      double Spread;  // beam spread (KeV)
      double Lumi;    // luminosity (pb^-1)
   };

   // List Runs for phasejpsiscan production
   static const runInfo ListRuns[] =  {
      { 28312, 28346,     3050.213,    0,     14.919 },
      { 28347, 28381,     3059.257,    0,     15.060 },
      { 28241, 28266,     3080.195,    0,     17.528 },
      { 28382, 28387,     3083.060,    0,      4.769 },
      { 28466, 28469,     3083.060,    0,      4.769 },
      { 28388, 28416,     3089.418,    0,     15.558 },
      { 28472, 28475,     3089.418,    0,     15.558 },
      { 28417, 28453,     3092.324,    0,     14.910 },
      { 28476, 28478,     3092.324,    0,     14.910 },
      { 28479, 28482,     3095.261,    0,      2.143 },
      { 28487, 28489,     3095.994,    0,      1.816 },
      { 28490, 28492,     3096.390,    0,      2.135 },
      { 28493, 28495,     3097.777,    0,      2.069 },
      { 28496, 28498,     3098.904,    0,      2.203 },
      { 28499, 28501,     3099.606,    0,      0.756 },
      { 28504, 28505,     3101.923,    0,      1.612 },
      { 28506, 28509,     3106.144,    0,      2.106 },
      { 28510, 28511,     3112.615,    0,      1.720 },
      { 28512, 28513,     3120.442,    0,      1.264 },
   };
   static const int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

   // List Runs for R-scan 2015 + RscanDQ
   static const runInfo ListRscan[] =  {
      { 39355, 39618,     3080.,    0,     126.21  },
      { 39619, 39650,     2950.,    0,      15.96  },
      { 39651, 39679,     2981.,    0,      16.046 },
      { 39680, 39710,     3000.,    0,      15.849 },
      { 39711, 39738,     3020.,    0,      17.315 },
      { 39775, 40069,     2900.,    0,     105.53  },
   };
   static const int Np2 = sizeof(ListRscan)/sizeof(ListRscan[0]);

   // List Runs for J/Psi scan 15th — 18th April, 2018
   //               + 3080 data 2019: 7th-11th Feb 2019
   static const runInfo ListJ2018[] =  {
      { 55060, 55065,     3087.659, 1.312, 2.50166 },
      { 55066, 55073,     3095.726, 1.058, 2.96535 },
      { 55074, 55083,     3096.203, 1.132, 5.10588 },
      { 55084, 55088,     3096.986, 1.009, 3.07306 },
      { 55089, 55091,     3097.226, 1.270, 1.70679 },
      { 55092, 55097,     3097.654, 1.047, 4.78731 },
      { 55098, 55103,     3098.728, 1.183, 5.75102 },
      { 55104, 55109,     3104.000, 1.010, 5.8398  },
      { 59016, 59141,     3080.,    0.,    83.935  },
   };
   static const int Np3 = sizeof(ListJ2018)/sizeof(ListJ2018[0]);

   // cache value for one run
   static int save_runNo = 0;
   static double Ecms = 3.097;
   if ( runNo == save_runNo ) {
      return Ecms;
   }

   int absrun = abs(runNo);
   bool found = false;

   // search in phasejpsiscan
   for(int i = 0; i < Np; i++) {
      if ( absrun >= ListRuns[i].runS &&
            absrun <= ListRuns[i].runE ) {
         Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV

         // correct energy according BAM-00268: -0.55 +/- 0.03 MeV
         Ecms -= 0.55e-3;

//          Lumi   = ListRuns[i].Lumi;
         found  = true;
         break;
      }
   }

   // search in R-scan
   if ( !found ) {
      for(int i = 0; i < Np2; i++) {
         if ( absrun >= ListRscan[i].runS &&
               absrun <= ListRscan[i].runE ) {
            Ecms   = ListRscan[i].Ebeam  * 1.e-3;  // MeV -> GeV
//             Lumi   = ListRscan[i].Lumi;
            found  = true;
            if ( runNo > 0 ) { // data
               RscanDQ rdq = RscanDQ(runNo);
               double Ebeam = rdq.getEbeam(); // beam energy in GeV
               int status = rdq.getStatus();

               if ( fabs(2*Ebeam-Ecms) > 1e-4 ) {
                  cout << " GetEcms::WARNING "
                     << "RscanDQ problem for run= " << runNo
                     << " Ebeam= " << Ebeam
                     << " Ecms= " << Ecms << endl;
                  Warning("RscanDQ problem");
               }

               if ( status != 1 ) {
                  cout << " GetEcms::RscanDQ bad run="
                     << runNo << endl;
                  Ecms = -Ecms;
               }
            }
            break;
         }
      }
   }

   // search in J/Psi scan 2018
   if ( !found ) {
      for(int i = 0; i < Np3; i++) {
         if ( absrun >= ListJ2018[i].runS &&
               absrun <= ListJ2018[i].runE ) {
            Ecms   = ListJ2018[i].Ebeam  * 1.e-3;  // MeV -> GeV
            found  = true;
            break;
         }
      }
   }

   if ( !found ) {
      // run not in the list
      if ( absrun >= 9947 && absrun <= 10878 ) { // J/Psi 2009
         Ecms =  3.097;
      } else {
         cout << " GetEcms::WARNING unknown run# " << runNo
              << " use Ecms= " << Ecms << endl;
         Warning("Unknown run");
      }
   }

   if( verbose ) {
      cout << " GetEcms: run= " << runNo << " -> "
           << " Ecms= " << Ecms
//            << " Spread= " << Spread
//            << " Lumi= " << Lumi
           << endl;
   }

   save_runNo = runNo;
   return Ecms;
}
*/

// {{{1 FillHistoMC() && ISRhistoMC()
//--------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, Select& Slct) {
//--------------------------------------------------------------------
   if ( !isMC ) {
      return;
   }

   const TMcEvent*  m_TMcEvent  = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   // EventTag
   m_EventTagSvc->setMcEvent(const_cast<TMcEvent*>(m_TMcEvent));
   unsigned int evTag = m_EventTagSvc->getEventTag();
   // check general event type
   if ( (evTag & 0xF) != 4 ) { // 4 is J/Psi event
      static int nprt = 0;
      if ( nprt < 1 ) {
         printf(" WARNING: MC is not J/Psi evTag= 0x%08x\n",evTag);
         ++nprt;
      }
      Warning("MC is not J/Psi: check that it is MCGPJ!");
      evTag = 0;
   }

   // decay code of J/Psi:
   int decJpsi  = (evTag >> 8) & 0xFF;
   if ( evTag == 0 ) { //MCGPJ
      decJpsi = 261; // K+ K- eta (+ISR gammas)
   }
   Slct.decJpsi = decJpsi;
   hst[100]->Fill(decJpsi);

   int idx_jpsi=-1;
   if ( evTag == 0 ) { //MCGPJ
      idx_jpsi = -99;  // == primary vertex
   }
   JpsiTbl.Reset();

   // momentum eta & phi -> K+,K-
   int idx_eta=-1;
   int idx_phi=-1;
   HepLorentzVector LVeta, LVphi;
   vector<HepLorentzVector> LVKp, LVKm;

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      if ( part_pdg == 443 ) { // J/Psi
         idx_jpsi = part->getTrackIndex();
      }

      if( part->getMother() == idx_jpsi ) {  // decays of J/Psi
         // collect decays of J/psi
         JpsiTbl.vecdec.push_back(part_pdg);

         if ( !(decJpsi == 261 || decJpsi == 68) ) {
            continue; // skip if it is not Jpsi -> phi eta
         }
         // MCGPJ generates "K+ K- eta" in the primary vertex
         if ( part_pdg == 221 ) {               // eta
            hst[131]->Fill(Vp.mag());
            hst[132]->Fill(Vp.rho());
            hst[133]->Fill(Vp.cosTheta());
            idx_eta = part->getTrackIndex();
            HepLorentzVector LV( Vp, sqrt(Vp.mag2()+SQ(meta)) );
            LVeta = LV;
         } else if ( part_pdg == 333 ) {        // phi
            hst[135]->Fill(Vp.mag());
            hst[136]->Fill(Vp.rho());
            hst[137]->Fill(Vp.cosTheta());
            idx_phi = part->getTrackIndex();

            HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(mphi)) );
            LVphi = LV;
         } else if ( abs(part_pdg) == 321 ) {   // K+ or K-
            HepLorentzVector LV( Vp, sqrt(Vp.mag2()+SQ(mk)) );
            if ( part_pdg > 0 ) {
               hst[141]->Fill( Vp.mag() );
               hst[142]->Fill(Vp.cosTheta());
               LVKp.push_back(LV);
            } else {
               hst[143]->Fill( Vp.mag() );
               hst[144]->Fill(Vp.cosTheta());
               LVKm.push_back(LV);
            }
         }
      }

      if ( part->getMother() == idx_phi ) { // J/Psi -> phi eta
                                            //          |-> K+ K-
         if ( abs(part_pdg) == 321 ) {      // K+ or K-
            HepLorentzVector LV( Vp, sqrt(Vp.mag2()+SQ(mk)) );
            if ( part_pdg > 0 ) {
               hst[141]->Fill( Vp.mag() );
               hst[142]->Fill(Vp.cosTheta());
               LVKp.push_back(LV);
            } else {
               hst[143]->Fill( Vp.mag() );
               hst[144]->Fill(Vp.cosTheta());
               LVKm.push_back(LV);
            }
         }
      }

      if ( part->getMother() == idx_eta ) { // J/Psi -> phi eta
                                            //              |-> 2gamma
         if ( part_pdg != 22 ) {  // it's not 2gamma decay
            Slct.dec_eta = -1;
         }
      }
   } // end of while

   if ( decJpsi == 261 || decJpsi == 68 ) { // J/Psi -> phi eta
      // save invariant masses: M(K+K-), M(K+eta)
      if ( LVKp.size() == 1 && LVKm.size() == 1  ) {
         Slct.mc_mkk = (LVKp[0]+LVKm[0]).m();
         hst[146]->Fill( Slct.mc_mkk );
         Slct.mc_mkpeta = (LVKp[0]+LVeta).m();
         hst[147]->Fill( Slct.mc_mkpeta );
         hst[148]->Fill( SQ(Slct.mc_mkk), SQ(Slct.mc_mkpeta) );
      }

      // save Slct.dec_eta
      if ( Slct.dec_eta == 0 ) {
         Slct.dec_eta = 1;  // eta -> 2 gamma decay
      } else {
         Slct.dec_eta = 0;
      }
      hst[149]->Fill( Slct.dec_eta );
   }
}

//--------------------------------------------------------------------
static void ISRhistoMC(const ReadDst* selector, Select& Slct) {
//--------------------------------------------------------------------
   // summary for ISR according of MCGPJ model: 1 or 2 ISR photons
   // calculate Xisr = s'/s,
   // where s is center of mass system energy squared
   // and s' is energy after ISR

   if ( !isMC ) {
      return;
   }

   const TMcEvent*  m_TMcEvent  = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   HepLorentzVector Pr;        // sum of primary particles

   double Eg_isr[2] = {0.,0.}; // energy of first two gammas
   int    Ng_isr    = 0;       // number of gammas
   double Etot_isr  = 0.;      // total energy
   HepLorentzVector Pisr;      // sum of ISR photons


   TIter mcParticlesIter(mcParticles);
   while ( TMcParticle* part =
              static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      hst[111]->Fill(part_pdg);
      if ( part->getMother() == -99 ) {
         hst[112]->Fill(part_pdg);
         HepLorentzVector Pp = HepLorentzVector(
                                  part->getInitialMomentumX(),
                                  part->getInitialMomentumY(),
                                  part->getInitialMomentumZ(),
                                  part->getInitialMomentumE()
                               ) * 1e3; // ->MeV
         Pr += Pp;

//       cout << " PDG= " << part_pdg << " Pp= " << Pp << endl;

         if ( part_pdg == 22 ) { // gamma from primary vertex => ISR
            double Eg = part->getInitialMomentumE() * 1e3; // ->MeV
            Etot_isr += Eg;
            if( Ng_isr < 2 ) {
               Eg_isr[Ng_isr] = Eg;
            }
            Ng_isr++;

            hst[113]->Fill(Eg);

            Pisr += Pp;
         }
      }
   } // end of while()

   hst[114]->Fill(Ng_isr);
   hst[115]->Fill(Etot_isr);
   if ( Ng_isr == 2 ) {
      hst[119]->Fill(Eg_isr[0],Eg_isr[1]);
   }

   double S = Pr*Pr;
   HepLorentzVector Prp = Pr - Pisr;
   double Sp = Prp*Prp;
   Slct.Xisr = Sp/S;
   hst[116]->Fill(Slct.Xisr);
   hst[117]->Fill(1.-Slct.Xisr);
   hst[118]->Fill( sqrt(S) - Slct.LVcms.e() * 1e3);

   // initial number of events for efficiency: Mkk(MC) for Xisr>0.9
   if ( Slct.Xisr > 0.9 ) {
      hst[121]->Fill( Slct.mc_mkk );
   }
}


// {{{1 select K+K- candidates, no other tracks.
//--------------------------------------------------------------------
static bool ChargedTracksKK(ReadDst* selector, Select& Slct) {
//--------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
//    static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   int Nother = 0; // other good tracks (not K)

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.); // initial point for MDC rec.
      HepPoint3D IP(Slct.xorig[0],Slct.xorig[1],Slct.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      // the nearest distance to IP
      double Rvxy0 = vecipa[0]; // in xy plane
      double Rvz0  = vecipa[3]; // in z direction

      hst[11]->Fill(Rvxy0);
      hst[12]->Fill(Rvz0);
      hst[13]->Fill(cosTheta);
      hst[14]->Fill(RtoD(theta));
      hst[15]->Fill( mdcTrk->charge()*mdcTrk->p() );

      if( fabs(Rvxy0) >= Rvxy0_max ) {
         continue;
      }
      if( fabs(Rvz0) >= Rvz0_max ) {
         continue;
      }
      if ( fabs(cosTheta) >= cosTheta_max ) {
         continue;
      }

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      if ( !mdcKalTrk ) {
         continue;
      }

      // PID information
      pid->init();
      pid->setMethod(pid->methodProbability());
      pid->setChiMinCut(4);
      pid->setRecTrack(itTrk);

      // list of systems used for PID
      pid->usePidSys( pid->useDedx() |
                      pid->useTof1() | pid->useTof2() |
                      pid->useTofE() );

      pid->identify(
         pid->onlyPion()    |
         pid->onlyKaon()    |
         pid->onlyProton()  |
         pid->onlyMuon()    |
         pid->onlyElectron()
      );

      pid->calculate(Slct.runNo);
      int check_kaon = 0;
      if ( pid->IsPidInfoValid() ) {
         if( pid->probKaon() > 0 ) {
            hst[21]->Fill( log10(pid->probKaon()) );
         } else {
            hst[21]->Fill(-5.);
         }

         // check that track is kaon:
         if ( (pid->probKaon() > 0.001)             &&
               (pid->probKaon() > pid->probPion())   &&
               (pid->probKaon() > pid->probProton())
            ) {
            check_kaon = 1;
         }
      } // end IsPidInfoValid
      hst[22]->Fill( double(check_kaon) );
      if( check_kaon == 0 ) {
         Nother++;
         continue;
      }

      mdcKalTrk->setPidType(RecMdcKalTrack::kaon);
      hst[23]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      if( mdcKalTrk->charge() > 0 ) {
         Slct.trk_Kp.push_back(mdcKalTrk);
      } else {
         Slct.trk_Km.push_back(mdcKalTrk);
      }
   }  //-------------------------------------------------End for (i)

   int np = Slct.trk_Kp.size();
   int nm = Slct.trk_Km.size();
   hst[24]->Fill(np,nm);
   hst[25]->Fill(np+nm);
   hst[26]->Fill(Nother);

   // check angle between tracks of one charge
   for ( int i = 1; i < np; ++i ) {
      const auto& t1 = Slct.trk_Kp[i];
      Hep3Vector Vp1(t1->px(), t1->py(), t1->pz());
      for ( int ii = 0; ii < i; ++ii ) {
         const auto& t2 = Slct.trk_Kp[ii];
         Hep3Vector Vp2(t2->px(), t2->py(), t2->pz());
         double ang = Vp1.angle(Vp2);
         hst[6]->Fill(cos(ang));
      }
   }
   for ( int i = 1; i < nm; ++i ) {
      const auto& t1 = Slct.trk_Km[i];
      Hep3Vector Vp1(t1->px(), t1->py(), t1->pz());
      for ( int ii = 0; ii < i; ++ii ) {
         const auto& t2 = Slct.trk_Km[ii];
         Hep3Vector Vp2(t2->px(), t2->py(), t2->pz());
         double ang = Vp1.angle(Vp2);
         hst[7]->Fill(cos(ang));
      }
   }

   // require exactly one "+" and one "-"
   if ( np != 1 || nm != 1 ) {
      return false;
   }
   hst[27]->Fill(Nother);
   // no other good tracks from primary vertex
   if ( Nother > 0 ) {
      return false;
   }

   hst[1]->Fill(1); // "cuts"
   if ( isMC ) {
      hst[101]->Fill( Slct.decJpsi );
   }

   return true;
}

// {{{1 select gammas candidates
//--------------------------------------------------------------------
static bool NeutralTracks(ReadDst* selector, Select& Slct) {
//--------------------------------------------------------------------
   // parameters of reconstruction
   static const double min_angle = 10 * M_PI/180; // 10 grad

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   for ( int i = evtRecEvent->totalCharged();
         i < evtRecEvent->totalTracks(); i++ ) {
      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if ( !itTrk->isEmcShowerValid() ) {
         continue;
      }
      RecEmcShower* emcTrk = itTrk->emcShower();

      // *) good EMC time:
      if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) {
         continue;
      }

      // *) good EMC energy deposited in the barrel (endcap)
      double eraw = emcTrk->energy();
      double absCosTheta = fabs(  cos(emcTrk->theta()) );

      bool GoodCluster=false;
      if ( absCosTheta < 0.8 ) {  //barrel
         GoodCluster = eraw > 25E-3;
      } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) {//endcap
         GoodCluster = eraw > 50E-3;
      }
      if ( !GoodCluster ) {
         continue;
      }

      // *) the nearest charged track is far from cluster
      Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

      double tang = 200.; // min angle between cluster and track
      for ( int j = 0; j < evtRecEvent->totalCharged(); j++ ) {
         DstEvtRecTracks* jtTrk =
            static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(j));
         if ( !jtTrk->isExtTrackValid() ) {
            continue;
         }
         RecExtTrack* extTrk = jtTrk->extTrack();
         if ( extTrk->emcVolumeNumber() == -1 ) {
            continue;   //track does not hit EMC
         }
         Hep3Vector extpos = extTrk->emcPosition();
         double angd = fabs(extpos.angle(emcpos)); // [0,pi]
         if ( angd < tang ) {
            tang = angd;
         }
      } //--------------------------------------------End for(j)

      hst[31]->Fill(RtoD(tang));
      if ( tang < min_angle ) {
         continue;
      }

      Slct.gtrk.push_back(emcTrk);
      Slct.angt.push_back(tang);

      // (E,Px,Py,Pz) for good gammas
      Hep3Vector p3 = emcpos - Slct.xorig;
      p3 *= eraw / p3.mag();
      Slct.Pg.push_back( HepLorentzVector(p3,eraw) );
   } //-----------------------------------------------End for (i)

   int ng = Slct.gtrk.size();
   hst[32]->Fill(ng);

   if ( ng < 2 ) {
      return false;
   }

   hst[1]->Fill(2); // "cuts"
   if ( isMC ) {
      hst[102]->Fill( Slct.decJpsi );
   }

   // momentum of selected pair of Kaons and Gammas
   hst[33]->Fill( Slct.trk_Kp[0]->p() );
   hst[33]->Fill( -Slct.trk_Km[0]->p() );
   for(int i = 0; i < ng; i++) {
      hst[34]->Fill(Slct.Pg[i].e());
   }

   return true;
}

// {{{1 Vertex & Kinematic Fit
//--------------------------------------------------------------------
static bool VertKinFit(Select& Slct) {
//--------------------------------------------------------------------

   // search for a good vertex:
   WTrackParameter wp[2] = {
      WTrackParameter(mk, Slct.trk_Kp[0]->getZHelix(),
                      Slct.trk_Kp[0]->getZError() ),
      WTrackParameter(mk, Slct.trk_Km[0]->getZHelix(),
                      Slct.trk_Km[0]->getZError() )
   };

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

   VertexFit* vtxfit = VertexFit::instance();
   vtxfit->init();
   for(int k = 0; k < 2; k++) {
      vtxfit->AddTrack(k,wp[k]);
   }
   vtxfit->AddVertex(0, vxpar,0, 1);
   bool ok = vtxfit->Fit(0);

   hst[41]->Fill( double(ok) );
   if ( !ok ) {
      return false;
   }
   hst[1]->Fill(3);
   hst[42]->Fill(vtxfit->chisq(0));

   // 4C kinematic fit: select two gammas with the best chisq
   vtxfit->BuildVirtualParticle(0);
   vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0
   for(int k = 0; k < 2; k++) {
      wp[k] = vtxfit->wtrk(k);
   }

   KinematicFit* kmfit = KinematicFit::instance();

   double chisq = 9999.;
   int ng = Slct.gtrk.size();
   for(int i = 0; i < ng-1; i++) {
      RecEmcShower* g1Trk = Slct.gtrk[i];
      for(int j = i+1; j < ng; j++) {
         RecEmcShower* g2Trk = Slct.gtrk[j];

         kmfit->init();
         for(int k = 0; k < 2; k++) {
            kmfit->AddTrack(k, wp[k]);
         }
         kmfit->AddTrack(2, 0.0, g1Trk);
         kmfit->AddTrack(3, 0.0, g2Trk);
         kmfit->AddFourMomentum(0, Slct.LVcms);
         bool oksq = kmfit->Fit();
         if ( oksq ) {
            double chi2 = kmfit->chisq();
            if ( chi2 < chisq && chi2 > 0. ) {
               chisq = chi2;
               Slct.g4f[0] = g1Trk;
               Slct.g4f[1] = g2Trk;
            }
         }
      } //-----------------------------------------------End for(j)
   } //-----------------------------------------------End for(i)

   hst[51]->Fill(chisq);
   if ( !Slct.g4f[0] ) {
      return false;   // can not find good pair of gammas
   }

   hst[1]->Fill(4); // "cuts"
   if ( isMC ) {
      hst[104]->Fill( Slct.decJpsi );
   }

   // ++ check that one more gamma does not give better fit ++
   double chisq_3g = 9999.;
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] || gt == Slct.g4f[1] ) {
         continue;
      }
      kmfit->init();
      for(int k = 0; k < 2; k++) {
         kmfit->AddTrack(k, wp[k]);
      }
      kmfit->AddTrack(2, 0.0, Slct.g4f[0]);
      kmfit->AddTrack(3, 0.0, Slct.g4f[1]);
      kmfit->AddTrack(4, 0.0, gt); // additional gamma
      kmfit->AddFourMomentum(0, Slct.LVcms);
      bool oksq = kmfit->Fit();
      if ( oksq ) {
         double chi2 = kmfit->chisq();
         if ( chi2 < chisq_3g && chi2 > 0. ) {
            chisq_3g = chi2;
         }
      }
   } // end for( additional gamma )

   if ( chisq_3g < 200. ) {
      hst[52]->Fill( chisq - chisq_3g );
      if ( chisq_3g < chisq ) {
         hst[53]->Fill( chisq_3g );
         hst[54]->Fill( chisq );
      }
   }
   // ++ repeat fit for two best photons ++
   kmfit->init();
   for ( int k = 0; k < 2; k++ ) {
      kmfit->AddTrack(k, wp[k]);
   }
   kmfit->AddTrack(2, 0.0, Slct.g4f[0]);
   kmfit->AddTrack(3, 0.0, Slct.g4f[1]);
   kmfit->AddFourMomentum(0, Slct.LVcms);
   bool oksq = kmfit->Fit();
   if ( !oksq || fabs(kmfit->chisq()-chisq) > 0.1 ) {
      cout << " WARNING: Bad second 4C fit ";
      if( oksq ) {
         cout << "chisq= "<<chisq<<" new fit chisq= "<<kmfit->chisq();
      }
      cout << endl;
      Warning("Bad second 4C fit");
      return false;
   }
   chisq=kmfit->chisq();

   //-----------------------------------------------------------------
   // Final histograms and ntuples
   //-----------------------------------------------------------------
   double Eg_max = 0.; // max momentum of gammas not selected by fit
   double mang = 200.; // min angle wrt of charge trk for sel. gammas
   HepLorentzVector Pgg0; // sum of momentums of good photons
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] || gt == Slct.g4f[1] ) {
         // two selected photons
         Pgg0 += Slct.Pg[i];
         double deg = RtoD(Slct.angt[i]);
         hst[61]->Fill(deg);
         if( deg < mang ) {
            mang = deg;
         }
         continue;
      }
      // rejected photons
      double Eg = Slct.Pg[i].e();
      hst[62]->Fill(Eg);
      if( Eg > Eg_max ) {
         Eg_max = Eg;
      }
   }
   if ( Eg_max > 0.01 ) {
      hst[63]->Fill(Eg_max);
   }

   // Momentums after kinematic corrections:
   // K+ and K-
   HepLorentzVector Pkp = kmfit->pfit(0);
   HepLorentzVector Pkm = kmfit->pfit(1);
   HepLorentzVector Pkk = Pkp + Pkm;
   double Mkk  = Pkk.m();

   // two gammas
   HepLorentzVector Pg1 = kmfit->pfit(2);
   HepLorentzVector Pg2 = kmfit->pfit(3);
   HepLorentzVector Pgg = Pg1 + Pg2;
   double Mgg  = Pgg.m();
   hst[71]->Fill(Mgg*Mgg,Mkk*Mkk);

   // J/Psi
   HepLorentzVector Pjpsi = Pkk + Pgg;
   Hep3Vector beta_jpsi = Pjpsi.boostVector();
   HepLorentzVector Pphi = Pkk;
   Pphi.boost(-beta_jpsi); // momentum phi in J/Psi rest system
   if ( chisq < 80 &&
         ( Mkk > 2*mk && Mkk < 1.08 ) &&
         fabs(Mgg-meta) < 0.024
      ) {
      hst[72]->Fill( Pphi.cosTheta() );
   }

   // M^2(Kg)
//    double M2kpg1 = (Pkp + Pg1).m2();
//    double M2kpg2 = (Pkp + Pg2).m2();
//    double M2kmg1 = (Pkm + Pg1).m2();
//    double M2kmg2 = (Pkm + Pg2).m2();
   // M^2(K eta)
   double M2kpeta = (Pkp + Pgg).m2();
   double M2kmeta = (Pkm + Pgg).m2();

   // HepLorentzVector.perp == HepLorentzVector.vect.rho
   // HepLorentzVector.rho == HepLorentzVector.vect.mag
   Double_t xfill[] = {
      chisq, chisq_3g,
      Pkp.rho(), Pkp.cosTheta(), Pkp.phi(),
      Pkm.rho(), Pkm.cosTheta(), Pkm.phi(),
      Pg1.e(),   Pg1.cosTheta(), Pg1.phi(),
      Pg2.e(),   Pg2.cosTheta(), Pg2.phi(),
      Pgg.rho(), Pgg.cosTheta(), Pgg.phi(),
      Mkk, Mgg,
      M2kpeta, M2kmeta,
      double(Slct.decJpsi),
      Slct.Xisr,
      Slct.mc_mkk,
      Slct.mc_mkpeta
   };
//       M2kpg1,M2kpg2, M2kmg1,M2kmg2,
   m_tuple[0]->Fill( xfill );

   if ( isMC ) {
      // fill decay table of Jpsi
      static const double seta = 0.008;
      if ( chisq < 80
         && Mkk > 2*mk && Mkk < 1.08
         && abs(Mgg-meta) < 3*seta ) {
         JpsiTbl.Add();
      }
   }

   return true;
}

// {{{1 MAIN: Event
//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool SelectKKggEvent( ReadDst*       selector,
                      TEvtHeader*    m_TEvtHeader,
                      TDstEvent*     m_TDstEvent,
                      TEvtRecObject* m_TEvtRecObject,
                      TMcEvent*      m_TMcEvent,
                      TTrigEvent*    m_TTrigEvent,
                      TDigiEvent*    m_TDigiEvent,
                      THltEvent*     m_THltEvent        ) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   m_abscor->AbsorptionCorrection(selector);

   Select Slct; // information for the current event

   //-----------------------------------------------------------------
   //-- Get event information --
   //-----------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   int eventNo = m_TEvtHeader->getEventId();
   Slct.runNo  = runNo;
   Slct.event  = eventNo;
   hst[0]->Fill(abs(runNo));
   hst[1]->Fill(0); // "cuts"

   double Ecms =  GetEcms(runNo);
   if ( Ecms < 0 ) { // skip bad runs
      return false;
   }
   hst[9]->Fill(Ecms);
   Slct.LVcms = HepLorentzVector(Ecms*sin(beam_angle), 0, 0, Ecms);

   isMC = (runNo < 0);
   if ( (isMC && m_TMcEvent->getMcParticleCol()->IsEmpty() ) ||
         ( !isMC &&
           m_TMcEvent != 0 &&
           !m_TMcEvent->getMcParticleCol()->IsEmpty() )
      ) {
      cout << " WARNING: something wrong: isMC= " << isMC
         << " McParticles= "
         << m_TMcEvent->getMcParticleCol()->GetEntries() << endl;
      Warning("Incorrect number of MC particles");
      return false;
   }

   if ( isMC ) {
      FillHistoMC(selector,Slct);       // MC histo
      ISRhistoMC( selector,Slct);       // MC ISR histo
   }

   Hep3Vector xorigin = getVertexOrigin(runNo);
   Slct.xorig = xorigin;

   //-----------------------------------------------------------------
   if ( !ChargedTracksKK(selector, Slct) ) {
      return false;
   }
   //-----------------------------------------------------------------
   if ( !NeutralTracks(selector, Slct) ) {
      return false;
   }
   //-----------------------------------------------------------------
   if ( !VertKinFit(Slct) ) {
      return false;
   }

   return true;
}

// {{{1 EndJob
//--------------------------------------------------------------------
BeanUserShared_EXPORT
void SelectKKggEndJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
   }

   if ( isMC ) {
      // print tables of decays
      cout << string(65,'#') << endl;
      cout << "Decays of J/psi -> K+ K- 2gammas" << endl
           << "       final: " << JpsiTbl.ntot << " decays" << endl
           << "       size of table is " << JpsiTbl.Size() << endl;
      JpsiTbl.Print(0.01); // do not print decays wich < 0.01% of all
      cout << "Enddecay" << endl << endl;
   }

   string module = string(__func__);
   module = module.substr(0,module.size()-6);
   if ( warning_msg.empty() ) {
      cout << " There are no warnings in " << module << endl;
   } else {
      cout << " Check output for WARNINGS in " << module << endl;
      for(auto it = warning_msg.begin();
            it != warning_msg.end(); ++it) {
         cout << it->first << " : " << it->second << endl;
      }
   }
}

//-------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
