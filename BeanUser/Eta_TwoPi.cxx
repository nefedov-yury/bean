//======================================================================//
//                                                                      //
// Eta_TwoPi - search for e+ e- -> Eta Pi+ Pi-                          //
//                                                                      //
//======================================================================//
//
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

//-----------------------------------------------------------------------------
// Structure to save variables for a single event
//-----------------------------------------------------------------------------
struct SelectEtaTwoPi {
   // Run-info
   int runNo;              // run-number
   int event;              // event-number
   HepLorentzVector LVcms; // Momentum in center of mass system
   Hep3Vector xorig;       // interaction point from DB

   // MC information:
   int decJpsi;            // decay codes for J/Psi MC events
   double Xisr;            // Xisr = s'/s
   int decEta;             // eta decay:1=2gamma, 2=3pi0, 3=pi+pi-pi0
   vector<HepLorentzVector> mcPg; // 4-mom. of gammas from eta decay
   double M2_pip_pim;       //squared invariant mass of (pi+, pi-)
   double M2_pip_eta;       //squared invariant mass of (pi+, eta)

   // Pion candidate
   vector<RecMdcKalTrack*> trk_Pip;  // positive tracks
   vector<RecMdcKalTrack*> trk_Pim;  // negative tracks
   vector<double> pid_Pip; //pid of positive tracks
   vector<double> pid_Pim; //pid of negative tracks
   vector<double> E_p_Pip; //"E over p" of positive tracks
   vector<double> E_p_Pim; //"E over p" of negative tracks

   // gamma tracks
   vector<RecEmcShower*> gtrk;
   vector<double> angt;             // angle with closest charged track
   vector<HepLorentzVector> Pg;     // 4-momentum of gammas
   vector<RecEmcShower*> g4f;       // six best photons after 4C-fit

   SelectEtaTwoPi() {
      decJpsi = -1;
      Xisr = 0.;
      decEta = -1;

      g4f.resize(6,nullptr);
      mcPg.reserve(6);
   }
};
//-----------------------------------------------------------------------------
// Global variables
//-----------------------------------------------------------------------------
const static bool init_selection=false;

const static double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
const static double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
const static double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
const static double meta   = 0.547862; // 547.862   +/- 0.017   MeV
const static double momega = 0.78265;  // 782.65    +/- 0.12    MeV
const static double mk     = 0.493677; // 493.677   +/- 0.016   MeV
const static double mphi   = 1.019461; //1019.461   +/- 0.019   MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1*> hst;
static std::vector<TH2*> hst2;
static std::vector<TNtupleD*> m_tuple;

// container for warnings
static map<string,int> warning_msg;

static bool isMC = true;

//Dalitz plots limits
static map<double,tuple<double,double,int,int>> Dalitz_Parameters {
  {2.000,  make_tuple(2.5,4.0,50,50)},
  {2.050,  make_tuple(2.5,4.0,50,50)},
  {2.100,  make_tuple(2.5,4.0,50,50)},
  {2.12655,  make_tuple(3.0,4.0,50,50)},
  {2.150,  make_tuple(3.0,4.5,50,50)},
  {2.175,  make_tuple(3.0,4.5,50,50)},
  {2.200,  make_tuple(3.0,4.5,50,50)},
  {2.2324, make_tuple(3.0,5.0,50,50)},
  {2.3094, make_tuple(3.5,5.5,50,50)},
  {2.3864, make_tuple(3.5,5.5,50,50)},
  {2.396,  make_tuple(3.5,5.5,50,50)},
  {2.500,  make_tuple(4.5,6.5,50,50)},
  {2.6444, make_tuple(4.5,6.5,50,50)},
  {2.6464, make_tuple(4.5,6.5,50,50)},
  {2.700,  make_tuple(5.0,7.0,50,50)},
  {2.8,    make_tuple(5.5,8.0,50,50)},
  {2.900,  make_tuple(6.0,8.0,50,50)},
  {2.950,  make_tuple(6.0,8.0,50,50)},
  {2.981,  make_tuple(6.5,8.5,50,50)},
  {3.000,  make_tuple(6.5,8.5,50,50)},
  {3.020,  make_tuple(6.5,9.0,50,50)},
  {3.080,  make_tuple(7.0,9.0,50,50)},
};

// Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
static inline void Warning(const char* msg) {
   warning_msg[string(msg)] += 1;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
static inline Double_t RtoD(Double_t ang) {
   return ang*180/M_PI;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
static inline double SQ(double x) {
   return x*x;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void Eta_TwoPiStartJob(ReadDst* selector) {
//-----------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }
   if ( init_selection ) {
      cout << " ===== INIT SELECTION =====" << endl;
   }

   hst.resize(500,nullptr);
   hst2.resize(10, nullptr);
   m_tuple.resize(10,nullptr);

   // init Absorption Correction
   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

   //initialize EventTag
   m_EventTagSvc = EventTagSvc::instance();

   // set paths to pdg & decayCodes files:
   m_EventTagSvc->setPdtFile(
      selector->AbsPath("Analysis/EventTag/share/pdt_bean.table")
                            );
   m_EventTagSvc->setDecayTabsFile(
      selector->AbsPath(
         "Analysis/EventTag/share/DecayCodes/dcode_charmonium.txt"
      )                           );
   m_EventTagSvc->setIgnorePhotons(false);

   if ( selector->Verbose() ) {
      m_EventTagSvc->setVerbose(1);
   }

   m_EventTagSvc->initialize();

   // We have to initialize DatabaseSvc ---------------------------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   // We have to initialize Magnetic field ------------------------------------
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if ( !(mf->GetPath()).empty() ) {
      cout << " WARNING:"
           << " MagneticFieldSvc has already been initialized" << endl
           << "                         path = " << mf->GetPath() << endl;
      Warning("MagneticFieldSvc has already been initialized");
   }
   // set path to directory with magnetic fields tables
   mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
   mf->UseDBFlag(false); // like in the boss program
   mf->RunMode(3); // like in the boss program

   // set path for ParticleID algorithm
   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   //--------- Book histograms ------------------------------------------------

   hst[0] = new TH1D("Runs", "run numbers", 1000,0.,50000);
   hst[1] = new TH1D("cuts_0","preselections cuts", 20,-0.5,19.5);
   hst[9] = new TH1D("Ecms","E_{cms} GeV", 1100, 2., 3.1);

   //  ChargedTracksPiPi:
   hst[11] = new TH1D("Rxy","R_{xy}", 200,-2.,2.);
   hst[12] = new TH1D("Rz","R_{z}", 200,-20.,20.);
   hst[13] = new TH1D("cos_theta","cos(#theta)", 200,-1.,1.);
   hst[14] = new TH1D("theta","#theta", 180,0.,180.);
   hst[15] = new TH1D("PQ","Charged momentum (Q*P)", 600,-3.,3.);

   hst[16] = new TH1D("EoP","E/P for all charged tracks", 150,0.,1.5);
   hst[17] = new TH1D("Pid_clPi","lg(CL_{#pi})", 100,-4.,0.);

   hst[23] = new TH1D("Pi_Q","Charged momentum (Q*P)for Pi", 300,-1.5,1.5);
   hst[24] = new TH2D("Pi_pm","N(Pi+) vs N(Pi-)", 10,-0.5,9.5, 10,-0.5,9.5);
   hst[25] = new TH1D("Pi_N","number of good pions", 10,-0.5,9.5);

   hst[26] = new TH1D("Mrec2","recoil mass square of pi+ pi-", 100,0.,1.);

   // NeutralTracks:
   hst[31] = new TH1D("G_ang","angle with closest chg.trk", 180,0.,180.);
   hst[32] = new TH1D("G_n","N_{#gamma} in event", 15,-0.5,14.5);
   hst[33] = new TH1D("Pi_PM","momentum Pi+ and Pi-", 240,-1.2,1.2);
   hst[34] = new TH1D("G_P","Momentum of gammas", 200,0.,2.);

   hst[35] = new TH1D("Ng_mc2g","N_{#gamma} for eta->2g", 11,-0.5,10.5);
   hst[36] = new TH1D("Ng_mc6g","N_{#gamma} for eta->3pi0", 15,-0.5,14.5);

   // VertKinFit_2:
   hst[41] = new TH1D("vtx_fit_2", "vertex fit 2g: 0/1 - bad/good", 2,-0.5,1.5);
   hst[42] = new TH1D("vtx_chi2_2", "vertex fit 2g #chi^2", 100,0.,100.);
   hst[51] = new TH1D("f4c_chi2_2","4C-fit 2g: min #chi^{2}", 300,0.,300.);

   hst[61] = new TH1D("f4c_ang_2","4C-fit 2g: min angle gamma trk", 180,0.,180.);
   hst[62] = new TH1D("fit_eg_2","fit 2g: E(#gamma) rejected", 500,0.,0.5);
   hst[63] = new TH1D("f4c_egm_2","4C-fit 2g: E_{max}(#gamma) dropped", 100,0.,2.);

   hst[71] = new TH2D("Mgg_Mpipi_2", "M^{2}(Pi^{+}Pi^{-}) vs M^{2}(2 #gamma)",
                       220,0.,1.1, 200,0.,10.);
   hst[72] = new TH1D("M2g", "M^{2}(2 #gamma)", 200,0.,1.); // expected eta

   // VertKinFit_6:
   hst[241] = new TH1D("vtx_fit_6", "vertex fit 6g: 0/1 - bad/good", 2,-0.5,1.5);
   hst[242] = new TH1D("vtx_chi2_6", "vertex fit 6g #chi^2", 100,0.,100.);
   hst[251] = new TH1D("f4c_chi2_6","4C-fit 6g: min #chi^{2}", 300,0.,300.);

   hst[261] = new TH1D("f4c_ang_6","4C-fit 6g: min angle gamma trk", 180,0.,180.);
   hst[262] = new TH1D("fit_eg_6","fit 6g: E(#gamma) rejected", 500,0.,0.5);
   hst[263] = new TH1D("f4c_egm_6","4C-fit 6g: E_{max}(#gamma) dropped", 100,0.,2.);

   hst[271] = new TH2D("Mgg_Mpipi_6", "M^{2}(Pi^{+}Pi^{-}) vs M^{2}(6 #gamma)",
                       220,0.,1.1, 200,0.,10.);
   hst[272] = new TH1D("M6g", "M^{2}(6 #gamma)", 200,0.,1.); // expected eta


   // Monte Carlo histograms:
   hst[100] = new TH1D("mc_deccode0", "decCode nocut",258,-1.5,256.5);
   hst[101] = new TH1D("mc_deccode1", "decCode cut#1",258,-1.5,256.5);
   hst[102] = new TH1D("mc_deccode2", "decCode cut#2",258,-1.5,256.5);
   hst[103] = new TH1D("mc_deccode3", "decCode cut#3",258,-1.5,256.5);
   hst[105] = new TH1D("mc_deccodeF", "decCode final",258,-1.5,256.5);

   hst[111] = new TH1D("mc_pdg", "PDG for all", 2001,-1000.5,1000.5);
   hst[112] = new TH1D("mc_pdg0","PDG from primary vtx.",
                        2001,-1000.5,1000.5);
   hst[113] = new TH1D("mc_Eg", "ISR: E(#gamma) (MeV)", 100,0.,1000.);
   hst[114] = new TH1D("mc_Ng", "ISR: N_{#gamma}", 10,-0.5,9.5);
   hst[115] = new TH1D("mc_Etot","ISR: Eisr (MeV)",100,0.,2000.);
   hst[116] = new TH1D("mc_xisr_2g","ISR: s'/s for 2g",1000,0.,1.);
   hst[117] = new TH1D("mc_xisr_6g","ISR: s'/s for 6g", 1000, 0., 1.);
   hst[118] = new TH1D("mc_bspr","MC: energy spread (MeV)",100,-5.,5.);
   hst[119] = new TH1D("mc_dec_Eta", "Decay codes for eta (1 - 2g, 2 - 6g)",
                       5, -1.5, 3.5);

   hst[121] = new TH1D("mc_PipP", "Momentum of #pi^{+}", 100,0.,2.);
   hst[122] = new TH1D("mc_PipPt","Pt of #pi^{+}", 100,0.,2.);
   hst[123] = new TH1D("mc_PipC", "cos(#Theta) of #pi^{+}", 100,-1.,1.);
   hst[124] = new TH1D("mc_PimP", "Momentum of #pi^{-}", 100,0.,2.);
   hst[125] = new TH1D("mc_PimPt","Pt of #pi^{-}", 100,0.,2.);
   hst[126] = new TH1D("mc_PimC", "cos(#Theta) of #pi^{-}", 100,-1.,1.);
   hst[131] = new TH1D("mc_EtaP", "Momentum of #eta", 100,0.,2.);
   hst[132] = new TH1D("mc_EtaPt","Pt of #eta", 100,0.,2.);
   hst[133] = new TH1D("mc_EtaC", "cos(#Theta) of #eta", 100,-1.,1.);
   hst[134] = new TH1D("mc_2gE", "E(g) from eta->2g", 100,0.,2.);
   hst[135] = new TH1D("mc_2gC", "cos(g) from eta->2g", 100,-1.,1.);
   hst[136] = new TH1D("mc_3pi0P", "P(pi0) from eta->3pi0", 100,0.,2.);
   hst[137] = new TH1D("mc_3pi0C", "cos(pi0) from eta->3pi0", 100,-1.,1.);
   hst[138] = new TH1D("mc_6gE", "E(g) from eta->3pi0->6g", 100,0.,2.);
   hst[139] = new TH1D("mc_6gC", "cos(g) from eta->3pi0->6g", 100,-1.,1.);
   hst[141] = new TH1D("mc_3pi1P", "P(pi+) eta->pi+pi-pi0", 100,0.,2.);
   hst[142] = new TH1D("mc_3pi1C", "cos(pi+) eta->pi+pi-pi0", 100,-1.,1.);
   hst[143] = new TH1D("mc_3pi2P", "P(pi-) eta->pi+pi-pi0", 100,0.,2.);
   hst[144] = new TH1D("mc_3pi2C", "cos(pi-) eta->pi+pi-pi0", 100,-1.,1.);
   hst[145] = new TH1D("mc_3pi3P", "P(pi0) eta->pi+pi-pi0", 100,0.,2.);
   hst[146] = new TH1D("mc_3pi3C", "cos(pi0) eta->pi+pi-pi0", 100,-1.,1.);
   hst[147] = new TH1D("mc_3pi2gE", "E(g) eta->pi+pi-pi0", 100,0.,2.);
   hst[148] = new TH1D("mc_3pi2gC", "cos(g) eta->pi+pi-pi0", 100,-1.,1.);

   hst[151] = new TH1D("mc_cos_eta_2g", "cos(#Theta) of eta for eta->2g", 100, -1., 1.);
   hst[152] = new TH1D("mc_cos_pip_2g", "cos(#Theta) of pi+ for eta->2g", 100, -1., 1.);
   hst[153] = new TH1D("mc_cos_pim_2g", "cos(#Theta) of pi- for eta->2g", 100, -1., 1.);
   hst[154] = new TH1D("mc_cos_eta_6g", "cos(#Theta) of eta for eta->3pi0->6g", 100, -1., 1.);
   hst[155] = new TH1D("mc_cos_pip_6g", "cos(#Theta) of pi+ for eta->3pi0->6g", 100, -1., 1.);
   hst[156] = new TH1D("mc_cos_pim_6g", "cos(#Theta) of pi- for eta->3pi0->6g", 100, -1., 1.);

   hst2[2] = new TH2D("Dalitz_ini_2g", "Dalitz Plot for initial momenta of particles, eta->2g", 50,0.,10.,80,0.,10.);
   hst2[6] = new TH2D("Dalitz_ini_6g", "Dalitz Plot for initial momenta of particles, eta->3pi0->6g", 40,0.,10.,64,0.,10.);
   hst2[3] = new TH2D("Dalitz_ini_2g_0", "ini Dalitz Plot eta->2g for x>0.9", 50,0.,10.,80,0.,10.);

   // final ntuple for e+ e- -> Pi+ Pi- 2 gammas
   m_tuple[2] = new TNtupleD("a4c2","after 4C kinematic fit (2 gammas)",
                "ch2:"            // chi^2 of 4C fit
                "Ppip:Cpip:phipip:pidpip:E_ppip:"  // P, cos(Theta), phi, PID and "E over p" of Pi+
                "Ppim:Cpim:phipim:pidpim:E_ppim:"  // P, cos(Theta), phi, PID and "E over p" of Pi-
                "Eg1:Cg1:phig1:"   // E, cos(Theta) and phi of gamma-1
                "Eg2:Cg2:phig2:"   // E, cos(Theta) and phi of gamma-2
                "Pgg:Cgg:phigg:"   // P, cos(Theta) and phi of eta (2gamma)
                "M2pi:Mgg:"        // invariant masses of Pi+Pi- and 2gammas
                "M2pp_true:M2pe_true:" //squared invariant masses of (pi+,pi-) and (pi+,eta) with momenta taken before kinematic fit for 2 gamma
                "dec:"            // MC: decay codes of J/Psi or eta
                "xisr"            // MC: Xisr = s'/s
                            );

   // final ntuple for e+ e- -> Pi+ Pi- 6 gammas
   m_tuple[6] = new TNtupleD("a4c6","after 4C kinematic fit (6 gammas)",
                "ch2:"            // chi^2 of 4C fit
                "Ppip:Cpip:phipip:pidpip:E_ppip:"  // P, cos(Theta), phi, PID and "E over p" of Pi+
                "Ppim:Cpim:phipim:pidpim:E_ppim:"  // P, cos(Theta), phi, PID and "E over p" of Pi-
                "Eg1:Cg1:phig1:"   // E, cos(Theta) and phi of gamma-1
                "Eg2:Cg2:phig2:"   // E, cos(Theta) and phi of gamma-2
                "Eg3:Cg3:phig3:"   // E, cos(Theta) and phi of gamma-3
                "Eg4:Cg4:phig4:"   // E, cos(Theta) and phi of gamma-4
                "Eg5:Cg5:phig5:"   // E, cos(Theta) and phi of gamma-5
                "Eg6:Cg6:phig6:"   // E, cos(Theta) and phi of gamma-6
                "Pgg:Cgg:phigg:"   // P, cos(Theta) and phi of eta (6gamma)
                "M2pi:Mgg:"        // invariant masses of Pi+Pi- and 6gammas
                "M2pp_true:M2pe_true:" //squared invariant masses of (pi+,pi-) and (pi+,eta) with momenta taken before kinematic fit for 6 gamma
                "dec:"            // MC: decay codes of J/Psi or eta
                "xisr"            // MC: Xisr = s'/s
                            );

   // register in selector to save in given directory
   const char* SaveDir = "Eta_TwoPi";
   VecObj his1o(hst.begin(),hst.end());
   selector->RegInDir(his1o,SaveDir);
   VecObj his2o(hst2.begin(),hst2.end());
   selector->RegInDir(his2o,SaveDir);
   VecObj ntuples(m_tuple.begin(),m_tuple.end());
   selector->RegInDir(ntuples,SaveDir);
}

//-----------------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//-----------------------------------------------------------------------------
   // cache value for one run
   static int save_runNo = 0;
   static Hep3Vector xorigin;
   if ( runNo == save_runNo ) {
      return xorigin;
   }

   // update vertex for new run
   xorigin.set(0.,0.,0.);
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();

   // for R-scan (2015)
   int run = abs(runNo);
   if (run >=  9947 && run <= 10878 ) { // J/Psi 2009 no 6.6.4 !
      vtxsvc->SetBossVer("6.6.3");
   } else if( run >=39355 && run <= 43253 ) { // R-scan(2015)+Y(2175)
     #if (BOSS_VER < 700)
          vtxsvc->SetBossVer("6.6.5.p01");
     #else
          vtxsvc->SetBossVer("7.0.3");
     #endif
   } else {
      cout << " WARNING:" << __func__ << " :" << " run=" << runNo
         << " We will use default Boss Version: " << vtxsvc->GetBossVer()
         << endl;
      Warning( (string("getVertexOrigin: default Boss Version ")
            + vtxsvc->GetBossVer()).c_str() );
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
           " Cannot obtain vertex information for run#" << runNo << endl;
      exit(1);
   }

   save_runNo = runNo;
   return xorigin;
}

//-----------------------------------------------------------------------------
// GetEcms: return negative value for bad runs
//-----------------------------------------------------------------------------
static double GetEcms(int runNo, bool verbose = false) {
//-----------------------------------------------------------------------------
   struct runInfo {
      int runS, runE; // first and last runs in period
      double Ebeam;   // energy of beam (MeV)
      double Spread;  // beam spread (KeV)
      double Lumi;    // luminosity (pb^-1)
   };

   // List Runs for R-Scan Phase 1 (2014/2015)
   // https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples
   static const runInfo ListRscan[] =  {
       { 39355, 39618,      3080.,      0,      126.185 },
       { 39711, 39738,      3020.,      0,      17.290 },
       { 39680, 39710,      3000.,      0,      15.881 },
       { 39651, 39679,      2981.,      0,      16.071 },
       { 39619, 39650,      2950.,      0,      15.942 },
       { 39775, 40069,      2900.,      0,      105.253 },
       { 40128, 40296,      2644.4,     0,      34.003 },
       { 40300, 40435,      2646.4,     0,      33.722 },
       { 40436, 40439,      2700.,      0,      1.034 },
       { 40440, 40443,      2800.,      0,      1.008 },
       { 40459, 40769,      2396.,      0,      66.869 },
       { 40771, 40776,      2500.,      0,      1.098 },
       { 40806, 40951,      2386.4,     0,      22.549 },
       { 40989, 41121,      2200.,      0,      13.699 }, //site -- 2646.4
       { 41122, 41239,      2232.4,     0,      11.856 },
       { 41240, 41411,      2309.4,     0,      21.089 },
       { 41416, 41532,      2175.,      0,      10.625 },
       { 41533, 41570,      2150.,      0,      2.841 },
       { 41588, 41727,      2100.,      0,      12.167 },
       { 41729, 41909,      2000.,      0,      10.074 },
       { 41911, 41958,      2050.,      0,      3.343 },
       { 42004, 43253,      2126.55,    0,      108.49},
       // Last 2 runs are with separated beams =>
       // should not be used in the analysis
//        { 40777, 40804,      2644.4,     0,      0. },
//        { 41959, 41999,      2232.4,     0,      0. },
   };
   static const int Np = sizeof(ListRscan)/sizeof(ListRscan[0]);

   // cache value for one run
   static int save_runNo = 0;
   static double Ecms = -3.; // bad run if not found
   if ( runNo == save_runNo ) {
      return Ecms;
   }

   int absrun = abs(runNo);
   bool found = false;

   // search in R-scan
   for(int i = 0; i < Np; i++) {
      if ( absrun >= ListRscan[i].runS && absrun <= ListRscan[i].runE ) {
         Ecms   = ListRscan[i].Ebeam  * 1.e-3;  // MeV -> GeV
//             Lumi   = ListRscan[i].Lumi;
         found  = true;
         if ( runNo > 0 ) { // data
            RscanDQ rdq = RscanDQ(runNo);
            double Ebeam = rdq.getEbeam(); // beam energy in GeV
            int status = rdq.getStatus();

            if ( fabs(2*Ebeam-Ecms) > 1e-4 ) {
               cout << " GetEcms::WARNING RscanDQ problem for run= "
                  << runNo << " Ebeam= " << Ebeam
                  << " Ecms= " << Ecms << endl;
               Warning("RscanDQ problem");
            }

            if ( status != 1 ) {
               cout << " GetEcms::RscanDQ bad run=" << runNo << endl;
               Ecms = -Ecms;
            }
         }
         break;
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
      cout << " GetEcms: run= " << runNo << " -> " << " Ecms= " << Ecms
//            << " Spread= " << Spread
//            << " Lumi= " << Lumi
           << endl;
   }

   save_runNo = runNo;
   return Ecms;
}

//-----------------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, SelectEtaTwoPi& Slct) {
//-----------------------------------------------------------------------------
   if ( !isMC ) {
      return;
   }

   const TMcEvent*  m_TMcEvent  = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   // EventTag: only for J/Psi MC:
   int absrun = abs(Slct.runNo);
   if ( absrun >= 9947 && absrun <= 10878 ) { // J/Psi 2009
      m_EventTagSvc->setMcEvent(const_cast<TMcEvent*>(m_TMcEvent));
      unsigned int evTag = m_EventTagSvc->getEventTag();
      // check general event type
      if ( (evTag & 0xF) != 4 ) { // 4 is J/Psi event
         static int nprt = 0;
         if ( nprt < 1 ) {
            printf(" WARNING: MC is not J/Psi evTag= 0x%08x\n",evTag);
            ++nprt;
         }
         Warning("MC is not J/Psi");
      } else {
         // decay code of J/Psi:
         Slct.decJpsi = (evTag >> 8) & 0xFF;
      }

      hst[100]->Fill( Slct.decJpsi );
      return;
   }

   // study ISR
   int    Ng_isr    = 0;       // number of gammas
   double Etot_isr  = 0.;      // total energy
   HepLorentzVector Pisr;      // sum of ISR photons
   HepLorentzVector Pr;        // sum of primary particles
   HepLorentzVector Ppip;      // momentum of Pi+ for Dalitz plot
   HepLorentzVector Ppim;      // momentum of Pi- for Dalitz plot
   HepLorentzVector Peta;      // momentum of Eta for Dalitz plot
   double M2_pp, M2_pe;        // squared invariant masses for Dalitz plots

   // e+ e- -> g(ISR) + virt_ptcl -> g(ISR) + pi+ pi- eta
   static const int Vpdg = 70022; // 'pdg' of virtual particle
   int idx_virt = -99; // primary vertex

   // search for decay eta->2gamma or eta->3pi0 or eta->pi+pi-pi0
   int idx_eta = -1;
   vector<int> Pdg_eta;

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );
      double Ep = part->getInitialMomentumE();

      hst[111]->Fill(part_pdg);
      if ( part->getMother() == idx_virt ) {
         hst[112]->Fill(part_pdg);
         Pr += HepLorentzVector( Vp, Ep );
         if ( abs(part_pdg) == 22 ) { // ISR gamma
            Ng_isr++;
            Etot_isr += Ep;
            Pisr += HepLorentzVector( Vp, Ep );
            hst[113]->Fill(Ep * 1e3); // => MeV
         } else if ( part_pdg ==  Vpdg ) { // virtual particle
            idx_virt = part->getTrackIndex();
         } else if ( part_pdg ==  211 ) { // pi+
//             idx_pip = part->getTrackIndex();
            hst[121]->Fill(Vp.mag());
            hst[122]->Fill(Vp.rho());
            hst[123]->Fill(Vp.cosTheta());
            Ppip.setE(Ep);
            Ppip.setVect(Vp);
         } else if ( part_pdg == -211 ) { // pi-
//             idx_pim = part->getTrackIndex();
            hst[124]->Fill(Vp.mag());
            hst[125]->Fill(Vp.rho());
            hst[126]->Fill(Vp.cosTheta());
            Ppim.setE(Ep);
            Ppim.setVect(Vp);
         } else if ( part_pdg ==  221 ) { // eta
            idx_eta = part->getTrackIndex();
            hst[131]->Fill(Vp.mag());
            hst[132]->Fill(Vp.rho());
            hst[133]->Fill(Vp.cosTheta());
            Peta.setE(Ep);
            Peta.setVect(Vp);
         }
      } // end of primary vertex

      if ( part->getMother() == idx_eta ) {
         Pdg_eta.push_back(int(part_pdg));  // eta -> ...
      }
   } // end of while

   // ISR histograms
   hst[114]->Fill(Ng_isr);
   hst[115]->Fill(Etot_isr * 1e3); // => MeV
   double S = Pr*Pr;
   HepLorentzVector Prp = Pr - Pisr;
   double Sp = Prp*Prp;
   Slct.Xisr = Sp/S;
   hst[118]->Fill( (sqrt(S) - Slct.LVcms.e()) * 1e3 ); // => MeV

   // paranoid check
   if ( idx_eta == -1 ) {
      cout << "WARNING: idx_eta= " << idx_eta << endl;
      selector->PrintMcDecayTree(-99,0); // What is it?
      Warning("MC: bad eta index");
      return;
   }

   // set Slct.decEta
   Slct.decEta = 0;
   if ( Pdg_eta.size() == 2 && Pdg_eta[0]==22 && Pdg_eta[1]==22 ) {
      Slct.decEta = 1; // eta -> 2 gamma
      hst[116]->Fill(Slct.Xisr);
      M2_pp = (Ppip+Ppim).m2();
      M2_pe = (Ppip+Peta).m2();
      Slct.M2_pip_pim = M2_pp;
      Slct.M2_pip_eta = M2_pe;
      hst2[2]->Fill(M2_pp, M2_pe);
      if ( Slct.Xisr > 0.9 ) {
         hst2[3]->Fill(M2_pp, M2_pe);
      }
      hst[151]->Fill(Peta.cosTheta());
      hst[152]->Fill(Ppip.cosTheta());
      hst[153]->Fill(Ppim.cosTheta());
   } else if ( Pdg_eta.size() == 3 ) {
      if ( Pdg_eta[0]==111 && Pdg_eta[1]==111 && Pdg_eta[2]==111 ) {
         Slct.decEta = 2; // eta -> 3pi0
         hst[117]->Fill(Slct.Xisr);
         M2_pp = (Ppip+Ppim).m2();
         M2_pe = (Ppip+Peta).m2();
         Slct.M2_pip_pim = M2_pp;
         Slct.M2_pip_eta = M2_pe;
         hst2[6]->Fill(M2_pp, M2_pe);
         hst[154]->Fill(Peta.cosTheta());
         hst[155]->Fill(Ppip.cosTheta());
         hst[156]->Fill(Ppim.cosTheta());
      } else {
         sort(Pdg_eta.begin(),Pdg_eta.end());
         if ( Pdg_eta[0]==-211 && Pdg_eta[1] ==111 && Pdg_eta[2]==211 ) {
            Slct.decEta = 3; // eta -> pi- pi0 pi+
         }
      }
   }
   //Fill histo with decay codes of eta
   hst[119]->Fill(Slct.decEta);

   // fill MC-histograms for eta decay
   vector<int> idx_pi0; // indexes of pi0 from eta->pi0+...
   TIter mcIter2(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter2.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );
      double Ep = part->getInitialMomentumE();

      if ( part->getMother() == idx_eta ) {
         if ( Slct.decEta == 1 ) { // eta -> 2gamma
            if ( part_pdg == 22 ) {
               Slct.mcPg.push_back( HepLorentzVector(Vp,Ep) );
               hst[134]->Fill(Ep);
               hst[135]->Fill(Vp.cosTheta());
            }
         } else if ( Slct.decEta == 2 ) { // eta -> 3pi0
            if ( part_pdg == 111 ) {
               hst[136]->Fill(Vp.mag());
               hst[137]->Fill(Vp.cosTheta());
               // fill gammas from pi0 decays
               idx_pi0.push_back( int(part->getTrackIndex()) );
            }
         } else if ( Slct.decEta == 3 ) { // eta -> pi+ pi- pi0
            if ( part_pdg == 211 ) {         // pi+
               hst[141]->Fill(Vp.mag());
               hst[142]->Fill(Vp.cosTheta());
            } else if ( part_pdg == -211 ) { // pi-
               hst[143]->Fill(Vp.mag());
               hst[144]->Fill(Vp.cosTheta());
            } else if ( part_pdg == 111 ) {  // pi0
               hst[145]->Fill(Vp.mag());
               hst[146]->Fill(Vp.cosTheta());
               idx_pi0.push_back( int(part->getTrackIndex()) );
            }
         }
         continue;
      }

      // search for gammas from pi0 decays
      if ( part_pdg == 22 ) {
         for (auto idx : idx_pi0) {
            if ( part->getMother() == idx ) {
               if ( Slct.decEta == 2 ) { // eta -> 3pi0
                  hst[138]->Fill(Ep);
                  hst[139]->Fill(Vp.cosTheta());
                  Slct.mcPg.push_back( HepLorentzVector(Vp,Ep) );
               } else if ( Slct.decEta == 3 ) { // eta -> pi+ pi- pi0
                  hst[147]->Fill(Ep);
                  hst[148]->Fill(Vp.cosTheta());
                  Slct.mcPg.push_back( HepLorentzVector(Vp,Ep) );
               }
            }
         }
      } // end if gamma
   } // end of while-2

   return;
}

//-----------------------------------------------------------------------------
// select Pi+ Pi- candidates, no other tracks.
//-----------------------------------------------------------------------------
static bool ChargedTracksTwoPi(ReadDst* selector, SelectEtaTwoPi& Slct) {
//-----------------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
//    static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }

      RecMdcTrack *mdcTrk = itTrk->mdcTrack();

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.); // initial point for MDC reconstruction
      HepPoint3D IP(Slct.xorig[0],Slct.xorig[1],Slct.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      double Rvxy0 = vecipa[0]; // the nearest distance to IP in xy plane
      double Rvz0  = vecipa[3]; // ... in z direction

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

      // E/p check (no cut here, E_p could be used to add cut while making histograms)
      double E_p = -1.;
      if ( itTrk->isEmcShowerValid() ) {
         RecEmcShower *emcTrk = itTrk->emcShower();
         double EnergyEMC = emcTrk->energy();
         double pMDC =  mdcTrk->p();
         E_p =  EnergyEMC / pMDC;
      }
      hst[16]->Fill( E_p );

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
      double check_pion = 0.0;
      if ( pid->IsPidInfoValid() ) {
         if( pid->probPion() > 0 ) {
            hst[17]->Fill( log10(pid->probPion()) );
         } else {
            hst[17]->Fill(-5.);
         }
         // check that track is pion:
         if ( (pid->probPion() > 0.001)             &&
               (pid->probPion() > pid->probKaon())   &&
               (pid->probPion() > pid->probProton())
            ) {
            check_pion = pid->probPion();
         }
      } // end IsPidInfoValid

      //Check_pion cut (could be used during makinng histograms or right here)
      /*if( check_pion == 0 ) {
         continue;
      }*/

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      if ( !mdcKalTrk ) {
         continue;
      }

      mdcKalTrk->setPidType(RecMdcKalTrack::pion);
      hst[23]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      if( mdcKalTrk->charge() > 0 ) {
         Slct.trk_Pip.push_back(mdcKalTrk);
         Slct.pid_Pip.push_back(check_pion);
         Slct.E_p_Pip.push_back(E_p);
      } else {
         Slct.trk_Pim.push_back(mdcKalTrk);
         Slct.pid_Pim.push_back(check_pion);
         Slct.E_p_Pim.push_back(E_p);
      }
   }  //-------------------------------------------------End for (i)

   int np = Slct.trk_Pip.size();
   int nm = Slct.trk_Pim.size();
   hst[24]->Fill(np,nm);
   hst[25]->Fill(np+nm);

   // require exactly one "+" and one "-"
   if ( np != 1 || nm != 1 ) {
      return false;
   }

   hst[1]->Fill(1);
   if ( isMC ) {
      if ( Slct.decJpsi > -1 ) {
         hst[101]->Fill( Slct.decJpsi );
      } else {
         hst[101]->Fill( Slct.decEta );
      }
   }

   // recoil mass of pi+ pi-
   auto trp = Slct.trk_Pip[0];
   Hep3Vector Vp(trp->px(), trp->py(), trp->pz());;
   HepLorentzVector LVp( Vp, sqrt( Vp.mag2() + SQ(mpi) ) );
   auto trm = Slct.trk_Pim[0];
   Hep3Vector Vm(trm->px(), trm->py(), trm->pz());;
   HepLorentzVector LVm( Vm, sqrt( Vm.mag2() + SQ(mpi) ) );
   double Mrec2 = (Slct.LVcms - LVp - LVm).m2();
   hst[26]->Fill( Mrec2 );

   return true;
}

// select gammas candidates
//-----------------------------------------------------------------------------
static int NeutralTracks(ReadDst* selector, SelectEtaTwoPi& Slct) {
//-----------------------------------------------------------------------------
   // parameters of reconstruction
   static const double min_angle = 10 * M_PI/180; // 10 grad

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
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

      // *) good EMC energy deposited in the barrel (endcap) part of EMC
      double eraw = emcTrk->energy();
      double absCosTheta = fabs(  cos(emcTrk->theta()) );

      bool GoodCluster=false;
      if ( absCosTheta < 0.8 ) {  //barrel
         GoodCluster = eraw > 25E-3;
      } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) { //endcap
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

   if ( isMC ) {
      if ( Slct.decEta == 1 ) { // eta->2gamma
         hst[35]->Fill(ng);
      } else if ( Slct.decEta == 2 ) { // eta->3pi0
         hst[36]->Fill(ng);
      }
   }

   hst[1]->Fill(2); // "cuts"
   if ( isMC ) {
      if ( Slct.decJpsi > -1 ) {
         hst[102]->Fill( Slct.decJpsi );
      } else {
         hst[102]->Fill( Slct.decEta );
      }
   }

   // momentum of selected pair of pions and Gammas
   hst[33]->Fill( Slct.trk_Pip[0]->p() );
   hst[33]->Fill( -Slct.trk_Pim[0]->p() );
   for(int i = 0; i < ng; i++) {
      hst[34]->Fill(Slct.Pg[i].e());
   }

   return ng;
}

//-----------------------------------------------------------------------------
// Vertex & Kinematic Fit
//-----------------------------------------------------------------------------
static bool VertKinFit_2(SelectEtaTwoPi& Slct) {
//-----------------------------------------------------------------------------
   //Clearing Slct.g4f
   Slct.g4f.clear();
   Slct.g4f.resize(6, nullptr);

   // search for a good vertex:
   WTrackParameter wp[2] = {
      WTrackParameter(mpi, Slct.trk_Pip[0]->getZHelix(),
                      Slct.trk_Pip[0]->getZError() ),
      WTrackParameter(mpi, Slct.trk_Pim[0]->getZHelix(),
                      Slct.trk_Pim[0]->getZError() )
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

   // 4C kinematic fit
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
      return false;   // can not find good set of gammas
   }
   hst[1]->Fill(4);
   if ( isMC ) {
      if ( Slct.decJpsi > -1 ) {
         hst[103]->Fill( Slct.decJpsi );
      } else {
         hst[103]->Fill( Slct.decEta );
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
      cout << " WARNING: Bad second 4C fit for 2gammas";
      if( oksq ) {
         cout << "2g: chisq= "<<chisq<<" new fit chisq= "<<kmfit->chisq();
      }
      cout << endl;
      Warning("Bad second 4C fit for 2gammas");
      return false;
   }
   chisq=kmfit->chisq();

   //--------------------------------------------------------------------------
   // Final histograms and ntuples
   //--------------------------------------------------------------------------

   double Eg_max = 0.; // max momentum of gammas not selected by fit
   double mang = 200.; // min angle wrt of charge trk for selected gammas
   HepLorentzVector Pgg0; // sum of momentums of good photons
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] || gt == Slct.g4f[1] ) {
         // one of selected photons
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
   // Pi+ and Pi-
   HepLorentzVector Ppp = kmfit->pfit(0);
   HepLorentzVector Ppm = kmfit->pfit(1);
   HepLorentzVector P2pi = Ppp + Ppm;
   double M2pi  = P2pi.m();

   // two gammas
   HepLorentzVector Pg1 = kmfit->pfit(2);
   HepLorentzVector Pg2 = kmfit->pfit(3);
   HepLorentzVector Pgg = Pg1 + Pg2;
   double Mgg  = Pgg.m();

   hst[71]->Fill(Mgg*Mgg,M2pi*M2pi);
   hst[72]->Fill(Mgg*Mgg);

   // decay code to save in n-tuple
   double dec = double(Slct.decJpsi);
   if ( Slct.decJpsi < 0 ) {
      dec = double(Slct.decEta);
      hst[105]->Fill(dec);
   }

   // HepLorentzVector.perp == HepLorentzVector.vect.rho
   // HepLorentzVector.rho == HepLorentzVector.vect.mag
   Double_t xfill[] = {
      chisq,
      Ppp.rho(), Ppp.cosTheta(), Ppp.phi(), Slct.pid_Pip[0], Slct.E_p_Pip[0],
      Ppm.rho(), Ppm.cosTheta(), Ppm.phi(), Slct.pid_Pim[0], Slct.E_p_Pim[0],
      Pg1.e(), Pg1.cosTheta(), Pg1.phi(),
      Pg2.e(), Pg2.cosTheta(), Pg2.phi(),
      Pgg.rho(), Pgg.cosTheta(), Pgg.phi(),
      M2pi, Mgg,
      Slct.M2_pip_pim, Slct.M2_pip_eta,
      dec, Slct.Xisr
   };
   m_tuple[2]->Fill( xfill );

   return true;
}

//-----------------------------------------------------------------------------
// Vertex & Kinematic Fit for 6 gamma
//-----------------------------------------------------------------------------
static bool VertKinFit_6(SelectEtaTwoPi& Slct) {
//-----------------------------------------------------------------------------
   //Clearing Slct.g4f
   Slct.g4f.clear();
   Slct.g4f.resize(6, nullptr);

   // search for a good vertex:
   WTrackParameter wp[2] = {
      WTrackParameter(mpi, Slct.trk_Pip[0]->getZHelix(),
                      Slct.trk_Pip[0]->getZError() ),
      WTrackParameter(mpi, Slct.trk_Pim[0]->getZHelix(),
                      Slct.trk_Pim[0]->getZError() )
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

   hst[241]->Fill( double(ok) );
   if ( !ok ) {
      return false;
   }
   hst[1]->Fill(3);
   hst[242]->Fill(vtxfit->chisq(0));

   // 4C kinematic fit
   vtxfit->BuildVirtualParticle(0);
   vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0
   for(int k = 0; k < 2; k++) {
      wp[k] = vtxfit->wtrk(k);
   }

   KinematicFit* kmfit = KinematicFit::instance();

   double chisq = 9999.;
   int ng = Slct.gtrk.size();
   for(int i1 = 0; i1 < ng-5; i1++) {
      RecEmcShower* g1Trk = Slct.gtrk[i1];
      for(int i2 = i1+1; i2 < ng-4; i2++) {
         RecEmcShower* g2Trk = Slct.gtrk[i2];
         for(int i3 = i2+1; i3 < ng-3; i3++){
           RecEmcShower* g3Trk = Slct.gtrk[i3];
           for(int i4 = i3+1; i4 < ng-2; i4++){
             RecEmcShower* g4Trk = Slct.gtrk[i4];
             for(int i5 = i4+1; i5 < ng-1; i5++){
               RecEmcShower* g5Trk = Slct.gtrk[i5];
               for(int i6 = i5+1; i6 < ng; i6++){
                 RecEmcShower* g6Trk = Slct.gtrk[i6];

                 kmfit->init();
                 for(int k = 0; k < 2; k++) {
                    kmfit->AddTrack(k, wp[k]);
                 }
                 kmfit->AddTrack(2, 0.0, g1Trk);
                 kmfit->AddTrack(3, 0.0, g2Trk);
                 kmfit->AddTrack(4, 0.0, g3Trk);
                 kmfit->AddTrack(5, 0.0, g4Trk);
                 kmfit->AddTrack(6, 0.0, g5Trk);
                 kmfit->AddTrack(7, 0.0, g6Trk);
                 kmfit->AddFourMomentum(0, Slct.LVcms);
                 bool oksq = kmfit->Fit();
                 if ( oksq ) {
                    double chi2 = kmfit->chisq();
                    if ( chi2 < chisq && chi2 > 0. ) {
                       chisq = chi2;
                       Slct.g4f[0] = g1Trk;
                       Slct.g4f[1] = g2Trk;
                       Slct.g4f[2] = g3Trk;
                       Slct.g4f[3] = g4Trk;
                       Slct.g4f[4] = g5Trk;
                       Slct.g4f[5] = g6Trk;
                    }
                 }
               } //-----------------------------------------------End for(i6)
             } //-----------------------------------------------End for(i5)
           } //-----------------------------------------------End for(i4)
         } //-----------------------------------------------End for(i3)
      } //-----------------------------------------------End for(i2)
   } //-----------------------------------------------End for(i1)

   hst[251]->Fill(chisq);
   if ( !Slct.g4f[0] ) {
      return false;   // can not find good set of gammas
   }
   hst[1]->Fill(4);
   if ( isMC ) {
      if ( Slct.decJpsi > -1 ) {
         hst[103]->Fill( Slct.decJpsi );
      } else {
         hst[103]->Fill( Slct.decEta );
      }
   }

   // ++ repeat fit for six best photons ++
   kmfit->init();
   for ( int k = 0; k < 2; k++ ) {
      kmfit->AddTrack(k, wp[k]);
   }
   kmfit->AddTrack(2, 0.0, Slct.g4f[0]);
   kmfit->AddTrack(3, 0.0, Slct.g4f[1]);
   kmfit->AddTrack(4, 0.0, Slct.g4f[2]);
   kmfit->AddTrack(5, 0.0, Slct.g4f[3]);
   kmfit->AddTrack(6, 0.0, Slct.g4f[4]);
   kmfit->AddTrack(7, 0.0, Slct.g4f[5]);
   kmfit->AddFourMomentum(0, Slct.LVcms);
   bool oksq = kmfit->Fit();
   if ( !oksq || fabs(kmfit->chisq()-chisq) > 0.1 ) {
      cout << " WARNING: Bad second 4C fit for 6 gammas";
      if( oksq ) {
         cout << "6g: chisq= "<<chisq<<" new fit chisq= "<<kmfit->chisq();
      }
      cout << endl;
      Warning("Bad second 4C fit for 6gammas");
      return false;
   }
   chisq=kmfit->chisq();

   //--------------------------------------------------------------------------
   // Final histograms and ntuples
   //--------------------------------------------------------------------------

   double Eg_max = 0.; // max momentum of gammas not selected by fit
   double mang = 200.; // min angle wrt of charge trk for selected gammas
   HepLorentzVector Pgg0; // sum of momentums of good photons
   for(int i = 0; i < ng; i++ ) {
      auto gt = Slct.gtrk[i];
      if( gt == Slct.g4f[0] ||
          gt == Slct.g4f[1] ||
          gt == Slct.g4f[2] ||
          gt == Slct.g4f[3] ||
          gt == Slct.g4f[4] ||
          gt == Slct.g4f[5]) {
         // one of selected photons
         Pgg0 += Slct.Pg[i];
         double deg = RtoD(Slct.angt[i]);
         hst[261]->Fill(deg);
         if( deg < mang ) {
            mang = deg;
         }
         continue;
      }
      // rejected photons
      double Eg = Slct.Pg[i].e();
      hst[262]->Fill(Eg);
      if( Eg > Eg_max ) {
         Eg_max = Eg;
      }
   }
   if ( Eg_max > 0.01 ) {
      hst[263]->Fill(Eg_max);
   }

   // Momentums after kinematic corrections:
   // Pi+ and Pi-
   HepLorentzVector Ppp = kmfit->pfit(0);
   HepLorentzVector Ppm = kmfit->pfit(1);
   HepLorentzVector P2pi = Ppp + Ppm;
   double M2pi  = P2pi.m();

   // two gammas
   HepLorentzVector Pg1 = kmfit->pfit(2);
   HepLorentzVector Pg2 = kmfit->pfit(3);
   HepLorentzVector Pg3 = kmfit->pfit(4);
   HepLorentzVector Pg4 = kmfit->pfit(5);
   HepLorentzVector Pg5 = kmfit->pfit(6);
   HepLorentzVector Pg6 = kmfit->pfit(7);
   HepLorentzVector Pgg = Pg1 + Pg2 + Pg3 + Pg4 + Pg5 + Pg6;
   double Mgg  = Pgg.m();

   hst[271]->Fill(Mgg*Mgg,M2pi*M2pi);
   hst[272]->Fill(Mgg*Mgg);

   // decay code to save in n-tuple
   double dec = double(Slct.decJpsi);
   if ( Slct.decJpsi < 0 ) {
      dec = double(Slct.decEta);
      hst[105]->Fill(dec);
   }

   // HepLorentzVector.perp == HepLorentzVector.vect.rho
   // HepLorentzVector.rho == HepLorentzVector.vect.mag
   Double_t xfill[] = {
      chisq,
      Ppp.rho(), Ppp.cosTheta(), Ppp.phi(), Slct.pid_Pip[0], Slct.E_p_Pip[0],
      Ppm.rho(), Ppm.cosTheta(), Ppm.phi(), Slct.pid_Pim[0], Slct.E_p_Pim[0],
      Pg1.e(), Pg1.cosTheta(), Pg1.phi(),
      Pg2.e(), Pg2.cosTheta(), Pg2.phi(),
      Pg3.e(), Pg3.cosTheta(), Pg3.phi(),
      Pg4.e(), Pg4.cosTheta(), Pg4.phi(),
      Pg5.e(), Pg5.cosTheta(), Pg5.phi(),
      Pg6.e(), Pg6.cosTheta(), Pg6.phi(),
      Pgg.rho(), Pgg.cosTheta(), Pgg.phi(),
      M2pi, Mgg,
      Slct.M2_pip_pim, Slct.M2_pip_eta,
      dec, Slct.Xisr
   };
   m_tuple[6]->Fill( xfill );

   return true;
}

//-----------------------------------------------------------------------------
bool Eta_TwoPiEvent( ReadDst*       selector,
                      TEvtHeader*    m_TEvtHeader,
                      TDstEvent*     m_TDstEvent,
                      TEvtRecObject* m_TEvtRecObject,
                      TMcEvent*      m_TMcEvent,
                      TTrigEvent*    m_TTrigEvent,
                      TDigiEvent*    m_TDigiEvent,
                      THltEvent*     m_THltEvent        ) {
//-----------------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   m_abscor->AbsorptionCorrection(selector);

   SelectEtaTwoPi Slct; // information for the current event

   //--------------------------------------------------------------------------
   //-- Get event information --
   //--------------------------------------------------------------------------
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
           << " McParticles= " << m_TMcEvent->getMcParticleCol()->GetEntries()
           << endl;
      Warning("Incorrect number of MC particles");
      return false;
   }

   if ( isMC ) {
     double Ecms_local = Ecms;
     try{
       Dalitz_Parameters.at(Ecms_local);
     }
     catch(std::out_of_range const&){
       Ecms_local = round(Ecms_local*100000.)/100000.;
     }
     if (!( Ecms_local == 2.12655 || Ecms_local == 2.396 || Ecms_local == 2.900 || Ecms_local == 3.080) && ( hst[151]->GetNbinsX() != 50 )){
       hst[151]->Rebin(2); hst[152]->Rebin(2); hst[153]->Rebin(2);
       hst[154]->Rebin(2); hst[155]->Rebin(2); hst[156]->Rebin(2);
     }
     double dp_lim_x = get<0>(Dalitz_Parameters.at(Ecms_local));
     double dp_lim_y = get<1>(Dalitz_Parameters.at(Ecms_local));
     if(hst2[2]->GetXaxis()->GetXmax() != dp_lim_x){
       hst2[2]->GetXaxis()->SetLimits(0., dp_lim_x);
     }
     if(hst2[2]->GetYaxis()->GetXmax() != dp_lim_y){
       hst2[2]->GetYaxis()->SetLimits(0., dp_lim_y);
     }
     if(hst2[6]->GetXaxis()->GetXmax() != dp_lim_x){
       hst2[6]->GetXaxis()->SetLimits(0., dp_lim_x);
     }
     if(hst2[6]->GetYaxis()->GetXmax() != dp_lim_y){
       hst2[6]->GetYaxis()->SetLimits(0., dp_lim_y);
     }
   }

   if ( isMC ) {
      FillHistoMC(selector,Slct);       // MC histo
      //ISRhistoMC( selector, Slct);      // MC ISR histo
   }

   Hep3Vector xorigin = getVertexOrigin(runNo);
   Slct.xorig = xorigin;
   //--------------------------------------------------------------------------
   if ( !ChargedTracksTwoPi(selector, Slct) ) {
      return false;
   }
   //--------------------------------------------------------------------------
   if ( init_selection ) {
      return true;
   }
   //--------------------------------------------------------------------------
   int n_gammas = NeutralTracks(selector, Slct);
   if ( n_gammas < 2 ) {
      return false;
   }
   if ( VertKinFit_2(Slct) ) {
     return true;
   }
   if ( n_gammas < 6 ) {
     return false;
   }
   if ( VertKinFit_6(Slct) ) {
     return true;
   }
   return false;
}

//-----------------------------------------------------------------------------
void Eta_TwoPiEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
   }

   string module = string(__func__);
   module = module.substr(0,module.size()-6);
   if ( warning_msg.empty() ) {
      cout << " There are no warnings in " << module << endl;
   } else {
      cout << " Check output for WARNINGS in " << module << endl;
      for(auto it = warning_msg.begin(); it != warning_msg.end(); ++it) {
         cout << it->first << " : " << it->second << endl;
      }
   }
}

#ifdef __cplusplus
}
#endif
