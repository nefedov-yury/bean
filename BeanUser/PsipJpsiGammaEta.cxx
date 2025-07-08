//==================================================================//
//                                                                  //
// PsipJpsiGammaEta:                                                //
// Study of the eta->2gammas detection efficiency in process:       //
// e+ e- -> Psi(2S) -> pi+ pi- J/Psi                                //
//                             |-> gamma  eta                       //
//                                         |-> 2gamma               //
//                                                                  //
//==================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
// #include <TNtupleD.h>
#include <TTree.h>
#include <TDatabasePDG.h>

// #include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>

typedef HepGeom::Point3D<double> HepPoint3D;
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

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

//--------------------------------------------------------------------
// {{{1 Structure to save variables for a single event
//--------------------------------------------------------------------
struct PimPipGammas {
   // Run-info:
   int runNo;              // run-number
   int event;              // event-number
   HepLorentzVector LVcms; // Momentum in center of mass system
   Hep3Vector xorig;       // interaction point from DB

   // MC information:
   int decPsip, decJpsi;      // decay codes for MC events
   int dec_eta;               // 1 if eta->2gamma
   HepLorentzVector mcLVjpsi; // true 4-momentum of J/Psi
   HepLorentzVector mcLVg;    // true 4-mom. of gamma from J/Psi decay
   HepLorentzVector mcLVeta;  // true 4-mom. of eta from J/Psi decay
   int mcidx_g;               // index in LVg which corresponds mcLVg

   // soft pions:
   RecMdcKalTrack* trk_Pip; // pi+
   RecMdcKalTrack* trk_Pim; // pi-
   HepLorentzVector LVjpsi; // estimated 4-momentum J/Psi
   double Mrec;             // recoil mass
   double cosPM, invPM;     // cos(Theta pi+ pi-) & inv.mass

   // gamma tracks
   vector<RecEmcShower*> trk_g;
   vector<HepLorentzVector> LVg; // 4-momentum of gammas

   // candidates for gamma in decay J/Psi -> gamma eta
   vector<int> idx_g;           // indexes in LVg vector
   vector<double> M2rec_g;      // vec of recoil mass square of gamma

   PimPipGammas() {
      decPsip = decJpsi = -1;
      dec_eta = 0;
      mcidx_g = -1;
      trk_Pip = trk_Pim = nullptr;
   }
};

//--------------------------------------------------------------------
// {{{1 Structure for root-trees
//--------------------------------------------------------------------
// eta reconstruction efficiency
struct Xnt_gammaeta {
   double Ptp,Ptm,Mrec; // Pt(pi+), Pt(pi-), Mrec(pi+pi-)
   double Eg0,Cg0;      // E,cos(Theta) of gamma in J/Psi -> g0 eta
   double Peta,Ceta;    // P,cos(Theta) of eta
   double Eg1,Cg1;      // E,cos(Theta) of rec gamma in eta -> g1 g2
   double Eg2,Cg2;      // E,cos(Theta) of predicted gamma (g2)
   double Egr,Cgr;      // E and cos(Theta) of found gamma (gr)
   double rE,dTh;       // relations btw predicted(g2) and found(gr)
                        // * variables used in selection:
   double M2gr;         // Mrec^2(pi+pi-g0)
   double M2gg;         // Minv^2(g1 g2)
   double m2fr;         // Mrec^2(pi+pi-g0g1) or g2 missing mass^2
   double mggf;         // Minv(g1,gr) after 4C kinematic constraints
   double ch2;          // chi^2 of 4C kinematic constraints
                        // * MC-variables
   int    decj;         // MC: decay codes of J/Psi
   int    dgam;         // MC: 1 if g0 found correctly
   int    deta;         // MC: 1 if eta decays into two gammas
};

//--------------------------------------------------------------------
// {{{1 Global variables
//--------------------------------------------------------------------
static const double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
static const double mpsip  = 3.68610;  // 3686.10   +/- 0.06    MeV
static const double mjpsi  = 3.096900; // 3096.900  +/- 0.006   MeV
static const double mpi    = 0.13957;  // 139.57039 +/- 0.00018 MeV
static const double mpi0   = 0.13498;  // 134.9768  +/- 0.0005  MeV
static const double meta   = 0.547862; // 547.862   +/- 0.017   MeV
static const double momega = 0.78265;  // 782.65    +/- 0.12    MeV
static const double mk     = 0.493677; // 493.677   +/- 0.016   MeV
static const double mk0    = 0.497611; // 497.611   +/- 0.013   MeV
static const double mphi   = 1.019461; //1019.461   +/- 0.016   MeV
static const double Mn     = 0.939565; // 939.5654133+/-0.0000058 MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

// histograms
static vector<TH1*> hst;
// static vector<TNtupleD*> m_tuple;
static TTree* m_nt_geta = nullptr;
static Xnt_gammaeta xnt_geta;

// decays classification
static DecayTable JpsiTbl;

// container for warnings
static map<string,int> warning_msg;

// expect one data taking period for all events
static int DataPeriod = 0;
static bool isMC = false;

//--------------------------------------------------------------------
// {{{1 Functions: use C-linkage names
//--------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

//--------------------------------------------------------------------
static inline void Warning(const string& msg) {
   warning_msg[msg] += 1;
}
//--------------------------------------------------------------------

//--------------------------------------------------------------------
static inline double RtoD(double ang) {
   return ang*180/M_PI;
}
//--------------------------------------------------------------------

//--------------------------------------------------------------------
static inline double SQ(double x) {
   return x*x;
}
//--------------------------------------------------------------------

//--------------------------------------------------------------------
static bool SelectPM(double cosPM, double invPM) {
//--------------------------------------------------------------------
   // criteria for pair pions (+/-) for selection good Mrec
   bool ret = true;
   if ( (cosPM > 0.90) ||         // flying in one direction
         (abs(invPM-mk0) < 0.008 ) // pions from K^0_s decay
      ) {
      ret = false;
   }
   return ret;
}

// {{{1 StartJob, book histograms
//--------------------------------------------------------------------
void PsipJpsiGammaEtaStartJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
      printf("Masses of particles from PDG\n");
      printf("M_Psi(2S) = %f MeV\n", mpsip*1e3);
      printf("M_J/Psi   = %f MeV\n", mjpsi*1e3);
      printf("M_pi^\\pm  = %f MeV\n", mpi*1e3);
      printf("M_pi^0    = %f MeV\n", mpi0*1e3);
      printf("M_eta     = %f MeV\n", meta*1e3);
      printf("Momega    = %f MeV\n", momega*1e3);
      printf("M_K^\\pm   = %f MeV\n", mk*1e3);
      printf("M_K^0     = %f MeV\n", mk0*1e3);
      printf("M_phi     = %f MeV\n", mphi*1e3);
      printf("M_n       = %f MeV\n", Mn*1e3);
   }

   hst.resize(300,nullptr);
   // m_tuple.resize(5,nullptr);

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
      // Warning("EventTagSvc has already been initialized");
   }

   // initialize DatabaseSvc -----------------------------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      dbs->SetDBFilePath(
            selector->AbsPath("Analysis/DatabaseSvc/dat"));
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
      // Warning("MagneticFieldSvc has already been initialized");
   }

   // set path for ParticleID algorithm ------------------------------
   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   // Book histograms ------------------------------------------------

   hst[1] = new TH1D("All_cuts","selections cuts", 20,-0.5,19.5);

   // ChargedTracks:
   hst[11] = new TH1D("Rxy","R_{xy}", 200,-2.,2.);
   hst[12] = new TH1D("Rz","R_{z}", 200,-20.,20.);
   hst[13] = new TH1D("cos_theta","cos(#theta)", 200,-1.,1.);
   hst[14] = new TH1D("theta","#theta", 180,0.,180.);
   hst[15] = new TH1D("Pid_clpi","lg(CL_{#pi})", 100,-4.,0.);
   hst[16] = new TH1D("Pid_ispi","1 - #pi, 0 - another particle",
         2,-0.5,1.5);
   hst[18] = new TH1D("pi_QP","Charged momentum (Q*P)", 200,-1.,1.);
   hst[19] = new TH1D("pi_ct","soft pions cos(#theta)", 200,-1.,1.);

   hst[20] = new TH1D("Nch","number of good trk", 20,-0.5,19.5);
   hst[21] = new TH2D("pi_pm","Npi+ vs Npi-", 5,-0.5,4.5, 5,-0.5,4.5);

   hst[22] = new TH1D("Mrec","recol mass of pi+ pi-", 100,3.,3.2);
   hst[23] = new TH1D("cosPM", "cos(Theta_pi+_pi-)", 100,-1.0,1.0);
   hst[24] = new TH1D("invPM", "inv.mass(pi+ pi-)", 500,0.25,0.75);

   hst[31] = new TH1D("SMrec",  "Mrec(pi+ pi-)", 100,3.05,3.15);
   hst[32] = new TH1D("ScosPM", "cos(Theta_pi+_pi-)", 100,-1.0,1.0);
   hst[33] = new TH1D("SinvPM", "inv.mass(pi+ pi-)", 500,0.25,0.75);

   hst[41] = new TH1D("Pip_P","P(#pi+)", 100,0.,0.5);
   hst[42] = new TH1D("Pip_T","Pt(#pi+)", 100,0.,0.5);
   hst[43] = new TH1D("Pip_C","cos #Theta(#pi+)", 200,-1.0,1.0);
   hst[44] = new TH1D("Pim_P","P(#pi-)", 100,0.,0.5);
   hst[45] = new TH1D("Pim_T","Pt(#pi-)", 100,0.,0.5);
   hst[46] = new TH1D("Pim_C","cos #Theta(#pi-)", 200,-1.0,1.0);

   // NeutralTracks:
   hst[51] = new TH1D("n_ang","angle with closest chg.trk",
         180,0.,180.);
   hst[52] = new TH1D("n_Ng","N_{#gamma} in event", 11,-0.5,10.5);
   hst[53] = new TH1D("n_Eg","E_{#gamma}", 200,0.,2.);

   // search for "hot" chanels in EMC
   // hst[55] = new TH2D("mapB","3) barrel Z vs phi",
         // 270,-135.,+135., 360,-180.,180.);
   // hst[56] = new TH2D("mapC1","3) endcap Z<0 rho vs phi",
         // 40, 50.,90., 360,-180.,180.);
   // hst[57] = new TH2D("mapC2","3) endcap Z>0 rho vs phi",
         // 40, 50.,90., 360,-180.,180.);

   hst[61] = new TH1D("Mrg2","M^{2}rec(#pi+,#pi-,#gamma)",
         200,-0.2,0.8);
   hst[62] = new TH1D("Eg_jpsi","E#gamma in J/#Psi rest.sys.",
         100,1.,2.);
   hst[63] = new TH1D("Ngj","N#gamma selected", 5,-0.5,4.5);

   hst[71] = new TH1D("Mg2_pi0","M^2(#gamma#gamma)", 200,0.,0.04);
   hst[72] = new TH1D("Mg2_all","M^2(#gamma#gamma)", 500,0.,10.);
   hst[73] = new TH1D("Npi0","N(#pi^{0})", 5,-0.5,4.5);

   // EtaEff:
   hst[81] = new TH1D("Mgg2_a","M^{2}(2#gamma) all", 200,0.,0.6);
   hst[82] = new TH1D("M2fr_min","M^2(fr) minimal", 100,0.,0.02);
   hst[85] = new TH1D("rE","Eg(pred)/Eg(rec)", 250,0.,2.5);
   hst[86] = new TH1D("dTh","#delta#Theta (pred-rec)", 100,0.,20.);
   hst[87] = new TH2D("mch2D","ratE vs dTheta;#Theta",
         50,0.,20.,50,0.5,2.0);

   hst[91] = new TH1D("M2rl_min","M^2(real) min", 200,-0.02,0.02);
   hst[92] = new TH1D("Mgg2_rl","M^{2}(2#gamma) real", 200,0.,0.6);

   hst[95] = new TH1D("vtx_chi2", "vertex fit #chi^2", 100,0.,100.);
   hst[96] = new TH1D("fit_ch2","kin-fit: #chi^{2}", 250,0.,250.);
   hst[97] = new TH1D("Mgg_fit","M(2#gamma) fit", 100,0.5,0.6);

   // Monte Carlo histograms -----------------------------------------
   hst[100] = new TH1D("mc_dec0", "decPsi(2S) nocut",256,-0.5,255.5);
   hst[101] = new TH1D("mc_dec2", "decPsi(2S) cut#2",256,-0.5,255.5);
   hst[102] = new TH1D("mc_dec3", "decPsi(2S) cut#3",256,-0.5,255.5);
   hst[103] = new TH1D("mc_dec4", "decPsi(2S) cut#4",256,-0.5,255.5);
   hst[104] = new TH1D("mc_dec5", "decPsi(2S) cut#5",256,-0.5,255.5);

   hst[105] = new TH1D("mc_dcj0", "dec J/psi nocut",300,-0.5,299.5);
   hst[106] = new TH1D("mc_dcj2", "dec J/psi cut#2",300,-0.5,299.5);
   hst[107] = new TH1D("mc_dcj3", "dec J/psi cut#3",300,-0.5,299.5);
   hst[108] = new TH1D("mc_dcj4", "dec J/psi cut#4",300,-0.5,299.5);
   hst[109] = new TH1D("mc_dcj5", "dec J/psi cut#5",300,-0.5,299.5);

   // FillHistoMC:
   hst[111] = new TH1D("mc_pdg",
         "PDG codes of all particles", 2001,-1000.5,1000.5);
   hst[112] = new TH1D("mc_pdg0",
         "PDG of particles from primary vertex", 2001,-1000.5,1000.5);

   hst[121] = new TH1D("mc_PsipPt", "Pt of #Psi(2S)", 100,0.,1.);
   hst[122] = new TH1D("mc_PsipC", "cos(#Theta) of Psi(2S)",
         100,-1.,1.);

   hst[123] = new TH1D("mc_JpsiP", "Momentum of J/#Psi", 100,0.,1.);
   hst[124] = new TH1D("mc_JpsiPt","Pt of J/#Psi", 100,0.,1.);
   hst[125] = new TH1D("mc_JpsiC", "cos(#Theta) of J/#Psi",
         100,-1.,1.);

   hst[126] = new TH1D("mc_PipP", "Momentum of #pi^{+}", 100,0.,1.);
   hst[127] = new TH1D("mc_PipC", "cos(#Theta) of #pi^{+}",
         100,-1.,1.);
   hst[128] = new TH1D("mc_PimP", "Momentum of #pi^{-}", 100,0.,1.);
   hst[129] = new TH1D("mc_PimC", "cos(#Theta) of #pi^{-}",
         100,-1.,1.);

   hst[131] = new TH1D("mc_EtaP", "Momentum of #eta", 100,1.,2.);
   hst[132] = new TH1D("mc_EtaC", "cos(#Theta) of #eta", 100,-1.,1.);
   hst[133] = new TH1D("mc_GamE", "E(#gamma)", 100,1.,2.);
   hst[134] = new TH1D("mc_GamC", "cos(#Theta) of #gamma",
         100,-1.,1.);

   hst[135] = new TH1D("mc_EtaPrc", "P(#eta) in J/#Psi RS",
         100,1.4,1.6);
   hst[136] = new TH1D("mc_EtaCrc", "cos(#eta) in J/#Psi RS",
         100,-1.,1.);
   hst[137] = new TH1D("mc_GamPrc", "P(#gamma) in J/#Psi RS",
         100,1.4,1.6);
   hst[138] = new TH1D("mc_GamCrc", "cos(#gamma) in J/#Psi RS",
         100,-1.,1.);
   hst[139] = new TH1D("mc_A_EtaG",
         "angle(#eta,#gamma) in J/#Psi RS", 100,179.005,181.005);
   hst[140] = new TH1D("mc_GamEtaE","E(#gamma) from #eta decay",
         200,0.,2.);
   hst[141] = new TH1D("mc_GamEtaEirc","E(g) from #eta dec. RS",
         200,0.,2.);
   hst[142] = new TH1D("mc_GamEtaAng","Angle(g1,g2)", 180,0.,180.);

   // MC ChargedTracks:
   hst[151] = new TH1D("mc_dPjpsi","dP(J/#Psi) mc-rec",
         100,-0.1,0.1);
   hst[152] = new TH1D("mc_dAjpsi","dAngle(J/#Psi) mc-rec",
         100,0.,0.2);

   // MatchGammaMC:
   hst[161] = new TH1D("mc_dPx_g","MC-REC dPx #gamma",
         150,-0.075,0.075);
   hst[162] = new TH1D("mc_dPy_g","MC-REC dPy #gamma",
         150,-0.075,0.075);
   hst[163] = new TH1D("mc_dPz_g","MC-REC dPz #gamma",
         150,-0.075,0.075);
   hst[164] = new TH1D("mc_dPa_g","MC-REC dAngle #gamma",
         100,0.,0.1);
   // hst[165] = new TH1D("mc_dEm_g","MC-REC min dE #gamma",
         // 200,-0.1,0.1);
   hst[166] = new TH1D("mc_dEg","dE(#gamma) mc-rec", 100,-0.2,0.2);
   hst[167] = new TH1D("mc_dAg","dAngle(#gamma) mc-rec",100,0.,0.1);
   hst[168] = new TH1D("mc_EgF","no MC-REC E(#gamma) ", 100,1.,2.);
   hst[169] = new TH1D("mc_CgF","no MC-REC cos(#gamma)",100,-1.,1.);

   // MC NeutralTracks:
   hst[171] = new TH1D("mcMrg2T","M^{2}rec(#pi+,#pi-,#gamma) T",
         200,-.2,.8);
   hst[172] = new TH1D("mcMrg2F","M^{2}rec(#pi+,#pi-,#gamma) F",
         200,-.2,.8);
   hst[173] = new TH1D("mcMrg2TT","M^{2}rec(#pi+,#pi-,#gamma) TT",
         200,-.2,.8);
   hst[174] = new TH1D("mcMrg2TE","M^{2}rec(#pi+,#pi-,#gamma) TE",
         200,-.2,.8);
   hst[176] = new TH1D("mcEg_jpsiT","E#gamma in J/#Psi RS T",
         100,1.,2.);
   hst[177] = new TH1D("mcEg_jpsiF","E#gamma in J/#Psi RS F",
         100,1.,2.);
   hst[178] = new TH1D("mcEg_jpsiTT","E#gamma in J/#Psi RS TT",
         100,1.,2.);
   hst[179] = new TH1D("mcEg_jpsiTF","E#gamma in J/#Psi RS TF",
         100,1.,2.);

   hst[181] = new TH1D("mcEg_jpsiST","E#gamma in J/#Psi RS ST",
         100,1.3,1.7);
   hst[182] = new TH1D("mcEg_jpsiSF","E#gamma in J/#Psi RS SF",
         100,1.3,1.7);
   hst[185] = new TH1D("mcNgjT","N#gamma selected T", 5,-0.5,4.5);
   hst[186] = new TH1D("mcNgjF","N#gamma selected F", 5,-0.5,4.5);
   hst[187] = new TH1D("mcNgjTT","N#gamma selected TT", 5,-0.5,4.5);

   hst[191] = new TH1D("mcMg2_pi0T","M^2(#gamma#gamma) T",
         200,0.,0.04);
   hst[192] = new TH1D("mcMg2_pi0F","M^2(#gamma#gamma) F",
         200,0.,0.04);
   hst[193] = new TH1D("mcMg2_pi0TE","M^2(#gamma#gamma) TE",
         200,0.,0.04);

   hst[195] = new TH1D("mc_dcjPi0", "dec J/psi Npi0>0",
         300,-0.5,299.5);
   hst[196] = new TH1D("mcNpi0T","N(#pi^{0}) T", 5,-0.5,4.5);
   hst[197] = new TH1D("mcNpi0F","N(#pi^{0}) F", 5,-0.5,4.5);
   hst[198] = new TH1D("mcNpi0TE","N(#pi^{0}) TE", 5,-0.5,4.5);

   // MC EtaEff:
   hst[201] = new TH1D("mcMgg2_aT","M^{2}(2#gamma) all T",
         200,0.,0.6);
   hst[202] = new TH1D("mcMgg2_aTT","M^{2}(2#gamma) all TT",
         200,0.,0.6);
   hst[203] = new TH1D("mcMgg2_aF","M^{2}(2#gamma) all F",
         200,0.,0.6);
   hst[204] = new TH1D("mcMgg2_aTE","M^{2}(2#gamma) all TE",
         200,0.,0.6);
   hst[205] = new TH1D("mcM2fr_minT","M^2(fr) minimal T",
         100,0.,0.02);
   hst[206] = new TH1D("mcM2fr_minTT","M^2(fr) minimal TT",
         100,0.,0.02);
   hst[207] = new TH1D("mcM2fr_minF","M^2(fr) minimal F",
         100,0.,0.02);
   hst[208] = new TH1D("mcM2fr_minTE","M^2(fr) minimal TE",
         100,0.,0.02);

   hst[211] = new TH1D("mcrET","Eg(pred)/Eg(rec) T", 250,0.,2.5);
   hst[212] = new TH1D("mcdThT","#delta#Theta (pred-rec) T",
         100,0.,20.);
   hst[213] = new TH1D("mcrEF","Eg(pred)/Eg(rec) F", 250,0.,2.5);
   hst[214] = new TH1D("mcdThF","#delta#Theta (pred-rec) F",
         100,0.,20.);
   hst[215] = new TH2D("mc_2DT","ratE vs dTheta T;#Theta",
         50,0.,20.,50,0.5,2.0);
   hst[216] = new TH2D("mc_2DF","ratE vs dTheta F;#Theta",
         50,0.,20.,50,0.5,2.0);

   hst[221] = new TH1D("mcM2rl_minT","M^2(real) T", 200,-0.02,0.02);
   hst[222] = new TH1D("mcM2rl_minTT","M^2(real) TT", 200,-0.02,0.02);
   hst[223] = new TH1D("mcM2rl_minF","M^2(real) F", 200,-0.02,0.02);
   hst[224] = new TH1D("mcMgg2_rlT","M^{2}(2#gamma) real T",
         200,0.,0.6);
   hst[225] = new TH1D("mcMgg2_rlTT","M^{2}(2#gamma) real TT",
         200,0.,0.6);
   hst[226] = new TH1D("mcMgg2_rlF","M^{2}(2#gamma) real F",
         200,0.,0.6);
   hst[228] = new TH1D("mcMgg_fitT","M(2#gamma) fit T", 100,0.5,0.6);
   hst[229] = new TH1D("mcMgg_fitF","M(2#gamma) fit F", 100,0.5,0.6);


   // Tree for eta efficiency study in gamma eta
   m_nt_geta = new TTree("eff_eta", "eta efficiency in gamma eta");
   m_nt_geta->Branch("Ptp"   , &xnt_geta.Ptp   );
   m_nt_geta->Branch("Ptm"   , &xnt_geta.Ptm   );
   m_nt_geta->Branch("Mrec"  , &xnt_geta.Mrec  );
   m_nt_geta->Branch("Eg0"   , &xnt_geta.Eg0   );
   m_nt_geta->Branch("Cg0"   , &xnt_geta.Cg0   );
   m_nt_geta->Branch("Peta"  , &xnt_geta.Peta  );
   m_nt_geta->Branch("Ceta"  , &xnt_geta.Ceta  );
   m_nt_geta->Branch("Eg1"   , &xnt_geta.Eg1   );
   m_nt_geta->Branch("Cg1"   , &xnt_geta.Cg1   );
   m_nt_geta->Branch("Eg2"   , &xnt_geta.Eg2   );
   m_nt_geta->Branch("Cg2"   , &xnt_geta.Cg2   );
   m_nt_geta->Branch("Egr"   , &xnt_geta.Egr   );
   m_nt_geta->Branch("Cgr"   , &xnt_geta.Cgr   );
   m_nt_geta->Branch("rE"    , &xnt_geta.rE    );
   m_nt_geta->Branch("dTh"   , &xnt_geta.dTh   );

   m_nt_geta->Branch("M2gr"  , &xnt_geta.M2gr  );
   m_nt_geta->Branch("M2gg"  , &xnt_geta.M2gg  );
   m_nt_geta->Branch("m2fr"  , &xnt_geta.m2fr  );
   m_nt_geta->Branch("mggf"  , &xnt_geta.mggf  );
   m_nt_geta->Branch("ch2"   , &xnt_geta.ch2   );

   m_nt_geta->Branch("decj"  , &xnt_geta.decj  );
   m_nt_geta->Branch("dgam"  , &xnt_geta.dgam  );
   m_nt_geta->Branch("deta"  , &xnt_geta.deta  );


   // register in selector to save in given directory
   const char* SaveDir = "PsipJpsiGammaEta";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
   // VecObj ntuples(m_tuple.begin(),m_tuple.end());
   // selector->RegInDir(ntuples,SaveDir);

   VecObj Vreg;
   Vreg.push_back(m_nt_geta);
   selector->RegInDir(Vreg,SaveDir);
}

// {{{1 getVertexOrigin()
//--------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//--------------------------------------------------------------------
   static int save_runNo = 0;
   static Hep3Vector xorigin;

   if ( runNo == save_runNo ) {
      return xorigin;
   }

   // update vertex for new run
   xorigin.set(0.,0.,0.);
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();

#if (BOSS_VER > 700)
   string BossVer("7.0.9");
#else
   string BossVer("6.6.4");
   int run = abs(runNo);
   if ( (run >= 25338 && run <= 27090)) {  // 2012 Psi(2S)
      BossVer = "6.6.4.p03";
   }
#endif
   vtxsvc->SetBossVer(BossVer);
   vtxsvc->handle(runNo);

   if ( vtxsvc->isVertexValid() ) {
      double* dbv = vtxsvc->PrimaryVertex();
      xorigin.set(dbv[0],dbv[1],dbv[2]);
      if( verbose ) {
         cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
      }
   } else {
      cout << " FATAL ERROR:"
         " Cannot obtain vertex information for run#" << runNo
         << endl;
      exit(1);
   }

   save_runNo = runNo;
   return xorigin;
}

// {{{1 FillHistoMC()
//--------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, PimPipGammas& ppg) {
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
   if ( (evTag & 0xF) != 5 ) { // 5 is Psi(2S) event
      printf(" WARNING: MC is not Psi(2s) evTag= 0x%08x\n",evTag);
      Warning("MC is not Psi(2s)");
      evTag = 0;
   }

   // decay code of Psi(2S):
   int decPsip  = (evTag >> 8) & 0xFF;
   ppg.decPsip = decPsip;
   hst[100]->Fill(decPsip);

   // decay code of J/Psi from Psi(2S) -> pi+ pi- J/Psi
   int decJpsi = 0;
   if( decPsip == 64 ) {
      decJpsi = ( evTag >> 16 ) & 0xFF;
      ppg.decJpsi = decJpsi;
      hst[105]->Fill(decJpsi);
   }

   int idx_psip=-1;
   int idx_jpsi=-1;
   int idx_eta=-1;
   JpsiTbl.Reset();

   // momentum gamma & eta in rest frame of J/Psi
   Hep3Vector beta_jpsi;
   HepLorentzVector LVeta, LVgamma;
   vector<Hep3Vector> gammas_from_eta;

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      hst[111]->Fill(part_pdg);
      if ( part->getMother() == -99 ) { // primary vertex
         hst[112]->Fill(part_pdg);
      }

      if ( part_pdg == 100443 ) { // Psi(2S)
         hst[121]->Fill(Vp.rho());
         hst[122]->Fill(Vp.cosTheta());
         idx_psip = part->getTrackIndex();
      }

      if ( part->getMother() == idx_psip ) {      // decays of Psi(2S)
         if ( decPsip == 64 ) {                    // -> pi+ pi- J/Psi
            if ( part_pdg == 443 ) {                // J/Psi
               hst[123]->Fill(Vp.mag());
               hst[124]->Fill(Vp.rho());
               hst[125]->Fill(Vp.cosTheta());
               idx_jpsi = part->getTrackIndex();
               HepLorentzVector LVjpsi( Vp,
                     sqrt(Vp.mag2() + SQ(mjpsi)) );
               ppg.mcLVjpsi = LVjpsi;
               beta_jpsi = LVjpsi.boostVector();
            } else if ( part_pdg == 211 ) {         // pi+
               hst[126]->Fill(Vp.mag());
               hst[127]->Fill(Vp.cosTheta());
            } else if ( part_pdg == -211 ) {        // pi-
               hst[128]->Fill(Vp.mag());
               hst[129]->Fill(Vp.cosTheta());
            }
         }
      } // end Psi(2S) decays

      if( part->getMother() == idx_jpsi ) {  // decays of J/Psi
         // collect decays of J/psi
         JpsiTbl.vecdec.push_back(part_pdg);
      }

      if ( decJpsi == 22 ) {                     // J/Psi -> eta gamma
         if ( part->getMother() == idx_jpsi ) {
            if ( part_pdg == 221 ) {               // eta
               idx_eta = part->getTrackIndex();
               hst[131]->Fill(Vp.mag());
               hst[132]->Fill(Vp.cosTheta());

               HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(meta)) );
               ppg.mcLVeta = LV;
               LV.boost(-beta_jpsi);
               LVeta=LV;
               hst[135]->Fill( LV.vect().mag() );
               hst[136]->Fill( LV.cosTheta() );

            } else if ( part_pdg == 22 ) {         // gamma
               hst[133]->Fill(Vp.mag());
               hst[134]->Fill(Vp.cosTheta());

               HepLorentzVector LV( Vp, Vp.mag() );
               ppg.mcLVg = LV;
               LV.boost(-beta_jpsi);
               LVgamma=LV;
               hst[137]->Fill( LV.vect().mag() );
               hst[138]->Fill( LV.cosTheta() );
            }
         } // end J/Psi decay

         if ( part->getMother() == idx_eta ) {     // eta -> 2 gamma
            if ( part_pdg == 22 ) {
               gammas_from_eta.push_back(Vp);
            } else {
               ppg.dec_eta = -1;       // it's not 2gamma decay
            }
         }
      }
   } // end of while

   if ( decJpsi == 22 ) {     // J/Psi -> eta gamma
      // check angle(eta gamma) in the rest frame of J/Psi
      double alpha= LVeta.angle(LVgamma.vect());
      hst[139]->Fill( RtoD(alpha) );

      if ( ppg.dec_eta == 0 ) {
         ppg.dec_eta = 1;  // eta -> 2 gamma decay
         for( const auto& Vp : gammas_from_eta ) {
            HepLorentzVector LV( Vp, Vp.mag() );
            hst[140]->Fill(LV.e());

            LV.boost(-beta_jpsi);
            hst[141]->Fill(LV.e());
         }
         if ( gammas_from_eta.size() == 2 ) {
            double ang = gammas_from_eta[0].angle(gammas_from_eta[1]);
            hst[142]->Fill(RtoD(ang));
         }
      } else {
         ppg.dec_eta = 0;
      }
   }
}

// {{{1 MatchGammaMC()
// find the reconstructed gamma that best matches the gamma
// in decay of J/Psi -> gamma eta (decJpsi=22)
//--------------------------------------------------------------------
static void MatchGammaMC(PimPipGammas& ppg) {
//--------------------------------------------------------------------
   if ( !isMC ) {
      return;
   }
   if ( ppg.decJpsi != 22 ) {
      return;
   }
   const static double maxdp = 0.05;   // 50MeV

   // MC parameters of gamma (see FillHistoMC())
   Hep3Vector mcg = ppg.mcLVg.vect();
   double mcEg = ppg.mcLVg.e();

   double min_de = 2*maxdp; // may be asymmetric
   int Ng = ppg.LVg.size();
   for ( int i = 0; i < Ng; ++i ) {
      Hep3Vector gi = ppg.LVg[i].vect();
      double Egi = ppg.LVg[i].e();
      hst[161]->Fill( gi[0]-mcg[0] );
      hst[162]->Fill( gi[1]-mcg[1] );
      hst[163]->Fill( gi[2]-mcg[2] );
      if (     fabs(gi[0]-mcg[0]) < maxdp
            && fabs(gi[1]-mcg[1]) < maxdp
            && fabs(gi[2]-mcg[2]) < maxdp ) {
         hst[164]->Fill( mcg.angle(gi) );
         if ( fabs(Egi - mcEg) < fabs(min_de) ) {
            min_de = Egi - mcEg;
            ppg.mcidx_g = i;
         }
      }
   }

   if ( ppg.mcidx_g >= 0 ) { // matching for gamma
      // hst[165]->Fill( min_de );

      // difference in momentum of gamma: true and rec
      int ig = ppg.mcidx_g;
      double dE = ppg.mcLVg.e() - ppg.LVg[ig].e();
      double dA = ppg.mcLVg.angle(ppg.LVg[ig].vect());
      hst[166]->Fill(dE);
      hst[167]->Fill(dA);
   } else {
      // save E and cos(theta) no matching gammas
      hst[168]->Fill(mcEg);
      hst[169]->Fill(ppg.mcLVg.cosTheta());
   }
}

// {{{1 Charged tracks
//    - two soft pions with recoil mass of M(J/psi)
//    - no other charged tracks
//--------------------------------------------------------------------
static bool ChargedTracks(ReadDst* selector, PimPipGammas& ppg) {
//--------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
   // static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   int Ncharged = 0; // counter for all good charged tracks

   // vectors for soft pions
   vector<RecMdcKalTrack*> trk_p; // plus
   vector<RecMdcKalTrack*> trk_m; // minus

   for ( int i = 0; i < evtRecEvent->totalCharged(); i++ ) {

      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }
      if( !itTrk->isMdcKalTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();
      if ( mdcTrk->stat() == -222 ) { // skip cloned track
         continue;
      }

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC rec.
      HepPoint3D IP(ppg.xorig[0],ppg.xorig[1],ppg.xorig[2]);
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

      if( fabs(Rvxy0) >= Rvxy0_max ) {
         continue;
      }
      if( fabs(Rvz0) >= Rvz0_max ) {
         continue;
      }
      if ( fabs(cosTheta) >= cosTheta_max ) {
         continue;
      }

      Ncharged++;

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

      pid->calculate(ppg.runNo);
      int check_pion = 0;
      if ( pid->IsPidInfoValid() ) {
         // check that track is pion:
         if( pid->probPion() > 0 ) {
            hst[15]->Fill( log10(pid->probPion()) );
         } else {
            hst[15]->Fill(-5.);
         }

         if ( (pid->probPion() > 0.001)             &&
               (pid->probPion() > pid->probKaon())   &&
               (pid->probPion() > pid->probProton())
            ) {
            check_pion = 1;
         }
      } // end IsPidInfoValid
      hst[16]->Fill( double(check_pion) );
      if( check_pion == 0 ) {
         continue;
      }

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      if ( !mdcKalTrk ) {
         continue;
      }
      if ( std::isnan(mdcKalTrk->px())
            || std::isnan(mdcKalTrk->py())
            || std::isnan(mdcKalTrk->pz()) ) {
         continue;
      }
      mdcKalTrk->setPidType(RecMdcKalTrack::pion);

      hst[18]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      hst[19]->Fill( cos(mdcKalTrk->theta()) );
      if( mdcKalTrk->p() > 0.45 ) {
         continue;
      }

      if( mdcKalTrk->charge() > 0 ) {
         trk_p.push_back(mdcKalTrk);
      } else {
         trk_m.push_back(mdcKalTrk);
      }
   } // end for charged tracks loop

   hst[20]->Fill( double(Ncharged) );
   int Ntrkp = trk_p.size();
   int Ntrkm = trk_m.size();
   hst[21]->Fill(Ntrkp,Ntrkm);
   if ( Ncharged !=2 ||
        Ntrkp != 1 || Ntrkm != 1 ) {
      return false;
   }
   hst[1]->Fill(1); // "cuts"

   // calculate recoil mass of pi+ pi-
   auto& tp = trk_p[0];
   Hep3Vector Vp(tp->px(), tp->py(), tp->pz());
   HepLorentzVector LVp( Vp, sqrt( Vp.mag2() + SQ(mpi) ) );
   auto& tm = trk_m[0];
   Hep3Vector Vm(tm->px(), tm->py(), tm->pz());
   HepLorentzVector LVm( Vm, sqrt( Vm.mag2() + SQ(mpi) ) );

   // Mrec(pi+ pi-)
   ppg.LVjpsi = ppg.LVcms - LVp - LVm;
   double Mrec2 = ppg.LVjpsi.m2();
   if( Mrec2 < 0 ) {
      return false;
   }
   double Mrec = sqrt(Mrec2);
   hst[22]->Fill( Mrec );
   if( Mrec <= 3.092 || Mrec >= 3.102 ) { // narrow CUT on M(J/Psi)
      return false;
   }

   // Theta(pi+ pi-)
   double cosPM = Vm.cosTheta(Vp);
   hst[23]->Fill( cosPM );

   // Minv( pi+ pi-)
   double invPM = (LVp+LVm).m();
   hst[24]->Fill( invPM );

   // check Mrec candidate
   if ( !SelectPM( cosPM, invPM ) ) {
      return false;
   }
   hst[1]->Fill(2); // "cuts"
   if ( isMC ) {
      hst[101]->Fill(ppg.decPsip);
      hst[106]->Fill(ppg.decJpsi);
   }

   // save selected pi+ pi-
   ppg.trk_Pip = tp;
   ppg.trk_Pim = tm;
   ppg.Mrec = Mrec;
   ppg.cosPM = cosPM;
   ppg.invPM = invPM;

   hst[31]->Fill( Mrec );
   hst[32]->Fill( cosPM );
   hst[33]->Fill( invPM );

   hst[41]->Fill(Vp.mag());
   hst[42]->Fill(Vp.perp());
   hst[43]->Fill(Vp.cosTheta());
   hst[44]->Fill(Vm.mag());
   hst[45]->Fill(Vm.perp());
   hst[46]->Fill(Vm.cosTheta());

   if ( isMC ) {
      if ( ppg.decPsip == 64 ) {
         // diff. J/Psi momentum true and rec.
         double dP = ppg.mcLVjpsi.vect().mag()
            - ppg.LVjpsi.vect().mag();
         double dA = ppg.mcLVjpsi.angle(ppg.LVjpsi.vect());
         hst[151]->Fill(dP);
         hst[152]->Fill(dA);
      }
   }

   return true;
}

// {{{1 select neutral tracks
//    - collect good photons
//    - search for photons from J/Psi->gamma X decay with recoil mass
//      close to mass of eta-particle
//    - remove events with gammas from pi0 decays
//--------------------------------------------------------------------
static bool NeutralTracks(ReadDst* selector, PimPipGammas& ppg) {
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
      // map of gammas
      // if ( absCosTheta < 0.8 ) {  //barrel
         // hst[55]->Fill(emcpos.z(), RtoD(emcpos.phi()) );
      // } else { // endcap
         // int hz = (emcpos.z() > 0);
         // hst[56+hz]->Fill(emcpos.rho(), RtoD(emcpos.phi()) );
      // }

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
      hst[51]->Fill(RtoD(tang));
      if ( tang < min_angle ) {
         continue;
      }

      ppg.trk_g.push_back( emcTrk );
      Hep3Vector p3 = emcpos - ppg.xorig;
      p3 *= eraw / p3.mag();
      ppg.LVg.push_back( HepLorentzVector(p3,eraw) );

   } // end of neutrals loop

   int Ng = ppg.trk_g.size(); // number of good photons
   hst[52]->Fill(Ng);

   if ( Ng <= 1 ) {
      return false;
   }
   hst[1]->Fill(3); // "cuts"
   for ( int i = 0; i < Ng; ++i ) {
      hst[53]->Fill(ppg.LVg[i].e());
   }
   if ( isMC ) {
      hst[102]->Fill(ppg.decPsip);
      hst[107]->Fill(ppg.decJpsi);
   }

   // find the reconstructed gamma that best suts the
   // gamma in decay of J/Psi -> gamma eta (decJpsi=22)
   MatchGammaMC(ppg);

   // =========================================================
   // Recoil mass square of gamma in J/Psi->gamma+X decay
   // should be in region: Meta^2 +/- 0.2
   // this is corresponds E_gamma(in J/Psi rest system) [1.45,1.55]GeV
   // M^2rec = M^2(J/Psi)-2*M(J/Psi)*E_gamma(in J/Psi RS)

   Hep3Vector beta_jpsi = ppg.LVjpsi.boostVector();
   for ( int i = 0; i < Ng; ++i ) {
      double Mrg2 = (ppg.LVjpsi - ppg.LVg[i]).m2();
      hst[61]->Fill(Mrg2);
      if ( isMC ) {
         if ( ppg.decJpsi == 22 ) {
            hst[171]->Fill(Mrg2);
            if ( i == ppg.mcidx_g ) {
               hst[173]->Fill(Mrg2);
            }
            if ( ppg.dec_eta == 1 ) {
               hst[174]->Fill(Mrg2);
            }
         } else {
            hst[172]->Fill(Mrg2);
         }
      }

      HepLorentzVector gj(ppg.LVg[i]);
      gj.boost(-beta_jpsi);
      double Egj = gj.e();
      hst[62]->Fill(Egj);
      if ( isMC ) {
         if ( ppg.decJpsi == 22 ) {
            hst[176]->Fill(Egj);
            if ( i == ppg.mcidx_g ) {
               hst[178]->Fill(Egj);
            } else {
               hst[179]->Fill(Egj);
            }
         } else {
            hst[177]->Fill(Egj);
         }
      }

      if ( fabs(Mrg2-SQ(meta)) < 0.2 ) { // CUT for separated photon
         ppg.idx_g.push_back(i);
         ppg.M2rec_g.push_back(Mrg2); // to save in ntuple
         if ( isMC ) {
            if ( ppg.decJpsi == 22 ) {
               hst[181]->Fill(Egj);
            } else {
               hst[182]->Fill(Egj);
            }
         }
      }

   } // end of gammas loop
   int Ngj = ppg.idx_g.size();
   hst[63]->Fill(Ngj);
   if ( isMC ) {
      if ( ppg.decJpsi == 22 ) {
         hst[185]->Fill(Ngj);
         if ( ppg.mcidx_g >= 0 ) {
            hst[187]->Fill(Ngj);
         }
      } else {
         hst[186]->Fill(Ngj);
      }
   }

   if ( Ngj == 0 ) {
      return false;
   }
   hst[1]->Fill(4); // "cuts"
   if ( isMC ) {
      hst[103]->Fill(ppg.decPsip);
      hst[108]->Fill(ppg.decJpsi);
   }

   // =========================================================
   // analyse gg pairs: reject pi0->2gamma
   int Npi0 = 0;
   for ( int i = 0; i < Ng-1; ++i ) {
      auto& LVgi = ppg.LVg[i];
      for ( int j = i+1; j < Ng; ++j ) {
         auto& LVgj = ppg.LVg[j];

         double Mgg2 = (LVgi+LVgj).m2();
         hst[71]->Fill(Mgg2);
         hst[72]->Fill(Mgg2);

         if ( isMC ) {
            if ( ppg.decJpsi == 22 ) {
               hst[191]->Fill(Mgg2);
               if ( ppg.dec_eta == 1 ) {
                  hst[193]->Fill(Mgg2);
               }
            } else {
               hst[192]->Fill(Mgg2);
            }
         }

         if ( 0.013 < Mgg2 && Mgg2 < 0.022 ) { // CUT pi0
            Npi0 += 1;
         }
      } // end of for(j)
   } // end of for(i)
   hst[73]->Fill(Npi0);

   if ( isMC ) {
      if ( Npi0 > 0 ) {
         hst[195]->Fill(ppg.decJpsi);
      }
      if ( ppg.decJpsi == 22 ) {
         hst[196]->Fill(Npi0);
         if ( ppg.dec_eta == 1 ) {
            hst[198]->Fill(Npi0);
         }
      } else {
         hst[197]->Fill(Npi0);
      }
   }

   if ( Npi0 > 0 ) {
      return false;
   }
   hst[1]->Fill(5); // "cuts"
   if ( isMC ) {
      hst[104]->Fill(ppg.decPsip);
      hst[109]->Fill(ppg.decJpsi);
   }

   return true;
}

// {{{1 Eta reconstruction efficiency
//--------------------------------------------------------------------
static void EtaEff( ReadDst* selector, PimPipGammas& ppg ) {
//--------------------------------------------------------------------
   // window for selection of eta: see cuts.h (in PsipJpsiPhiEta)
   static const double seta = 0.008;
   static const double weta = 3*seta; // standard

   // const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   int Ng   = ppg.LVg.size();

   // search for << J/Psi -> gamma+eta  >> such that the
   // recoil mass of full system is the smallest
   double M2fr_min = 10;
   int k_min = -1;
   int i_min = -1;
   vector<HepLorentzVector> LVg_min(4); // save gammas here
   double M2gr_min = 0; // recoil M^2 of gamma in J/Psi -> gamma eta
   double M2gg_min = 0;

   for ( size_t ig0 = 0; ig0 < ppg.idx_g.size(); ++ig0 ) {
      int k = ppg.idx_g[ig0]; // candidate to "separated" gamma
      const auto& LVg0 = ppg.LVg[k];

      for ( int ig = 0; ig < Ng; ++ig ) {
         if ( ig == k ) {
            continue;
         }
         const auto& LVg1 = ppg.LVg[ig];

         // calculate recoil moment
         HepLorentzVector LVrec = ppg.LVjpsi - LVg0 - LVg1;
         // it must be gamma (from eta decay)
         HepLorentzVector LVg2;
         LVg2.setVectM(LVrec.vect(), 0.); // set zero mass

         // check that these gammas satisfy the decay of eta
         HepLorentzVector LVgg = LVg1 + LVg2;
         double Mgg2 = LVgg.m2();
         hst[81]->Fill(Mgg2);
         if ( isMC ) {
            if ( ppg.decJpsi == 22 ) {
               hst[201]->Fill(Mgg2);
               if ( k == ppg.mcidx_g ) {
                  hst[202]->Fill(Mgg2);
               }
               if ( ppg.dec_eta == 1 ) {
                  hst[204]->Fill(Mgg2);
               }
            } else {
               hst[203]->Fill(Mgg2);
            }
         }

         // Mgg^2 ~ 0.295 +/- 0.021
         if ( fabs(Mgg2-0.295) < 0.045 ) { // Select eta CUT
            // total recoil mass
            HepLorentzVector LVfr = LVrec - LVg2;
            double M2fr = LVfr.m2();
            if ( fabs(M2fr) < fabs(M2fr_min) ) {
               M2fr_min = M2fr;
               k_min = k;
               i_min = ig;
               LVg_min[0] = LVg0;
               LVg_min[1] = LVg1;
               LVg_min[2] = LVg2;
               M2gr_min = ppg.M2rec_g[ig0];
               M2gg_min = Mgg2;
            }
         }
      } // end of gammas loop
   } // end of "separated" gamma loop

   hst[82]->Fill(M2fr_min);
   if ( isMC ) {
      if ( ppg.decJpsi == 22 ) {
         hst[205]->Fill(M2fr_min);
         if ( k_min == ppg.mcidx_g ) {
            hst[206]->Fill(M2fr_min);
         }
         if ( ppg.dec_eta == 1 ) {
            hst[208]->Fill(M2fr_min);
         }
      } else {
         hst[207]->Fill(M2fr_min);
      }
   }

   if ( fabs(M2fr_min) > 0.01 ) { // CUT minimal total recoil mass
      return;
   }

   // search near the predicted photon for a real photon
   // OLD: with minimal recoil mass of full system
   // NEW: with a minimum angle to the predicted one

   const auto& LVgp = LVg_min[2]; // predicted photon

   // HepLorentzVector LVrec = ppg.LVjpsi - LVg_min[0] - LVg_min[1];
   double M2real_min = 10; // old
   int j_min = -1;
   double rE_min = 0, dTh_min = 100.;
   for ( int jg = 0; jg < Ng; ++jg ) {
      if ( jg == k_min || jg == i_min ) {
         continue;
      }
      const auto& LVgj = ppg.LVg[jg];

      double rE = LVgp.e()/LVgj.e();
      double cosTh = LVgp.vect().cosTheta( LVgj.vect() );
      double dTh = RtoD( acos(cosTh) ); // [0,180]
      hst[85]->Fill(rE);
      hst[86]->Fill(dTh);
      if ( isMC ) {
         if ( ppg.decJpsi == 22 ) {
            hst[211]->Fill(rE);
            hst[212]->Fill(dTh);
         } else {
            hst[213]->Fill(rE);
            hst[214]->Fill(dTh);
         }
      }

      // if ( 0.4 < rE && rE < 2. && dTh < 20 ) { // CUT
         // HepLorentzVector LVfr = LVrec - LVgj;
         // double M2fr = LVfr.m2();
         // if ( fabs(M2fr) < fabs(M2real_min) ) {

      if ( dTh < dTh_min ) {
         // M2real_min = M2fr;
         j_min = jg;
         rE_min = rE;
         dTh_min = dTh;
      }
   } // end of for(jg)
   int found = int( j_min != -1 );

   double chisq = 9999.;
   double Mgg_found = 0.;
   if ( found ) {
      // save found real gamma in LVg_min
      const auto& LVgj = ppg.LVg[j_min];
      LVg_min[3] = LVgj;

      hst[87]->Fill(dTh_min,rE_min);
      hst[91]->Fill(M2real_min);
      if ( isMC ) {
         if ( ppg.decJpsi == 22 ) {
            hst[215]->Fill(dTh_min,rE_min);
         } else {
            hst[216]->Fill(dTh_min,rE_min);
         }

         if ( ppg.decJpsi == 22 ) {
            hst[221]->Fill(M2real_min);
            if ( k_min == ppg.mcidx_g ) {
               hst[222]->Fill(M2real_min);
            }
         } else {
            hst[223]->Fill(M2real_min);
         }
      }

      // check eta-mass
      HepLorentzVector LVgg = LVg_min[1] + LVgj;
      double Mgg2 = LVgg.m2();
      hst[92]->Fill(Mgg2);
      if ( isMC ) {
         if ( ppg.decJpsi == 22 ) {
            hst[224]->Fill(Mgg2);
            if ( k_min == ppg.mcidx_g ) {
               hst[225]->Fill(Mgg2);
            }
         } else {
            hst[226]->Fill(Mgg2);
         }
      }

      // here we apply 4C kinamatic constraints and select "eta"
      // as we do in the main analysis

      // Vertex fit
      WTrackParameter wp[2] = {
         WTrackParameter(mpi,ppg.trk_Pip->getZHelix(),
                         ppg.trk_Pip->getZError() ),
         WTrackParameter(mpi,ppg.trk_Pim->getZHelix(),
                         ppg.trk_Pim->getZError() ),
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
      if ( ok ) {
         hst[95]->Fill(vtxfit->chisq(0));

         // Kinematik fit
         vtxfit->BuildVirtualParticle(0);
         // translate parameters of tracks to the vertex# 0
         vtxfit->Swim(0);
         for(int k = 0; k < 2; k++) {
            wp[k] = vtxfit->wtrk(k);
         }
         KinematicFit* kmfit = KinematicFit::instance();
         kmfit->init();
         for ( int k = 0; k < 2; k++ ) {
            kmfit->AddTrack(k, wp[k]);
         }
         kmfit->AddTrack(2, 0.0, ppg.trk_g[k_min]);
         kmfit->AddTrack(3, 0.0, ppg.trk_g[i_min]);
         kmfit->AddTrack(4, 0.0, ppg.trk_g[j_min]);
         kmfit->AddFourMomentum(0, ppg.LVcms);
         bool oksq = kmfit->Fit();
         if ( oksq ) {
            chisq = kmfit->chisq();
            hst[96]->Fill(chisq);
            // Momentums after kinematic corrections:
            HepLorentzVector Pg1 = kmfit->pfit(3);
            HepLorentzVector Pg2 = kmfit->pfit(4);
            HepLorentzVector Pgg = Pg1 + Pg2;
            Mgg_found  = Pgg.m();

            hst[97]->Fill(Mgg_found);
            if ( isMC ) {
               if ( ppg.decJpsi == 22 ) {
                  hst[228]->Fill(Mgg_found);
               } else {
                  hst[229]->Fill(Mgg_found);
               }
            }

            // set found flag to 2 if we found good eta
            if ( fabs(Mgg_found-meta) < weta ) { // see cuts.h
               found = 2;
            }
         } // end kinematic fit
      } // end vertex fit
   } // end if( found )

   // OLD:
   // double Mgg_pred = LVeta.m();
   // double(found),
   // M2real_min,

   // fill and save Ttree
   xnt_geta.Ptp  = ppg.trk_Pip->pxy();
   xnt_geta.Ptm  = ppg.trk_Pim->pxy();
   xnt_geta.Mrec = ppg.Mrec;
   xnt_geta.Eg0  = LVg_min[0].e();
   xnt_geta.Cg0  = LVg_min[0].cosTheta();
   HepLorentzVector LVeta = LVg_min[1] + LVg_min[2];
   xnt_geta.Peta = LVeta.vect().mag();
   xnt_geta.Ceta = LVeta.cosTheta();
   xnt_geta.Eg1  = LVg_min[1].e();
   xnt_geta.Cg1  = LVg_min[1].cosTheta();
   xnt_geta.Eg2  = LVg_min[2].e();
   xnt_geta.Cg2  = LVg_min[2].cosTheta();
   xnt_geta.Egr  = LVg_min[3].e();
   xnt_geta.Cgr  = LVg_min[3].cosTheta();
   xnt_geta.rE   = rE_min;
   xnt_geta.dTh  = dTh_min;

   xnt_geta.M2gr = M2gr_min;
   xnt_geta.M2gg = M2gg_min;
   xnt_geta.m2fr = M2fr_min;
   xnt_geta.mggf = Mgg_found;
   xnt_geta.ch2  = chisq;

   if ( isMC ) {
      xnt_geta.decj = ppg.decJpsi;
      xnt_geta.dgam = int(k_min == ppg.mcidx_g);
      xnt_geta.deta = ppg.dec_eta;
   } else {
      xnt_geta.decj = 0;
      xnt_geta.dgam = 0;
      xnt_geta.deta = 0;
   }

   m_nt_geta->Fill();

   if ( isMC ) {
      // more or less final cuts:
      if (
            abs(xnt_geta.Cg0) < 0.8 &&
            abs(xnt_geta.Cg1) < 0.8 &&
            xnt_geta.m2fr < 0.002 &&
            xnt_geta.Eg1 > 0.25 && xnt_geta.Eg1 < 1.7
         )
      {
         JpsiTbl.Add(); // save in table
      }
   }

}

// {{{1 MAIN: Loop for each event
//--------------------------------------------------------------------
bool PsipJpsiGammaEtaEvent( ReadDst*       selector,
                            TEvtHeader*    m_TEvtHeader,
                            TDstEvent*     m_TDstEvent,
                            TEvtRecObject* m_TEvtRecObject,
                            TMcEvent*      m_TMcEvent,
                            TTrigEvent*    m_TTrigEvent,
                            TDigiEvent*    m_TDigiEvent,
                            THltEvent*     m_THltEvent      ) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   PimPipGammas ppg; // information for the current event

   //-----------------------------------------------------------------
   //-- Get event information --
   //-----------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   int eventNo = m_TEvtHeader->getEventId();
   ppg.runNo = runNo;
   ppg.event  = eventNo;
   hst[1]->Fill(0); // "cuts"

   // define data taking period --------------------------------------
   // ATTENTION: this part is calculated only once: we get constants
   // that are the same in all events of one job
   if ( DataPeriod == 0 ) {
      int idx_ac = -1; // index for abscor
      int run = abs(runNo);
      if ( (run >= 8093 && run <= 9025)  ) {         // 2009 Psi(2S)
         DataPeriod = 2009;
         idx_ac = 0;
      } else if ( (run >= 25338 && run <= 27090) ) { // 2012 Psi(2S)
         DataPeriod = 2012;
         idx_ac = 1;
      } else if ( (run >= 66257 && run <= 69292) ) { // 2021 Psi(2S)
         DataPeriod = 2021;
         idx_ac = 2;
      } else if ( (run >= 9613 && run <= 9779) ) { // 3650 2009-data
         DataPeriod = 3650;
         idx_ac = 0;
      } else if ( (run >= 33725 && run <= 33772) ) { // 3650 2012-data
         DataPeriod = 3650;
         idx_ac = 1;
      } else if ( (run >= 69612 && run <= 70132) ) { // 3650 2021-data
         DataPeriod = 3650;
         idx_ac = 2;
      } else {
         cout << " FATAL: Data taking period undefined for runNo= "
              << runNo << endl;
         exit(EXIT_FAILURE);
      }
#if (BOSS_VER == 709)
      string dir("Analysis/AbsCor/dat/00-00-41/");
      vector<string> c3p {
         "c3ptof2009Jpsi.txt",
         "c3ptof2012Jpsi.txt",
         "c3ptof2021psip.txt"
      };
      vector<string> evsetToF {
         "evsetTofCorFunctionPar2009Jpsi.txt",
         "evsetTofCorFunctionPar2012Jpsi.txt",
         "evsetTofCorFunctionPar2021psip.txt"
      };
      string data_c3p = selector->AbsPath( dir+c3p[idx_ac] );
      string cor_evsetTof = selector->AbsPath( dir+evsetToF[idx_ac] );
      m_abscor->ReadDatac3p( data_c3p, true );
      m_abscor->ReadCorFunpara( cor_evsetTof, true );
#else
      (void)idx_ac; // suppres warning about unused var
#endif

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
         exit(EXIT_FAILURE);
      }
   } // end of DataPeriod definition ---------------------------------

   // call AbsCorr after reading c3p and cor_evsetTof files:
   // second call AbsCorr() is possible, the result does not change
   m_abscor->AbsorptionCorrection(selector);

   double Ecms = 3.686; // (GeV) energy in center of mass system
   if ( DataPeriod == 3650 ) {
      Ecms = 3.650;
   }
   ppg.LVcms = HepLorentzVector(Ecms*sin(beam_angle), 0, 0, Ecms);

   FillHistoMC(selector,ppg); // MC histo

   Hep3Vector xorigin = getVertexOrigin(runNo);
   ppg.xorig = xorigin;

   // Study of the photon reconstruction efficiency
   if ( !ChargedTracks(selector, ppg) ) {
      return false;
   }
   if ( !NeutralTracks(selector, ppg) ) {
      return false;
   }

   EtaEff(selector,ppg);

   return true;
}

// {{{1 EndJob
//--------------------------------------------------------------------
void PsipJpsiGammaEtaEndJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
   }

   if ( isMC ) {
      // print tables of decays
      cout << string(65,'#') << endl;
      cout << "Decays of J/psi in PsipJpsiGammaEta" << endl
           << "   search for Eta-Gamma: "
           << JpsiTbl.ntot << " events" << endl
           << "       size of table is " << JpsiTbl.Size() << endl;
      JpsiTbl.Print(0.1); // do not print decays with P<0.1% of all
      cout << "Enddecay" << endl << endl;
   }

   string module = string(__func__);
   module = module.substr(0,module.size()-6);
   int nw = warning_msg.size();
   string tw = (nw > 0) ? to_string(nw)+" WARNINGS" : "no warnings";
   printf(" There are %s in %s\n",tw.c_str(),module.c_str());
   for(auto it = warning_msg.begin(); it != warning_msg.end(); ++it) {
      cout << it->first << " : " << it->second << endl;
   }
}

//--------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
