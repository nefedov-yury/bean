//======================================================================//
//                                                                      //
// SelectOmegaEta - search for e+ e- -> pi+ pi- pi0 Eta                 //
//                                                   |-> 2 gammas       //
//                                                                      //
//======================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <vector>
// #include <map>

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/TwoVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using CLHEP::pi;
using CLHEP::twopi;


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

#include "DstEvtRecTracks.h"
#include "TofHitStatus.h"
#include "ReadDst.h"


//add ISR photon to 4C fit
//#define add_ISR OK


using namespace std;

struct Select_OmegaEta {
  double Ecms;                  // energy in center of mass system

  // charge tracks:
  vector<int> qtrks;            // trk# in EvtRecTrkCol
  int ipip, ipim;               // trk# for charge +/-
  HepLorentzVector ppip, ppim;  // 4-momentum of trks
  int ipip_mc, ipim_mc;         // correspondent mcParticle
  int pip_pdg, pim_pdg;         // PDG code of these mc-particles (or 0)

  HepLorentzVector P_missing;   // missing 4-momentum
  double SoverS0;               // for check MCGPJ

  // neutral (gammas) tracks
  vector<int> gammas;           // trk# in EvtRecTrkCol
  vector<HepLorentzVector> Pg;  // 4-momentum of gammas
  int ig1, ig2, ig3, ig4;       // the best gammas after 4C kin.fit
  HepLorentzVector ppi01, ppi02;// 4-momentum of pairs (g+g) ???

  HepLorentzVector Ptot;        // 4-momentum of the overall system

  // momentums after kinematic fit
  HepLorentzVector Ppip, Ppim, Pgg1,Pgg2, Pomega, Pgg;

  Select_OmegaEta() {
    Ecms = -1;
    ipip = ipim = -1;
    ipip_mc = ipim_mc = -1;
    pip_pdg = pim_pdg = 0;
    ig1 = ig2 = ig3 = ig4 = -1;
    SoverS0 = -1;
  }

  int nGood() { return qtrks.size(); }
  int mc_trk(int i) { return (qtrks[i] == ipip) ? ipip_mc : ipim_mc; }
  int isPiPlus(int i) { return qtrks[i] == ipip; }
  int isPiMinus(int i) { return qtrks[i] == ipim; }
  int nGam() { return gammas.size(); }
  int gamma1() { return (ig1 >= 0) ? gammas[ig1] : -1; }
  int gamma2() { return (ig2 >= 0) ? gammas[ig2] : -1; }
  int gamma3() { return (ig3 >= 0) ? gammas[ig3] : -1; }
  int gamma4() { return (ig4 >= 0) ? gammas[ig4] : -1; }
  bool GoodGammas()
        {return (ig1 != -1 && ig2 != -1 && ig3 != -1 && ig4 != -1); }
};

typedef Select_OmegaEta Select;

//-----------------------------------------------------------------------------
// Global file variables
//-----------------------------------------------------------------------------
const static bool init_selection=true;

const static double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
const static double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
const static double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
const static double meta   = 0.547862; // 547.862   +/- 0.017   MeV
const static double momega = 0.78265;  // 782.65    +/- 0.12    MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;
static std::vector<TNtuple* > m_tuple;

static int ERROR_WARNING = 0;
static bool isMC = false;

// static bool omega_to_pis=false;
static bool omega_best_mass=false;
static bool one_ISR=false;
// static bool one_ISR_less50=false;
// static bool one_ISR_more50=false;
static bool two_ISR=false;
// static bool two_ISR_less50=false;
// static bool two_ISR_more50=false;
static bool SoverS0 = false;
static bool EtaTo2Gammas = false;

// Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
inline Double_t RtoD(Double_t ang){return ang*180/M_PI;}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SelectOmegaEtaStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << " Start: " << __func__ << "()" << endl;
  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }

  his1.resize(500,(TH1D*)0);
  his2.resize(500,(TH2D*)0);
  m_tuple.resize(10,0);

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

  // We have to initialize DatabaseSvc ----------------------------------------
  DatabaseSvc* dbs = DatabaseSvc::instance();
  if ( (dbs->GetDBFilePath()).empty() ) {
    // set path to directory with databases:
    dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
  }

  // We have to initialize Magnetic field--------------------------------------
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

  //--------- Book histograms--------------------------------------------------

  his1[0] = new TH1D("Runs", "run numbers", 1000,0.,50000);
  his1[1] = new TH1D("cuts_0","preselections cuts", 20,-0.5,19.5);
  his1[9] = new TH1D("Ecms","E_{cms}", 1500, 3.0, 3.15);

  his1[11] = new TH1D("Rxy","R_{xy}", 100,-5.,5.);
  his1[12] = new TH1D("Rz","R_{z}", 100,-50.,50.);
  his1[13] = new TH1D("cos_theta","cos(#theta)", 100,0.,1.);
  his1[14] = new TH1D("QPtrack","Charged momentum (Q*P)", 600,-3.,3.);
  his1[15] = new TH1D("Ntrks","N_{charged tracks} in event", 11,-0.5,10.5);
  his1[16] = new TH1D("Ncharge","N_{charge} in event", 11,-5.5,5.5);

  his1[21] = new TH1D("PidInfo", "0 - no pid, 1 - pid", 2,-0.5,1.5);
  his1[22] = new TH1D("LogCLpi","lg(CL_{#pi})", 100,-5.,0.);
  his1[23] = new TH1D("LogCLK","lg(CL_{K})", 100,-5.,0.);
  his1[24] = new TH1D("LogCLp","lg(CL_{p})", 100,-5.,0.);
  his1[25] = new TH1D("KalTrk", "0,1 - no/yes mdcKalTrk", 2,-0.5,1.5);
  his2[21] = new TH2D("SignTrks","sign of trks (+) vs (-)",
                3,-0.5,2.5, 3,-0.5,2.5);

  his1[31] = new TH1D("AngleTrack","The angle between the photon"
                " and the nearest charged track, degrees", 180,0.,180.);
  his1[32] = new TH1D("Ngamma","N_{gamma} in event", 11,-0.5,10.5);
  his1[33] = new TH1D("ini_Ppi","Momentum of selected pions (Q*P)",
                600,-3.,3.);
  his1[34] = new TH1D("ini_Pgam","Momentum of selected gammas", 300,0.,3.);

/*

  his1[11] = new TH1D("E_p","#frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);
  his1[12] = new TH1D("E_p_new","NEW  #frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);
  his1[13] = new TH1D("E_p_fin","FINAL  #frac{Energy_{EMC}}{p_{MDC}}", 100, 0, 1.1);

  his1[14] = new TH1D("vtx_fit", "0 - no fit, 1 - good fit", 2, -0.5, 1.5);
  his1[15] = new TH1D("vtx_chi2", "vertex fit #chi^2", 100, 0., 100.);
  his1[16] = new TH1D("Mgg_omega_5C_Mggcut", "Final  M_{#gamma#gamma} omega", 400, 0.06, 0.2);
  his1[180] = new TH1D("Mgg_pi0_5C_Mggcut", "Final  M_{#gamma#gamma} pi0", 400, 0.06, 0.2);
  his1[181] = new TH1D("Momega_5C_Mggcut", "Final  M_{omega}", 400, 0.55, 0.95);
  his1[182] = new TH1D("chi2_5C", "min #chi^2   5C", 100, 0., 200.);
  his1[183] = new TH1D("Mgg_omega_in5C", "5C-fit: M_{#gamma#gamma} without min chi2", 400, 0.06, 0.2);
  his2[184] = new TH2D("chi2_5C_all", "all #chi^2   5C vs. g1g2g3g4-combinations", 200, 0., 200., 6,-0.5,5.5);
  his1[185] = new TH1D("Mgg_omega_5C", "5C-fit: M_{#gamma#gamma} omega", 400, 0.06, 0.2);
  his1[186] = new TH1D("Mgg_pi0_5C", "5C-fit: M_{#gamma#gamma} pi0", 400, 0.06, 0.2);
  his2[180] = new TH2D("Momega_Meta_c", "Final M_{#omega} vs M_{#eta} c", 400, 0.4, 3., 125, 0.3, 0.8);
  his2[181] = new TH2D("Momega_Meta_xy", "Final M_{#omega} vs M_{#eta} xy", 400, 0.4, 3., 125, 0.3, 0.8);
  his2[182] = new TH2D("Momega_Meta_d", "Final M_{#omega} vs M_{#eta} d", 400, 0.4, 3., 125, 0.3, 0.8);
  his2[183] = new TH2D("Momega_Meta_all", "Final M_{#omega} vs M_{#eta}", 400, 0.4, 3., 125, 0.3, 0.8);
  his2[185] = new TH2D("Momega_Meta", "Final M_{#omega} vs M_{#eta}", 400, 0.4, 3., 125, 0.3, 0.8);


  his1[187] = new TH1D("chi2_5C_69", "min #chi^2   5C  decCode=69", 100, 0., 200.);
  his1[188] = new TH1D("chi2_5C_not69", "min #chi^2   5C  decCode!=69", 100, 0., 200.);
  his1[189] = new TH1D("chi2_4C_69", "min #chi^{2} for 4C-Fit   decCode=69", 200, 0, 200. );
  his1[190] = new TH1D("chi2_4C_not69", "min #chi^{2} for 4C-Fit   decCode!=69", 200, 0, 200. );
  his1[191] = new TH1D("Momega_fin_not69", "Final  M_{omega}  decCode!=69", 400, 0.55, 0.95);
  his1[192] = new TH1D("chi2_4C_SBc", "min #chi^{2}  4C   SB-center", 200, 0, 200. );
  his1[193] = new TH1D("chi2_5C_SBc", "min #chi^{2}   5C  SB-center", 200, 0., 200.);
  his1[194] = new TH1D("chi2_4C_SBxy", "min #chi^{2}  4C  SB-xy", 200, 0, 200. );
  his1[195] = new TH1D("chi2_5C_SBxy", "min #chi^{2}  5C  SB-xy", 200, 0., 200.);
  his1[196] = new TH1D("chi2_4C_SBd", "min #chi^{2}  4C  SB-diag.", 200, 0, 200. );
  his1[197] = new TH1D("chi2_5C_SBd", "min #chi^{2}  5C  SB-diag.", 200, 0., 200.);

  his1[17] = new TH1D("Etot_4C","4C-fit: E_{tot}", 100, 2.5, 3.5);
  his1[18] = new TH1D("Ptot_4C","4C-fit: P_{tot}", 200, 0, 1.0);


  his2[1] = new TH2D("depth", "Depth_{-} vs Depth_{+}", 100, -10., 60., 100, -10., 60.);

  his1[31] = new TH1D("chi2_4C", "min #chi^{2} for 4C-Fit", 300, 0, 300 );
  his1[32] = new TH1D("Etot","E_{3#pi} rec.", 100, 0., 3.5);
  his1[33] = new TH1D("Ptot","P_{3#pi} rec.", 100, 0., 2.);

  his1[40] = new TH1D("Rxy_MCpart","R_{xy} MCpart for omega->3pi", 100, 0., 20.);
  his1[41] = new TH1D("Rz_MCpart","R_{z} MCpart for omega->3pi", 100, -50., 50.);

  his2[51] = new TH2D("Mpm_Mp0", "M^{2}(#pi^{+}#pi^{-}) vs M^{2}(#pi^{+}#pi^{0})", 100,0.,9., 100,0.,9.);
  his2[52] = new TH2D("Mpm_Mm0", "M^{2}(#pi^{+}#pi^{-}) vs M^{2}(#pi^{-}#pi^{0})", 100,0.,9., 100,0.,9.);
  his2[53] = new TH2D("Mp0_Mm0", "M^{2}(#pi^{+}#pi^{0}) vs M^{2}(#pi^{-}#pi^{0})", 100,0.,9., 100,0.,9.);


  // Monte Carlo histograms
  his1[101] = new TH1D("mc_type", "0 - bkg, 1 - 3pi, 2 - pi+rho", 3, -0.5, 2.5);
  his1[102] = new TH1D("mc_deccode", "decCode",258,-1.5,256.5);
  his1[103] = new TH1D("mc_deccode_final", "decCode after select",258,-1.5,256.5);
  his1[106] = new TH1D("mc_deccode_P", "decCode for Ptot_rec > 0.2",258,-1.5,256.5);


  his1[107] = new TH1D("cuts","cuts: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-3#sigma",9,-0.5,8.5);
  his1[108] = new TH1D("fEtaToGammas","Final #eta -> #gamma#gamma: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-3#sigma",9,-0.5,8.5);
  his1[280] = new TH1D("cuts_omega_best_mass","cuts for omega best mass (MCGPJ): 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-3#sigma",9,-0.5,8.5);
  his1[209] = new TH1D("cuts_67","cuts deccode=67: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-3#sigma",9,-0.5,8.5);
  his1[300] = new TH1D("cuts_1ISR","cuts: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-5Cfit, 6-M_{gg}/c, 7-M-missing",9,-0.5,8.5);
  his1[301] = new TH1D("cuts_2ISR","cuts: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-5Cfit, 6-M_{gg}/c, 7-M-missing",9,-0.5,8.5);
  his1[302] = new TH1D("cuts_1ISR_bestmass","cuts: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-5Cfit, 6-M_{gg}/c, 7-M-missing",9,-0.5,8.5);
  his1[303] = new TH1D("cuts_2ISR_bestmass","cuts: 0-all, 1-good charged, 2-good gamma, 3-4Cfit, 4-chi2, 5-5Cfit, 6-M_{gg}/c, 7-M-missing",9,-0.5,8.5);
  //ISR
  his1[111] = new TH1D("mc_Eg", "energy of one gamma ISR (MeV)",1000,0.,1000.);
  his1[112] = new TH1D("mc_Egs", "E(ISR) (MeV)",1000,0.,3000.);
  his1[115] = new TH1D("mc_Egmax", "E max gamma ISR (MeV)",1000,0.,1000.);
  his1[116] = new TH1D("mc_Ng", "N(gammas ISR)",10,-0.5,9.5);
  his1[117] = new TH1D("mc_Ng1","E(ISR) one gamma",1000,0,1000.);
  his1[118] = new TH1D("mc_Ng2","E(ISR) two gamma",1000,0,1000.);
  his1[119] = new TH1D("mc_Ng2m","E(ISR max) two gamma",1000,0,1000.);
  his2[252] = new TH2D("mc_Eisr1_Eisr2","Eisr 1 vs Eisr 2", 1500, 0., 1500., 1500, 0., 1500.);
  //new
  his1[113] = new TH1D("mc_E4pi", "E(#omega #pi^{0}) (MeV)",1000,2200.,3200.);
  his1[114] = new TH1D("mc_Etot", "E(#omega #pi^{0})+E(ISR) (MeV)", 2000,3000.,3200.);


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

  //new
  his1[151] = new TH1D("mc_ini_pdg", "PDG codes of particles",2001,-1000.5,1000.5);
  his1[152] = new TH1D("mc_ini_pdg0", "PDG of particles from primary vertex",2001,-1000.5,1000.5);
  his1[153] = new TH1D("mc_ini_pdgEta", "PDG codes of particles from #eta-decay",2001,-1000.5,1000.5);
  his1[155] = new TH1D("mc_ini_Ppi0", "Momentum of #pi0 MC", 150,0.,1.5);
  his1[156] = new TH1D("mc_ini_PpiZ", "Momentum of #pi^{#pm} MC", 300,-1.5,1.5);





  his2[114] = new TH2D("mc_depth_pi", "Depth_{#pi^{-}} vs Depth_{#pi^{+}} when #pi^{-} && #pi^{+} ", 100, -10., 60., 100, -10., 60.);
  his2[115] = new TH2D("mc_depth_mu", "Depth_{#mu^{-}} vs Depth_{#mu^{+}} when #mu^{-} || #mu^{+}", 100, -10., 60., 100, -10., 60.);
  his2[116] = new TH2D("mc_depth_others", "Depth_{-} vs Depth_{+} without #pi & #mu", 100, -10., 60., 100, -10., 60.);

  his2[117] = new TH2D("ch2_Mgg","ch2: 30-100 vs  Mgg +- (1-5)sigma",8,-0.5,7.5,5,-0.5,4.5);
  his2[118] = new TH2D("ch2_Mgg_152","ch2: 30-100 vs  Mgg +- (1-5)sigma for decCode=152",8,-0.5,7.5,5,-0.5,4.5);
  his2[119] = new TH2D("mc_chi2_Mgg_ratio", "ratio 152/all", 8,-0.5,7.5,5,-0.5,4.5);
  his2[119]->Sumw2();

  his1[140] = new TH1D("sys_chi2", "#chi^{2} for 4C-Fit", 100, 0, 100 );
  his1[141] = new TH1D("sys_Mgg", "M_{#gamma#gamma}", 100, 0.109, 0.159);

  his1[35] = new TH1D("final_Etot","E_{3#pi} final",50,2.5,3.5);
  his1[36] = new TH1D("final_ Ptot","P_{3#pi} final", 50,0.,1.);
  his1[142] = new TH1D("final_cos_theta_ppi","cos(#theta_{#pi^{+}})", 40, -1., 1.);
  his1[143] = new TH1D("final_cos_theta_mpi","cos(#theta_{#pi^{-}})", 40, -1., 1.);
  his1[144] = new TH1D("final_cos_theta_pi0","cos(#theta_{#pi^{0}})", 40, -1., 1.);
  his1[145] = new TH1D("final_Ppip","Final P_{#pi^{+}} rec.", 40, 0., 2.);
  his1[146] = new TH1D("final_Ppim","Final P_{#pi^{-}} rec.", 40, 0., 2.);
  his1[147] = new TH1D("final_Ppi0","Final P_{#pi^{0}} rec.", 40, 0., 2.);

  his1[160] = new TH1D("mc_pdg_Ppini", "Momentum of #pi^{+} init", 40,0.,2.);
  his1[161] = new TH1D("mc_pdg_Ppfin", "Momentum of #pi^{+} final", 40,0.,2.);
  his2[160] = new TH2D("mc_pdg_Pp_vs_Pm_ini","#pi^{+} vs  #pi^{-} init", 40,0.,2.,40,0.,2.);
  his1[162] = new TH1D("mc_pdg_Pmini", "Momentum of #pi^{-} init", 40,0.,2.);
  his1[163] = new TH1D("mc_pdg_Pmfin", "Momentum of #pi^{-} final", 40,0.,2.);
  his1[164] = new TH1D("mc_pdg_P0ini", "Momentum of #pi^{0} init", 40,0.,2.);
  his1[165] = new TH1D("mc_pdg_P0fin", "Momentum of #pi^{0} final", 40,0.,2.);


  his1[166] = new TH1D("mc_pdg_Ppip","Final P_{#pi^{+}} rec.", 40, 0., 2.);
  his1[167] = new TH1D("mc_pdg_Ppim","Final P_{#pi^{-}} rec.", 40, 0., 2.);
  his1[168] = new TH1D("mc_pdg_Pmup","Final P_{#mu^{+}} rec.", 40, 0., 2.);
  his1[169] = new TH1D("mc_pdg_Pmum","Final P_{#mu^{-}} rec.", 40, 0., 2.);
  his1[170] = new TH1D("mc_pdg_Pep","Final P_{#e^{+}} rec.", 40, 0., 2.);
  his1[171] = new TH1D("mc_pdg_Pem","Final P_{#e^{-}} rec.", 40, 0., 2.);
  his1[172] = new TH1D("mc_pdg_Potherp","Final P_other_{+} rec.", 40, 0., 2.);
  his1[173] = new TH1D("mc_pdg_Potherm","Final P_other_{-} rec.", 40, 0., 2.);

  his1[201] = new TH1D("Mw_4C", "M_{omega} after 4C-Fit", 400, 0.55, 0.95);
  his1[202] = new TH1D("Mw_4Cbig", "M_{omega} after 4C-Fit", 400, 0.4, 3.);
  his1[203] = new TH1D("Mgg_eta_5C_Mggcut", "Final  M_{#gamma#gamma} eta", 500, 0.3, 0.8);
  his1[204] = new TH1D("Mgg_all", "all combinations of  M_{#gamma#gamma} after 4C", 1000, 0., 1.);
  his1[205] = new TH1D("Mgg_all_3sigma", "all combinations of  M_{#gamma#gamma} after 4C (3 #sigma)", 420, 0.06, 0.2);
  his1[206] = new TH1D("Mw_4Cbig_3sigma", "M_{omega} after 4C-Fit (M_{#gamma#gamma} 3 #sigma)", 400, 0.4, 3.);

  his2[200] = new TH2D("Dalitz_Momega_Mgg", "Dalitz Momega - y,  Mgg - x", 1500, 0., 1.5, 1500, 0.2,1.7);
  his1[225] = new TH1D("Mgg_eta", "M_{#gamma#gamma} after 4C (3#sigma eta)", 500, 0.3, 0.8);
  his1[226] = new TH1D("Mw_3sigmaEta", "M_{omega} after 4C-Fit (M_{#eta} 3#sigma)", 400, 0.4, 3.);
  his1[227] = new TH1D("Mgg_eta_pi0", "M_{#gamma#gamma} after 4C (3#sigma eta, pi0)", 500, 0.3, 0.8);
  his1[228] = new TH1D("Mw_3sigmaEta_pi0", "M_{omega} after 4C-Fit (3#sigma #eta & #pi^{0})", 400, 0.4, 3.);
  his1[229] = new TH1D("Mpi0_3sigmaEta", "M_{#gamma#gamma} for pi0 after 4C (3#sigma Eta)", 420, 0.06, 0.2);
  his1[230] = new TH1D("Meta_3sigmapi0", "M_{#gamma#gamma} for eta after 4C (3#sigma pi0)", 500, 0.3, 0.8);
  his1[231] = new TH1D("Mpi0_3sigmaW", "M_{#gamma#gamma} for pi0 after 4C (3#sigma W)", 420, 0.06, 0.2);
  his1[232] = new TH1D("Meta_3sigmaW", "M_{#gamma#gamma} for eta after 4C (3#sigma W)", 500, 0.3, 0.8);
  his1[233] = new TH1D("mc_deccode_wEta", "decCode after select Omega-Eta",258,-1.5,256.5);
  //-------------------check MCGPJ
  his1[211] = new TH1D("Mw_4C_mcgpj", "M_{omega} after 4C-Fit", 400, 0.55, 0.95);
  his1[212] = new TH1D("Mw_4Cbig_mcgpj", "M_{omega} after 4C-Fit", 400, 0.4, 3.);
  his1[214] = new TH1D("Mgg_all_mcgpj", "all combinations of  M_{#gamma#gamma} after 4C", 1000, 0., 1.);
  his1[215] = new TH1D("Mgg_all_3sigma_mcgpj", "all combinations of  M_{#gamma#gamma} after 4C (3 #sigma)", 420, 0.06, 0.2);
  his1[216] = new TH1D("Mw_4Cbig_3sigma_mcgpj", "M_{omega} after 4C-Fit (M_{#gamma#gamma} 3 #sigma)", 400, 0.4, 3.);
  his1[217] = new TH1D("M_2pi_mcgpj", "M_{2#pi} after 4C-Fit", 3000, 0., 3.);
  his1[218] = new TH1D("M_2pi_3sigma_mcgpj", "M_{2#pi} after 4C-Fit (M_{#gamma#gamma} 3 #sigma)", 3000, 0., 3.);
  his1[219] = new TH1D("M_2pi_3sigmaW_mcgpj", "M_{2#pi} after 4C-Fit (M_{#gamma#gamma}+W 3 #sigma)", 3000, 0., 3.);
  his2[210] = new TH2D("Dalitz_Momega_Mgg_mcgpj", "Dalitz Momega - y,  Mgg - x", 1500, 0., 1.5, 1500, 0.2,1.7);
  his1[220] = new TH1D("mc_Mw", "M_{#pi^{+}#pi^{-}#gamma#gamma}", 3000, 0.4, 3.4);
  his1[37] = new TH1D("SoverS0_ini","MCGPJ: S'/S init. (for best mass)", 120, 0., 1.2);
  his1[38] = new TH1D("SoverS0_fin","MCGPJ: S'/S fin. (for best mass)",  120, 0., 1.2);


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
  his1[262] = new TH1D("M_mis", "M-missing", 2000, -1., 1.);
  his1[263] = new TH1D("mc_Pz_mis_69", "Pz-missing, iMC signal", 300,-1.5, 1.5);
  his1[264] = new TH1D("mc_Pt_mis_69", "Pt-missing, iMC signal", 150, 0., 1.5);
  his1[265] = new TH1D("mc_M_mis_69", "M-missing, iMC signal", 2000, -1., 1.);
  his1[266] = new TH1D("mc_Pz_mis_bkgr", "Pz-missing, iMC bkgr", 300,-1.5, 1.5);
  his1[267] = new TH1D("mc_Pt_mis_bkgr", "Pt-missing, iMC bkgr", 150, 0., 1.5);
  his1[268] = new TH1D("mc_M_mis_bkgr", "M-missing, iMC bkgr", 2000, -1., 1.);
  his1[269] = new TH1D("mc_deccode_Mmis008", "decCode, M-missing < 0.08",258,-1.5,256.5);
  his1[270] = new TH1D("mc_deccode_Mmis010", "decCode, M-missing < 0.1",258,-1.5,256.5);
  his1[271] = new TH1D("mc_deccode_Mmis012", "decCode, M-missing < 0.12",258,-1.5,256.5);
  his1[272] = new TH1D("mc_deccode_Mmis006", "decCode, M-missing < 0.06",258,-1.5,256.5);
  his1[273] = new TH1D("mc_deccode_Mmis005", "decCode, M-missing < 0.05",258,-1.5,256.5);
  his1[274] = new TH1D("mc_deccode_Mmis004", "decCode, M-missing < 0.04",258,-1.5,256.5);
  his1[275] = new TH1D("Nevents_Mmis","N_{events} for M-missing cuts 40-120 MeV", 11, 2.5, 13.5);
  his1[276] = new TH1D("M2_mis", "M^{2}-missing", 2000, -1., 1.);
  his1[277] = new TH1D("mc_M2_mis_69", "M^{2}-missing, iMC signal", 2000, -1., 1.);
  his1[278] = new TH1D("mc_M2_mis_bkgr", "M^{2}-missing, iMC bkgr", 2000, -1., 1.);
  his1[279] = new TH1D("mc_deccode_Mmis_108", "decCode, -0.1 < M-missing < 0.08",258,-1.5,256.5);
  his1[281] = new TH1D("mc_deccode_Mmis_107", "decCode, -0.1 < M-missing < 0.07",258,-1.5,256.5);
  his1[282] = new TH1D("mc_deccode_Mmis_106", "decCode, -0.1 < M-missing < 0.06",258,-1.5,256.5);
  his1[283] = new TH1D("mc_deccode_Mmis_105", "decCode, -0.1 < M-missing < 0.05",258,-1.5,256.5);
  his1[284] = new TH1D("M_mis_1ISR_less50", "M-missing 1 ISR, E_isr < 50. MeV", 2000, -1., 1.);
  his1[285] = new TH1D("M_mis_1ISR_more50", "M-missing 1 ISR, E_isr > 50. MeV", 2000, -1., 1.);
  his1[286] = new TH1D("M_mis_2ISR_less50", "M-missing 2 ISR, E_isr < 50. MeV", 2000, -1., 1.);
  his1[287] = new TH1D("M_mis_2ISR_more50", "M-missing 2 ISR, E_isr > 50. MeV", 2000, -1., 1.);

  //-------------weta + add ISR in 4C
  his1[360] = new TH1D("Pz_mis_weta", "Pz-missing #omega#eta", 300,-1.5, 1.5);
  his1[361] = new TH1D("Pt_mis_weta", "Pt-missing #omega#eta", 150, 0., 1.5);
  his1[362] = new TH1D("M_mis_weta", "M-missing #omega#eta", 2000, -1., 1.);
  his1[363] = new TH1D("mc_Pz_mis_67", "Pz-missing #omega#eta, iMC signal", 300,-1.5, 1.5);
  his1[364] = new TH1D("mc_Pt_mis_67", "Pt-missing #omega#eta, iMC signal", 150, 0., 1.5);
  his1[365] = new TH1D("mc_M_mis_67", "M-missing #omega#eta, iMC signal", 2000, -1., 1.);
  his1[366] = new TH1D("mc_Pz_mis_weta_bkgr", "Pz-missing #omega#eta, iMC bkgr", 300,-1.5, 1.5);
  his1[367] = new TH1D("mc_Pt_mis_weta_bkgr", "Pt-missing #omega#eta, iMC bkgr", 150, 0., 1.5);
  his1[368] = new TH1D("mc_M_mis_weta_bkgr", "M-missing #omega#eta, iMC bkgr", 2000, -1., 1.);
  his1[369] = new TH1D("mc_deccode_Mmis_weta_008", "decCode #omega#eta, M-missing < 0.08",258,-1.5,256.5);
  his1[370] = new TH1D("mc_deccode_Mmis_weta_010", "decCode #omega#eta, M-missing < 0.1",258,-1.5,256.5);
  his1[371] = new TH1D("mc_deccode_Mmis_weta_012", "decCode #omega#eta, M-missing < 0.12",258,-1.5,256.5);
  his1[372] = new TH1D("mc_deccode_Mmis_weta_006", "decCode #omega#eta, M-missing < 0.06",258,-1.5,256.5);
  his1[373] = new TH1D("mc_deccode_Mmis_weta_005", "decCode #omega#eta, M-missing < 0.05",258,-1.5,256.5);
  his1[374] = new TH1D("mc_deccode_Mmis_weta_007", "decCode #omega#eta, M-missing < 0.07",258,-1.5,256.5);
  his1[375] = new TH1D("Nevents_Mmis_weta","N_{events} #omega#eta  for M-missing cuts 5-12 (*10 MeV)", 11, 2.5, 13.5);
  his1[376] = new TH1D("M2_mis_weta", "M^{2}-missing #omega#eta", 2000, -1., 1.);
  his1[377] = new TH1D("mc_M2_mis_67", "M^{2}-missing #omega#eta, iMC signal", 2000, -1., 1.);
  his1[378] = new TH1D("mc_M2_mis_weta_bkgr", "M^{2}-missing #omega#eta, iMC bkgr", 2000, -1., 1.);

  his1[379] = new TH1D("mc_Emax_othergammas_weta67", "Emax for other good #gamma (MeV), #omega#eta, 67",3000,0.,3000.);
  his1[380] = new TH1D("mc_N_othergammas_weta67","N_{other good  #gamma} in event  #omega#eta, 67", 11, -0.5, 10.5);
  his1[381] = new TH1D("mc_Emax_othergammas_weta_bkgr", "Emax for other good #gamma (MeV), #omega#eta, bkgr",3000,0.,3000.);
  his1[382] = new TH1D("mc_N_othergammas_weta_bkgr","N_{other good  #gamma} in event  #omega#eta, bkgr", 11, -0.5, 10.5);
  his1[383] = new TH1D("Emax_othergammas_weta", "Emax for other good #gamma (MeV), #omega#eta",3000,0.,3000.);
  his1[384] = new TH1D("N_othergammas_weta","N_{other good  #gamma} in event  #omega#eta", 11, -0.5, 10.5);
  his1[385] = new TH1D("mc_Emax_othergammas_weta_w2pi0", "Emax for other good #gamma (MeV), #omega#eta, #omega#pi0#pi0",3000,0.,3000.);
  his1[386] = new TH1D("mc_N_othergammas_weta_w2pi0","N_{other good  #gamma} in event  #omega#eta, #omega#pi0#pi0", 11, -0.5, 10.5);

  m_tuple[0] = new TNtuple("a4C","after 4C fit",
                                "ch2:Eisr:PZisr:PTisr:"
                                "Etot:Ptot:Esum:PZsum:PTsum:"
                                "mcNg:mcEisr:mcEmisr:mcThisr:"
                                "dc:PZmis:PTmis:Mmis"

*/

  // register in selector to save in given directory
  const char* SaveDir = "SelectOmegaEta";
  VecObj his1o(his1.begin(),his1.end());
  selector->RegInDir(his1o,SaveDir);
  VecObj his2o(his2.begin(),his2.end());
  selector->RegInDir(his2o,SaveDir);
  VecObj ntuples(m_tuple.begin(),m_tuple.end());
  selector->RegInDir(ntuples,SaveDir);

}

//-----------------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo)
//-----------------------------------------------------------------------------
{
  static int save_runNo = 0;
  static Hep3Vector xorigin;

  if ( runNo == save_runNo ) return xorigin;

  // update vertex for new run
  xorigin.set(0.,0.,0.);
  VertexDbSvc* vtxsvc = VertexDbSvc::instance();

  int run = abs(runNo);
  if (
          (run >= 28241 && run <= 28266)  // 3080-old
       || (run >=  9947 && run <= 10878)  // J/Psi 2009
     ) {
    vtxsvc->SetBossVer("6.6.3");
  } else if ( run < 30000 ) {   // here: J/Psi scan
    vtxsvc->SetBossVer("6.6.4");
  } else {
    vtxsvc->SetBossVer("6.6.5");  // here: 3080 from R-scan (2015)
  }

  vtxsvc->handle(runNo);

  if ( vtxsvc->isVertexValid() ) {
    double* dbv = vtxsvc->PrimaryVertex();
    xorigin.set(dbv[0],dbv[1],dbv[2]);
    // cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
  } else {
    cout << "Cannot obtain vertex information for run#" << runNo << endl;
    ERROR_WARNING++;
    exit(1);
  }

  save_runNo = runNo;
  return xorigin;
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

  static runInfo ListRuns[] =  {
        { 28312, 28346,     3050.254,    915,     14.918 },
        { 28347, 28381,     3059.273,    842,     15.059 },
        { 28241, 28266,     3080.224,    839,     17.513 },
        { 28382, 28387,     3082.491,    883,      2.711 },
        { 28466, 28469,     3084.351,   1019,      2.057 },
        { 28388, 28416,     3089.213,    905,     13.934 },
        { 28472, 28475,     3091.284,    777,      1.624 },
        { 28417, 28449,     3092.104,    864,     12.354 },
        { 28476, 28478,     3093.716,    837,      1.843 },
        { 28479, 28482,     3095.299,    819,      2.143 },
        { 28487, 28489,     3096.017,   1077,      1.816 },
        { 28490, 28492,     3096.432,    894,      2.134 },
        { 28493, 28495,     3097.783,    710,      2.069 },
        { 28496, 28498,     3098.943,    833,      2.203 },
        { 28499, 28501,     3099.572,   1170,      0.756 },
        { 28504, 28505,     3101.914,    960,      1.612 },
        { 28506, 28509,     3106.146,    981,      2.106 },
        { 28510, 28511,     3112.609,    667,      1.720 },
        { 28512, 28513,     3120.294,   1116,      1.264 },
        { 33490, 33556,     3810.951 + 0.55,   1088,     50.54 },
        { 33572, 33657,     3899.608 + 0.55,   1025,     52.61 },
        { 33659, 33719,     4090.238 + 0.55,    722,     52.63 },
        { 30372, 30437,     4190.557 + 0.55,    644,     43.09 },
        { 32046, 32140,     4218.961 + 0.55,   1497,     54.13 },
        { 30438, 30491,     4229.893 + 0.55,   1025,     44.40 },
        { 32141, 32226,     4245.102 + 0.55,   1583,     55.59 },
        { 30492, 30557,     4310.388 + 0.55,   1071,     44.90 },
        { 31281, 31325,     4389.833 + 0.55,   1577,     55.18 },
                        };

  int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

  int absrun = abs(run);

  // we need some values any way
  double Ecms   = 3.097;
//   double Spread = 1000.;
  Lumi          = 0;

  int i = 0;
  for(; i < Np; i++) {
     if ( absrun >= ListRuns[i].runS && absrun <= ListRuns[i].runE ) {
       Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV
//        Spread = ListRuns[i].Spread * 1.e-6;  // KeV -> GeV
       Lumi   = ListRuns[i].Lumi;
       break;
     }
  }

  if ( i != Np ) {
    // correct energy according Yadi Wang fit of J/Psi peak"
    Ecms -= 0.55e-3; // -0.55MeV
  } else {
    // run not in the list
    if ( absrun >= 39355 && absrun <= 39618 ) { // 3080 from R-scan (2015)
      Ecms =  3.080; // no BEMS information
      Lumi = 120.11483;
    } else  if ( absrun >= 9947 && absrun <= 10878 ) { // J/Psi 2009
      Ecms =  3.097; // for inclusive MC
      Lumi = 0.;
    } else  if ( absrun >= 31327 && absrun <= 31390 ) {
      Ecms =  4.420; // no BEMS information
      Lumi =  44.67;
    } else if ( absrun >= 31983 && absrun <= 32045 ) {
      Ecms =  4.210; // no BEMS information
      Lumi =  54.55;
    } else if ( absrun >= 27147 && absrun <= 27233 ) {
      Ecms =  3.0827;
      Lumi = 13.5158;
      cout << " GetEcms::WARNING run# " << run
           << " this is old not recommended run" << endl;
    } else {
      cout << " GetEcms::WARNING unknown run# " << run
           << " use Ecms= " << Ecms << endl;
      ERROR_WARNING++;
    }
  }

  return Ecms;
}

/*
//----------------------------------------------------------------------------
static void IsPionsFromOmegaMC(const ReadDst* selector,
                TEvtHeader* m_TEvtHeader, TEvtRecObject* m_TEvtRecObject)
//----------------------------------------------------------------------------
{
 //-----------------------for tracks
  int runNo = m_TEvtHeader->getRunId();
  int eventNo = m_TEvtHeader->getEventId();
  Hep3Vector xorigin = getVertexOrigin(runNo);
  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  //-----------------------for Kalman tracks
  ParticleID *pid = ParticleID::instance();
  pid->set_path(selector->AbsPath("Analysis/ParticleID"));

  //---------------------for MC particles
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    if ( part_pdg==223) {

      vector<Int_t> dref = part->getDaughters ();
      for (unsigned int i = 0; i < dref.size(); i++) {
        const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
        int daughter_pdg = daughter->getParticleID ();

        if (daughter_pdg==211 || daughter_pdg==-211){


          // Match charged reconstructed tracks with MC particles
          Hep3Vector p3mc( daughter->getInitialMomentumX(),
                     daughter->getInitialMomentumY(),
                     daughter->getInitialMomentumZ() );

          double maxdp = 0.1;
          double maxangle = 0.035; // ~ 2deg.
          double angle = 3.14;
          double dP = 1.;
          int ipi_mc = -1;

          for(int i = 0; i < evtRecEvent->totalCharged(); i++){
            DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
            if( !itTrk->isMdcTrackValid() ) continue;


            //----------------------Kalman tracks info
            pid->init();
            pid->setMethod(pid->methodProbability());
            pid->setChiMinCut(4);
            pid->setRecTrack(itTrk);

            // list of systems used for PID
            pid->usePidSys( pid->useDedx() |
                      pid->useTof1() | pid->useTof2());
            //                 pid->useEmc()  |
            //                 pid->useMuc()


            pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton()| pid->onlyMuon() | pid->onlyElectron());

            pid->calculate(runNo);
            if ( !pid->IsPidInfoValid() ) {
        continue;
            }



            //-----------------------------------------------------------RecMdcKalTrac
            RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
            if ( !mdcKalTrk ) {
        cout << " No Kalman track ! EventNo= "<<eventNo << endl;
        continue;
            }
            //----------------------------------------------------
             if ( pid->probPion() >  0.001 ) {

        mdcKalTrk->setPidType(RecMdcKalTrack::pion);
        Hep3Vector p3pi(mdcKalTrk->px(), mdcKalTrk->py(), mdcKalTrk->pz());
        angle = p3pi.angle(p3mc);
        dP = abs( p3mc.mag() - p3pi.mag() );
            } else {
        continue;
            }

            if ( dP < maxdp && angle < maxangle ) {
        maxdp = dP;
        maxangle  = angle;
        ipi_mc = i;
            }

          }//----------------End for(i)
          if(ipi_mc!=-1){
            DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(ipi_mc);
            RecMdcTrack *mdcTrk = itTrk->mdcTrack();
            HepVector a = mdcTrk->helix();
            HepSymMatrix Ea = mdcTrk->err();
            HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
            HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
            VFHelix helixip(point0,a,Ea);
            helixip.pivot(IP);
            HepVector vecipa = helixip.a();
            double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
            double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction


            his1[40]->Fill(Rvxy0);
            his1[41]->Fill(Rvz0);

          }

        }//------------------End if(211)
      }//--------------------End for(daughters)
      return;
    }//----------------------End if(223)
  }//------------------------End while(parts)
  return;
}
*/

/*
//----------------------------------------------------------------------------
static  void PionsMomentumHistMC(const ReadDst* selector, bool ini_fin_flag)                  // flag  0 - ini, 1 - fin
//----------------------------------------------------------------------------
{
  if ( !isMC ) return;
  double Pp=0,Pm=0;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return;

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part = static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    if ( part_pdg==111 || part_pdg==211 || part_pdg==-211 ) {

      Hep3Vector pi_3p( part->getInitialMomentumX(),
                        part->getInitialMomentumY(),
                        part->getInitialMomentumZ() );

      if( part_pdg==211 ) {
        if ( !ini_fin_flag ){ Pp=pi_3p.mag(); his1[160]->Fill( Pp );}
        else his1[161]->Fill( pi_3p.mag());
      }

      if( part_pdg==-211 ) {
        if ( !ini_fin_flag ){Pm=pi_3p.mag();  his1[162]->Fill( Pm );}
        else his1[163]->Fill( pi_3p.mag());
      }

      if( part_pdg==111 ) {
        if ( !ini_fin_flag ) his1[164]->Fill( pi_3p.mag());
        else his1[165]->Fill( pi_3p.mag());
      }
    }
  }

  if (!ini_fin_flag )  his2[160]->Fill( Pp,Pm );
  return;
}
*/

/*
//-----------------------------------------------------------------------------
static int ParticleTypesMC(const ReadDst* selector, Select& Slct)
//-----------------------------------------------------------------------------
{
  omega_to_pis=false;
  omega_best_mass=false;
  one_ISR=false;
  two_ISR=false;
  one_ISR_less50=false;
  one_ISR_more50=false;
  two_ISR_less50=false;
  two_ISR_more50=false;
  SoverS0 = false;
  EtaTo2Gammas = false;
  if ( !isMC ) return -1;

  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol ();
  if ( !mcParticles ) return -1;

  int pi_count=0, omega_count=0, PiFromOmega_count=0, other_count=0;

  HepLorentzVector mc_pi[4]; // pi+ [0], pi- [1], pi0  [2,3]
  int n_mcpi[] = { 0, 0, 0, 0, 0}; // [4] for all pi0

  // ISR: initial state radiation
  int Ng_isr = 0;       // number of gammas
  int Eta_Ngammas = 0; // number of gammas for eta decay
  double E_isr = 0.;    // total energy
  double Emax_isr = 0.; // gammas with maximal energy
  double E_ISR[2] = {0.,0.};

  // TIterator *mcParticlesIter = mcParticles->MakeIterator ();
  TIter mcParticlesIter(mcParticles);
  while ( TMcParticle* part =
                static_cast<TMcParticle*>(mcParticlesIter.Next()) ) {
    long int part_pdg = part->getParticleID ();

    his1[151]->Fill(part_pdg);
    if ( part->getMother() == -99 ) his1[152]->Fill(part_pdg);

    if ( part_pdg == 22 ) { // gamma
      if ( part->getMother() == -99 ) { // from primary vertex => ISR gammas

        double Eg = part->getInitialMomentumE() * 1e3; // MeV
        E_ISR[Ng_isr] = Eg;
        Ng_isr++;
        his1[111]->Fill(Eg);
        E_isr += Eg;
        if ( Eg > Emax_isr ) Emax_isr = Eg;
      }
      continue; // ignore gammas
    }

    //eta
     if ( part_pdg==221) {

       vector<Int_t> dref = part->getDaughters ();
       for (unsigned int i = 0; i < dref.size(); i++) {
         const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
         int daughter_pdg = daughter->getParticleID ();
         his1[153]->Fill(daughter_pdg);
         if(daughter_pdg==22) Eta_Ngammas++;
       }//end for(Daughters)
       if(Eta_Ngammas==2) EtaTo2Gammas = true;
     }//end if(221)


     //all particles
     if ( part->getMother() == -99 ) { // from primary vertex

       if ( part_pdg==111 || part_pdg==211 || part_pdg==-211  || part_pdg==221) { // pi 0/+/-, eta
         pi_count++;

         Hep3Vector pi_3p( part->getInitialMomentumX(),
                           part->getInitialMomentumY(),
                           part->getInitialMomentumZ() );
         double pi_E = part->getInitialMomentumE();
         double pi_P = pi_3p.mag();

         if ( part_pdg == 111 ){
           his1[155]->Fill(pi_P);

           mc_pi[2] = HepLorentzVector( pi_3p, pi_E );
           n_mcpi[2]++;

         }else if(part_pdg==221) {

           mc_pi[3] = HepLorentzVector( pi_3p, pi_E );
           n_mcpi[3]++;


         } else {
           double sign = part_pdg/211.;
           his1[156]->Fill(pi_P*sign);

           if ( part_pdg==211 ) {
             mc_pi[0] = HepLorentzVector( pi_3p, pi_E );
             n_mcpi[0]++;
           }
           if ( part_pdg==-211 ) {
             mc_pi[1] = HepLorentzVector( pi_3p, pi_E );
             n_mcpi[1]++;
           }

         }
         continue;
       } // end of if(pions)
     }// end 'from primary vertex'





       if ( part_pdg==223) {
         omega_count++;

         vector<Int_t> dref = part->getDaughters ();
         for (unsigned int i = 0; i < dref.size(); i++) {
           const TMcParticle *daughter = m_TMcEvent->getMcParticle (dref[i]);
           int daughter_pdg = daughter->getParticleID ();
           if (daughter_pdg==111 ||
               daughter_pdg==211 || daughter_pdg==-211) PiFromOmega_count++;
         }
         continue;
       }

       other_count++;

     } //-----------------------------------------------End  while()

     // summary for ISR:
  Slct.SoverS0 = 4.*(0.5*Slct.E_cms - E_ISR[0])*(0.5*Slct.E_cms - E_ISR[1])/(Slct.E_cms*Slct.E_cms);
  if( Slct.SoverS0  >=  0.8 ) SoverS0 = true;

  his1[112]->Fill(E_isr);
  his1[115]->Fill(Emax_isr);
  his1[116]->Fill(Ng_isr);
  if ( Ng_isr == 1 ) his1[117]->Fill(E_isr);
  if ( Ng_isr == 2 ) {
    his1[118]->Fill(E_isr);
    his1[119]->Fill(Emax_isr);
  }

  if ( Ng_isr == 1 ){
    one_ISR=true;
    if( E_isr < 50.) one_ISR_less50=true;
    else one_ISR_more50=true;
  }
  if ( Ng_isr == 2 ){
    his2[252]->Fill(E_ISR[0],E_ISR[1]);

    two_ISR=true;
    if( E_isr < 50.) two_ISR_less50=true;
    else two_ISR_more50=true;
  }

  if ( selector->Verbose() ) {
    selector->PrintMCParticle(mcParticles);
    cout << " pi_count= " << pi_count << " omega_count= " << omega_count
         << " PiFromOmega= " << PiFromOmega_count << " other_count= "
         << other_count << endl;
  }

  // fill histograms for  pi+, pi-, pi0, eta
  if ( n_mcpi[0] == 1 && n_mcpi[1] == 1 && n_mcpi[2] == 1 && n_mcpi[3] == 1 ) {
    HepLorentzVector POmegapi0 = mc_pi[0] + mc_pi[1] + mc_pi[2] + mc_pi[3];
    double EOmegapi0 = POmegapi0.e() * 1e3;
    his1[113]->Fill( EOmegapi0 );
    his1[114]->Fill( EOmegapi0 + E_isr ); // Etot

    HepLorentzVector Pw = mc_pi[0] + mc_pi[1] + mc_pi[2];
    //  HepLorentzVector Pw2 = mc_pi[0] + mc_pi[1] + mc_pi[3];
     his1[220]->Fill( Pw.m());

    if( (Pw.m() < 0.852  &&  Pw.m() > 0.712) ) {omega_best_mass=true;   his1[37]->Fill(Slct.SoverS0); }

    if ( selector->Verbose() ) {
      cout << " pi+, pi-, pi0, pi0:  POmegapi0= " <<  POmegapi0  << endl;
      cout << " E_isr= " << E_isr << " Etot= " << EOmegapi0+E_isr << endl;
    }

    if (  omega_count==1 && PiFromOmega_count == 3 ) omega_to_pis=true;
  }

//   if ( other_count > 0 ) return 0; // background
  return 0;
}
*/

/*
//-----------------------------------------------------------------------------
static void MatchRecMcTrks(ReadDst* selector, Select& Slct)
//-----------------------------------------------------------------------------
{
  const static double maxdp = 0.1;
  const static double maxangle = 0.035; // ~ 2deg.

  double maxdp_p = maxdp,maxdp_m = maxdp;
  double maxa_p = maxangle, maxa_m = maxangle;

  // Match charged reconstructed tracks with MC particles
  if ( !isMC ) return;
  const TMcEvent* m_TMcEvent = selector->GetMcEvent();
  const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();

  Hep3Vector p3pip = Slct.ppip.vect();
  Hep3Vector p3pim = Slct.ppim.vect();

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

    if ( dP_p < maxdp && angle_p < maxa_p ) {
      maxdp_p = dP_p;
      maxa_p  = angle_p;
      Slct.ipip_mc = i;
      Slct.pip_pdg = pId;
    }

    if ( dP_m < maxdp && angle_m < maxa_m ) {
      maxdp_m = dP_m;
      maxa_m  = angle_m;
      Slct.ipim_mc = i;
      Slct.pim_pdg = pId;
    }
  }//----------------------------------------------------End for(i)


  // check matching:
  his1[121]->Fill( maxa_p );
  his1[121]->Fill( -maxa_m );
  his1[122]->Fill( maxdp_p );
  his1[122]->Fill( -maxdp_m );
  his1[123]->Fill( ( Slct.pip_pdg == 211 ) ? 1. : 2. );
  his1[123]->Fill( ( Slct.pim_pdg == -211 ) ? -1. : -2. );

  return;
}
*/



//-----------------------------------------------------------------------------
bool SelectOmegaEtaEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << " start " << __func__ << "()" << endl;

  m_abscor->AbsorptionCorrection(selector);

  Select Slct; // information about selected pions and gammas

  //---------------------------------------------------------------------------
  //-- Get event information --
  //---------------------------------------------------------------------------
  int runNo   = m_TEvtHeader->getRunId();
  int eventNo = m_TEvtHeader->getEventId();
  his1[0]->Fill(abs(runNo));
  his1[1]->Fill(0); // "cuts"

  double lum0 = 0;
  Slct.Ecms =  GetEcms(abs(runNo), lum0);
  his1[9]->Fill(Slct.Ecms);

  isMC = (runNo < 0);
  if ( (isMC && m_TMcEvent->getMcParticleCol()->IsEmpty() ) ||
       ( !isMC &&
         m_TMcEvent != 0 &&
         !m_TMcEvent->getMcParticleCol()->IsEmpty() )
     ) {
    cout << " ERROR: something wrong: isMC= " << isMC
         << " McParticles= " << m_TMcEvent->getMcParticleCol()->GetEntries()
         << endl;
    return false;
  }

  int mc_type = -1;
  int decCode = -1;
  if ( isMC ) {
    m_EventTagSvc->setMcEvent(m_TMcEvent);
//nef    mc_type = ParticleTypesMC(selector, Slct);
    decCode = ((m_EventTagSvc->getEventTag())&0xFF00) >> 8;
    his1[101]->Fill(mc_type);
    his1[102]->Fill(decCode);

    //for MCGPJ
    if(!SoverS0) return false;

//nef    if(omega_to_pis) IsPionsFromOmegaMC(selector, m_TEvtHeader, m_TEvtRecObject);

    if(EtaTo2Gammas) his1[108]->Fill(0);
    if(omega_best_mass) his1[280]->Fill(0);
    if(one_ISR) his1[300]->Fill(0);
    if(two_ISR) his1[301]->Fill(0);
    if(one_ISR && omega_best_mass) his1[302]->Fill(0);
    if(two_ISR && omega_best_mass) his1[303]->Fill(0);

//nef  PionsMomentumHistMC( selector, 0);  // init MC  momentum  hists
  }

  Hep3Vector xorigin = getVertexOrigin(runNo);

  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  //---------------------------------------------------------------------------
  // Charged tracks
  //---------------------------------------------------------------------------
  // 1) select two good tracks with opposite charges from primary vertex
  int Ncharge = 0;
  for(int i = 0; i < evtRecEvent->totalCharged(); i++) {

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

    his1[11]->Fill(Rvxy0);
    his1[12]->Fill(Rvz0);
    his1[13]->Fill(fabs(cosTheta));

    if(fabs(Rvz0) >= 10.0) continue;
    if(fabs(Rvxy0) >= 1.0) continue;
    if(fabs(cosTheta) >= 0.93) continue;

    Slct.qtrks.push_back(i);
    Ncharge += mdcTrk->charge();
    his1[14]->Fill( mdcTrk->charge()*mdcTrk->p() );
  }
  //-------------------------------------------------End for(i)

  int nGood = Slct.nGood();
  his1[15]->Fill(nGood);
  his1[16]->Fill(Ncharge);

  if ( nGood != 2 || Ncharge != 0 ) return false; // skip event
  his1[1]->Fill(1); // "cuts"

  // 2) check that these tracks have PID information
  //    and fitted as Kalman tracks
  ParticleID *pid = ParticleID::instance();
#if (BOSS_VER < 700)
  pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
  pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

  int npp = 0, npm = 0; // number of +/- in event
  for(int i = 0; i < nGood; i++) {
    DstEvtRecTracks* itTrk=(DstEvtRecTracks*) evtRecTrkCol->At(Slct.qtrks[i]);

    pid->init();
    pid->setMethod(pid->methodProbability());
    pid->setChiMinCut(4);
    pid->setRecTrack(itTrk);

    // list of systems used for PID
    pid->usePidSys( pid->useDedx() |
                    pid->useTof1() | pid->useTof2() |
                    pid->useTofE() );
    //                 pid->useEmc()  |
    //                 pid->useMuc()

    pid->identify(
                pid->onlyPion()    |
                pid->onlyKaon()    |
                pid->onlyProton()  |
                pid->onlyMuon()    |
                pid->onlyElectron()
        );

    pid->calculate(runNo);
    his1[21]->Fill( (double)pid->IsPidInfoValid() );
    if ( !pid->IsPidInfoValid() ) return false;

    if( pid->probPion() > 0 )   his1[22]->Fill( log10(pid->probPion()) );
    else                        his1[22]->Fill(-5.);
    if( pid->probKaon() > 0 )   his1[23]->Fill( log10(pid->probKaon()) );
    else                        his1[23]->Fill(-5.);
    if( pid->probProton() > 0 ) his1[24]->Fill( log10(pid->probProton()) );
    else                        his1[24]->Fill(-5.);

    if ( !init_selection ) {
      // select Pions only
      if ( pid->probPion() <  0.001 ) return false;
    }

    RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
    his1[25]->Fill( (double)(mdcKalTrk != 0) );
    if ( !mdcKalTrk ) {
      cout << " No Kalman track ! EventNo= "<<eventNo << endl;
      return false;
    }

    // define charge by Kalman fit
    mdcKalTrk->setPidType(RecMdcKalTrack::pion);

    //-- E,Px,Py,Pz for good charged tracks
    Hep3Vector p3pi(mdcKalTrk->px(), mdcKalTrk->py(), mdcKalTrk->pz());
    double epi = sqrt( p3pi.mag2() + mpi*mpi );
    if ( mdcKalTrk->charge() > 0 ) {
      npp++;
      Slct.ipip = Slct.qtrks[i];
      Slct.ppip = HepLorentzVector( p3pi, epi );
    } else {
      npm++;
      Slct.ipim = Slct.qtrks[i];
      Slct.ppim = HepLorentzVector( p3pi, epi );
    }
  }  //-------------------------------------------------End for(i)

  his2[21]->Fill( npp, npm );
  if ( !init_selection ) {
    // require exactly one "+" and one "-"
    if ( npp != 1 || npm != 1 ) return false;
  }
  his1[1]->Fill(2);

  if(isMC){
    if(omega_best_mass) his1[280]->Fill(1);
    if(one_ISR) his1[300]->Fill(1);
    if(two_ISR) his1[301]->Fill(1);
    if(one_ISR && omega_best_mass) his1[302]->Fill(1);
    if(two_ISR && omega_best_mass) his1[303]->Fill(1);
  }
  //---------------------------------- Finish Good Charged Track Selection

  // Match charged reconstructed tracks with MC particles
//nef  if ( isMC ) MatchRecMcTrks(selector, Slct);

  //---------------------------------------------------------------------------
  // Neutral tracks
  //---------------------------------------------------------------------------
  for(int i = evtRecEvent->totalCharged();
          i < evtRecEvent->totalTracks(); i++) {
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if(!itTrk->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = itTrk->emcShower();

    // 1) good EMC time:
    if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) continue;

    // 2) good EMC energy deposited in the barrel (endcap) part of EMC
    double eraw = emcTrk->energy();
    double absCosTheta = fabs(  cos(emcTrk->theta()) );

    bool GoodCluster=false;
    if ( absCosTheta < 0.8 ) {  //barrel
      GoodCluster = eraw > 25E-3;
    } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) { //endcap
      GoodCluster = eraw > 50E-3;
    }
    if ( !GoodCluster ) continue;

    // 3) the nearest charged track is far from cluster
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

    double tang = 200.; // min angle between cluster and track
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      DstEvtRecTracks* jtTrk = (DstEvtRecTracks*) evtRecTrkCol->At(j);
      if ( !jtTrk->isExtTrackValid() ) continue;
      RecExtTrack *extTrk = jtTrk->extTrack();
      if ( extTrk->emcVolumeNumber() == -1 )
        continue; //track does not hit EMC
      Hep3Vector extpos = extTrk->emcPosition();
      double angd = fabs(extpos.angle(emcpos)); // [0,pi]
      if ( angd < tang ) tang = angd;
    }
    //--------------------------------------------End for(j)

    his1[31]->Fill(RtoD(tang));
    static const double min_angle = 10 * M_PI/180; // 10 grad
    if ( tang < min_angle ) continue;

    Slct.gammas.push_back(i);

    // (E,Px,Py,Pz) for good gammas
    Hep3Vector p3 = emcpos - xorigin; // position emc cluster wrt of vertex
    p3 *= eraw / p3.mag();            // normalization on energy
    Slct.Pg.push_back( HepLorentzVector(p3,eraw) );
  }
  //-----------------------------------------------End for(i)

  int nGam = Slct.nGam();
  his1[32]->Fill(nGam);

  int min_Gammas = 4;
  if ( init_selection ) min_Gammas = 2;
  if ( nGam < min_Gammas ) return false; // skip event
  his1[1]->Fill(3);
  //------------------------------------- Finish Good Photon Selection

  his1[33]->Fill( Slct.ppip.vect().mag());
  his1[33]->Fill(-Slct.ppim.vect().mag());
  for(int i = 0; i < nGam; i++) {
     his1[34]->Fill(Slct.Pg[i].e());
  }

  //---------------------------------------------------------------------------
  if (init_selection) return true;
  //---------------------------------------------------------------------------

/*

  if(isMC){
    if(omega_best_mass) his1[280]->Fill(2);
    if(one_ISR) his1[300]->Fill(2);
    if(two_ISR) his1[301]->Fill(2);
    if(one_ISR && omega_best_mass) his1[302]->Fill(2);
    if(two_ISR && omega_best_mass) his1[303]->Fill(2);

    his1[125]->Fill(Slct.pip_pdg);
    his1[126]->Fill(Slct.pim_pdg);
  }

  // Vertex & Kinematic Fit
  RecMdcKalTrack *pipTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.ipip))->mdcKalTrack();
  RecMdcKalTrack *pimTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.ipim))->mdcKalTrack();

  WTrackParameter wvpipTrk, wvpimTrk;
  wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
  wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());

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
  vtxfit->AddTrack(0,  wvpipTrk);
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);
  bool ok = vtxfit->Fit(0);
  his1[14]->Fill( double(ok) );
  if ( !ok )  return false; //skip event

  vtxfit->BuildVirtualParticle(0);
  his1[15]->Fill(vtxfit->chisq(0));
  vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0

  WTrackParameter wpip = vtxfit->wtrk(0);
  WTrackParameter wpim = vtxfit->wtrk(1);

  KinematicFit * kmfit = KinematicFit::instance();
  //--------------------------------------------------------------------------------------------------------------------------------Ecms =  GetEcms(runNo, lum);
//   double lum = 0;
//   double Ecms =  GetEcms(abs(runNo), lum);

  //  double Ecms =  3.097;  // for inclusive MC


  HepLorentzVector ecms(Ecms*sin(beam_angle), 0, 0, Ecms);
  //------------------------------------------------------4C Fit
  // make fake ISR photon: Pt component must be small non zero
  // because of the error (bug) in the WTrackParameter
  double tmp_xy=1.e-6;
  double tmp_z=0.01;
  double tmp_e = sqrt(2*tmp_xy*tmp_xy+tmp_z*tmp_z);
  HepLorentzVector pisr(tmp_xy,tmp_xy,tmp_z,tmp_e);
  double dphi = 1.e-3;  // 1mrad
  double dthe = 1.e-3;  // 1mrad
  double dE = 3.;        // 3 GeV
  WTrackParameter wISR(xorigin,pisr,dphi,dthe,dE);


  int ig[4] = {-1,-1,-1,-1};
  double chisq = 9999.;
  double chisq4C = 9999., chisq5C = 9999.;
  for(int i = 0; i < nGam-3; i++) {
    RecEmcShower *g1Trk =
         ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[i]))->emcShower();
    for(int j = i+1; j < nGam-2; j++) {
      RecEmcShower *g2Trk =
         ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[j]))->emcShower();
      for(int k = j+1; k < nGam-1; k++){
        RecEmcShower *g3Trk =
          ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[k]))->emcShower();
        for(int l = k+1; l < nGam; l++){
          RecEmcShower *g4Trk =
            ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[l]))->emcShower();

          kmfit->init();
          kmfit->AddTrack(0, wpip);
          kmfit->AddTrack(1, wpim);
          kmfit->AddTrack(2, 0.0, g1Trk);
          kmfit->AddTrack(3, 0.0, g2Trk);
          kmfit->AddTrack(4, 0.0, g3Trk);
          kmfit->AddTrack(5, 0.0, g4Trk);
#ifdef  add_ISR
          kmfit->AddTrack(6, wISR);
#endif
          kmfit->AddFourMomentum(0, ecms);
          bool oksq = kmfit->Fit();
          //           cout << " i,j,k,l="<< i << j << k << l << " oksq= " << oksq
          //                << " chi2= " << kmfit->chisq() << endl
          //                << "    P_isr= " << kmfit->pfit(6) << endl;
          if ( oksq ) {
            double chi2 = kmfit->chisq();
            // his1[12]->Fill(chi2);
            if ( chi2 < chisq ) {
              chisq = chi2;
              ig[0] = i;
              ig[1] = j;
              ig[2] = k;
              ig[3] = l;
            }
          }
        }//------------------------------------------End for(l)
      }//--------------------------------------------End for(k)
    }//-----------------------------------------------End for(j)
  } //-----------------------------------------------End for(i)

  his1[31]->Fill(chisq);
  //  if ( !Slct.GoodGammas() ) return false; // can not find 4  good gammas
 if ( ig[0]==-1 || ig[1]==-1 ||  ig[2]==-1 || ig[3]==-1) return false; // can not find 4  good gammas
 his1[107]->Fill(3);
 if(isMC && omega_best_mass) his1[280]->Fill(3);
 if(isMC){
    if(one_ISR) his1[300]->Fill(3);
    if(two_ISR) his1[301]->Fill(3);
    if(one_ISR && omega_best_mass) his1[302]->Fill(3);
    if(two_ISR && omega_best_mass) his1[303]->Fill(3);
  }
 chisq4C = chisq;
 if (chisq > 100) return false;
 his1[107]->Fill(4);
 if(isMC && omega_best_mass) his1[280]->Fill(4);
 if(isMC){
    if(one_ISR) his1[300]->Fill(4);
    if(two_ISR) his1[301]->Fill(4);
    if(one_ISR && omega_best_mass) his1[302]->Fill(4);
    if(two_ISR && omega_best_mass) his1[303]->Fill(4);
  }


 //P_missing for E_mis & M_mis
 Slct.P_missing = ecms - Slct.ppip - Slct.ppim - Slct.Pg[ig[0]] - Slct.Pg[ig[1]] - Slct.Pg[ig[2]] - Slct.Pg[ig[3]];


  float xfill[] = {
                kmfit->chisq(),P_isr.e(),P_isr.z(),P_isr.perp(),
                P_tot.e(),P_tot.rho(),P_sum.e(),P_sum.z(),P_sum.perp(),
                Slct.Ng,Slct.Eisr,Slct.Emax_isr,Slct.Th_isr,
                decCode,Slct.Pmiss.z(),Slct.Pmiss.perp(),Slct.Pmiss.m()
                  };
  m_tuple[0]->Fill( xfill );


 //============================================================================================================    Omega-Eta   ===================================================================
 RecEmcShower *g11Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[ig[0]]))->emcShower();
 RecEmcShower *g21Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[ig[1]]))->emcShower();
 RecEmcShower *g31Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[ig[2]]))->emcShower();
 RecEmcShower *g41Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(Slct.gammas[ig[3]]))->emcShower();

 kmfit->init();
 kmfit->AddTrack(0, wpip);
 kmfit->AddTrack(1, wpim);
 kmfit->AddTrack(2, 0.0, g11Trk);
 kmfit->AddTrack(3, 0.0, g21Trk);
 kmfit->AddTrack(4, 0.0, g31Trk);
 kmfit->AddTrack(5, 0.0, g41Trk);
#ifdef  add_ISR
 kmfit->AddTrack(6, wISR);
#endif
 kmfit->AddFourMomentum(0, ecms);
 bool oksq1 = kmfit->Fit();
 if ( !oksq1 ) {
   cerr << " Bad fit: somthing wrong!" << endl;
   ERROR_WARNING++;
   return false;
 }

  HepLorentzVector P_pip = kmfit->pfit(0);
  HepLorentzVector P_pim = kmfit->pfit(1);
  HepLorentzVector P_tot = P_pip + P_pim +
    kmfit->pfit(2)+kmfit->pfit(3)+kmfit->pfit(4)+kmfit->pfit(5);

  his1[251]->Fill(kmfit->chisq());
  his1[255]->Fill(P_tot.e());
  his1[256]->Fill(P_tot.rho());

#ifdef  add_ISR
  HepLorentzVector P_isr = kmfit->pfit(6);
  HepLorentzVector P_sum = P_tot + P_isr;
  his1[252]->Fill(P_isr.e());
  his1[253]->Fill(P_isr.z());
  his1[254]->Fill(P_isr.perp());
  his1[257]->Fill(P_sum.e());
  his1[258]->Fill(P_sum.rho());
  his2[251]->Fill(P_isr.e(),P_tot.e()+P_isr.e());
#endif

*/

  /*
  //------------------------------------------------------Missing Mass


  his1[260]->Fill(Slct.P_missing.z());
  his1[261]->Fill(Slct.P_missing.perp());
  his1[262]->Fill(Slct.P_missing.m());
  his1[276]->Fill(Slct.P_missing.m2());

  if(isMC){
    if(one_ISR_less50) his1[284]->Fill(Slct.P_missing.m());
    if(one_ISR_more50) his1[285]->Fill(Slct.P_missing.m());
    if(two_ISR_less50) his1[286]->Fill(Slct.P_missing.m());
    if(two_ISR_more50) his1[287]->Fill(Slct.P_missing.m());
  }

  if( isMC && decCode==67){
    his1[363]->Fill(Slct.P_missing.z());
    his1[364]->Fill(Slct.P_missing.perp());
    his1[365]->Fill(Slct.P_missing.m());
    his1[377]->Fill(Slct.P_missing.m2());
  }

  if( isMC && decCode!=67){
    his1[366]->Fill(Slct.P_missing.z());
    his1[367]->Fill(Slct.P_missing.perp());
    his1[368]->Fill(Slct.P_missing.m());
    his1[378]->Fill(Slct.P_missing.m2());
  }

  if( isMC && decCode==69){
    his1[263]->Fill(Slct.P_missing.z());
    his1[264]->Fill(Slct.P_missing.perp());
    his1[265]->Fill(Slct.P_missing.m());
    his1[277]->Fill(Slct.P_missing.m2());
  }

  if( isMC && decCode!=69){
    his1[266]->Fill(Slct.P_missing.z());
    his1[267]->Fill(Slct.P_missing.perp());
    his1[268]->Fill(Slct.P_missing.m());
    his1[278]->Fill(Slct.P_missing.m2());
  }
  //------------------------------------------------------------------
  */
/*

  HepLorentzVector P_wgg[6],P_gg[6], P_w, P_2pi, P_2pigg;
  P_wgg[0] = kmfit->pfit(2) + kmfit->pfit(3);
  P_gg[0] = kmfit->pfit(4) + kmfit->pfit(5);

  P_wgg[1] = kmfit->pfit(2) + kmfit->pfit(4);
  P_gg[1] = kmfit->pfit(3) + kmfit->pfit(5);

  P_wgg[2] = kmfit->pfit(2) + kmfit->pfit(5);
  P_gg[2] = kmfit->pfit(3) + kmfit->pfit(4);

  P_wgg[3] = kmfit->pfit(3) + kmfit->pfit(4);
  P_gg[3] = kmfit->pfit(2) + kmfit->pfit(5);

  P_wgg[4] = kmfit->pfit(3) + kmfit->pfit(5);
  P_gg[4] = kmfit->pfit(2) + kmfit->pfit(4);

  P_wgg[5] = kmfit->pfit(4) + kmfit->pfit(5);
  P_gg[5] = kmfit->pfit(2) + kmfit->pfit(3);
  double M_w = 0, M_gg = 0., M_tmp = 10., M_2pi, M_2pigg, Mw_gg;
  double sigma_eta = 0.008;//old 0.009
  double sigma_w = 0.010;//old 0.011
  double Llimit = 4.;
  double Ulimit = 7.;

  int Nbest=6;
  for(int i=0; i<6; i++){
    if(fabs(P_gg[i].m() - meta) < M_tmp) {
      M_tmp = fabs(P_gg[i].m() - meta);
      M_gg = P_gg[i].m();
      Nbest = i;
    }
  }

  if(Nbest!=6){

    P_w  = P_pip + P_pim + P_wgg[Nbest];
    M_w = P_w.m();
    Mw_gg = P_wgg[Nbest].m();


     if((Mw_gg  > 0.120) && (Mw_gg < 0.150)){
      his1[230]->Fill(M_gg);
    }


    if( (M_gg  > (meta - 3.*sigma_eta)) && (M_gg < (meta + 3.*sigma_eta)) ){
        his1[225]->Fill(M_gg);
        his1[226]->Fill(M_w);
        his1[229]->Fill(Mw_gg);

      if((Mw_gg  > 0.120) && (Mw_gg < 0.150)){
        his1[227]->Fill(M_gg);
        his1[228]->Fill(M_w);

        if( ((0.78265 - 3.*0.011) < M_w)  &&  (M_w < (0.78265 + 3.*0.011)) ) {//-----main sigma_W = 0.01

#ifdef  add_ISR
          if(fabs( Slct.P_missing.m() ) < 0.1) {
#endif

          if( isMC ){

            his1[103]->Fill(decCode);
            if(EtaTo2Gammas) his1[108]->Fill(5);
            if(omega_best_mass) {his1[280]->Fill(5);   his1[38]->Fill(Slct.SoverS0); }
            if(decCode==67)  his1[209]->Fill(5);
          }


           his1[232]->Fill(M_gg);
           his1[231]->Fill(Mw_gg);
           his1[107]->Fill(5);


#ifdef  add_ISR
      }
#endif



           //------------------------------Missing Mass-cut (for add ISR)
#ifdef  add_ISR
      if(fabs( Slct.P_missing.m() ) < 0.05) {  his1[373]->Fill(decCode);  his1[375]->Fill(5); }
      if(fabs( Slct.P_missing.m() ) < 0.06) {  his1[372]->Fill(decCode);  his1[375]->Fill(6); }
      if(fabs( Slct.P_missing.m() ) < 0.07) {  his1[374]->Fill(decCode);  his1[375]->Fill(7); }
      if(fabs( Slct.P_missing.m() ) < 0.08) {  his1[369]->Fill(decCode);  his1[375]->Fill(8); }
      if(fabs( Slct.P_missing.m() ) < 0.1)  {  his1[370]->Fill(decCode);  his1[375]->Fill(10);}
      if(fabs( Slct.P_missing.m() ) < 0.12) {  his1[371]->Fill(decCode);  his1[375]->Fill(12); }
#endif
*/

           /*

           //other photons
           double Emax = 0;
           int nOtherGammas = 0;
           for(int i = 0; i < nGam; i++) {
             if((i != ig[0])  && (i != ig[1]) && (i != ig[2]) && (i != ig[3]) ){
               nOtherGammas++;
               if(Emax < Slct.Pg[i].e()) Emax = Slct.Pg[i].e();
             }
           }

           his1[383]->Fill(Emax*1000.);
           his1[384]->Fill(nOtherGammas);

           if( isMC && (decCode==81 || decCode==135) ){
             his1[385]->Fill(Emax*1000.);
             his1[386]->Fill(nOtherGammas);
           }

           if( isMC && decCode==67){
             his1[379]->Fill(Emax*1000.);
             his1[380]->Fill(nOtherGammas);
           }

           if( isMC && decCode!=67){
             his1[381]->Fill(Emax*1000.);
             his1[382]->Fill(nOtherGammas);
           }
           //-------------------------------
           */

/*
        }//end if(w 3 sigma)
      }//end if(pi0 3 sigma)
    }//end if(eta 3 sigma)
  }//end if(Nbest)

 //==================================================================================================================================================================================================

 //------------------------------------------------------------SB w & eta-------------------------------

#ifdef  add_ISR
  if(fabs( Slct.P_missing.m() ) < 0.1) {
#endif


  if(Nbest!=6 && (Mw_gg  > 0.120) && (Mw_gg < 0.150) ){

    his2[185]->Fill(M_w, M_gg);

    //c -  3 Sigma-cut
    if( ((meta - 3.*sigma_eta) < M_gg) && (M_gg < (meta + 3.*sigma_eta))  &&  ((0.78265 - 3.*sigma_w) < M_w)  &&  (M_w < (0.78265 + 3.*sigma_w)) ){

      his1[107]->Fill(6);
      if(isMC &&  decCode==67)  his1[209]->Fill(6);
      if(isMC && omega_best_mass) his1[280]->Fill(6);
      his2[180]->Fill(M_w, M_gg);
      his2[183]->Fill(M_w, M_gg);

    }//end c - 3 Sigma-cut


    //xy-y
    if(  ((0.78265 - 3.*sigma_w) < M_w)  &&   (M_w < (0.78265 + 3.*sigma_w)) ){
      if(( ((meta - Ulimit*sigma_eta) < M_gg)  &&  (M_gg < (meta - Llimit*sigma_eta)) ) || (  ((meta + Llimit*sigma_eta) < M_gg)  &&  (M_gg < (meta + Ulimit*sigma_eta)) )){

        his2[181]->Fill(M_w, M_gg);
        his2[183]->Fill(M_w, M_gg);
        his1[107]->Fill(7);
        if(isMC &&  decCode==67)  his1[209]->Fill(7);
        if(isMC && omega_best_mass) his1[280]->Fill(7);

      }
    }
    //xy-x
    if( ((meta - 3.*sigma_eta) < M_gg)  && (M_gg < (meta + 3.*sigma_eta))  ){
      if(( ((0.78265 - Ulimit*sigma_w) < M_w)   &&  (M_w < (0.78265 - Llimit*sigma_w))  ) || ( ((0.78265 + Llimit*sigma_w) < M_w)   &&  (M_w < (0.78265 + Ulimit*sigma_w))  )){

        his2[181]->Fill(M_w, M_gg);
        his2[183]->Fill(M_w, M_gg);
        his1[107]->Fill(7);
        if(isMC &&  decCode==67)  his1[209]->Fill(7);
        if(isMC && omega_best_mass) his1[280]->Fill(7);

      }
    }

    //d-up
    if( ((0.78265 - Ulimit*sigma_w) < M_w)   &&  (M_w < (0.78265 - Llimit*sigma_w))   &&  ((meta + Llimit*sigma_eta) < M_gg)   &&  (M_gg < (meta + Ulimit*sigma_eta))  ){

      his2[182]->Fill(M_w, M_gg);
      his2[183]->Fill(M_w, M_gg);
      his1[107]->Fill(8);
      if(isMC &&  decCode==67)  his1[209]->Fill(8);
      if(isMC && omega_best_mass) his1[280]->Fill(8);

    }
    //d-down
    if( ((0.78265 + Llimit*sigma_w) < M_w) &&  (M_w < (0.78265 + Ulimit*sigma_w))   && ((meta - Ulimit*sigma_eta) < M_gg)  && (M_gg < (meta - Llimit*sigma_eta))  ){

      his2[182]->Fill(M_w, M_gg);
      his2[183]->Fill(M_w, M_gg);
      his1[107]->Fill(8);
      if(isMC &&  decCode==67)  his1[209]->Fill(8);
      if(isMC && omega_best_mass) his1[280]->Fill(8);

    }
    //d`-up
    if( ((0.78265 + Llimit*sigma_w) < M_w)   &&  (M_w < (0.78265 + Ulimit*sigma_w))   &&  ((meta + Llimit*sigma_eta) < M_gg)   &&   (M_gg < (meta + Ulimit*sigma_eta))  ){

      his2[182]->Fill(M_w, M_gg);
      his2[183]->Fill(M_w, M_gg);
      his1[107]->Fill(8);
      if(isMC &&  decCode==67)  his1[209]->Fill(8);
      if(isMC && omega_best_mass) his1[280]->Fill(8);

    }
    //d`-down
    if( ((0.78265 - Ulimit*sigma_w) < M_w)   &&  (M_w < (0.78265 - Llimit*sigma_w))  &&  ((meta - Ulimit*sigma_eta) < M_gg)   &&   (M_gg < (meta - Llimit*sigma_eta))   ){

      his2[182]->Fill(M_w, M_gg);
      his2[183]->Fill(M_w, M_gg);
      his1[107]->Fill(8);
      if(isMC &&  decCode==67)  his1[209]->Fill(8);
      if(isMC && omega_best_mass) his1[280]->Fill(8);

    }


  }//end if(pi0 3 sigma)


#ifdef  add_ISR
  }
#endif

*/

  return true;
}

//-----------------------------------------------------------------------------
void SelectOmegaEtaEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
  if ( selector->Verbose() ) cout << __func__ << "()" << endl;

//nef   his2[119]->Divide(his2[118],his2[117],1.,1.,"B");

  if ( ERROR_WARNING != 0 ) {
    cout << " ERROR_WARNING= " << ERROR_WARNING << " Check output!" << endl;
  }

  if ( init_selection ) {
    cout << " ===== INIT SELECTION =====" << endl;
  }
}
#ifdef __cplusplus
}
#endif
