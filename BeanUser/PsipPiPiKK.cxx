//==================================================================//
//                                                                  //
// PsipPiPiKK:                                                      //
// Study of the track reconstruction efficiency for pi and K        //
// in process:                                                      //
//              e+ e- -> Psi(2S) -> pi+ pi- J/Psi                   //
//                                           |-> pi+ pi- K+ K-      //
//                                                                  //
//==================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
// #include <TNtupleD.h>
#include <TTree.h>
#include <TDatabasePDG.h>

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
// #include "VertexFit/KinematicFit.h"
// #include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "EventTag/EventTagSvc.h"

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

// {{{1 Structure with job options
//--------------------------------------------------------------------
struct JobOption {
   int Hshift;                  // shift for histograms
   double Mrec_min, Mrec_max;   // Mrec(pi+pi-) range
   double M2Kwin;               // M^2(K) range
   double M2Ksb;                // M^2(K) side-band shift
   double Minv_min, Minv_max;   // Minv(pi+pi-K+K-) range
   double M2Piwin;              // M^2(pi) range
   double M2Pisb;               // M^2(pi) side-band shift

   JobOption( int hshift,
         double mrmin, double mrmax, double m2kw, double m2ksb,
         double mimin, double mimax, double m2piw, double m2pisb ) :
      Hshift{hshift},
      Mrec_min{mrmin}, Mrec_max{mrmax}, M2Kwin{m2kw},   M2Ksb{m2ksb},
      Minv_min{mimin}, Minv_max{mimax}, M2Piwin{m2piw}, M2Pisb{m2pisb}
   {}
};
//--------------------------------------------------------------------

// {{{1 Structure to save variables for a single event
//--------------------------------------------------------------------
struct PipPimKpKm {
   // Run-info
   int runNo;              // run-number
   int event;              // event-number
   HepLorentzVector LVcms; // Momentum in center of mass system
   Hep3Vector xorig;       // interaction point from DB

   // MC information:
   int decPsip, decJpsi; // decay codes for MC events
   int good;             // == 1 for good ppKK and 0 otherwise
   int tag_gam;          // tag gamma in decay chain
   int tag_pi0;          // tag pi0 in decay chain
   int tag_ppppKK;       // tag for Psi(2S) -> pi+ pi- pi+ pi- K+ K-

   // tracks
   vector<RecMdcKalTrack*> trk_Pi;  // pion tracks
   vector<RecMdcKalTrack*> trk_K;   // kaon tracks
   vector<RecMdcKalTrack*> trk_o;   // other tracks
   vector<RecEmcShower*>   trk_g;   // gamma tracks

   vector<HepLorentzVector> LVpi;   // pion momentums
   vector<HepLorentzVector> LVk;    // kaon momentums
   vector<Hep3Vector>        Vo;    // momentums of oth.

   vector<HepLorentzVector> LVg;    // photons momentum
   double Eg_sum;                   // sum of energy of photons
   double Eg_max;                   // maximal energy

   // recoil masses of pions closest to M(J/Psi)
   double Mrec;

   // invariant masses of pi+pi-K+K- closest to M(J/Psi)
   double Minv;
   HepLorentzVector LVjpsi; // momentum of the best combination
   vector<int> fpi;    // indexes of "free" pions (not used in LVjpsi

   PipPimKpKm() :
      decPsip{0}, decJpsi{0},
      good{0}, tag_gam{0}, tag_pi0{0}, tag_ppppKK{0},
      Eg_sum{0}, Eg_max{0},
      Mrec{0.}, Minv{0.} {}
};
//--------------------------------------------------------------------

// {{{1 Structure for root-trees
//--------------------------------------------------------------------
// Kaon reconstruction efficiency
struct Xnt_effK {
   int    Zk;           // charge of predicted Kaon
   double Ptk,Ck;       // Pt,cos(Theta) of predicted Kaon
   double dP,dTh;       // relations btw predicted and found Kaon
   int    Nk;           // number of Kaons in an event (1 or 2)
   int    fl;           // 0/1/2 -> predicted/found/pid(K) (remove?)
                        // * variables used in selection:
   double Mrec;         // Mrec(pi+pi-) to search for J/Psi
   double Mrk2;         // M^2_recoil( 2(pi+pi-) K+/- )
                        // * MC-variable
   int    good;         // 1 for 'true' K+K- 2(pi+pi-) event
};

// Pion reconstruction efficiency
struct Xnt_effPi {
   int    Zpi;          // charge of predicted Pion
   double Ptpi,Cpi;     // Pt,cos(Theta) of predicted Pion
   double dP,dTh;       // relations btw predicted and found Pion
   int    Npi;          // number of Pions in an event (4 or 3)
   int    fl;           // 0/1/2 -> predicted/found/pid(pi) (remove?)
                        // * variables used in selection:
   double MppKK;        // Minv(pi+pi-K+K-) to search for J/Psi
   double Mrpi2;        // M^2_recoil( (pi+pi-K+K-) pi+/- )
                        // * MC-variable
   int    good;         // 1 for 'true' K+K- 2(pi+pi-) event
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

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

// job options
static vector<JobOption*> Jobs;

// histograms
static vector<TH1*> hst;
// static vector<TNtupleD*> m_tuple;
static TTree* m_nt_effK = nullptr;
static Xnt_effK xnt_effK;
static TTree* m_nt_effPi = nullptr;
static Xnt_effPi xnt_effPi;

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

//--------------------------------------------------------------------
static void CloneHistograms(int Hmin, int Hmax) {
//--------------------------------------------------------------------
// clone histograms in range [Hmin, Hmax] for all Jobs

   char ch = 'A'; // add underscore letter at end of cloned histogram
   for ( const auto& job : Jobs ) {
      int sh = job->Hshift;
      if ( sh == 0 ) {
         continue;
      }
      for ( int ih = Hmin; ih <= Hmax; ih++ ) {
         int nh = ih + sh;
         if ( nh > int(hst.size()) ) {
            cout << " ERROR: CloneHistograms: nh= " << nh
               << " out of histo size: " << hst.size() << endl;
            cout << "        shift= " << sh << endl;
            cout << "        Hrange= [" << Hmin << ", " << Hmax
               << "]\n";
            exit(EXIT_FAILURE);
         }

         if ( hst[ih] ) {
            string hname = string(hst[ih]->GetName())
               + string("_") + ch;
            hst[nh] = static_cast<TH1*>(
                  hst[ih]->Clone(hname.c_str()) );
         }
      }
      ch += 1; // next letter
   }
}

// {{{1 StartJob, book histograms
//--------------------------------------------------------------------
void PsipPiPiKKStartJob(ReadDst* selector) {
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
   }

   // set job options:
   //   Hshift,                                 :for CloneHistograms()
   //   Mrec_min, Mrec_max, M2Kwin, M2Ksb       :see KTrkRecEff()
   //   Minv_min, Minv_max, M2Piwin, M2Pisb     :see PiTrkRecEff()
   Jobs.push_back( new JobOption( 0, // must be 0 for first !
                                  3.095,3.099, 0.025, 0.06,
                                  3.087,3.107, 0.01,  0.025 ) );
   Jobs.push_back( new JobOption( 400,
                                  3.092,3.102, 0.06,  0.07,
                                  3.080,3.114, 0.02,  0.025 ) );
   // Note:
   // M2Ksb and M2Pisb are not used anywhere
   // mk^2  = 0.244  sigma = 0.015 (18)
   // mpi^2 = 0.0195 sigma = 0.0045 (6)
   // 0   -> mk2 in [.219;.269] mpi2 in [ .0095; .0295]
   // 400 -> mk2 in [.194;.294] mpi2 in [-.0005; .0395]


   hst.resize(1000,nullptr);
   // m_tuple.resize(10,nullptr);

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
   // nef: psi' decays
   m_EventTagSvc->setUserDecayTabsFile(
      selector->AbsPath( "BeanUser/mypsip_dec.codes" ) );

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

   // ChargedTracks:
   hst[1] = new TH1D("S1a_Npi","1A) N(#pi)", 10,-0.5,9.5);
   hst[2] = new TH1D("S1a_Nk", "1A) N(K)", 10,-0.5,9.5);
   hst[3] = new TH1D("S1a_No", "1A) N(oth)", 10,-0.5,9.5);
   hst[4] = new TH1D("S1a_Zpi","1A) Z(#pi)", 9,-4.5,4.5);
   hst[5] = new TH1D("S1a_Zk", "1A) Z(K)", 9,-4.5,4.5);
   hst[6] = new TH1D("S1a_Zo", "1A) Z(oth)", 9,-4.5,4.5);
   hst[8] = new TH1D("S1a_Ncl","1A) N cloned", 10,0.5,+10.5);
   hst[9] = new TH1D("S1a_Nt", "1A) Ntot", 10,-0.5,9.5);

   // NeutralTracks:
   hst[11] = new TH1D("S1b_Ng", "1B) N(#gamma)", 10,-0.5,9.5);
   hst[12] = new TH1D("S1b_Egm","1B) Eg(max)", 100,0.,1.);
   hst[13] = new TH1D("S1b_Egs","1B) Sum E(#gamma)", 100,0.,1.);
   hst[14] = new TH1D("S1b_Mg2","1B) M^2(#gamma#gamma)", 200,0.,0.04);
   hst[19] = new TH1D("S1b_Mg2a","1B) M^2(#gamma#gamma) all",
         500,0.,2.);
   hst[15] = new TH1D("S1b_Npi0","1B) Npi0", 10,-0.5,9.5);

//    hst[16] = new TH2D("S1b_mapB","1B)barrel Z vs phi",
//                       270,-135.,+135., 360,-180.,180.);
//    hst[17] = new TH2D("S1b_mapC1","1B)endcap Z<0 rho vs phi",
//                       40, 50.,90., 360,-180.,180.);
//    hst[18] = new TH2D("S1b_mapC2","1B)endcap Z>0 rho vs phi",
//                       40, 50.,90., 360,-180.,180.);

   // MrecPiPi:
   hst[21] = new TH1D("S2_Mrec","2) Mrec(#pi^{+}#pi^{-})",
         100,3.085,3.11);
   hst[22] = new TH1D("S2_MrecZ","2) Mrec(#pi^{+}#pi^{-})",
         100,3.,3.2);
   hst[23] = new TH1D("S2_cos", "2) cos(Theta_pi+_pi-)",
         100,-1.0,1.0);
   hst[24] = new TH1D("S2_invM","2) Minv(pi+ pi-)",
         500,0.25,0.75);
   hst[27] = new TH1D("S2_Mrb", "2) best Mrec(pi+pi-)",
         100,3.085,3.11);

   // KTrkRecEff:
   hst[31] = new TH1D("S3_M2K", "3) M^2rec(K)", 200,0.125,0.365);
   hst[32] = new TH1D("S3_dPK", "3) #deltaP(K)", 200,-0.2,0.2);
   hst[33] = new TH1D("S3_dThK","3) #delta#Theta(K)", 100,0.,20.);
   hst[34] = new TH2D("S3_2DK", "3) dP vs dTh (K);#Theta",
                       50,0.,20., 50,-0.1,0.1);
   hst[35] = new TH1D("S3_Pkp", "3) P(K^{+})", 200,0.,2.);
   hst[36] = new TH1D("S3_Pkm", "3) P(K^{-})", 200,0.,2.);
   hst[37] = new TH1D("S3_Ptkp","3) Pt(K^{+})", 200,0.,2.);
   hst[38] = new TH1D("S3_Ptkm","3) Pt(K^{-})", 200,0.,2.);
   hst[39] = new TH1D("S3_Ckp", "3) cos(#Theta) of K+", 100,-1.,1.);
   hst[40] = new TH1D("S3_Ckm", "3) cos(#Theta) of K-", 100,-1.,1.);

   hst[41] = new TH2D("S3_Kp5","Pt vs cos(Theta) for K+ 5tracks",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[42] = new TH2D("S3_Km5","Pt vs cos(Theta) for K- 5tracks",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[43] = new TH2D("S3_Kp6","Pt vs cos(Theta) for K+ 6tracks",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[44] = new TH2D("S3_Km6","Pt vs cos(Theta) for K- 6tracks",
         18,-0.9,0.9, 13,0.1,1.4);

   hst[45] = new TH2D("S6_Kp0","Pt vs cos(Theta) for K+ PID no",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[46] = new TH2D("S6_Km0","Pt vs cos(Theta) for K- PID no",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[47] = new TH2D("S6_Kp1","Pt vs cos(Theta) for K+ PID ok",
         18,-0.9,0.9, 13,0.1,1.4);
   hst[48] = new TH2D("S6_Km1","Pt vs cos(Theta) for K- PID ok",
         18,-0.9,0.9, 13,0.1,1.4);
   CloneHistograms(31, 49);

   // MinvPiPiKK:
   hst[51] = new TH1D("S4_MppKK", "4) M(pi+pi-K+K-)", 100,3.05,3.14);
   hst[52] = new TH1D("S4_MppKKZ","4) M(pi+pi-K+K-)", 100,3.,3.2);
   hst[54] = new TH1D("S4_MppKKb","4) best M(pi+pi-K+K-)",
         100,3.05,3.14);
   hst[55] = new TH1D("S4_Ppi", "4) P(#pi) non J/Psi decay",
         100,0.,1.0);

   // PiTrkRecEff:
   hst[61] = new TH1D("S5_M2Pi", "5) M^2rec(#pi)", 200,-0.01,0.05);
   hst[62] = new TH1D("S5_Ppl",  "5) P(#pi) lost", 100,0.,1.0);
   hst[63] = new TH1D("S5_dPpi", "5) #deltaP(pi)", 200,-0.2,0.2);
   hst[64] = new TH1D("S5_dThPi","5) #delta#Theta(pi)", 100,0.,20.);
   hst[65] = new TH2D("S5_2DPi", "5) dP vs dTh (pi);#Theta",
                       50,0.,20., 50,-0.1,0.1);
   hst[66] = new TH1D("S5_Ptpip","5) Pt(#pi^{+})", 100,0.,0.5);
   hst[67] = new TH1D("S5_Ptpim","5) Pt(#pi^{-})", 100,0.,0.5);
   hst[68] = new TH1D("S5_Cpip", "5) cos(#Theta) of #pi+",
         100,-1.,1.);
   hst[69] = new TH1D("S5_Cpim", "5) cos(#Theta) of #pi-",
         100,-1.,1.);

   hst[71] = new TH2D("S5_Pip5","Pt vs cos(Theta) for pi+ 5tracks",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[72] = new TH2D("S5_Pim5","Pt vs cos(Theta) for pi- 5tracks",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[73] = new TH2D("S5_Pip6","Pt vs cos(Theta) for pi+ 6tracks",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[74] = new TH2D("S5_Pim6","Pt vs cos(Theta) for pi- 6tracks",
         18,-0.9,0.9, 7,0.05,0.4);

   hst[75] = new TH2D("S6_Pip0","Pt vs cos(Theta) for pi+ PID no",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[76] = new TH2D("S6_Pim0","Pt vs cos(Theta) for pi- PID no",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[77] = new TH2D("S6_Pip1","Pt vs cos(Theta) for pi+ PID ok",
         18,-0.9,0.9, 7,0.05,0.4);
   hst[78] = new TH2D("S6_Pim1","Pt vs cos(Theta) for pi- PID ok",
         18,-0.9,0.9, 7,0.05,0.4);
   CloneHistograms(61, 79);

   // Monte Carlo histograms:
   // FillHistoMC:
   hst[101] = new TH1D("MC_decPsi", "dec Psi(2S) all",
         256,-0.5,255.5);
   hst[105] = new TH1D("MC_decJpsi","dec J/psi all",
         256,-0.5,255.5);

   // Psi(2S) -> J/Psi pi+ pi-
   hst[111] = new TH1D("MC_PsipPt", "Pt(#Psi(2S))", 100,0.,1.);
   hst[112] = new TH1D("MC_PsipC", "cos(#Theta) of Psi(2S)",
         100,-1.,1.);
   hst[116] = new TH1D("MC_JpsiP", "P(J/#Psi)", 100,0.,1.);
   hst[117] = new TH1D("MC_JpsiPt","Pt(J/#Psi)", 100,0.,1.);
   hst[118] = new TH1D("MC_JpsiC", "cos(#Theta) of J/#Psi",
         100,-1.,1.);
   hst[119] = new TH1D("MC_PipP", "P(#pi^{+})", 200,0.,1.);
   hst[120] = new TH1D("MC_PipC", "cos(#Theta) of #pi^{+}",
         100,-1.,1.);
   hst[121] = new TH1D("MC_PimP", "P(#pi^{-})", 200,0.,1.);
   hst[122] = new TH1D("MC_PimC", "cos(#Theta) of #pi^{-}",
         100,-1.,1.);
   hst[124] = new TH1D("MC_PgE", "E(#gamma) [-22]", 100,0.,0.2);
   hst[125] = new TH1D("MC_PgC", "cos(#Theta) of #gamma [-22]",
         100,-1.,1.);

   // ISR photons
   hst[126] = new TH1D("MC_ISRE", "E(#gamma) ISR", 100,0.,0.05);
   hst[127] = new TH1D("MC_ISRC", "cos(#Theta) of #gamma ISR",
         100,-1.,1.);

   // pi0 in decay chain
   hst[131] = new TH1D("MC_pi0N",  "N(#pi^{0})", 10,-0.5,9.5);
   hst[132] = new TH1D("MC_pi0dc1","pi0 dec Psi(2S)",256,-0.5,255.5);
   hst[133] = new TH1D("MC_pi0dc2","pi0 dec J/psi",256,-0.5,255.5);
   hst[134] = new TH1D("MC_nopi0", "nopi0 dec Psi(2S)",
         256,-0.5,255.5);

   // gamma in decay chain
   hst[141] = new TH1D("MC_gN", "N(#gamma)", 10,0.5,10.5);
   hst[142] = new TH1D("MC_gE", "E(#gamma) all", 100,0.,0.5);
   hst[143] = new TH1D("MC_gE1","E(#gamma) for 1#gamma", 100,0.,0.2);
   hst[144] = new TH1D("MC_gC1","cos(#Theta) for 1#gamma",
         100,-1.,1.);

   // Psi(2S) -> pi+ pi- pi+ pi- K+ K-
   hst[151] = new TH1D("MC_dPx","tag_ppppKK dPx", 100,-0.005,0.005);
   hst[152] = new TH1D("MC_dPy","tag_ppppKK dPy", 100,-0.005,0.005);
   hst[153] = new TH1D("MC_dPz","tag_ppppKK dPz", 100,-0.005,0.005);
   hst[154] = new TH1D("MC_4p2kdc1", "4p2k good dec Psi(2S)",
         256,-0.5,255.5);
   hst[155] = new TH1D("MC_4p2kdc2", "4p2k good dec J/psi",
         256,-0.5,255.5);
   hst[156] = new TH1D("MC_4p2kdc3", "4p2k bad dec Psi(2S)",
         256,-0.5,255.5);
   hst[157] = new TH1D("MC_4p2kdc4", "4p2k bad dec J/psi",
         256,-0.5,255.5);

   // 2(pi+pi-) K+K-
   hst[161] = new TH1D("MCG_Ppip","Pt of #pi^{+}", 200,0.,2.);
   hst[162] = new TH1D("MCG_Ppim","Pt of #pi^{-}", 200,0.,2.);
   hst[163] = new TH1D("MCG_Pkp", "Pt of K^{+}", 200,0.,2.);
   hst[164] = new TH1D("MCG_Pkm", "Pt of K^{-}", 200,0.,2.);
   hst[165] = new TH1D("MCG_Cpip","cos(#Theta) of #pi^{+}",
         100,-1.,1.);
   hst[166] = new TH1D("MCG_Cpim","cos(#Theta) of #pi^{-}",
         100,-1.,1.);
   hst[167] = new TH1D("MCG_Ckp", "cos(#Theta) of K^{+}",
         100,-1.,1.);
   hst[168] = new TH1D("MCG_Ckm", "cos(#Theta) of K^{-}",
         100,-1.,1.);

   // MC ChargedTracks:
   hst[201] = new TH1D("MC_1a1","1A)fin bad(0)/good(1)", 2,-0.5,1.5);
   hst[202] = new TH1D("MC_1a2","1A) dec Psi(2S)",256,-0.5,255.5);
   hst[203] = new TH1D("MC_1a3","1A) dec J/Psi in pipiJ/Psi",
         256,-0.5,255.5);

   // MC NeutralTracks:
   hst[211] = new TH1D("MC_NgF", "1B) Ng F", 10,-0.5,9.5);
   hst[212] = new TH1D("MC_NgT", "1B) Ng T", 10,-0.5,9.5);
   hst[213] = new TH1D("MC_EgmF","1B) Eg(max) F", 100,0.,1.);
   hst[214] = new TH1D("MC_EgmT","1B) Eg(max) T", 100,0.,1.);
   hst[215] = new TH1D("MC_EgsF","1B) Eg(sum) F", 100,0.,1.);
   hst[216] = new TH1D("MC_EgsT","1B) Eg(sum) T", 100,0.,1.);
   hst[217] = new TH1D("MC_M2ggF","1B) M^2(gg) F", 200,0.,0.04);
   hst[218] = new TH1D("MC_M2ggT","1B) M^2(gg) T", 200,0.,0.04);
   hst[219] = new TH1D("MC_Npi0F","1B) Npi0 F", 10,-0.5,9.5);
   hst[220] = new TH1D("MC_Npi0T","1B) Npi0 T", 10,-0.5,9.5);
   hst[221] = new TH1D("MC_1b1","1B)pi0 rejected bad/good",
         2,-0.5,1.5);
   hst[222] = new TH1D("MC_1b2","1B) bad(0)/good(1)", 2,-0.5,1.5);
   hst[223] = new TH1D("MC_1b3","1B) dec Psi(2S)",256,-0.5,255.5);
   hst[224] = new TH1D("MC_1b4","1B)dec J/Psi in pipiJ/Psi",
         256,-0.5,255.5);

   // MC MrecPiPi:
   hst[231] = new TH1D("MC_MrecF","2)nonpipiJ/Psi Mrec",
         100,3.085,3.11);
   hst[232] = new TH1D("MC_MrecT","2)   pipiJ/Psi Mrec",
         100,3.085,3.11);

   // MC KTrkRecEff:
   hst[241] = new TH1D("MC_M2KF", "3) M^2rec(K) F", 200,0.125,0.365);
   hst[242] = new TH1D("MC_M2KT", "3) M^2rec(K) T", 200,0.125,0.365);
   hst[243] = new TH1D("MC_dPK_F","3) #deltaP(K) F", 200,-0.2,0.2);
   hst[244] = new TH1D("MC_dPK_T","3) #deltaP(K) T", 200,-0.2,0.2);
   hst[245] = new TH1D("MC_dThK_F","3) #delta#Theta(K) F",
         100,0.,20.);
   hst[246] = new TH1D("MC_dThK_T","3) #delta#Theta(K) T",
         100,0.,20.);
   hst[247] = new TH2D("MC_2DK_F", "3) dP vs dTh (K) F;#Theta",
         50,0.,20., 50,-0.1,0.1);
   hst[248] = new TH2D("MC_2DK_T", "3) dP vs dTh (K) T;#Theta",
         50,0.,20., 50,-0.1,0.1);
   hst[249] = new TH1D("MC_3_nm6","3) No mach for 6 tracks",
         2,-0.5,1.5);

   hst[251] = new TH1D("MC_3_1","3) bad(0)/good(1)", 2,-0.5,1.5);
   hst[252] = new TH1D("MC_3_2","3) dec Psi(2S)",256,-0.5,255.5);
   hst[253] = new TH1D("MC_3_3F","3) dec J/psi F",256,-0.5,255.5);
   hst[254] = new TH1D("MC_3_3T","3) dec J/psi T",256,-0.5,255.5);

   // re-weighted histograms for MC
   // hst[255] = new TH1D("MC_3_WK","all weights for kaons",
         // 200, 0.9,1.1);
   // hst[256] = new TH2D("S6_Kp1W","Pt vs cos(Theta) for K+ PID ok",
         // 18,-0.9,0.9, 13,0.1,1.4);
   // hst[257] = new TH2D("S6_Km1W","Pt vs cos(Theta) for K- PID ok",
         // 18,-0.9,0.9, 13,0.1,1.4);
   CloneHistograms(241, 259);

   // MC MinvPiPiKK:
   hst[261] = new TH1D("MC_MppKK_F","4) M(pi+pi-K+K-) F",
         100,3.05,3.14);
   hst[262] = new TH1D("MC_MppKK_T","4) M(pi+pi-K+K-) T",
         100,3.05,3.14);
   hst[263] = new TH1D("MC_PpiF", "4) P(#pi) non J/Psi decay F",
         100,0.,1.0);
   hst[264] = new TH1D("MC_PpiT", "4) P(#pi) non J/Psi decay T",
         100,0.,1.0);

   // MC PiTrkRecEff:
   hst[271] = new TH1D("MC_M2piF", "5) M^2rec(pi) F",
         200,-0.01,0.05);
   hst[272] = new TH1D("MC_M2piT", "5) M^2rec(pi) T",
         200,-0.01,0.05);
   hst[273] = new TH1D("MC_PplF", "5) P(#pi) lost F", 100,0.,1.0);
   hst[274] = new TH1D("MC_PplT", "5) P(#pi) lost T", 100,0.,1.0);

   hst[275] = new TH1D("MC_dPpi_F","5) #deltaP(pi) F", 200,-0.1,0.1);
   hst[276] = new TH1D("MC_dPpi_T","5) #deltaP(pi) T", 200,-0.1,0.1);
   hst[277] = new TH1D("MC_dThPi_F","5) #delta#Theta(pi) F",
         100,0.,20.);
   hst[278] = new TH1D("MC_dThPi_T","5) #delta#Theta(pi) T",
         100,0.,20.);
   hst[279] = new TH2D("MC_2DPi_F", "5) dP vs dTh (pi) F;#Theta",
         50,0.,20., 50,-0.1,0.1);
   hst[280] = new TH2D("MC_2DPi_T", "5) dP vs dTh (pi) T;#Theta",
         50,0.,20., 50,-0.1,0.1);
   hst[281] = new TH1D("MC_5_nm6","5) No mach for 6 tracks",
         2,-0.5,1.5);

   hst[282] = new TH1D("MC_5_1","5) bad(0)/good(1)", 2,-0.5,1.5);
   hst[283] = new TH1D("MC_5_2","5) dec Psi(2S)",256,-0.5,255.5);
   hst[284] = new TH1D("MC_5_3F","5) dec J/Psi F",256,-0.5,255.5);
   hst[285] = new TH1D("MC_5_3T","5) dec J/Psi T",256,-0.5,255.5);

   // re-weighted histograms for MC
   // hst[286] = new TH1D("MC_5_WPi","all weights for pions",
         // 200, 0.9,1.1);
   // hst[287] = new TH2D("S6_Pip1W","Pt vs cos(Theta) for pi+ PID ok",
         // 18,-0.9,0.9, 7,0.05,0.4);
   // hst[288] = new TH2D("S6_Pim1W","Pt vs cos(Theta) for pi- PID ok",
         // 18,-0.9,0.9, 7,0.05,0.4);
   CloneHistograms(271, 289);

   // Tree for Kaon efficiency study
   m_nt_effK = new TTree("eff_K","K reconstruction efficiency");
   m_nt_effK->Branch("Zk"   , &xnt_effK.Zk   );
   m_nt_effK->Branch("Ptk"  , &xnt_effK.Ptk  );
   m_nt_effK->Branch("Ck"   , &xnt_effK.Ck   );
   m_nt_effK->Branch("dP"   , &xnt_effK.dP   );
   m_nt_effK->Branch("dTh"  , &xnt_effK.dTh  );
   m_nt_effK->Branch("Nk"   , &xnt_effK.Nk   );
   m_nt_effK->Branch("fl"   , &xnt_effK.fl   );

   m_nt_effK->Branch("Mrec" , &xnt_effK.Mrec );
   m_nt_effK->Branch("Mrk2" , &xnt_effK.Mrk2 );

   m_nt_effK->Branch("good" , &xnt_effK.good );


   // Tree for Pion efficiency study
   m_nt_effPi = new TTree("eff_Pi","Pi reconstruction efficiency");
   m_nt_effPi->Branch("Zpi"   , &xnt_effPi.Zpi   );
   m_nt_effPi->Branch("Ptpi"  , &xnt_effPi.Ptpi  );
   m_nt_effPi->Branch("Cpi"   , &xnt_effPi.Cpi   );
   m_nt_effPi->Branch("dP"    , &xnt_effPi.dP    );
   m_nt_effPi->Branch("dTh"   , &xnt_effPi.dTh   );
   m_nt_effPi->Branch("Npi"   , &xnt_effPi.Npi   );
   m_nt_effPi->Branch("fl"    , &xnt_effPi.fl    );

   m_nt_effPi->Branch("MppKK" , &xnt_effPi.MppKK );
   m_nt_effPi->Branch("Mrpi2" , &xnt_effPi.Mrpi2 );

   m_nt_effPi->Branch("good"  , &xnt_effPi.good  );


   // register in selector to save in given directory
   const char* SaveDir = "PipPimKpKm";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
   // VecObj ntuples(m_tuple.begin(),m_tuple.end());
   // selector->RegInDir(ntuples,SaveDir);

   VecObj Vreg;
   Vreg.push_back(m_nt_effK);
   Vreg.push_back(m_nt_effPi);
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
      exit(EXIT_FAILURE);
   }

   save_runNo = runNo;
   return xorigin;
}

// {{{1 FillHistoMC()
//--------------------------------------------------------------------
static void FillHistoMC(const ReadDst* selector, PipPimKpKm& ppKK) {
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
   ppKK.decPsip = decPsip;
   hst[101]->Fill(ppKK.decPsip);

   // decay code of J/Psi from Psi(2S) -> pi+ pi- J/Psi
   ppKK.decJpsi = 0;
   if( decPsip == 64 ) {
      int decJpsi = ( evTag >> 16 ) & 0xFF;
      ppKK.decJpsi = decJpsi;
      hst[105]->Fill(ppKK.decJpsi);
   }

   int idx_psip =-1;
   set<int> idx_ign;
   Hep3Vector Ppsip;
   Hep3Vector Pisr;              // ISR photons
   Hep3Vector Pgs;               // sum other gammas
   vector<Hep3Vector> Psave(6);  // save particles of interest

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      int p_pdg = part->getParticleID();
      int p_idx = part->getTrackIndex();
      int p_mother = part->getMother();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      // skip decay chains from "final" particles
      if ( idx_ign.find(p_mother) != idx_ign.end() ) {
         idx_ign.insert( p_idx );
         continue; // go to next particle
      }

      if ( p_pdg == 100443 ) { // Psi(2S)
         idx_psip = p_idx;
         Ppsip = Vp;
         hst[111]->Fill(Vp.rho());
         hst[112]->Fill(Vp.cosTheta());
         continue;
      }

      // fill tag_gam && tag_pi0
      if ( p_pdg == -22 || p_pdg == 22 ) {
         idx_ign.insert( p_idx );
         if ( idx_psip < 0 ) { // befor Psi(2S) -> ISR photons
            hst[126]->Fill(Vp.mag());
            hst[127]->Fill(Vp.cosTheta());
            Pisr += Vp;
         } else {
            Pgs += Vp;
            ppKK.tag_gam += 1;
         }
      }
      if ( p_pdg == 111 ) { // tag any pi0
         idx_ign.insert( p_idx );
         ppKK.tag_pi0 += 1;
         continue;
      }

      // fill histograms only
      if ( ppKK.decPsip == 64 ) {          // Psi(2S) -> pi+ pi- J/Psi
         if ( p_mother == idx_psip ) {     // decays of Psi(2S)
            if ( p_pdg == 443 ) {               // J/Psi
               hst[116]->Fill(Vp.mag());
               hst[117]->Fill(Vp.rho());
               hst[118]->Fill(Vp.cosTheta());
            } else if ( p_pdg == 211 ) {        // pi+
               hst[119]->Fill(Vp.mag());
               hst[120]->Fill(Vp.cosTheta());
            } else if ( p_pdg == -211 ) {       // pi-
               hst[121]->Fill(Vp.mag());
               hst[122]->Fill(Vp.cosTheta());
            } else if ( p_pdg == -22 ) {        // gamma
               hst[124]->Fill(Vp.mag());
               hst[125]->Fill(Vp.cosTheta());
            }
         }
      }

      // tag 2(pi+ pi-) K+ K- final state
      if ( p_pdg == 211 ) {             // pi+
         idx_ign.insert( p_idx );
         ppKK.tag_ppppKK += 1;
         Psave[ppKK.tag_ppppKK%2] = Vp; // 0,1
      } else if ( p_pdg == -211 ) {     // pi-
         idx_ign.insert( p_idx );
         ppKK.tag_ppppKK += 10;
         Psave[2+(ppKK.tag_ppppKK/10)%2] = Vp; // 2,3
      } else if ( p_pdg == 321 ) {      // K+
         idx_ign.insert( p_idx );
         ppKK.tag_ppppKK += 100;
         Psave[4] = Vp;
      } else if ( p_pdg == -321 ) {     // K-
         idx_ign.insert( p_idx );
         ppKK.tag_ppppKK += 1000;
         Psave[5] = Vp;
      }

   } // end of while

   // pi0 in final state:
   hst[131]->Fill(ppKK.tag_pi0);
   if ( ppKK.tag_pi0 > 0 ) {
      hst[132]->Fill(ppKK.decPsip);
      if ( ppKK.decPsip == 64 ) {
         hst[133]->Fill(ppKK.decJpsi);
      }
      return;
   }
   hst[134]->Fill(ppKK.decPsip);

   // gamma in final state:
   hst[141]->Fill(ppKK.tag_gam);
   double Eg = Pgs.mag();
   if ( ppKK.tag_gam > 0 ) {
      hst[142]->Fill(Eg);
      if ( ppKK.tag_gam == 1 ) {
         hst[143]->Fill(Eg);
         hst[144]->Fill(Pgs.cosTheta());
      }
   }

   // accept one gamma <25 MeV
   if ( ppKK.tag_gam > 1 || Eg > 0.025 ) {
      return;
   }

   // tag Psi(2S) -> pi+ pi- pi+ pi- K+ K-
   if ( ppKK.tag_ppppKK == 1122 ) {
      const double maxdp = 1e-4;
      Hep3Vector PppppKK = Ppsip;
      for ( auto p : Psave ) {
         PppppKK -= p;
      }
      PppppKK -= Pgs;

      hst[151]->Fill( PppppKK.x() );
      hst[152]->Fill( PppppKK.y() );
      hst[153]->Fill( PppppKK.z() );

      if (    fabs( PppppKK.x() ) < maxdp
              && fabs( PppppKK.y() ) < maxdp
              && fabs( PppppKK.z() ) < maxdp
         ) {
         ppKK.good = 1;
         hst[154]->Fill(ppKK.decPsip);
         if ( ppKK.decPsip == 64 ) {
            hst[155]->Fill(ppKK.decJpsi);
         }

         // fill histograms
         hst[161]->Fill(Psave[0].perp());
         hst[161]->Fill(Psave[1].perp());
         hst[162]->Fill(Psave[2].perp());
         hst[162]->Fill(Psave[3].perp());
         hst[163]->Fill(Psave[4].perp());
         hst[164]->Fill(Psave[5].perp());

         hst[165]->Fill(Psave[0].cosTheta());
         hst[165]->Fill(Psave[1].cosTheta());
         hst[166]->Fill(Psave[2].cosTheta());
         hst[166]->Fill(Psave[3].cosTheta());
         hst[167]->Fill(Psave[4].cosTheta());
         hst[168]->Fill(Psave[5].cosTheta());
      } else {
         hst[156]->Fill(ppKK.decPsip);
         if ( ppKK.decPsip == 64 ) {
            hst[157]->Fill(ppKK.decJpsi);
         }
         ppKK.tag_ppppKK *= 1;
      }
   } // end if tag

}

// {{{1 1A) select pions, kaons and other charged tracks
//--------------------------------------------------------------------
static bool ChargedTracks(ReadDst* selector, PipPimKpKm& ppKK) {
//--------------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
   // static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   ParticleID* pid = ParticleID::instance();

   int Ncl  = 0;   // number of cloned tracks
   int Nzpi = 0;   // charge of pions
   int Nzk  = 0;   // charge of kaons
   int Nzo  = 0;   // charge of other tracks

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
         Ncl++;
         continue;
      }

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC rec.
      HepPoint3D IP(ppKK.xorig[0],ppKK.xorig[1],ppKK.xorig[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      // the nearest distance to IP
      double Rvxy0 = vecipa[0]; // in xy plane
      double Rvz0  = vecipa[3]; // in z direction

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
      if ( std::isnan(mdcKalTrk->px())
            || std::isnan(mdcKalTrk->py())
            || std::isnan(mdcKalTrk->pz()) ) {
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

      pid->calculate(ppKK.runNo);
      int flag_pik = 0; // 1 - pi; 2 - K; 0 - other
      if ( pid->IsPidInfoValid() ) {
         // check that track is pion:
         if ( (pid->probPion() > 0.001)             &&
               (pid->probPion() > pid->probKaon())   &&
               (pid->probPion() > pid->probProton())
            ) {
            flag_pik = 1;
            mdcKalTrk->setPidType(RecMdcKalTrack::pion);

            ppKK.trk_Pi.push_back(mdcKalTrk);
            Nzpi += mdcKalTrk->charge();
            Hep3Vector tmp(
                  mdcKalTrk->px(),mdcKalTrk->py(),mdcKalTrk->pz()
                  );
            ppKK.LVpi.push_back(
                  HepLorentzVector( tmp, sqrt(tmp.mag2()+SQ(mpi)) )
                  );
         } else
            // check that track is kaon:
            if ( (pid->probKaon() > 0.001)             &&
                  (pid->probKaon() > pid->probPion())   &&
                  (pid->probKaon() > pid->probProton())
               ) {
               flag_pik = 2;
               mdcKalTrk->setPidType(RecMdcKalTrack::kaon);

               ppKK.trk_K.push_back(mdcKalTrk);
               Nzk += mdcKalTrk->charge();
               Hep3Vector tmp(
                     mdcKalTrk->px(),mdcKalTrk->py(),mdcKalTrk->pz()
                     );
               ppKK.LVk.push_back(
                     HepLorentzVector( tmp, sqrt(tmp.mag2()+SQ(mk)) )
                     );
            }
      } // end IsPidInfoValid

      // all other tracks
      if ( flag_pik == 0 ) {
         mdcKalTrk->setPidType(RecMdcKalTrack::pion);
         ppKK.trk_o.push_back(mdcKalTrk);
         Nzo += mdcKalTrk->charge();
         ppKK.Vo.push_back( Hep3Vector(
                  mdcKalTrk->px(),mdcKalTrk->py(),mdcKalTrk->pz()
                  ));
      }
   } // end for charged tracks loop

   int Npi = ppKK.trk_Pi.size();
   int Nk  = ppKK.trk_K.size();
   int No  = ppKK.trk_o.size();
   int Ntot= Npi + Nk + No;
   hst[1]->Fill(Npi);
   hst[2]->Fill(Nk);
   hst[3]->Fill(No);
   hst[8]->Fill(Ncl);

   // select (3pi + 2K) or (4pi + 1K) or (4pi + 2K)
   if ( Npi < 3 || Npi > 4 ) {
      return false;
   }
   if ( Nk < 1 || Nk > 2 ) {
      return false;
   }
   if ( Npi+Nk < 5 || Npi+Nk > 6 ) {
      return false;
   }
   if ( Ntot < 5 || Ntot > 6 ) {
      return false;
   }

   hst[4]->Fill(Nzpi);
   hst[5]->Fill(Nzk);
   hst[6]->Fill(Nzo);

   // select good charge: 0 for 6 tracks and +/-1 for 5 tracks
   if ( abs(Nzpi) > 1 ) {
      return false;
   }
   if ( abs(Nzk) > 1 ) {
      return false;
   }
   if ( abs(Nzpi+Nzk+Nzo) > 1 ) {
      return false;
   }

   hst[9]->Fill(Ntot);
   if ( isMC ) {
      hst[201]->Fill(ppKK.good);
      hst[202]->Fill(ppKK.decPsip);
      if ( ppKK.decPsip == 64 ) {
         hst[203]->Fill( ppKK.decJpsi );
      }
   }

   return true;
}

// {{{1 1B) analysis of the neutral tracks
//--------------------------------------------------------------------
static bool NeutralTracks(ReadDst* selector, PipPimKpKm& ppKK) {
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
//       if ( absCosTheta < 0.8 ) {  //barrel
//          hst[16]->Fill(emcpos.z(), RtoD(emcpos.phi()) );
//       } else { // endcap
//          int iz = (emcpos.z() > 0);
//          hst[17+iz]->Fill(emcpos.rho(), RtoD(emcpos.phi()) );
//       }

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
      if ( tang < min_angle ) {
         continue;
      }

      ppKK.trk_g.push_back( emcTrk );
      if ( ppKK.Eg_max < eraw ) {
         ppKK.Eg_max = eraw;
      }
      ppKK.Eg_sum += eraw;

      Hep3Vector p3 = emcpos - ppKK.xorig;
      p3 *= eraw / p3.mag();
      ppKK.LVg.push_back( HepLorentzVector(p3,eraw) );

   } // end of neutrals loop

   int Ng = ppKK.trk_g.size(); // number of good photons

   hst[11]->Fill(Ng);
   hst[12]->Fill( ppKK.Eg_max );
   hst[13]->Fill( ppKK.Eg_sum );

   if ( isMC ) {
      hst[211+ppKK.good]->Fill(Ng);
      hst[213+ppKK.good]->Fill( ppKK.Eg_max );
      hst[215+ppKK.good]->Fill( ppKK.Eg_sum );
   }

   // search for pi0
   int Npi0 = 0;
   for ( int i = 0; i < Ng-1; ++i ) {
      auto& LVgi = ppKK.LVg[i];
      for ( int j = i+1; j < Ng; ++j ) {
         auto& LVgj = ppKK.LVg[j];

         double Mgg2 = (LVgi+LVgj).m2(); // invariant mass of (gg)
         hst[14]->Fill(Mgg2);
         hst[19]->Fill(Mgg2); // in wide region
         if ( isMC ) {
            hst[217+ppKK.good]->Fill( Mgg2 );
         }

         // mpi0^2=0.0182, asymmetrical background
         if ( fabs(Mgg2-0.017) < 0.005 ) { // strong cut!
            Npi0 += 1;
         }
      }
   }

   hst[15]->Fill(Npi0);
   if ( isMC ) {
      hst[219+ppKK.good]->Fill(Npi0);
      if ( Npi0 > 0 ) {
         hst[221]->Fill(ppKK.good);
      }
   }
   if ( Npi0 > 0 ) {
      return false;
   }

   if ( isMC ) {
      hst[222]->Fill(ppKK.good);
      hst[223]->Fill(ppKK.decPsip);
      if ( ppKK.decPsip == 64 ) {
         hst[224]->Fill( ppKK.decJpsi );
      }
   }

   return true;
}

// {{{1 2) recoil masses of soft pi+ pi-
//--------------------------------------------------------------------
static void MrecPiPi(PipPimKpKm& ppKK) {
//--------------------------------------------------------------------
   int Npi = ppKK.trk_Pi.size();
   if ( Npi != 4 ) {
      return;
   }

   for ( int i = 1; i < Npi; ++i ) {
      auto trki = ppKK.trk_Pi[i];
      if ( trki->p() > 0.45 ) {
         continue;
      }
      int zi = trki->charge();
      const auto& LVi = ppKK.LVpi[i];

      for ( int j = 0; j < i; ++j ) {
         auto trkj = ppKK.trk_Pi[j];
         if ( trkj->p() > 0.45 || trkj->charge() != -zi ) {
            continue;
         }
         const auto& LVj = ppKK.LVpi[j];

         double Mrec2 = (ppKK.LVcms - LVi - LVj).m2();
         if( Mrec2 < 0 ) {
            continue;
         }
         double Mrec = sqrt(Mrec2);

         hst[21]->Fill(Mrec);
         hst[22]->Fill(Mrec);

         if( Mrec < 3.0 || Mrec > 3.2 ) { // preliminary selection
            continue;
         }

         // Theta(pi+ pi-)
         double cosPM = Hep3Vector(LVj).cosTheta(Hep3Vector(LVi));
         hst[23]->Fill( cosPM );

         // Minv( pi+ pi-)
         double invPM = (LVi+LVj).m();
         hst[24]->Fill( invPM );

         // check Mrec candidate
         if ( !SelectPM( cosPM, invPM ) ) {
            continue;
         }

         // the best candidate (closest to M(J/Psi))
         if ( fabs(Mrec-mjpsi) < fabs(ppKK.Mrec-mjpsi) ) {
            ppKK.Mrec = Mrec;
         }

      } // end for(j)
   } // end for(i)

   hst[27]->Fill(ppKK.Mrec);
   if ( isMC ) {
      int PiPiJpsi = ( ppKK.decPsip == 64 );
      hst[231+PiPiJpsi]->Fill( ppKK.Mrec );
   }
}

// {{{1 3) Kaons track reconstruction efficiency
//--------------------------------------------------------------------
static void KTrkRecEff(PipPimKpKm& ppKK, const JobOption* job) {
//--------------------------------------------------------------------
   if( ppKK.Mrec < job->Mrec_min || ppKK.Mrec > job->Mrec_max ) {
      return;
   }

   int hs = job->Hshift;

   int Npi = ppKK.trk_Pi.size();
   int Nk  = ppKK.trk_K.size();
   int No  = ppKK.trk_o.size();
   int Ntot= Npi + Nk + No;

   HepLorentzVector LVpi_sum;
   for( const auto& lv : ppKK.LVpi ) {
      LVpi_sum += lv;
   }

   for ( int k = 0; k < Nk; ++k ) {
      auto tr_k = ppKK.trk_K[k];     // reconstructed track
      const auto& LV_k = ppKK.LVk[k];

      // recoil momentum of lost kaon
      HepLorentzVector LV_kl = ppKK.LVcms - LVpi_sum - LV_k;
      double Mrk2 = LV_kl.m2();
      hst[hs+31]->Fill(Mrk2);

      if ( isMC ) {
         hst[hs+241+ppKK.good]->Fill( Mrk2 );
      }

      if ( fabs(Mrk2-SQ(mk)) > job->M2Kwin ) {
         continue;
      }

      int Zk     = -tr_k->charge();  // charge of lost kaon
      double Pk  = LV_kl.rho();      // total momentum
      double Ptk = LV_kl.perp();     // Pt
      double Ck  = LV_kl.cosTheta(); // polar angle

      // search for reconstructed track close to predicted
      int found = 0;
      double dP = 10, dTh = 100;
      if ( Ntot == 6 ) {
         Hep3Vector Prec;
         if ( Nk == 2 ) {
            Prec = ppKK.LVk[1-k].vect();
         } else {
            Prec = ppKK.Vo[0];
         }
         dP = Pk - Prec.mag();
         double cosTh = LV_kl.vect().cosTheta( Prec );
         dTh = RtoD( acos(cosTh) );
         hst[hs+32]->Fill( dP );
         hst[hs+33]->Fill( dTh );
         hst[hs+34]->Fill( dTh, dP );
         if ( isMC ) {
            hst[hs+243+ppKK.good]->Fill( dP );
            hst[hs+245+ppKK.good]->Fill( dTh );
            hst[hs+247+ppKK.good]->Fill( dTh,dP );
         }

         if ( -0.12 < dP && dP < 0.08 && dTh < 10 ) { // OK
            found = 1;
         } else {
            // there are 6 tracks and the predicted track does
            // not match any of them - this is a background event:
            if( isMC ) {
               hst[hs+249]->Fill(ppKK.good);
            }
            // still save it to the root-tree
            // continue;
         }
      } // end if(Ntot==6)

      int iz = (Zk < 0);
      hst[hs+35+iz]->Fill(Pk);
      hst[hs+37+iz]->Fill(Ptk);
      hst[hs+39+iz]->Fill(Ck);

      hst[hs+41+iz+2*found]->Fill(Ck,Ptk); // Eff trk K

      if( isMC ) {
         hst[hs+251]->Fill(ppKK.good);
         hst[hs+252]->Fill(ppKK.decPsip);
         if ( ppKK.decPsip == 64 ) {
            hst[hs+253+ppKK.good]->Fill(ppKK.decJpsi);
         }
      }

      // study PID of Kaons
      if ( found ) {
         int ipid = (Nk == 2 )*2;
         hst[hs+45+iz+ipid]->Fill(Ck,Ptk); // Eff pid K

         // if ( isMC ) {
            // // re-weighted: for eff(trk*PID)
            // if ( Nk == 2 ) {
               // double pt = ppKK.trk_K[1-k]->pxy();
               // double w  = ReWeightTrkPid(DataPeriod,1,pt);
               // hst[hs+255]->Fill(w);
               // static_cast<TH2D*>(hst[hs+256+iz])->Fill(Ck,Ptk,w);
            // }
         // }
      }

      if ( hs != 400 ) {
        continue;
      }

      // fill and save Ttree
      xnt_effK.Zk   = Zk;
      xnt_effK.Ptk  = Ptk;
      xnt_effK.Ck   = Ck;
      xnt_effK.dP   = dP;
      xnt_effK.dTh  = dTh;
      xnt_effK.Nk   = Nk;

      if ( found && ( Nk == 2 ) ) {
         found = 2;
      }
      xnt_effK.fl   = found;

      xnt_effK.Mrec = ppKK.Mrec;
      xnt_effK.Mrk2 = Mrk2;

      if ( isMC ) {
         xnt_effK.good = ppKK.good;
      } else {
         xnt_effK.good = 0;
      }

      m_nt_effK->Fill();

   } // end for Kaons
}

// {{{1 4) reconstruct invariant masses of pi+ pi- K+ K-
//--------------------------------------------------------------------
static void MinvPiPiKK(PipPimKpKm& ppKK) {
//--------------------------------------------------------------------
   int Npi = ppKK.trk_Pi.size();
   int Nk  = ppKK.trk_K.size();
   if ( Nk != 2 ) {
      return;
   }

   HepLorentzVector LVk_sum;
   for( const auto& lv : ppKK.LVk ) {
      LVk_sum += lv;
   }

   for ( int i = 0; i < Npi-1; i++  ) {
      auto trki = ppKK.trk_Pi[i];
      int zi = trki->charge();
      const auto& LVi = ppKK.LVpi[i];

      for ( int j = i+1; j < Npi; j++ ) {
         auto trkj = ppKK.trk_Pi[j];
         if ( trkj->charge() != -zi ) {
            continue;
         }
         const auto& LVj = ppKK.LVpi[j];

         HepLorentzVector LVppKK = LVi + LVj + LVk_sum;
         double Minv = LVppKK.m();
         hst[51]->Fill(Minv);
         hst[52]->Fill(Minv);

         if( Minv < 3.0 || Minv > 3.2 ) { // preliminary selection
            continue;
         }

         // save pions from non J/Psi decay
         vector<int> idxl;
         for ( int l = 0; l < Npi; ++l ) {
            if ( l == i || l == j ) {
               continue;
            }
            idxl.push_back(l);
         }

         // good candidate
         if ( fabs(Minv-mjpsi) < fabs(ppKK.Minv-mjpsi) ) {
            ppKK.Minv = Minv;
            ppKK.LVjpsi = LVppKK;
            ppKK.fpi = idxl;
         }
      } // end for(j)
   } // end for(i)

   hst[54]->Fill(ppKK.Minv);
   if ( isMC ) {
      hst[261+ppKK.good]->Fill( ppKK.Minv );
   }
   for ( int i : ppKK.fpi ) {   // pions from non J/Psi decay
      double Ppi = ppKK.LVpi[i].rho(); // total momentum
      hst[55]->Fill(Ppi);
      if ( isMC ) {
         hst[263+ppKK.good]->Fill( Ppi );
      }
   }

}

// {{{1 5) Pions track reconstruction efficiency
//--------------------------------------------------------------------
static void PiTrkRecEff(PipPimKpKm& ppKK, const JobOption* job) {
//--------------------------------------------------------------------
   if( ppKK.Minv < job->Minv_min || ppKK.Minv > job->Minv_max ) {
      return;
   }

   int hs = job->Hshift;

   int Npi = ppKK.trk_Pi.size();
   int Nk  = ppKK.trk_K.size();
   int No  = ppKK.trk_o.size();
   int Ntot= Npi + Nk + No;

   for ( int i : ppKK.fpi ) { // pions from non J/Psi decay
      auto tr_pi = ppKK.trk_Pi[i];
      const auto& LVi = ppKK.LVpi[i];
      if ( LVi.rho() > 0.45 ) { // total momentum
         continue;
      }

      // recoil momentum of lost pion
      HepLorentzVector LV_pil = ppKK.LVcms - ppKK.LVjpsi - LVi;
      double Ppi = LV_pil.rho();        // total momentum
      double MrPi2 = LV_pil.m2();
      hst[hs+61]->Fill(MrPi2);
      hst[hs+62]->Fill(Ppi);

      if ( isMC ) {
         hst[hs+271+ppKK.good]->Fill( MrPi2 );
         hst[hs+273+ppKK.good]->Fill( Ppi );
      }

      if ( Ppi > 0.45 ) {
         continue;
      }

      if ( fabs(MrPi2-SQ(mpi)) > job->M2Piwin ) {
         continue;
      }

      int Zpi     = -tr_pi->charge();  // charge of lost kaon
      double Ptpi = LV_pil.perp();     // Pt
      double Cpi  = LV_pil.cosTheta(); // polar angle

      // search for reconstructed track close to predicted
      int found = 0;
      double dP = 10, dTh = 100;
      if ( Ntot == 6 ) {
         Hep3Vector Prec;
         if ( Npi == 4 ) {
            for ( int j : ppKK.fpi ) {
               if ( j == i ) {
                  continue;
               }
               Prec = ppKK.LVpi[j].vect();
               break;
            }
         } else {
            Prec = ppKK.Vo[0];
         }
         dP = Ppi - Prec.mag();
         double cosTh = LV_pil.vect().cosTheta( Prec );
         dTh = RtoD( acos(cosTh) );
         hst[hs+63]->Fill( dP );
         hst[hs+64]->Fill( dTh );
         hst[hs+65]->Fill( dTh, dP );
         if ( isMC ) {
            hst[hs+275+ppKK.good]->Fill( dP );
            hst[hs+277+ppKK.good]->Fill( dTh );
            hst[hs+279+ppKK.good]->Fill( dTh,dP );
         }

         if ( -0.12 < dP && dP < 0.08 && dTh < 15 ) { // OK
            found = 1;
         } else {
            // there are 6 tracks and the predicted track does
            // not match any of them - this is a background event:
            if( isMC ) {
               hst[hs+281]->Fill(ppKK.good);
            }
            // still save it to the root-tree
            // continue;
         }
      } // end if(Ntot==6)

      int iz = (Zpi < 0);
      hst[hs+66+iz]->Fill(Ptpi);
      hst[hs+68+iz]->Fill(Cpi);

      hst[hs+71+iz+2*found]->Fill(Cpi,Ptpi); // Eff trk pi

      if ( isMC ) {
         hst[hs+282]->Fill(ppKK.good);
         hst[hs+283]->Fill(ppKK.decPsip);
         if ( ppKK.decPsip == 64 ) {
            hst[hs+284+ppKK.good]->Fill(ppKK.decJpsi);
         }
      }

      // study PID of pions
      if ( found ) {
         int ipid = (Npi == 4 )*2;
         hst[hs+75+iz+ipid]->Fill(Cpi,Ptpi);  // Eff pid pi

         // if ( isMC ) {
            // // re-weighted: for eff(trk*PID)
            // if ( Npi == 4 ) {
               // for ( int j : ppKK.fpi ) {
                  // if ( j == i ) {
                     // continue;
                  // }
                  // double pt = ppKK.trk_Pi[j]->pxy();
                  // double w  = ReWeightTrkPid(DataPeriod,0,pt);
                  // hst[hs+286]->Fill(w);
                  // static_cast<TH2D*>(hst[hs+287+iz])
                     // ->Fill(Cpi,Ptpi,w);
                  // break;
               // }
            // }
         // }
      }

      if ( hs != 400 ) {
        continue;
      }

      // fill and save Ttree
      xnt_effPi.Zpi   = Zpi;
      xnt_effPi.Ptpi  = Ptpi;
      xnt_effPi.Cpi   = Cpi;
      xnt_effPi.dP    = dP;
      xnt_effPi.dTh   = dTh;
      xnt_effPi.Npi   = Npi;

      if ( found && ( Npi == 4 ) ) {
         found = 2;
      }
      xnt_effPi.fl    = found;

      xnt_effPi.MppKK = ppKK.Minv;
      xnt_effPi.Mrpi2 = MrPi2;

      if ( isMC ) {
         xnt_effPi.good = ppKK.good;
      } else {
         xnt_effPi.good = 0;
      }

      m_nt_effPi->Fill();

   } // end loop for pions
}

// {{{1 MAIN: Loop for each event
//--------------------------------------------------------------------
bool PsipPiPiKKEvent( ReadDst*       selector,
                      TEvtHeader*    m_TEvtHeader,
                      TDstEvent*     m_TDstEvent,
                      TEvtRecObject* m_TEvtRecObject,
                      TMcEvent*      m_TMcEvent,
                      TTrigEvent*    m_TTrigEvent,
                      TDigiEvent*    m_TDigiEvent,
                      THltEvent*     m_THltEvent       ) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " start " << __func__ << "()" << endl;
   }

   PipPimKpKm ppKK;

   //-----------------------------------------------------------------
   //-- Get event information --
   //-----------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   int eventNo = m_TEvtHeader->getEventId();
   ppKK.runNo  = runNo;
   ppKK.event  = eventNo;

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
   ppKK.LVcms = HepLorentzVector(Ecms*sin(beam_angle), 0, 0, Ecms);

   FillHistoMC(selector,ppKK); // MC histo

   Hep3Vector xorigin = getVertexOrigin(runNo);
   ppKK.xorig = xorigin;

   // Study of the track reconstruction efficiency
   if ( !ChargedTracks(selector, ppKK) ) {
      return false;   // 1A)
   }
   if ( !NeutralTracks(selector, ppKK) ) {
      return false;   // 1B)
   }

   MrecPiPi(ppKK);                      // 2)
   for(const auto& job : Jobs) {
      KTrkRecEff(ppKK,job);             // 3)
   }

   MinvPiPiKK(ppKK);                    // 4)
   for(const auto& job : Jobs) {
      PiTrkRecEff(ppKK,job);            // 5)
   }

   return true;
}

// {{{1 EndJob
//--------------------------------------------------------------------
void PsipPiPiKKEndJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << __func__ << "()" << endl;
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
