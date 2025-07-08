//======================================================================//
//                                                                      //
// search for e+ e- -> pi+ pi- pi0 Eta                                  //
//                                  |-> 2 gammas                        //
//                                                                      //
//======================================================================//

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>
// #include <map>

#include <TH1.h>
#include <TH2.h>
#include <TNtupleD.h>
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
#include "EventTag/DecayTable.h"

#include "DstEvtRecTracks.h"
#include "TofHitStatus.h"
#include "ReadDst.h"

using namespace std;

struct Select_PipPimPi0Eta {
  double Ecms;                  // energy in center of mass system

  // MC information:
  int decJpsi;                  // decay codes for J/Psi MC events
  int dec_eta;                  // 1 if eta->2gamma
  int dec_phi;                  // 1 if phi->pi+pi-pi0

  // charge tracks:
  vector<int> qtrks;            // trk# in EvtRecTrkCol
  int ipip, ipim;               // trk# for charge +/-
  HepLorentzVector ppip, ppim;  // 4-momentum of trks

  // neutral (gammas) tracks
  vector<int> gammas;           // trk# in EvtRecTrkCol
  vector<HepLorentzVector> Pg;  // 4-momentum of gammas

  Select_PipPimPi0Eta() {
    Ecms = -1;
    decJpsi = -1;
    dec_eta = dec_phi = 0;
    ipip = ipim = -1;
  }

  int nGood() { return qtrks.size(); }
  int nGam() { return gammas.size(); }
};

typedef Select_PipPimPi0Eta Select;

//-----------------------------------------------------------------------------
// Global file variables
//-----------------------------------------------------------------------------

const static double beam_angle = 0.011; // 11 mrad

// masses of particles (GeV)           from PDG:
static const double mjpsi  = 3.096916; // 3096.916  +/- 0.011   MeV
const static double mpi    = 0.13957;  // 139.57018 +/- 0.00035 MeV
const static double mpi0   = 0.13498;  // 134.9766  +/- 0.0006  MeV
const static double meta   = 0.547862; // 547.862   +/- 0.017   MeV
const static double momega = 0.78265;  // 782.65    +/- 0.12    MeV
static const double mphi   = 1.019461; //1019.461   +/- 0.019   MeV

static AbsCor* m_abscor = 0;
static EventTagSvc* m_EventTagSvc = 0;

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;
static std::vector<TNtupleD* > m_tuple;

// decays classification
static DecayTable JpsiTbl;  // final selection for phi eta
static DecayTable JpsiTblW; // for omega-eta

// container for warnings
static map<string,int> warning_msg;

static bool isMC = false;

// Functions: use C-linkage names
#ifdef __cplusplus
extern "C" {
#endif

static inline void Warning(const char* msg) {
   warning_msg[string(msg)] += 1;
}

static inline Double_t RtoD(Double_t ang) {return ang*180/M_PI;}

static inline double SQ(double x) {
   return x*x;
}

void PipPimPi0EtaStartJob (ReadDst* selector)
{
  if ( selector->Verbose() ) cout << " Start: " << __func__ << "()" << endl;

  his1.resize(500,(TH1D*)0);
  his2.resize(500,(TH2D*)0);
  m_tuple.resize(10,0);

  // init Absorption Correction
  m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

  // initialize EventTag
  m_EventTagSvc = EventTagSvc::instance();
  if ( !m_EventTagSvc->IsInitialized() ) {
     // set paths to pdg & decayCodes files:
     m_EventTagSvc->setPdtFile(
        selector->AbsPath( "Analysis/EventTag/share/pdt_bean.table" )
                              );
     m_EventTagSvc->setDecayTabsFile(
        selector->AbsPath(
           "Analysis/EventTag/share/DecayCodes/dcode_charmonium.txt"
        )                           );
     m_EventTagSvc->setIgnorePhotons(false); // for "dcode_charmonium.txt"
     if ( selector->Verbose() ) {
        m_EventTagSvc->setVerbose(1);
     }
     m_EventTagSvc->initialize();
  } else {
     cout << " WARNING in " << __func__ << ": "
          << "EventTagSvc has already been initialized" << endl;
     Warning("EventTagSvc has already been initialized");
  }

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

  // set path for ParticleID algorithm
  ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
  pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
  pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif


  his1[0] = new TH1D("InvMBef","The invariant mass before kinematic fit", 500, 0, 5);
  his1[1] = new TH1D("JPsi1","The invariant mass for J/#psi one", 500, 0, 5);
  his1[2] = new TH1D("JPsi2","The invariant mass for J/#psi two", 500, 0, 5);
  his1[3] = new TH1D("Pi0","The invariant mass for #pi_{0}", 1000, 0, 1);
  his1[4] = new TH1D("Eta","The invariant mass for #eta", 500, 0, 5);
  his1[5] = new TH1D("Phi","The invariant mass for #phi", 500, 0, 5);
  his1[6] = new TH1D("RecoilMass","The recoil mass", 100, 0.98, 1.08);


  // Selection histograms
  his1[17] = new TH1D("Runs", "run numbers", 1000,0.,50000);
  his1[18] = new TH1D("cuts_0","selections cuts", 20,-0.5,19.5);

  his1[11] = new TH1D("Rxy","R_{xy}", 200,-2.,2.);
  his1[12] = new TH1D("Rz","R_{z}", 200,-20.,20.);
  his1[13] = new TH1D("cos_theta","cos(#theta)", 200,-1.,1.);
  his1[14] = new TH1D("QPtrack","Charged momentum (Q*P)", 400,-2.,2.);
  his1[15] = new TH1D("Ntrks","N_{charged tracks} in event", 11,-0.5,10.5);
  his1[16] = new TH1D("Ncharge","N_{charge} in event", 11,-5.5,5.5);

  his1[21] = new TH1D("PidInfo", "0 - no pid, 1 - pid", 2,-0.5,1.5);
  his1[22] = new TH1D("LogCLpi","lg(CL_{#pi})", 100,-5.,0.);
  his1[23] = new TH1D("LogCLK","lg(CL_{K})", 100,-5.,0.);
  his1[24] = new TH1D("LogCLp","lg(CL_{p})", 100,-5.,0.);

  his2[21] = new TH2D("SignTrks","sign of trks (+) vs (-)",
                3,-0.5,2.5, 3,-0.5,2.5);

  his1[19] = new TH1D("NGamma","N_{#gamma}",10 ,0 ,10);
  his1[20] = new TH1D("G_ang","angle with closest chg.trk", 180,0.,180.);

  his1[25] = new TH1D("Mpi0","M_{#pi^{0}}", 200, 0., 0.2);
  his1[26] = new TH1D("Mpi0phi","M_{#pi^{0}} phi-eta", 200, 0., 0.2);
  his1[27] = new TH1D("Mpi0w","M_{#pi^{0}} omega-eta", 200, 0., 0.2);

  his1[41] = new TH1D("vtx_fit", "vertex fit: 0/1 - bad/good", 2,-0.5,1.5);
  his1[42] = new TH1D("vtx_chi2", "vertex fit #chi^2", 100,0.,100.);

  his1[45] = new TH1D("f5c_loop","0,1=false Fit0,Fit1,5=oksq",10,-0.5,10.5);
  his1[51] = new TH1D("f5c_chisq","5C-fit: min #chi^{2}", 300,0.,300.);
  his1[52] = new TH1D("f5c_chi2","5C-fit: final min #chi^{2}", 200,0.,200.);

  // Monte Carlo histograms:
  his1[100] = new TH1D("mc_deccode0", "decCode nocut",258,-1.5,256.5);
  his1[102] = new TH1D("mc_deccode2", "decCode cut2",258,-1.5,256.5);
  his1[103] = new TH1D("mc_deccode3", "decCode cut3",258,-1.5,256.5);
  his1[105] = new TH1D("mc_deccode5", "decCode cut5",258,-1.5,256.5);

  his1[111] = new TH1D("mc_pdg", "PDG codes of all particles",
                      2001,-1000.5,1000.5);
  his1[112] = new TH1D("mc_pdg0", "PDG of particles from primary vertex",
                      2001,-1000.5,1000.5);
  his1[121] = new TH1D("mc_JpsiPt","Pt of J/#Psi", 100,0.,1.);
  his1[122] = new TH1D("mc_JpsiC", "cos(#Theta) of J/#Psi", 100,-1.,1.);
  his1[123] = new TH1D("mc_EtaP", "Momentum of #eta", 1000,0.,2.);
  his1[124] = new TH1D("mc_EtaPt","Pt of #eta", 1000,0.,2.);
  his1[125] = new TH1D("mc_EtaC", "cos(#Theta) of #eta", 100,-1.,1.);
  his1[126] = new TH1D("mc_PhiP", "Momentum of #phi", 1000,0.,2.);
  his1[127] = new TH1D("mc_PhiPt","Pt of #phi", 1000,0.,2.);
  his1[128] = new TH1D("mc_PhiC", "cos(#Theta) of #phi", 100,-1.,1.);
  his1[131] = new TH1D("mc_EtaRC",
           "cos(#Theta) of #eta in the rest frame of J/#Psi", 100,-1.,1.);
  his1[132] = new TH1D("mc_PhiRC",
           "cos(#Theta) of #phi in the rest frame of J/#Psi", 100,-1.,1.);
  his1[141] = new TH1D("mc_2gE", "E(g) from eta->2g", 100,0.,2.);
  his1[142] = new TH1D("mc_2gC", "cos(g) from eta->2g", 100,-1.,1.);
  his1[143] = new TH1D("mc_3pipP", "P(pi+) from phi", 100,0.,2.);
  his1[144] = new TH1D("mc_3pipC", "cos(pi+) from phi", 100,-1.,1.);
  his1[145] = new TH1D("mc_3pimP", "P(pi-) from phi", 100,0.,2.);
  his1[146] = new TH1D("mc_3pimC", "cos(pi-) from phi", 100,-1.,1.);
  his1[147] = new TH1D("mc_3pi0P", "P(pi0) from phi", 100,0.,2.);
  his1[148] = new TH1D("mc_3pi0C", "cos(pi0) from phi", 100,-1.,1.);
  his1[149] = new TH1D("mc_2g3pi", "1=eta2g,2=phi3pi and 7=1&2,8=all",
                       10,0.5,10.5);

//   m_tuple[0] = new TNtuple("a5C","after 5C fit",
//                                 "ch2:"
//                                 "Etot:Ptot"
//                                 );

  m_tuple[0] = new TNtupleD("a5C","after 5C fit",
               "ch2:"               // chi^2 of 5C fit
               "Ppip:Cpip:phipip:"  // P, cos(Theta), phi of pi+
               "Ppim:Cpim:phipim:"  // P, cos(Theta), phi of pi-
               "Ppi0:Cpi0:phipi0:"  // P, cos(Theta), phi of pi0
               "P2g:C2g:phi2g:"     // P, cos(Theta), phi of gamma-gamma
               "M2pi3:M2gg:M2rec:"  // squares of invariant mass
               "dec:deceta:decphi:" // MC decay codes of J/Psi, eta, phi
               "M2pi0:Npi0"
                                );

  const char* SaveDir = "one";
  VecObj his1o(his1.begin(),his1.end());
  selector->RegInDir(his1o,SaveDir);
  VecObj ntuples(m_tuple.begin(),m_tuple.end());
  selector->RegInDir(ntuples,SaveDir);


}

static Hep3Vector getVertexOrigin(int runNo, bool verbose = false)
{
  static int save_runNo = 0;
  static Hep3Vector xorigin;

  if ( runNo == save_runNo ) return xorigin;

  // update vertex for new run
  xorigin.set(0.,0.,0.);
  VertexDbSvc* vtxsvc = VertexDbSvc::instance();

  int run = abs(runNo);
  if (
       (run >=  9947 && run <= 10878)  // J/Psi 2009
     ) {
    vtxsvc->SetBossVer("6.6.3");
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

static void FillHistoMC(const ReadDst* selector, Select& Slct)
{
   if ( !isMC ) {
      return;
   }

   const TMcEvent*  m_TMcEvent  = selector->GetMcEvent();
   const TObjArray* mcParticles = m_TMcEvent->getMcParticleCol();
   if ( !mcParticles ) {
      return;
   }

   // EventTag: for J/Psi MC:
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
   his1[100]->Fill( Slct.decJpsi );

   JpsiTbl.Reset();

   int idx_jpsi=-1;
   int idx_eta=-1;
   int idx_phi=-1;
   int idx_rho=-1;
   // momentum eta & phi in rest frame of J/Psi
   Hep3Vector beta_jpsi;
   HepLorentzVector LVeta, LVphi;
   // search for decay eta -> 2gamma
   vector<int> Pdg_eta;
   // search for decay phi -> rho pi -> pi+ pi- pi0
   vector<int> Pdg_phi, Pdg_rho;

   TIter mcIter(mcParticles);
   while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
      long int part_pdg = part->getParticleID ();

      Hep3Vector Vp( part->getInitialMomentumX(),
                     part->getInitialMomentumY(),
                     part->getInitialMomentumZ() );

      his1[111]->Fill(part_pdg);
      if ( part->getMother() == -99 ) { // primary vertex
         his1[112]->Fill(part_pdg);
      }

      if ( part_pdg == 443 ) {                // J/Psi
         his1[121]->Fill(Vp.rho());
         his1[122]->Fill(Vp.cosTheta());
         idx_jpsi = part->getTrackIndex();
         HepLorentzVector LVjpsi( Vp, sqrt(Vp.mag2() + SQ(mjpsi)) );
         beta_jpsi = LVjpsi.boostVector();
      }

      if( part->getMother() == idx_jpsi ) {  // decays of J/Psi
         // collect decays of J/psi
         JpsiTbl.vecdec.push_back(part_pdg);
      }

      if ( part->getMother() == idx_jpsi
            && Slct.decJpsi == 68        ) {    // J/Psi -> phi eta
         if ( part_pdg == 221 ) {               // eta
            his1[123]->Fill(Vp.mag());
            his1[124]->Fill(Vp.rho());
            his1[125]->Fill(Vp.cosTheta());
            idx_eta = part->getTrackIndex();

            HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(meta)) );
            LV.boost(-beta_jpsi);
            LVeta=LV;
            his1[131]->Fill( LV.cosTheta() );

         } else if ( part_pdg == 333 ) {        // phi
            his1[126]->Fill(Vp.mag());
            his1[127]->Fill(Vp.rho());
            his1[128]->Fill(Vp.cosTheta());
            idx_phi = part->getTrackIndex();

            HepLorentzVector LV( Vp, sqrt(Vp.mag2() + SQ(mphi)) );
            LV.boost(-beta_jpsi);
            LVphi=LV;
            his1[132]->Fill( LV.cosTheta() );
         }
      }

      if ( part->getMother() == idx_eta ) { // J/Psi -> phi eta
         Pdg_eta.push_back(int(part_pdg));  // eta -> ...
      }
      if ( part->getMother() == idx_phi ) { // J/Psi -> phi eta
         Pdg_phi.push_back(int(part_pdg));  // phi -> ...
         if ( part_pdg == 113 || abs(part_pdg) == 213 ) { // rho
            idx_rho = part->getTrackIndex();
         }
      }
      if ( part->getMother() == idx_rho ) { // phi -> rho X
         Pdg_rho.push_back(int(part_pdg));  // rho -> ...
      }
   } // end of while

   // paranoid check
   if ( Slct.decJpsi == 68  && ( idx_eta == -1 || idx_phi == -1 ) ) {
      cout << "WARNING: Slct.decJpsi= " << Slct.decJpsi
           << " idx_eta,idx_phi= " << idx_eta << "," << idx_phi << endl;
      selector->PrintMcDecayTree(-99,0); // What is it?
      Warning("MC: bad eta or phi indexes");
   }

   if ( Slct.decJpsi == 68 ) {  // J/Psi -> phi eta

      // set Slct.dec_eta
      Slct.dec_eta = 0;
      if ( Pdg_eta.size() == 2 &&
           Pdg_eta[0] == 22 && Pdg_eta[1] == 22 ) {
         Slct.dec_eta = 1; // eta -> 2 gamma
         his1[149]->Fill(1);
      }

      // set Slct.dec_phi
      Slct.dec_phi = 0;
      if ( Pdg_phi.size() == 2 && idx_rho != -1 ) { // check phi->rho pi
         sort(Pdg_phi.begin(),Pdg_phi.end());
         if (    (Pdg_phi[0] ==-211 && Pdg_phi[1] == 213) // pi- rho+
              || (Pdg_phi[0] ==-213 && Pdg_phi[1] == 211) // rho- pi+
              || (Pdg_phi[0] == 111 && Pdg_phi[1] == 113) // pi0 rho0
            ) {
            if ( Pdg_rho.size() == 2 ) { // check rho->2pi
               sort(Pdg_rho.begin(),Pdg_rho.end());
               if (    (Pdg_rho[0] ==-211 && Pdg_rho[1] == 111) // pi-pi0
                    || (Pdg_rho[0] ==-211 && Pdg_rho[1] == 211) // pi-pi+
                    || (Pdg_rho[0] == 111 && Pdg_rho[1] == 211) // pi0pi+
                  ) {
                  Slct.dec_phi = 1;  // phi -> rho pi -> pi+pi-pi0
                  his1[149]->Fill(2);
               }
            }
         }
      }

      if ( Slct.dec_eta == 1 && Slct.dec_phi == 1 ) {
         his1[149]->Fill(7);
         TIter mcIter(mcParticles);
         while( auto part = static_cast<TMcParticle*>(mcIter.Next()) ) {
            long int part_pdg = part->getParticleID ();

            Hep3Vector Vp( part->getInitialMomentumX(),
                           part->getInitialMomentumY(),
                           part->getInitialMomentumZ() );
            double Ep = part->getInitialMomentumE();

            if ( part->getMother() == idx_eta ) {
               if ( part_pdg == 22 ) {
                  his1[141]->Fill(Ep);
                  his1[142]->Fill(Vp.cosTheta());
               }
            }

            if (   part->getMother() == idx_phi
                || part->getMother() == idx_rho ) {
               if ( part_pdg == 211 ) {         // pi+
                  his1[143]->Fill(Vp.mag());
                  his1[144]->Fill(Vp.cosTheta());
               } else if ( part_pdg == -211 ) { // pi-
                  his1[145]->Fill(Vp.mag());
                  his1[146]->Fill(Vp.cosTheta());
               } else if ( part_pdg == 111 ) {  // pi0
                  his1[147]->Fill(Vp.mag());
                  his1[148]->Fill(Vp.cosTheta());
               }
            }
         } // end of while-2
      } else {
         his1[149]->Fill(8);
      } // end if
   }
   return;
}


bool PipPimPi0EtaEvent (ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
{
  if ( selector->Verbose() ) cout << " start " << __func__ << "()" << endl;

  m_abscor->AbsorptionCorrection(selector);

  Select Slct;

  int runNo   = m_TEvtHeader->getRunId();
  int eventNo = m_TEvtHeader->getEventId();

  his1[17]->Fill(abs(runNo));
  his1[18]->Fill(0); // "cuts"

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
     FillHistoMC(selector,Slct);       // MC histo
  }

  Hep3Vector xorigin = getVertexOrigin(runNo);

  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

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
    double  Rvxy0=vecipa[0];  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];   //the nearest distance to IP in z direction

    his1[11]->Fill(Rvxy0);
    his1[12]->Fill(Rvz0);
    his1[13]->Fill(cosTheta);

    if(fabs(Rvz0) >= 10.0) continue;
    if(fabs(Rvxy0) >= 1.0) continue;
    if(fabs(cosTheta) >= 0.93) continue;

    Slct.qtrks.push_back(i);
    Ncharge += mdcTrk->charge();
    his1[14]->Fill( mdcTrk->charge()*mdcTrk->p() );
  }

  int nGood = Slct.nGood();
  his1[15]->Fill(nGood);
  his1[16]->Fill(Ncharge);

  if ( nGood != 2 || Ncharge != 0 ) return false;
  his1[18]->Fill(1); // "cuts"

  ParticleID* pid = ParticleID::instance();

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

//     if ( !pid->IsPidInfoValid() ) return false;
    if ( !pid->IsPidInfoValid() ) continue;

    if( pid->probPion() > 0 )   his1[22]->Fill( log10(pid->probPion()) );
    else                        his1[22]->Fill(-5.);
    if( pid->probKaon() > 0 )   his1[23]->Fill( log10(pid->probKaon()) );
    else                        his1[23]->Fill(-5.);
    if( pid->probProton() > 0 ) his1[24]->Fill( log10(pid->probProton()) );
    else                        his1[24]->Fill(-5.);

    // select Pions only
//     if ( pid->probPion() <  0.001 ) return false;
    if ( pid->probPion() <  0.001 || pid->probPion() < pid->probKaon() || pid->probPion() < pid->probProton() ) continue;

    RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
    //his1[25]->Fill( (double)(mdcKalTrk != 0) );
    if ( !mdcKalTrk ) {
      cout << " No Kalman track ! EventNo= "<< eventNo << endl;
//       return false;
      continue;
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
  }

  his2[21]->Fill( npp, npm );
  // require exactly one "+" and one "-"
  if ( npp != 1 || npm != 1 ) return false;
  his1[18]->Fill(2); // "cuts"
  if ( isMC ) {
    his1[102]->Fill( Slct.decJpsi );
  }

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

    his1[20]->Fill(RtoD(tang));
    static const double min_angle = 10 * M_PI/180; // 10 grad
    if ( tang < min_angle ) continue;

    Slct.gammas.push_back(i);

    // (E,Px,Py,Pz) for good gammas
    Hep3Vector p3 = emcpos - xorigin; // position emc cluster wrt of vertex
    p3 *= eraw / p3.mag();            // normalization on energy
    Slct.Pg.push_back( HepLorentzVector(p3,eraw) );
  }

  int nGam = Slct.nGam();
  his1[19]->Fill(Slct.nGam());
  if ( nGam < 4 ) return false;

  his1[18]->Fill(3); // "cuts"
  if ( isMC ) {
    his1[103]->Fill( Slct.decJpsi );
  }

  // makes no sense, because we can have more than 4 photons
  his1[0]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[0] + Slct.Pg[1] + Slct.Pg[2] + Slct.Pg[3] ).m() );

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
  his1[41]->Fill( double(ok) );
  if ( !ok )  return false; //skip event
  his1[18]->Fill(4); // "cuts"

  his1[42]->Fill(vtxfit->chisq(0));

  vtxfit->BuildVirtualParticle(0);
  vtxfit->Swim(0); // translate parameters of tracks to the vertex# 0

  WTrackParameter wpip = vtxfit->wtrk(0);
  WTrackParameter wpim = vtxfit->wtrk(1);

  KinematicFit * kmfit = KinematicFit::instance();
  //------------------------------------------------------

  double Ecms =  3.097; // (GeV) energy in center of mass system
  HepLorentzVector ecms(Ecms*sin(beam_angle), 0, 0, Ecms);

  //------------------------------------------------------5C Fit
  int ig[6] = {-1,-1,-1,-1, -1, -1};
  double chisq = 9999.;
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

          for ( int i1 = 2; i1 < 5; i1++ ) {
          for ( int j1 = i1+1; j1 <= 5; j1++ ) {

            kmfit->init();
            kmfit->AddTrack(0, wpip);
            kmfit->AddTrack(1, wpim);
            kmfit->AddTrack(2, 0.0, g1Trk);
            kmfit->AddTrack(3, 0.0, g2Trk);
            kmfit->AddTrack(4, 0.0, g3Trk);
            kmfit->AddTrack(5, 0.0, g4Trk);

            kmfit -> AddResonance(0, mpi0, i1, j1);
            kmfit -> AddFourMomentum(1, ecms);

            if ( !kmfit->Fit(0) ) {
               his1[45]->Fill(0);
               continue;
            }
            if ( !kmfit->Fit(1) ) {
               his1[45]->Fill(1);
               continue;
            }
            bool oksq = kmfit->Fit();
            if ( oksq ) {
              his1[45]->Fill(5);
              double chi2 = kmfit->chisq();
              // his1[12]->Fill(chi2);
              if ( chi2 < chisq ) {
                chisq = chi2;
                ig[0] = i;
                ig[1] = j;
                ig[2] = k;
                ig[3] = l;
                ig[4] = ig[i1-2];
                ig[5] = ig[j1-2];
              }
            } // end if oksq
          }} // end of for (i1,j1)
        }//------------------------------------------End for(l)
      }//--------------------------------------------End for(k)
    }//-----------------------------------------------End for(j)
  } //-----------------------------------------------End for(i)
  his1[51]->Fill(chisq);

  if ( ig[0]==-1 || ig[1]==-1 || ig[2]==-1 || ig[3]==-1 ) {
    return false; // can not find 4  good gammas
  }
  his1[18]->Fill(5); // "cuts"
  if ( isMC ) {
    his1[105]->Fill( Slct.decJpsi );
  }

  // move photons from pi0 decay to places (indices) 0 and 1
  for ( int i = 0; i < 4; i++ ) {
     if ( ig[i] == ig[4] )
       {
         int j = ig[0];
         ig[0] = ig[4];
         ig[i] = j;
       }
  }
  for ( int i = 0; i < 4; i++ ) {
     if ( ig[i] == ig[5] )
       {
         int j = ig[1];
         ig[1] = ig[5];
         ig[i] = j;
       }
  }


  // ++ repeat fit for best selected photons ++
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

  kmfit -> AddResonance(0, mpi0, 2, 3);
  kmfit -> AddFourMomentum(1, ecms);
  bool oksq = kmfit->Fit();

  double chi2 = (oksq) ? kmfit->chisq() : 9999.;
  // check correctness of second fit
  if ( !oksq || fabs(chisq - chi2) > 0.1 ) {
    cout << " WARNING: Bad second fit"
         << " chisq= "<<chisq<<" new fit chi2= "<< chi2
         << endl;
    Warning("Bad second fit");
    return false;
  }
  his1[52]->Fill(chi2);

  HepLorentzVector P_pip = kmfit->pfit(0);
  HepLorentzVector P_pim = kmfit->pfit(1);
  HepLorentzVector P_pi0 = kmfit->pfit(2) + kmfit->pfit(3);
  HepLorentzVector P_3pi = P_pip + P_pim + P_pi0;
  HepLorentzVector P_pi02 = Slct.Pg[ig[0]]+Slct.Pg[ig[1]]; //g11Trk & g21Trk;

  HepLorentzVector P_2g = kmfit->pfit(4) + kmfit->pfit(5);
  HepLorentzVector P_tot = P_3pi + P_2g;

  HepLorentzVector P_rec = ecms - P_2g; // recoil momentum of 'last' 2g

  double M2_3pi = P_3pi.m2();
  double M2_2g  = P_2g.m2();
  double M2_rec = P_rec.m2();
  double M_2g   = sqrt(M2_2g);
  double M_3pi  = sqrt(M2_3pi);

  int Npi0 = 0;
  for ( int i = 0; i < nGam-1; i++ ) {
     if ( i == ig[0] || i == ig[1] ) {
        continue;
     }
     for ( int j = i+1; j < nGam; j++ ) {
        if ( j == ig[0] || j == ig[1] ) {
            continue;
        }
        double Mgg0 = (Slct.Pg[i] + Slct.Pg[j]).m();
        his1[25]->Fill( Mgg0 );

        if ( Mgg0 > 0.115 && Mgg0 < 0.155 ) {
           Npi0++;
        }

        if ( chi2 < 100
               && fabs(M_2g-meta) < 0.016
               && M_3pi > 1.01 && M_3pi < 1.03
           ) {
           his1[26]->Fill( Mgg0 );
        }
        if ( chi2 < 100
               && fabs(M_2g-meta) < 0.024
               && M_3pi > 0.76 && M_3pi < 0.80
           ) {
           his1[27]->Fill( Mgg0 );
        }
     }
  }

  his1[1]->Fill( ( Slct.ppip + Slct.ppim + Slct.Pg[ig[0]] + Slct.Pg[ig[1]] + Slct.Pg[ig[2]] + Slct.Pg[ig[3]] ).m() );
  his1[2]->Fill( P_tot.m() );
  his1[3]->Fill( P_pi0.m() );
  his1[4]->Fill( (P_pip + P_pim + P_pi0).m() );

  his1[6]->Fill( P_rec.m() );

  double xfill[] = {
      chi2,
      P_pip.rho(), P_pip.cosTheta(), P_pip.phi(),
      P_pim.rho(), P_pim.cosTheta(), P_pim.phi(),
      P_pi0.rho(), P_pi0.cosTheta(), P_pi0.phi(),
      P_2g.rho(),  P_2g.cosTheta(),  P_2g.phi(),
      M2_3pi, M2_2g, M2_rec,
      double(Slct.decJpsi), double(Slct.dec_eta), double(Slct.dec_phi),
      P_pi02.m2(), double(Npi0)
  };
  m_tuple[0]->Fill( xfill );

  if ( isMC ) {
     // fill decay table of Jpsi for very narrow region of phi-eta
     if ( chi2 < 100
            && fabs(M_2g-meta) < 0.016
            && M_3pi > 1.01 && M_3pi < 1.03
        ) {
            JpsiTbl.Add();
          }
     // fill decay table of Jpsi for region of omega-eta
     if ( chi2 < 100
            && fabs(M_2g-meta) < 0.024
            && M_3pi > 0.76 && M_3pi < 0.80
        ) {
            JpsiTblW.vecdec = JpsiTbl.vecdec;
            JpsiTblW.Add();
          }
  }

  return true;
}

void PipPimPi0EtaEndJob (ReadDst* selector)
{
  if ( selector->Verbose() ) cout << __func__ << "()" << endl;

  if ( isMC ) {
     // print tables of decays for phi-eta
     string my_cuts(
        " List of cuts for J/Psi -> pi+ pi- pi0 eta:\n "
        "      chi2<100 && abs(Mgg-meta)<0.016 && 1.01<M3pi<1.03"
     );
     cout << my_cuts << endl << endl;

     cout << string(65,'#') << endl;
     cout << "Decays of J/psi" << endl
          << JpsiTbl.ntot << " decays" << endl
          << "       size of table is " << JpsiTbl.Size() << endl;
     JpsiTbl.Print(0.1); // do not print decays with P<0.1% of all
     cout << "Enddecay" << endl << endl;

     // print tables of decays for omega-eta
     string my_cutsW(
        " List of cuts for J/Psi -> pi+ pi- pi0 eta:\n "
        "      chi2<100 && abs(Mgg-meta)<0.024 && 0.76<M3pi<0.80"
     );
     cout << my_cutsW << endl << endl;

     cout << string(65,'#') << endl;
     cout << "Decays of J/psi W" << endl
          << JpsiTblW.ntot << " decays" << endl
          << "       size of table is " << JpsiTblW.Size() << endl;
     JpsiTblW.Print(0.1); // do not print decays with P<0.1% of all
     cout << "Enddecay" << endl << endl;
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
