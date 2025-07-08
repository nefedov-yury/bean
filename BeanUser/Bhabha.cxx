//#ifdef BHABHA
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Bhabha                                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
/*=====================================================

1. Select Bhabha events: 2 tracks, each with E>0.8*ebeam, p>0.8*pbeam
  a) Plot distribution of Theta1+Theta2,
     separately for 20<Theta*<40, 40<Theta*<90
     Theta*=90-|90-0.5*(180+Theta1-Theta2)|
  b) Plot distribution of |Phi1-Phi2|,
     separately for 20<Theta*<40, 40<Theta*<90
  c) Plot Theta1+Theta2 in small bins of Theta*
     Fit each distribution by Gaussian and plot:
     Mean versus Theta*   and  Sigma versus Theta*
  d) Do the same in bins of Phi*, plot Mean and Sigma versus Phi*
     Phi*=0.5*(Phi1+Phi2-180)
  e) Plot |Phi1-Phi2| in bins of Theta*
     Fit and plot Mean and Sigma versus Theta*
  f) Plot |Phi1-Phi2| in bins of Phi, separately
     for 20<Theta*<40, 40<Theta*<90.
     Fit and plot Meand and Sigma versus Phi*

2. Select Bhabha events:
  2 tracks, Theta1+Theta2 close to 180, Phi1-Phi2 close to 180,
  for *one track* p>0.8*pbeam, E>0.8*ebeam.
  For *another* track plot:
  a) <p(electron)> versus Theta, <p(positron)> versus Theta
  b) Resolution of p versus Theta (separately electrons/positrons).
  c) <E> versus Theta
  d) Sigma(E) versus Theta
  e) Same as a-d but versus Phi (separately barrel and End-Cap).

3. Select Bhabha events:
  2 tracks, Theta1+Theta2 close to 180, Phi1-Phi2 close to 180,
  for each track p>0.8*pbeam, E>0.8*ebeam. Plot:
  a) <dE/dx> versus Theta, <dE/dx> versus Phi
  b) Sigma(dE/dx) versus Theta, versus Phi
  c) Probability to have dE/dx<0.5 MIP (or no dE/dx) versus Theta, Phi
  d) <beta> versus Theta, versus Phi
  e) Sigma(beta) versus Theta, versus Phi
  f) Probability of beta<0.5 (or no beta) versus Theta, Phi
  g) Histogram of number of MUC layers
  h) Mean number of MUC layers versus Theta, versus Phi
  i) Number of MUC layers beyond 2nd layer versus Theta, versus Phi

4. In principle, the study 1-3 can be repeated
  for dimuon events, using MUC hits instead of EMC

5. Simple Tau selection.
  a) Select 2-track events. Plot MAp1,p2), MIN(p1,p2), p1+p2
  b) Make cuts: MAp1,p2)<0.95, MIN(p1+p2)<0.9
     Plot Theta1+Theta2, |Phi1-Phi2|
  c) Make cuts: |Theta1+Theta2-180|>something
                ||Phi1-Phi2|-180|>something
     Plot E1,E2,E1+E2. Reject events with E=Ebema (if necessary)
  d) Plot p1+p2, p1+p2+Eneu (energy of all neutral EMC clusters)
  e) Plot transverse component abs(PT) of the vector sum p1+p2
     Plot the same when the vector momenta of neutral clusters
     are added to p1+p2 vector.
  f) Plot the correlation between p1+p2 and PT.
     Decide what cut would gamma-gamma (low values)
     from tau-tau (large values). Make this cut.
  g) "Measure" the branching tau->e:
     Require that 1st tau is muon (MUC hits).
     For the second tau, plot E/p.
     Estimate the fraction of events in the peak E/p=1

===================================================== */

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include "TLorentzVector.h"

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "DstEvtRecTracks.h"
#include "BMdcKalTrack.h"
#include "TofHitStatus.h"

#include "ReadDst.h"


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

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"





#define speed_of_light 299792458
const static double beam_angle = 0.011; // 11 mrad
const static int minR=28241, maxR=28513;

using namespace std;
static int ERROR_WARNING = 0;


#ifdef __cplusplus
extern "C" {
#endif

VecObj his1;
static std::vector<TH1D*> hisEC, hisB;


static TH1D* h1a1;
static TH1D* h1a2;
static TH1D* h1b1;
static TH1D* h1b2;

static TH1D* haux1;
static TH1D* haux2;
static TH1D* haux3;

static TH2D* h1c;
static TH2D* h1d;
// static TH2D * h1d_boost;
static TH2D* h1e;
static TH2D* h1f1;
static TH2D* h1f2;

static TH2D* h3a1;
static TH2D* h3a2;
static TH2D* h3a3;

static TH2D* h3d1;
static TH2D* h3d2;
static TH2D* h3d2a;
static TH2D* h3h1;
static TH2D* h3h2;

static TH2D* h3j;
static TH1D* h3k;


static TH2D* h2a1;
static TH2D* h2a2;

static TH1D* h2a3;
static TH1D* h2a4;

static TH2D* h2c;

static TH1D* h2c2;
static TH1D* h2c3;

static TH2D* h22a1;
static TH2D* h22a2;
static TH2D* h22a3;
static TH2D* h22a4;

static TH2D* h22c;

static TH1D* h6a;
static TH2D* h6b;
static TH2D* haux4;

static TH2D* haux5;
static TH2D* haux6;

static TH1D* haux7;
static TH1D* haux8;
static TH1D* haux9;

static TH1D* haux11;
static TH1D* haux12;
static TH1D* haux13;

static TH1D* hEc_phi;
static TH1D* hEc_theta;
static TH1D* hEb_phi;
static TH1D* hEb_theta;
static TH1D* hEce;
static TH1D* hEbe;
static TH1D* hEc_sum;
static TH1D* hEb_sum;
static TH1D* hEc_trks;
static TH1D* hEb_trks;
static TH1D* hEmc_sum;

// static TH2D * haux10;

static TH1D* htheta1;


inline Double_t DtoR(Double_t ang) {
   return ang/180*M_PI;
}
inline Double_t RtoD(Double_t ang) {
   return ang*180/M_PI;
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
   };

   int Np = sizeof(ListRuns)/sizeof(ListRuns[0]);

   int absrun = abs(run);

   // we need some values any way
   double Ecms   = 3.097;
//    double Spread = 1000.;
   Lumi          = 0;

   int i = 0;
   for(; i < Np; i++) {
      if ( absrun >= ListRuns[i].runS && absrun <= ListRuns[i].runE ) {
         Ecms   = ListRuns[i].Ebeam  * 1.e-3;  // MeV -> GeV
//          Spread = ListRuns[i].Spread * 1.e-6;  // KeV -> GeV
         Lumi   = ListRuns[i].Lumi;
         break;
      }
   }

   if ( i != Np ) {
      // correct energy according Yadi Wang fit of J/Psi peak"
      Ecms -= 0.55e-3; // -0.55MeV
   } else {
      // run not in the list
      if ( absrun >= 39355 && absrun <= 39618 ) { // new 3080 !
         Ecms =  3.080; // no BEMS information
         Lumi = 120114.83;
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


//-----------------------------------------------------------------------------
void BhabhaStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " BhabhaStartJob() " << endl;
   }

   // reserve
   his1.resize(100,(TH1D*)0);
//    his1.resize(100,(TH2D*)0);
   hisEC.resize((maxR-minR+1),(TH1D*)0);
   hisB.resize((maxR-minR+1),(TH1D*)0);


   // book histograms

   htheta1 =  new TH1D("htheta1", "90 - (theta(1)-theta(2))/2", 360, 0,
         180);

   h1a1 =  new TH1D("h1a1_Theta1p2",
         "1. a.1) : #theta_{1}+#theta_{2} ( 20<#Theta*<40 ) ", 200, 178, 182);
   h1a2 =  new TH1D("h1a2_Theta1p2",
         "1. a.2)  : #theta_{1}+#theta_{2} ( 40<#Theta*<90 ) ", 200, 178, 182);
   h1b1 =  new TH1D("h1b1_Phi1m2",
         "1. b.1) : |#phi_1-#phi_2| ( 20<#Theta*<40 ) ", 200, 177, 183);
   h1b2 =  new TH1D("h1b2_Phi1m2",
         "1. b.2) : |#phi_1-#phi_2|  ( 40<#Theta*<90 ) ", 100, 178, 182);
   //1c theta*-bins
   h1c = new TH2D("h1c_Theta1p2_vs_ThetaSt",
         "1. c.1) ;#theta*;#theta_{1}+#theta_{2}",
         70, 20, 90, 100, 174, 186); //setOptionStat 1100011  !!!0.6
   //1d phi*-bins
   h1d = new TH2D("h1d_Theta1p2_vs_PhiSt",
         "1. d.1);#phi*;#theta_{1}+#theta_{2}",
         180, -90, 90, 100,175, 185);

   //1. e)  Plot |Phi1-Phi2| in bins of Theta*
   h1e = new TH2D("h1e_Phi1m2_vs_ThetaSt",
         "1. e.1) ;#Theta*;|#phi_1-#phi_2|",
         140, 20, 90,   200, 178,  182); //setOptionStat 1100011  !!!0.6

   // 1. f) Plot |Phi1-Phi2| in bins of Phi, separately
   h1f1 = new TH2D("h1f1_Phi1m2_vs_PhiSt",
         "1. f.1) ( 20<#Theta*<40 )  ;#phi*;#|#phi_1-#phi_2|",
         180, -90, 90, 200, 176, 184); //setOptionStat 1100011  !!!0.6

   h1f2 = new TH2D("h1f2_Phi1m2_vs_PhiSt",
         "1. f.2) ( 40<#Theta*<90 )  ;#phi*;#|#phi_1-#phi_2|",
         90,  -90, 90, 100, 179, 181); //setOptionStat 1100011  !!!0.6


   h3a1 = new TH2D("h3a1_DedxTheta", ";#theta;dE/dx", 200, 0, 180, 100,
         0.6,2.);
   h3a2 = new TH2D("h3a2_DedxPhi_fwd",
         "Endcap only:  |#Theta-90|>50;#varphi;dE/dx", 400, -180, 180, 100,0.6,
         2.);
   h3a3 = new TH2D("h3a3_DedxPhi",
         "Barrel only:  |#Theta-90|<50 ;#varphi;dE/dx", 400, -180, 180, 100,
         0.6,2.);

   h3d1 = new TH2D("h3d1_BetaTheta", ";#theta;#beta", 200, 0, 180, 400,
         -0.1,1.5);

   h3d2 = new TH2D("h3d2_BetaPhi",
         "Barrel only:  |#Theta-90|<50 ;#varphi;#beta", 400, -180, 180, 400,
         -0.1,1.5);
   h3d2a = new TH2D("h3d2a_BetaPhi",
         "Endcap only:  |#Theta-90|>50 ;#varphi;#beta", 400, -180, 180, 400,
         -0.1,1.5);

   h3h1 = new TH2D("h3h1_NMucLayersTheta", ";#theta;MUC Layers", 50, 0,
         180, 10,0,10.);
   h3h2 = new TH2D("h3h2_NMucLayersPhi", ";#varphi;MUC Layers", 100,
         -180, 180, 10,0,10.);

   h3j = new TH2D("h3j_dT_Theta",";#theta; T - T_{exp}",90, 0, 180,1200,
         -2,10);
   h3k = new TH1D("h3k_dT","T - T_{exp}",1200,-2,10);


   h2a1 = new TH2D("h2a1_pelTheta", ";#theta;p(electron)", 50, 0, 180,
         100,0,M_PI);
   h2a2 = new TH2D("h2a2_pposTheta", ";#theta;p(positron)", 50, 0, 180,
         100,0,M_PI);

   h2a3 = new TH1D("h2a3_p_fwd","(|#Theta-90| > 50 );p",800,-2.1,2.1);
   h2a4 = new TH1D("h2a4_p","(|#Theta-90|<50 );p",800,-2.1,2.1);


   h2c = new TH2D("h2c_ETheta", ";#theta;E;", 90, 0, 90, 400,0,2.1);

   h2c2 = new TH1D("h2c2_E_fwd","(|#Theta-90| > 50 );E",400,0,2.1);
   h2c3 = new TH1D("h2c3_E","(|#Theta-90|<50 );E",400,0,2.1);



   h22a1 = new TH2D("h22a1_pelPhi", ";#varphi;p(electron)", 360, -180,
         180, 400,0,2.2);
   h22a2 = new TH2D("h22a2_pposPhi", ";#varphi;p(positron)", 360, -180,
         180, 400,0,2.2);

   h22a3 = new TH2D("h22a3_pPhi_fwd", "(|#Theta-90| > 50 );#varphi;p",
         360, -180, 180, 800,-2.2,2.2);
   h22a4 = new TH2D("h22a4_pPhi", "(|#Theta-90| < 50 );#varphi;p", 360,
         -180, 180, 800,-2.2,2.2);


   h22c = new TH2D("h22c_EPhi", ";#varphi;E;", 360, -180, 180, 100,0,
         2.2);

//   haux1 = new TH2D("haux1_Phi1_Phi2", ";#varphi_1;#varphi_2;", 100,-M_PI,M_PI, 100,-M_PI,M_PI);

   haux1 = new TH1D("haux1_x", ";x_{poca}", 400, -1, 1 );
   haux2 = new TH1D("haux2_y", ";y_{poca}", 400, -1, 1 );
   haux3 = new TH1D("haux3_z", ";z", 400, -10, 10 );
   haux4 = new TH2D("haux4_x_y", ";y;x", 100, -1, 1, 100, -1,1 );

   haux5 = new TH2D("haux5_x_z", ";z;x", 100, -5, 5, 100, -1,1 );
   haux6 = new TH2D("haux6_y_z", ";z;x", 100, -5, 5, 100, -1,1 );

   haux7 = new TH1D("haux7", "raw px", 800, -0.08, 0.08 );
   haux8 = new TH1D("haux8", "raw py", 500, -0.05, 0.05 );
   haux9 = new TH1D("haux9", "raw pz", 800, -0.2, 0.2 );


   haux11 = new TH1D("haux11", "corrected px", 800, -0.08, 0.08 );
   haux12 = new TH1D("haux12", "corrected py", 500, -0.05, 0.05 );
   haux13 = new TH1D("haux13", "corrected pz", 800, -0.2, 0.2 );


   h6a = new TH1D("h6a_sigmaPt_inv","sigma(1/Pt)",1000,0,0.1);
   h6b = new TH2D("h6b_relsigma_Pt",";Pt;sigma(Pt)/Pt",400,-2,2,500,0,
         0.1);

   hEc_phi =  new TH1D("E_c_phi", "|#phi_1-#phi_2| EndCap", 200, 177,
         183);
   hEc_theta =  new TH1D("E_c_theta", "#theta_{1}+#theta_{2} EndCap",
         200, 178, 182);
   hEb_phi =  new TH1D("E_b_phi", "|#phi_1-#phi_2| Barell", 200, 177,
         183);
   hEb_theta =  new TH1D("E_b_theta", "#theta_{1}+#theta_{2} Barell",
         200, 178, 182);

   hEce = new TH1D("E_c_tracks","E of tracks EndCap",400,1.,2.);
   hEbe = new TH1D("E_b_tracks","E of tracks Barell",400,1.,2.);

   hEc_sum = new TH1D("E_c_sum","E total EndCap",400,2.5,3.5);
   hEb_sum = new TH1D("E_b_sum","E total Barell",400,2.5,3.5);
   hEc_trks = new TH1D("E_c_trks","E Kalman EndCap",400,2.5,3.5);
   hEb_trks = new TH1D("E_b_trks","E Kalman Barell",400,2.5,3.5);

   hEmc_sum = new TH1D("Emc_sum","EMC total energy",400,2.5,3.5);

   char buf1[100],buf2[100];
   for (int j=0; j<(maxR-minR+1); j++) {
      sprintf(buf1,"E_ec_%d",minR+j);
      sprintf(buf2,"E  total EndCap, run # %d",minR+j);
      hisEC[j] = new TH1D(buf1,buf2,400,2.5,3.5);
      sprintf(buf1,"E_b_%d",minR+j);
      sprintf(buf2,"E  total Barrel, run # %d",minR+j);
      hisB[j] =  new TH1D(buf1,buf2,400,2.5,3.5);
   }

   int i = 0;
   his1[i++] = h1a1;
   his1[i++] = h1a2;
   his1[i++] = h1b1;
   his1[i++] = h1b2;

   his1[i++] = haux1;
   his1[i++] = haux2;
   his1[i++] = haux3;
   his1[i++] = h2a3;
   his1[i++] = h2a4;
   his1[i++] = h2c2;
   his1[i++] = h2c3;
   his1[i++] = h6a;
   his1[i++] = h3k;

   his1[i++] = haux7;
   his1[i++] = haux8;
   his1[i++] = haux9;
   his1[i++] = haux11;
   his1[i++] = haux12;
   his1[i++] = haux13;



   his1[++i] = h1c;
   his1[++i] = h1d;
   his1[++i] = h1e;
   his1[++i] = h1f1;
   his1[++i] = h1f2;
   his1[++i] = h3a1;
   his1[++i] = h3a2;
   his1[++i] = h3a3;

   his1[++i] = h3d1;

   his1[++i] = h3d2;
   his1[++i] = h3d2a;

   his1[++i] = h3h1;
   his1[++i] = h3h2;

   his1[++i] = h2a1;
   his1[++i] = h2a2;
   his1[++i] = h2c;

   his1[++i] = h22a1;
   his1[++i] = h22a2;
   his1[++i] = h22c;
   his1[++i] = h22a3;
   his1[++i] = h22a4;
   his1[++i] = h6b;
   his1[++i] = haux4;
   his1[++i] = haux5;
   his1[++i] = haux6;
   his1[++i] = h3j;

   his1[++i] = htheta1;

   his1[++i] = hEc_phi;
   his1[++i] = hEb_phi;
   his1[++i] = hEc_theta;
   his1[++i] = hEb_theta;
   his1[++i] = hEce;
   his1[++i] = hEbe;
   his1[++i] = hEc_sum;
   his1[++i] = hEb_sum;
   his1[++i] = hEc_trks;
   his1[++i] = hEb_trks;
   his1[++i] = hEmc_sum;

   // register in selector to save in given directory
   selector->RegInDir(&his1,"Bhabha");
   VecObj hisECR(hisEC.begin(),hisEC.end());
   selector->RegInDir(&hisECR,"BhabhaS");
   VecObj hisBR(hisB.begin(),hisB.end());
   selector->RegInDir(&hisBR,"BhabhaS");
}

//-----------------------------------------------------------------------------
bool BhabhaEvent(ReadDst* selector,
      TEvtHeader* m_TEvtHeader,
      TDstEvent* m_TDstEvent,
      TEvtRecObject* m_TEvtRecObject,
      TMcEvent* m_TMcEvent,
      TTrigEvent* m_TTrigEvent,
      TDigiEvent* m_TDigiEvent,
      THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " UserTestEvent() " << endl;
   }

   int runNo = m_TEvtHeader->getRunId();
   // const double e_beam = 1.5462;
   double lum = 0;
   double e_beam =  GetEcms(runNo, lum) /2;


   // Double_t beta;
   int i;
   char flag;
   DstEvtRecTracks* trackInfo[2] ;
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol(); // nef

   TMdcTrack*  mdc  [2];
//nef     TMdcKalTrack* kal  [2];
   BMdcKalTrack* kal  [2];
   Double_t mdc_path[2];
   TEmcTrack*  emc  [2];
   TMdcDedx*   dedx [2];
   TMucTrack*  muc    [2];
   Double_t tof  [2];

   TLorentzVector p4 [2];
   Double_t phi[2];
   Double_t theta[2];
   Double_t energy[2];
   Double_t p[2];

   Double_t   kpx[2];
   Double_t   kpy[2];
   Double_t   kpz[2];
   // Double_t   kp[2];


   Double_t tof_sum;
   Double_t tof_path_sum;

   // Double_t pxy;

   int tof_cnt;

   /* 1. Select Bhabha events: 2 tracks, each with E>0.8*ebeam, p>0.8*pbeam
      a) Plot distribution of Theta1+Theta2, separately for  20<Theta*<40, 40<Theta*<90
                Theta*=90-|90-0.5*(180+Theta1-Theta2)| */

   /* 1. Select Bhabha events: 2 tracks */



   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif
   for (i=0; i< evtRecTrkCol->GetEntriesFast(); i++) {
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);


      pid->init();

      pid->setMethod(pid->methodProbability());
      pid->setChiMinCut(4);
      pid->setRecTrack(itTrk);
      pid->usePidSys(pid->useDedx() |
            pid->useTof1() |
            pid->useTof2() |
            pid->useTofE() |
//          pid->useEmc()  |
            pid->useMuc()  ); // PID sub-system

      pid->identify(pid->onlyElectron() |
            pid->onlyMuon());    // seperate Pion/Kaon
      pid->calculate(runNo);
      if(!(pid->IsPidInfoValid())) {
         continue;
      }

      //  cout << "Prob e = " <<  pid->probElectron() <<
      //    " Prob mu = " <<  pid->probMuon() << endl;
   }



   if  (m_TEvtRecObject->getEvtRecTrackCol()->GetEntriesFast()  == 2) {
      //cout << "begin" << endl;

      for (i=0; i<2; i++) {
// nef         trackInfo[i] = new DstEvtRecTracks((TEvtRecTrack*) m_TEvtRecObject->getEvtRecTrackCol()->At(i), m_TDstEvent);
         trackInfo[i] = (DstEvtRecTracks*) evtRecTrkCol->At(i);
         mdc[i]  = trackInfo[i]->mdcTrack(); //accessor
         emc[i]  = trackInfo[i]->emcShower(); //accessor
         dedx[i] = trackInfo[i]->mdcDedx();
         muc[i]  = trackInfo[i]->mucTrack();
         kal[i]  = trackInfo[i]->mdcKalTrack();
         if(kal[i]) {
            kal[i]->setPidType(RecMdcKalTrack::electron);
         }
         // nef: and use dr,fi0,kappa,dz,tanl instead of getZHelixE
         // nef: or kal[i]->getTMdcKalTrack()->getZHelixE()

         if ( emc[i] ) {

            if ( kal[i] ) {
               // pxy = kal[i]->pxy();

               p4[i].SetPxPyPzE(
                     kal[i]->px(),kal[i]->py(),kal[i]->pz(),
                     emc[i]->energy() );



            } else if (mdc[i]) {
               p4[i].SetPxPyPzE(
                     mdc[i]->px(),mdc[i]->py(),mdc[i]->pz(),
                     emc[i]->energy() );

            }
            if ((mdc[i])||(kal[i])) {
               // kp[i] = p4[i].P();
               kpx[i] = p4[i].Px();
               kpy[i] = p4[i].Py();
               kpz[i] = p4[i].Pz();



               // p4[i].Boost(-1.07E-2,0,-2E-4);
               p4[i].Boost(-2*sin(beam_angle/2),0.,0.);
               theta[i] = p4[i].Theta();
               phi[i] =  p4[i].Phi();
               energy[i] = p4[i].E();
               p[i]   =  p4[i].P();
            }
         }

         tof_cnt = 0;
         tof_sum = 0.0;
         tof_path_sum = 0.0;

         // calculate the average tof
         for(vector<RecTofTrack* >::const_iterator it =
                     trackInfo[i]->tofTrack().begin();
               it != trackInfo[i]->tofTrack().end(); ++it) {
            // do the mysterious check to reject wrong TOF events inspired by numerous BOS examples
            if  ( TofHitStatus::is_counter(  ((RecTofTrack*)(*it))->status() ) ) {
               tof_sum += 0.0 + ((RecTofTrack*)(*it))->tof();
               tof_path_sum += 0.0 + ((RecTofTrack*)(
                                 *it))->path(); // there is not a lot of sense in average path
               tof_cnt += 1;

            }
         }
         if ( tof_cnt ) {
            tof[i] = tof_sum / tof_cnt;
            mdc_path[i] = tof_path_sum / tof_cnt;
         } else {
            tof[i] = -100;
            mdc_path[i] =  0;
         }
      }

      flag = 0;
      char flags[2] = {0,0};
      // int j;

      for (i=0; i<2; i++) {
         if (((mdc[i])||(kal[i])))
            if (emc[i])
               if (( p[i] > 0.8*e_beam ) &&
                     ( energy[i] > 0.8*e_beam ) ) {
                  flag+=1;
                  flags[i] = 1;
               }
      }

      if (flag > 0 ) { //at least one with E>0.8*ebeam, p>0.8*pbeam
         // 2)
         if ( ((mdc[0])||(kal[0])) && ((mdc[1])||(kal[1])) && emc[0]
               && emc[1] ) {




            Double_t t1_p_t2 = theta[0] +  theta[1] ;
            Double_t p1_m_p2 = TMath::Abs(phi[0] -  phi[1]) ;

            if ( ( TMath::Abs(t1_p_t2 - M_PI) / M_PI  < 0.01 ) &&
                  ( TMath::Abs(p1_m_p2 - M_PI) / M_PI  <
                        0.01 ) ) {// Theta1+Theta2, Phi1-Phi2 close to 180
               for (i=0; i<2; i++) {
                  if (flags[i]) {  // for *one track* p>0.8*pbeam, E>0.8*ebeam.; [1-i] - another track
                     if (mdc[1-i]->charge() == -1 ) {
                        h2a1->Fill(RtoD(theta[1-i]),
                              p[1-i] );//2.a) <p(electron)> versus Theta,
                        h22a1->Fill(RtoD(phi[1-i]), p[1-i] );

                     } else {
                        h2a2->Fill(RtoD( theta[1-i]),
                              p[1-i] );//2.a) <p(positron)> versus Theta
                        h22a2->Fill(RtoD(phi[1-i]), p[1-i] );
                     }

                     if ( RtoD( TMath::Abs( theta[1-i] - M_PI/2 ))  > 50 ) { //forward
                        h2a3->Fill( mdc[1-i]->charge()* p[1-i] );
                        h2c2->Fill( energy[1-i] );
                        h22a3->Fill( RtoD(phi[1-i]),  p[1-i] * mdc[1-i]->charge() );
                     } else {
                        h2a4->Fill( mdc[1-i]->charge()* p[1-i] );
                        h2c3->Fill( energy[1-i] );
                        h22a4->Fill( RtoD(phi[1-i]),  p[1-i] * mdc[1-i]->charge() );
                     }


                     h2c -> Fill(RtoD(theta[1-i]), energy[1-i] ); // <E> versus Theta
                     h22c -> Fill(RtoD(phi[1-i]), energy[1-i] );

                  }
               }
            }
         }

         if (flag == 2) {     /* each with E>0.8*ebeam, p>0.8*pbeam */
            //cout << p4[0].P() <<" , " <<  mdc[0]->p() << endl;

            //------------------------------------------------------------------------------------E endcap/ E barrel for kalman only--------------------------------

            if ( kal[0] && kal[1] && emc[0] && emc[1]
                  && (kal[0]->charge() * kal[1]->charge() == -1) ) {


               double Esum = energy[0] + energy[1];
               double Phi = RtoD(TMath::Abs(phi[0] - phi[1]));
               double Theta =  RtoD(theta[0] + theta[1]);
               double Etrks = sqrt(kpx[0]*kpx[0]+kpy[0]*kpy[0]+kpz[0]*kpz[0]) +
                     sqrt(kpx[1]*kpx[1]+kpy[1]*kpy[1]+kpz[1]*kpz[1]);

               bool good_event = false;
               if ( RtoD( TMath::Abs( theta[0] - M_PI/2 ))  > 50 &&
                     RtoD( TMath::Abs( theta[1] - M_PI/2 ))  > 50     ) { //endcap
                  hEc_phi->Fill(Phi);
                  hEc_theta->Fill(Theta);
                  if ( fabs(Phi-180) < 0.5 && fabs(Theta-180) < 0.5 ) {
                     good_event = true;
                     hEce->Fill( energy[0] );
                     hEce->Fill( energy[1] );
                     hEc_sum->Fill( Esum );
                     hEc_trks->Fill( Etrks );
                     hisEC[abs(runNo)-minR]->Fill( Etrks );

                  }
               } else if (
                     RtoD( TMath::Abs( theta[0] - M_PI/2 ))  < 50 &&
                     RtoD( TMath::Abs( theta[1] - M_PI/2 ))  < 50    ) { // barell
                  hEb_phi->Fill(Phi);
                  hEb_theta->Fill(Theta);
                  if ( fabs(Phi-180) < 0.5 && fabs(Theta-180) < 0.5 ) {
                     good_event = true;
                     hEbe->Fill( energy[0] );
                     hEbe->Fill( energy[1] );
                     hEb_sum->Fill( Esum );
                     hEb_trks->Fill( Etrks );
                     hisB[abs(runNo)-minR]->Fill( Etrks );
                  }
               }

               if ( good_event) {
                  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
                  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

                  double Emc_sum = 0;
                  for(int i = 0 ; i< evtRecEvent->totalTracks(); i++) {
                     DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
                     if(!itTrk->isEmcShowerValid()) {
                        continue;
                     }
                     RecEmcShower* emcTrk = itTrk->emcShower();
                     // if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) continue;
                     double eraw = emcTrk->energy();
                     Emc_sum += eraw;
                  }
                  hEmc_sum->Fill(Emc_sum);
               } // end if

            }//------------end if
            //------------------------------------------------------------------------------------------------------------------------------------------------------



            double theta_mean = 90 - 0.5 * RtoD(theta[0]-theta[1]);
            htheta1->Fill(theta_mean);

            Double_t t1_p_t2 = theta[0] +  theta[1] ;
            Double_t p1_m_p2 = TMath::Abs(phi[0] -  phi[1]) ;

            Double_t th_star = M_PI/2 - TMath::Abs(theta[0]  -
                        theta[1])/2;  //  Theta*=90-|90-0.5*(180+Theta1-Theta2)|
            Double_t phi_star = 0.5* (phi[0] + phi[1]); //  Phi*=0.5*(Phi1+Phi2)


            // Plot distribution of Theta1+Theta2, |Phi1-Phi2|,
            // Plot |Phi1-Phi2| vs Phi*, all separately
            // for 20<Theta*<40, 40<Theta*<90.

            if ( (th_star > DtoR(20)) && (th_star <  DtoR(40)) ) {
               h1a1->Fill(  RtoD(theta[0] +
                           theta[1]) ); // 1. a.1) Theta1+Theta2, 20<Theta*<40
               h1b1->Fill(RtoD( TMath::Abs( phi[0] -
                                 phi[1]))); // 1. b.1) |Phi1-Phi2|, 20<Theta*<40
               h1f1->Fill(RtoD(phi_star),
                     RtoD(TMath::Abs( phi[0] - phi[1]))); // 1. f.1) Plot |Phi1-Phi2| in bins of Phi

            };
            if ( (th_star > DtoR(40)) && (th_star <  DtoR(90)) ) {
               h1a2->Fill(  RtoD(theta[0] +  theta[1] ));
               h1b2->Fill(RtoD( TMath::Abs( phi[0] - phi[1])));
               h1f2->Fill(RtoD(phi_star), RtoD(TMath::Abs( phi[0] - phi[1])));
            };

            //haux7->Fill(  sqrt(  4 + ( 1.0 / tan(phi[0]) + 1.0 / tan(phi[1]) ) * ( 1.0 / tan(phi[0]) + 1.0 / tan(phi[1]) ) ) ,
//                                1.0/ tan(phi[0]) - 1.0/ tan(phi[1])  );


            h1c->Fill(RtoD(th_star),
                  RtoD(theta[0] +  theta[1])); // 1. c) Plot Theta1+Theta2 in small bins of Theta*



            h1d->Fill(RtoD(phi_star),
                  RtoD( theta[0] +  theta[1])); // 1. d) Do the same in bins of Phi*

            h1e->Fill(RtoD(th_star),
                  RtoD(TMath::Abs(phi[0] -  phi[1]))); //1. e) Plot |Phi1-Phi2| in bins of Theta*



            if ( ( TMath::Abs(t1_p_t2 - M_PI) / M_PI  < 0.01 ) &&
                  ( TMath::Abs(p1_m_p2 - M_PI) / M_PI  <
                        0.01 ) ) { // Theta1+Theta2, Phi1-Phi2 close to 180
               if ( (theta_mean > 70) && (theta_mean < 110) ) {
                  return 1;
               }

               if (kal[0] && kal[1]) {
                  haux7->Fill( ( kpx[0] + kpx[1] )  /  2 / e_beam );
                  haux8->Fill( ( kpy[0] + kpy[1] )  /  2 / e_beam );
                  haux9->Fill( ( kpz[0] + kpz[1] )  /  2 / e_beam );

                  haux11->Fill( ( p4[0].Px() +  p4[1].Px() )  /  2 / e_beam );
                  haux12->Fill( ( p4[0].Py() +  p4[1].Py() )  /  2 / e_beam );
                  haux13->Fill( ( p4[0].Pz() +  p4[1].Pz() )  /  2 / e_beam );

               }

               for (i = 0; i < 2; i++) {
                  if( dedx[i] ) {
                     h3a1->Fill( RtoD(theta[i]),
                           dedx[i]->normPH() );   // 3.a-b) <dE/dx> versus Theta, <dE/dx> versus Phi
                  }

                  //**** TRecMdcKalTrack - where is pathSM ?!
                  Double_t dT = tof[i] - mdc_path[i]*1.0 /
                        (speed_of_light/1E7); // T - T_expected
                  h3j->Fill( RtoD(theta[i]),dT);
                  h3k->Fill( dT);


                  h3d1->Fill( RtoD(theta[i]),
                        mdc_path[i]*1.0 / tof[i] /
                        (speed_of_light/1E7) ); // def) beta versus Theta, versus Phi

                  if ( RtoD( TMath::Abs( theta[i] - M_PI/2 ))  < 50 ) { //barell
                     h3d2->Fill( RtoD(phi[i]),
                           mdc_path[i]*1.0 / tof[i] /(speed_of_light/1E7) );
                     if( dedx[i] ) {
                        h3a3->Fill(RtoD( phi[i]), dedx[i]->normPH() );
                     }
                  } else { //endcap
                     h3d2a->Fill( RtoD(phi[i]),
                           mdc_path[i]*1.0 / tof[i] /(speed_of_light/1E7) );
                     if( dedx[i] ) {
                        h3a2->Fill(RtoD( phi[i]), dedx[i]->normPH() );
                     }
                  }

                  if( muc[i] ) { // nef
                     h3h1->Fill( RtoD(theta[i]), muc[i]->numLayers() );
                     h3h2->Fill( RtoD(phi[i]),
                           muc[i]->numLayers() );     //ghi) number of MUC layers versus Theta, versus Phi
                  }
                  haux1 -> Fill( mdc[i]->x() );
                  haux2 -> Fill( mdc[i]->y() );
                  haux3 -> Fill( mdc[i]->z() );
                  haux4 -> Fill( mdc[i]->y(), mdc[i]->x() );
                  haux5 -> Fill( mdc[i]->z(), mdc[i]->x() );
                  haux6 -> Fill( mdc[i]->z(), mdc[i]->y() );

                  if (kal[i]) {
                     //                      h6a-> Fill( sqrt( kal[i]->getZErrorE(2,2) ) );
                     //                      h6b-> Fill( 1./kal[i]->getZHelixE(2) ,sqrt( kal[i]->getZErrorE(2,2) ) / fabs(kal[i]->getZHelixE(2)) );
                     HepSymMatrix Err = kal[i]->err();
                     h6a-> Fill( sqrt( Err[2][2] ) );
                     h6b-> Fill( 1./kal[i]->kappa(),
                           sqrt( Err[2][2] ) / fabs(kal[i]->kappa()) );
                  }

               }
            }


            //~ return 1;
         }

         return 0; // do no save event if at least one of track is bhabha
      }

// nef:     for (i=0; i<2; i++ )
//          if (trackInfo[i])
//             delete trackInfo[i];


   }




   return 0;
}

//-----------------------------------------------------------------------------
void BhabhaEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " UserTestEndJob() " << endl;
   }
}

#ifdef __cplusplus
}
#endif

//~ h1c->FitSlicesY();
//~ h1d->FitSlicesY();
//~ h1e->FitSlicesY();
//~ h1f1->FitSlicesY();
//~ h1f2->FitSlicesY();

//~ h3a1->ProfileX();
//~ h3a2->ProfileX();
//~ h3d1->ProfileX();
//~ h3d2->ProfileX();
//#endif
