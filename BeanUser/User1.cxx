//////////////////////////////////////////////////////////////////////
//
// User1
//
// This is example of user functions (here name is User1)
// IMPORTANT: It MUST contain "name"Event() function.
// Two functions "name"StartJob() and "name"EndJob() are optional.
//
// The return value of Event() function used only if
// option -o (define output ROOT tree file name) is defined:
//      true  -- save this event in output ROOT
//      false -- skip event
//
// NOTE: Please do not forget to use `extern "C"` directive for
//       functions.
//       The macros BeanUserShared_EXPORT is defined in file
//       "DLLDefines.h". This is important for WIN programming
//       but you may want to skip it in UNIX.
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>

#include <TH1.h>
#include <TH2.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

// #include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;

inline Double_t DtoR(Double_t ang){ return ang/180*M_PI; }

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void User1StartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " User1StartJob() " << endl;

   // reserve
   his1.resize(100,(TH1D*)0);
   his2.resize(100,(TH2D*)0);

   // book histograms:
   his1[1] = new TH1D("nmdctrk","Ntracks in MDC", 20,-0.5,19.5);
   his2[2] = new TH2D("dedx_nmdc","Ndedx vs Nmdc",
         20,-0.5,19.5, 20,-0.5,19.5);
   his2[3] = new TH2D("nmuc_nmdc","Nmuc vs Nmdc",
         20,-0.5,19.5, 20,-0.5,19.5);
   his2[4] = new TH2D("ntof_nmdc","Ntof vs Nmdc",
         20,-0.5,19.5, 20,-0.5,19.5);

   his1[2] = new TH1D("nmdcdedx","Nmdc Dedx", 20,-0.5,19.5);
   his1[3] = new TH1D("nmuctrk","Nmuc tracks", 20,-0.5,19.5);
   his1[4] = new TH1D("ntoftrk","Ntof tracks", 20,-0.5,19.5);

   his1[8] = new TH1D("th", "#theta*", 50, 0, M_PI/2*1.2);
   his1[9] = new TH1D("ph", "#varphi*", 60, -M_PI*1.2, M_PI*1.2);

// theta[0]+theta[1] - mTh,MTh hTheta12m
   //1a
   his1[10] = new TH1D("hTheta12m", "", 50, 3, 3.3);
   his1[11] = new TH1D("hTheta12M", "1:#theta_{1}+#theta_{2}",
         50, 3, 3.3);
// phi[0]+phi[1] - mTh, Mth
   //1b
   his1[12] = new TH1D("hPhi12m", "hPhi12m", 50, 3, 3.3);
   his1[13] = new TH1D("hPhi12M", "1:|#varphi_{1}-#varphi_{2}|",
         50, 3, 3.3);

// theta[0]+theta[1]:
   //1c theta*-bins
   his2[10] = new TH2D("Theta12Th", ";#theta*;#theta_{1}+#theta_{2}",
         50, 0.6, M_PI/2, 50, 3, 3.3); //setOptionStat 1100011  !!!0.6
   //1d phi*-bins
   his2[11] = new TH2D("Theta12Ph", ";#phi*;#theta_{1}+#theta_{2}",
         50, -M_PI/2, M_PI/2, 50, 3, 3.3);

   //1e theta*-bins
   his2[12] = new TH2D("Phi12Th",
         ";#theta*;|#varphi_{1}-#varphi_{2}|",
         50, 0.6, M_PI/2, 50, 3, 3.3);
   //1f phi*bins
   his2[13] = new TH2D("Phi12Phm",
         ";#varphi*;|#varphi_{1}-#varphi_{2}|",
         50, -M_PI/2, M_PI/2, 50, 3, 3.3);
   his2[14] = new TH2D("Phi12PhM",
         ";#varphi*;|#varphi_{1}-#varphi_{2}|",
         50, -M_PI/2, M_PI/2, 50, 3, 3.3);

//----------------------------p----------------------------

   his2[21] = new TH2D("h2pPTheta", ";#theta;p",
         50, 0, M_PI, 50, 0.5, 2.4);
   his2[22] = new TH2D("h2ePTheta", ";#theta;p",
         50, 0, M_PI, 50, 0.5, 2.4);

//----------------------------e----------------------------

   his2[31] = new TH2D("h2pETheta", ";#theta;E",
         50, 0, M_PI, 50, 0.5, 2.4);
   his2[32] = new TH2D("h2eETheta", ";#theta;E",
         50, 0, M_PI, 50, 0.5, 2.4);

//---------------------------dedx---------------------------

   his2[41] = new TH2D("hDedxP", ";p;dedx", 50, 1.4, 2, 100,0,2.);

   his2[42] = new TH2D("hDedxTheta", ";#theta;dedx",
         50, 0, M_PI, 100,0,2.);
   his2[43] = new TH2D("hDedxPhi", ";#varphi;dedx",
         50, -M_PI, M_PI, 100,0,2.);

//---------------------------beta---------------------------

   his2[51] = new TH2D("hBetaP", ";p;#beta", 50, 1.4, 2, 50, 0, 1.4);

   his2[52] = new TH2D("hBetaTheta", ";#theta;#beta",
         50, 0, M_PI, 50, 0, 1.4);
   his2[53] = new TH2D("hBetaPhi", ";#varphi;#beta",
         50, -M_PI, M_PI, 50, 0, 1.4);

//---------------------------MUC----------------------------

   his1[61] = new TH1D("hMucLayers", "MUC Layers", 10, 0, 10);
   his2[62] = new TH2D("hMucLayersTheta", ";#theta;N",
         50, 0, M_PI, 10, 0, 10);
   his2[63] = new TH2D("hMucLayersPhi", ";#phi;N",
         50, -M_PI, M_PI, 10, 0, 10);


   // register in selector to save in given directory
   VecObj his1o(his1.begin(),his1.end());
   selector->RegInDir( his1o     ,"User1");

   VecObj his2o(his2.begin(),his2.end());
   selector->RegInDir( his2o     ,"User1");

}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool User1Event(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " User1Event() " << endl;

   const double e_beam = 1.8495;

   const TObjArray* m_mdcTrackCol = m_TDstEvent->getMdcTrackCol();
   const TObjArray* m_mdcDedxCol  = m_TDstEvent->getMdcDedxCol();
   const TObjArray* m_mucTrackCol = m_TDstEvent->getMucTrackCol();
   const TObjArray* m_tofTrackCol = m_TDstEvent->getTofTrackCol();

   int NmdcTracks = m_mdcTrackCol->GetEntries();
   int NmdcDedx   = m_mdcDedxCol->GetEntries();
   int NmucTracks = m_mucTrackCol->GetEntries();
   int NtofTracks = m_tofTrackCol->GetEntries();

   his1[1]->Fill(NmdcTracks);
   his2[2]->Fill(NmdcTracks,NmdcDedx);
   his2[3]->Fill(NmdcTracks,NmucTracks);
   his2[4]->Fill(NmdcTracks,NtofTracks);
   if( NmdcTracks != 2 ) return false;

   his1[2]->Fill(NmdcDedx);
   his1[3]->Fill(NmucTracks);
   his1[4]->Fill(NtofTracks);

   vector<double> theta(2,0);
   vector<double> phi(2,0);
   vector<double> e(2,0);
   vector<double> p(2,0);
   vector<int>    charge(2,0);
   vector<double> dedx(2,0);
   vector<int>    muc_numLayers(2,0);
   vector<double> beta(2,0);

   TIter mdcTrackIter(m_mdcTrackCol);
   for(int i = 0; i < 2; i++) {
      TMdcTrack* mdcTrack = (TMdcTrack*)mdcTrackIter.Next();
      theta[i] = mdcTrack->theta();
      phi[i]   = mdcTrack->phi();
      e[i]     = mdcTrack->p(); // energy
      p[i]     = mdcTrack->p(); // total momentum
      charge[i]= mdcTrack->charge();
   }
   if( NmdcDedx == 2 ) {
     TIter mdcDedxIter(m_mdcDedxCol);
     for(int i = 0; i < 2; i++) {
        TMdcDedx* mdcDedx = (TMdcDedx*)mdcDedxIter.Next();
        dedx[i] = mdcDedx->normPH();
     }
   }
   if( NmucTracks == 1 || NmucTracks == 2 ) {
     TIter mucTrackIter(m_mucTrackCol);
     TMucTrack* mucTrack = 0;
     while ((mucTrack = (TMucTrack*)mucTrackIter.Next())) {
       int id = mucTrack->trackId();
       if( id < 0 || id > 1 ) {
         cout << "mucTrack id= " << id << endl;
         continue;
       }
       muc_numLayers[id] = mucTrack->numLayers();
     }
   }
   if( NtofTracks > 0 ) {
     TIter tofTrackIter(m_tofTrackCol);
     TTofTrack* tofTrack = 0;
     while ((tofTrack = (TTofTrack*)tofTrackIter.Next())) {
       int id = tofTrack->trackID();
       if( id < 0 || id > 1 ) {
         cout << "tofTrack id= " << id << endl;
         continue;
       }
       if( beta[id] > 0.01 ) continue;
//        Double_t beta = tof1path[i]/tof1[i]/29.9;
       beta[id] = tofTrack->path() / tofTrack->tof()/29.9;
     }
   }


   Double_t th = M_PI/2 - TMath::Abs(theta[0] - theta[1])/2; // 0,pi/2
   Double_t ph = phi[0] + phi[1];  // -pi/2 .. pi/2
   his1[8]->Fill(th);
   his1[9]->Fill(ph);

   Double_t minTh = 34, midTh = 50, maxTh = 90;
   //1a
   if ((e[0] > 0.8*e_beam) && (p[0] > 0.8*e_beam) &&
       (e[1] > 0.8*e_beam) && (p[1] > 0.8*e_beam)){
     if ((th > DtoR(minTh)) && (th < DtoR(midTh)) ){
       his1[10]->Fill(theta[0]+theta[1]);
       his1[12]->Fill(TMath::Abs(phi[0]-phi[1]));
       his2[13]->Fill(ph, TMath::Abs(phi[0]-phi[1]));
     }
     if ((th > DtoR(midTh)) && (th < DtoR(maxTh))){
       his1[11]->Fill(theta[0]+theta[1]);
       his1[13]->Fill(TMath::Abs(phi[0]-phi[1]));
       his2[14]->Fill(ph, TMath::Abs(phi[0]-phi[1]));
     }
     his2[10]->Fill(th, theta[0]+theta[1]);
     his2[11]->Fill(ph, theta[0]+theta[1]);
     his2[12]->Fill(th, TMath::Abs(phi[0]-phi[1]));
   }

   //2
   if ((TMath::Abs(theta[0]+theta[1] - M_PI) < 0.05) &&
       (TMath::Abs(TMath::Abs(phi[0]-phi[1]) - M_PI) < 0.1)){
     //
     if ((e[0] > 0.8*e_beam) && (p[0] > 0.8*e_beam)){
       if (charge[1] == 1){  //positron
         his2[21]->Fill(theta[1], p[1]);
         his2[31]->Fill(theta[1], e[1]);
       }
       if (charge[1] == -1){ //electron
         his2[22]->Fill(theta[1], p[1]);
         his2[32]->Fill(theta[1], e[1]);
       }
     }
     // the same
     if ((e[1] > 0.8*e_beam) && (p[1] > 0.8*e_beam)){
       if (charge[0] == 1){  //positron
         his2[21]->Fill(theta[0], p[0]);
         his2[31]->Fill(theta[0], e[0]);
       }
       if (charge[0] == -1){ //electron
         his2[22]->Fill(theta[0], p[0]);
         his2[32]->Fill(theta[0], e[0]);
       }
     }
   }

   //3
   if ((e[0] > 0.8*e_beam) && (p[0] > 0.8*e_beam)  &&
       (e[1] > 0.8*e_beam) && (p[1] > 0.8*e_beam)  &&
       (TMath::Abs(theta[0]+theta[1] - M_PI) < 0.05) &&
       (TMath::Abs(TMath::Abs(phi[0]-phi[1]) - M_PI) < 0.1))
     for (int i=0; i<2; i++){
        his2[41]->Fill(p[i], dedx[i]);
        his2[42]->Fill(theta[i], dedx[i]);
        his2[43]->Fill(phi[i], dedx[i]);

        his2[51]->Fill(p[i], beta[i]);
        his2[52]->Fill(theta[i], beta[i]);
        his2[53]->Fill(phi[i], beta[i]);

        his1[61]->Fill(muc_numLayers[i]);
        his2[62]->Fill(theta[i], muc_numLayers[i]);
        his2[63]->Fill(phi[i], muc_numLayers[i]);
     }

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void User1EndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " User1EndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
