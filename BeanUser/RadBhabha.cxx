#include "DLLDefines.h"         // mandatory!

#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TNtuple.h>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THashList.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "DstEvtRecTracks.h"
#include "ParticleID/ParticleID.h"
#include "EventTag/EventTagSvc.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "AbsCor/AbsCor.h"

#include "ReadDst.h"
#include <TObjString.h>
#include <assert.h>

using namespace std;


// static EventTagSvc* m_EventTagSvc;
static THashList  histograms_cache;
static TNtuple* ntp_electrons;
static AbsCor* m_abscor = 0;

#define hist(name) ( (TH1 *) histograms_cache.FindObject(name) )



#ifdef __cplusplus
extern "C" {
#endif



static void init_pid(ParticleID* pid, DstEvtRecTracks* dst_tracks, int run) {
   pid->init();
   pid->setMethod(pid->methodProbability());
   pid->setChiMinCut(4); //FIXME: 2013-12-06 which one should we use?

   pid->setRecTrack(dst_tracks);
   //~ pid->usePidSys(pid->useDedx() | pid->useTof()| pid->useTof1() | pid->useTof2() | pid->useEmc() | pid->useMuc());
   //~ pid->usePidSys(0xff);

   pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE()
                  | pid->useEmc());

}




//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void RadBhabhaStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   VecObj his1, his_aux, his1_init, his_dtag;

   //~ his_aux.push_back( new TH1D("", "miss momentum", 200, -1, 1));
   ntp_electrons = new TNtuple("ntp_electrons", "ntp_electrons",
                               "e_emc:p:cost:charge:evt_tot_emc:evt_dphi:pid_prob_e:pid_prob_sum:evt_phi_low:evt_phi_high:evt_theta_low:evt_theta_high"  );
   his1.push_back( ntp_electrons);
   his1.push_back( new TH1D("h_ncharged", "h_ncharged", 10,-0.5,9.5  ));
   his1.push_back( new TH1D("h_delta_phi", "h_delta_phi", 800,-1.9, 1.9  ));
   his1.push_back( new TH1D("h_delta_phi2", "h_delta_phi2", 800,-TMath::Pi(),
                            TMath::Pi()  ));


   his1.push_back( new TH1D("h_total_emc", "h_total_emc", 400,0, 4  ));
   his1.push_back( new TH1D("h_max_gamma_emc", "h_max_gamma_emc", 400,0, 2  ));

   his1.push_back( new TH1D("h_high_emc_e", "h_high_emc_e", 400,0, 4  ));
   his1.push_back( new TH1D("h_high_tof_e", "h_high_tof_e", 400,0, 4  ));
   his1.push_back( new TH1D("h_high_total_e", "h_high_total_e", 400,0, 4  ));

   his1.push_back( new TH1D("h_high_p_delta", "h_high_p_delta", 400,0, 2  ));
   his1.push_back( new TH1D("h_high_ep", "h_high_ep", 240, 0, 1.2));


   selector->RegInDir(his1,"RadBhabha");
   histograms_cache.AddAll(selector->GetOutputList());


   // We have to initialize DatabaseSvc
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if( (dbs->GetDBFilePath()).empty() ) {
      dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   // We have to initialize Magnetic field
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if( !(mf->GetPath()).empty() ) {
      cerr << " RhopiStartJob::WARNING:"
           << " MagneticFieldSvc has already been initialized" << endl
           << "                         path = " << mf->GetPath() << endl;
   }
   // set path to directory with magnetic fields tables
   mf->SetPath( selector->AbsPath("Analysis/MagneticField") );
   mf->UseDBFlag(false); // like in the boss program
   mf->RunMode(3); // like in the boss program


   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

}


static HepVector getPoca(BMdcKalTrack* mdcKalTrk, Hep3Vector& xorigin) {
   HepVector a = mdcKalTrk->getZHelix();
   HepSymMatrix Ea = mdcKalTrk->getZError();
   HepPoint3D point0(0.,0.,0.);
   HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
   VFHelix helixip3(point0,a,Ea);
   helixip3.pivot(IP);
   HepVector  vecipa = helixip3.a();

   return vecipa;
}

static  bool isGoodPoca(BMdcKalTrack* mdcKalTrk, Hep3Vector& xorigin) {
   HepVector  vecipa = getPoca(mdcKalTrk, xorigin);
   double dr=fabs(vecipa[0]);
   double dz=fabs(vecipa[3]);

   if (  dr>= 1.0 ) {
      return false;
   }
   if (  dz>= 10.0 ) {
      return false;
   }

   return true;
}


static Hep3Vector getVertexOrigin(int runNo) {
   Hep3Vector xorigin(0, 0, 0);

   VertexDbSvc* vtxsvc = VertexDbSvc::instance();
   vtxsvc->SetBossVer("6.6.2");
   vtxsvc->handle(runNo);

   if(vtxsvc->isVertexValid()) {
      double* dbv = vtxsvc->PrimaryVertex();
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
      //~ cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
   }

   return xorigin;

}



//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
bool RadBhabhaEvent( ReadDst* selector,
                     TEvtHeader* m_TEvtHeader,
                     TDstEvent* m_TDstEvent,
                     TEvtRecObject* m_TEvtRecObject,
                     TMcEvent* m_TMcEvent,
                     TTrigEvent* m_TTrigEvent,
                     TDigiEvent* m_TDigiEvent,
                     THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{

//    int event = m_TEvtHeader->getEventId();
   int run = m_TEvtHeader->getRunId();

   // evtRecTrkCol
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   const TEvtRecEvent*  m_evtRecEvent = m_TEvtRecObject->getEvtRecEvent();

   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   m_abscor->AbsorptionCorrection(selector);
   m_abscor->SuppressHotCrystals(selector);



   Hep3Vector xorigin = getVertexOrigin(run);



   hist("h_ncharged")->Fill(m_evtRecEvent->totalCharged());

   if (!(m_evtRecEvent->totalCharged() == 2)) {
      return false;
   };


   DstEvtRecTracks* track_1 = (DstEvtRecTracks*) evtRecTrkCol->At(0);
   DstEvtRecTracks* track_2 = (DstEvtRecTracks*) evtRecTrkCol->At(1);

   for (int i = 0; i < 2 ; ++ i) {
      DstEvtRecTracks* track = (DstEvtRecTracks*) evtRecTrkCol->At(i);
      if ( !track_1->isMdcKalTrackValid()) {
         return false;
      }
      track->mdcKalTrack()->setPidType(BMdcKalTrack::electron);

      if (track->mdcKalTrack()->charge() == 0) {
         return false;
      }

      if (!isGoodPoca(track->mdcKalTrack(), xorigin)) {
         return false;
      }

      if ( TMath::Abs(TMath::Cos(track->mdcKalTrack()->theta())) > 0.93) {
         return false;
      }
   }


   if (track_1->mdcKalTrack()->charge() != -1 * track_2->mdcKalTrack()->charge()) {
      return false;
   }

   HepVector poca_1 = getPoca(track_1->mdcKalTrack(), xorigin);
   HepVector poca_2 = getPoca(track_1->mdcKalTrack(), xorigin);
   if (TMath::Abs(poca_1[3] - poca_2[3]) >= 2.0) {
      return false;
   };

   double delta_phi = TMath::Pi() - TMath::Abs(track_1->mdcKalTrack()->phi() -
                      track_2->mdcKalTrack()->phi());
   //~ double delta_phi = track_1->mdcKalTrack()->phi() - track_2->mdcKalTrack()->phi();
   double delta_phi2 = track_1->mdcKalTrack()->phi() +
                       track_2->mdcKalTrack()->phi();

   // add delta_phi cut?
   hist("h_delta_phi")->Fill(delta_phi);
   hist("h_delta_phi2")->Fill(delta_phi2);

   double total_showers_energy = 0;
   double max_gamma_energy = 0;

   for (int i = 0; i < m_evtRecEvent->totalTracks(); ++i) {
      DstEvtRecTracks* track = (DstEvtRecTracks*) evtRecTrkCol->At(i);
      if (track->isEmcShowerValid()) {
         RecEmcShower* emcTrk = track->emcShower();
         double energy = emcTrk->energy();
         total_showers_energy += energy;

         if (i >= m_evtRecEvent->totalCharged()) {
            if (energy > max_gamma_energy) {
               max_gamma_energy = energy;
            }
         }
      }
   }



   //~ if (total_showers_energy < 3.25) { return false;}
   if (max_gamma_energy < 30E-3) {
      return false;
   }
   hist("h_total_emc")->Fill(total_showers_energy);
   hist("h_max_gamma_emc")->Fill(max_gamma_energy);







   DstEvtRecTracks* track_low;
   DstEvtRecTracks* track_high;

   if (track_1->mdcKalTrack()->p() < track_2->mdcKalTrack()->p()) {
      track_low = track_1;
      track_high = track_2;
   } else {
      track_low = track_2;
      track_high = track_1;
   }



   double high_tof_energy = 0;
   const std::vector<RecTofTrack* >& tofTracks = track_high->tofTrack();
   for (size_t i = 0; i < tofTracks.size(); ++i) {
      high_tof_energy += tofTracks[i]->energy();
   }

   if (!track_high->isEmcShowerValid()) {
      return false;
   }


   double high_emc_energy = track_high->emcShower()->energy();

   double high_total_energy = high_emc_energy + high_tof_energy;


   hist("h_high_emc_e")->Fill(high_emc_energy);
   hist("h_high_tof_e")->Fill(high_tof_energy);
   hist("h_high_total_e")->Fill(high_total_energy);

   //~ if (high_total_energy  < 1.5) { return false;}
   if (high_emc_energy  < 1.5) {
      return false;
   }

   double high_p_delta = TMath::Abs(track_high->mdcKalTrack()->p() - 3.77292/2);
   hist("h_high_p_delta")->Fill(high_p_delta);
   if (high_p_delta > 0.20) {
      return false;
   }

   //~ double high_ep_ratio = high_total_energy / track_high->mdcKalTrack()->p();
   double high_ep_ratio = high_emc_energy / track_high->mdcKalTrack()->p();
   hist("h_high_ep")->Fill(high_ep_ratio);


   if ( (high_ep_ratio > 1.05) || (high_ep_ratio < 0.85)) {
      return false;
   }

   if (track_low->mdcKalTrack()->p() > 1.3) {
      return false;
   }



   // pid

   init_pid(pid, track_low, run);

   pid->identify(
      pid->onlyPion()    | pid->onlyKaon() |
      pid->onlyElectron()
   );

   pid->calculate(run);

   double pid_prob_e = pid->probElectron() ;
   double pid_prob_sum = pid->probPion() + pid->probKaon() + pid->probElectron();






   double low_emc_energy = -1;
   if (track_low->isEmcShowerValid()) {
      low_emc_energy = track_low->emcShower()->energy();
   }


   //e_emc:p:cost:charge
   ntp_electrons->Fill(   low_emc_energy,
                          track_low->mdcKalTrack()->p(),
                          TMath::Cos(track_low->mdcKalTrack()->theta()),
                          track_low->mdcKalTrack()->charge(),
                          total_showers_energy,
                          delta_phi,
                          pid_prob_e,
                          pid_prob_sum,
                          track_low->mdcKalTrack()->phi(),
                          track_high->mdcKalTrack()->phi(),

                          track_low->mdcKalTrack()->theta(),
                          track_high->mdcKalTrack()->theta()


                      );


   return true;

}


//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void RadBhabhaEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " TestPIDEndJob() " << endl;
   }
}

#ifdef __cplusplus
}
#endif
