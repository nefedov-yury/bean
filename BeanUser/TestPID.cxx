//////////////////////////////////////////////////////////////////////
//
// TestPID
//
// This is a test example of using the ParticleID library functions
//
// The main differences from BOSS version:
// 1) You MUST point on directory with data-files.
//    In most cases it is:
//      pid->set_path( selector->AbsPath("Analysis/ParticleID") );
//    (AbsPath(dirname) returns "path_to_base_Bean_dir/" + dirname)
// 2) You MUST call pid->calculate(run_number) function because
//    there are no ways to know run_number inside of ParticleID
//
//  Example with pseudocode:
//  +++++++++++++++++++++++++++++++++++++++++
//  ParticleID *pid = ParticleID::instance();
//#if (BOSS_VER < 700)
//  pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
//#else
//  pid->set_path(selector->AbsPath("Analysis/ParticleID"));
//#endif
//  for (int i = 0; i < ncharg; i++) {
//      pid->init();
//      pid->setMethod(pid->methodProbability());
//      pid->setChiMinCut(4); // default=5.
//      pid->setRecTrack(dst_tracks[i]);
//      pid->usePidSys(pid->useDedx() | pid->useTof1() |
//                      pid->useTof2() | ...);
//      pid->identify(pid->onlyPion() | pid->onlyKaon() | ...);
//      pid->calculate(run_number); // run>0 for data; othervise MC
//      if(!(pid->IsPidInfoValid())) continue;
//      // user's selections
//  }
//  if(pid) delete pid;
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <cstdlib>

#include <TH1.h>
#include <TH2.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "DstEvtRecTracks.h"
#include "ParticleID/ParticleID.h"

#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static FILE* fout = 0;

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestPIDStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestPIDStartJob() " << endl;

   const char outfile[] = "TestPID.out";
   fout = fopen(outfile, "w+");
   if( !fout ) {
      cerr << " can not open: " << outfile << endl;
      exit(1);
   }
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool TestPIDEvent( ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestPIDEvent() " << endl;

   const TEvtRecEvent*
      m_evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   int event = m_TEvtHeader->getEventId();
   int run = m_TEvtHeader->getRunId();
   int ncharge = m_evtRecEvent->totalCharged();

   string starline(50,'*');
   fprintf(fout, "%s\n", starline.c_str() );
   fprintf(fout,"event: %i, run: %i, total charged: %i\n",
         event, run, ncharge);
   fprintf(fout, "%s\n", starline.c_str() );

   cout << "event: " << event <<", run: " << run
      << ", total charged: " << ncharge << endl;

   if( ncharge == 0 ) return false;

   // evtRecTrkCol
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   // initialization
   ParticleID *pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif
   // cout << " pid= " << pid << endl;

   //-------------------------------------------------------
   //--------------------standard-PID-----------------------
   //-------------------------------------------------------

   for(int i = 0; i < ncharge; i++){
      DstEvtRecTracks*
         dst_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(i);
      if( dst_tracks->isMdcTrackValid() ) {
         // TMdcTrack *mdcTrk = dst_tracks->mdcTrack();

         cout << " i= " << i << endl;

         pid->init();
         pid->setMethod(pid->methodProbability());
         pid->setChiMinCut(4);

         pid->setRecTrack(dst_tracks);
         pid->usePidSys( pid->useTof1() | pid->useTof2() |
               pid->useDedx() );
         // pid->usePidSys(pid->useDedx());
         // pid->usePidSys(pid->useTof1() | pid->useTof2());
         pid->identify(pid->onlyPion() | pid->onlyKaon() |
               pid->onlyProton() | pid->onlyMuon() |
               pid->onlyElectron());
         cout << " pid->calculate(" << run << ")" << endl;
         pid->calculate(run);

         if(pid->IsPidInfoValid()) {
            fprintf(fout,
                  "track #%i "
                  "chi2(electron, muon, pion, kaon, proton):",i);
            fprintf(fout,
                  "%15.10f %15.10f %15.10f %15.10f %15.10f\n",
                  pid->chi(0),pid->chi(1),pid->chi(2),pid->chi(3)
                  ,pid->chi(4));

            fprintf(fout,
                  "track #%i "
                  "prob(electron, muon, pion, kaon, proton):",i);
            fprintf(fout,
                  "%15.10f %15.10f %15.10f %15.10f %15.10f\n",
                  pid->probElectron(), pid->probMuon(),
                  pid->probPion(), pid->probKaon(),
                  pid->probProton());
         } else {
            fprintf(fout,"NOPID INFO\n");
         }
      } else {
         fprintf(fout,"BAD TRACK\n");
      }
   }

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestPIDEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestPIDEndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
