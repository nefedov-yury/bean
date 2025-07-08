//////////////////////////////////////////////////////////////////////
//
// PsipHelixParCor
//
// For Monte Carlo generated events, helix parameters corrections are
// applied to all pions and Kaons. The program is based on the
// "TrackCorrection" class from the "Analysis" utilities.
//
// Since there is no way to check that corrections have already been
// applied or not, this function should be called once before all
// other functions.
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <string>

#include <TH1.h>

#include "CLHEP/Random/Randomize.h"

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "TrackCorrection/TrackCorrection.h"

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

// {{{1 Global variables
//--------------------------------------------------------------------
// TrackCorrection class for the helix parameters corrections for MC
static TrackCorrection* helix_cor = nullptr;

// expect one data taking period for all events
static int DataPeriod = 0;
static bool isMC = false;

// histograms
static vector<TH1*> hst;

// container for warnings
static map<string,int> warning_msg;

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

// {{{1 StartJob, book histograms
//--------------------------------------------------------------------
void PsipHelixParCorStartJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   // The initialization of helix_cor must be based on the DataPeriod,
   // which is determined by the run number

   // Book histograms ------------------------------------------------
   hst.resize(10,nullptr);

   hst[1] = new TH1D("mc_K_QP0","K-momentim Q*P(K) before corr.",
         300,-1.5,1.5);
   hst[2] = new TH1D("mc_K_QP1","K-momentim Q*P(K) after corr.",
         300,-1.5,1.5);
   hst[5] = new TH1D("mc_pi_QP0","Pi-momentim Q*P(pi) before corr.",
         300,-1.5,1.5);
   hst[6] = new TH1D("mc_pi_QP1","Pi-momentim Q*P(pi) after corr.",
         300,-1.5,1.5);

   // register in selector to save in given directory
   const char* SaveDir = "PsipHelixParCor";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
}

// {{{1 MAIN: Event()
//--------------------------------------------------------------------
bool PsipHelixParCorEvent( ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   //-----------------------------------------------------------------
   //-- Get event information --
   //-----------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   isMC = (runNo < 0);
   if ( !isMC ) {
      return false;
   }

   // define data taking period --------------------------------------
   // ATTENTION: this part is calculated only once: we get constants
   // that are the same in all events of one job
   if ( DataPeriod == 0 ) {
      int run = abs(runNo);
      if ( (run >= 8093 && run <= 9025)  ) {         // 2009 Psi(2S)
         DataPeriod = 2009;
      } else if ( (run >= 25338 && run <= 27090) ) { // 2012 Psi(2S)
         DataPeriod = 2012;
      } else if ( (run >= 66257 && run <= 69292) ) { // 2021 Psi(2S)
         DataPeriod = 2021;
      } else if ( (run >= 9613 && run <= 9779) ) { // 3650 2009-data
         DataPeriod = 3650;
      } else if ( (run >= 33725 && run <= 33772) ) { // 3650 2012-data
         DataPeriod = 3650;
      } else if ( (run >= 69612 && run <= 70132) ) { // 3650 2021-data
         DataPeriod = 3650;
      } else {
         cout << " FATAL: Data taking period undefined for runNo= "
              << runNo << endl;
         exit(EXIT_FAILURE);
      }

      // Setting the SEED for the CLHEP::HepRandom generator
      // used in TrackCorrection class

      // Print the default random engine for the static generator
      cout << " INFO:: The name of the default CLHEP random engine"
         " is '" << CLHEP::HepRandom::getTheEngine()->name() << "'"
         << endl;

      //  MixMaxRng can be seeded by supplying a single integer in the
      //  range [1, 2 147 483 647]. Streams created from seeds
      //  differing by at least one bit somewhere are guaranteed
      //  absolutely to be independent and non-colliding for at least
      //  the next 10^100 random numbers.
      long seed(run);
      CLHEP::HepRandom::getTheEngine()->setSeed(seed,0);

      // set helix parameters corrections (helix_cor)
      if ( DataPeriod >= 2009 && DataPeriod <= 2021 ) {
         string period = to_string(DataPeriod) +
#if (BOSS_VER < 700)
            "_v6.6.4";
#else
            "_v7.0.9";
#endif
         helix_cor = new TrackCorrection(period);
         const auto& w = helix_cor->Warning;
         if ( !w.empty() ){
            cout << " WARNING: " << w << endl;
            Warning(w);
         }
         helix_cor->prt_table();
      } else {
         string W = "no helix parameters corrections"
            " for DataPeriod= " + to_string(DataPeriod);
         cout << " WARNING: " << W << endl;
         Warning(W);
      }
   } // end of DataPeriod definition ---------------------------------

   // no helix parameters corrections: nothing to do
   if ( helix_cor == nullptr ) {
      return false;
   }

   // we correct all charged tracks in the event
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   for(int i = 0; i < evtRecEvent->totalCharged(); i++) {
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

      // require Kalman fit
      RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();
      if ( !mdcKalTrk ) {
         continue;
      }
      if ( std::isnan(mdcKalTrk->px())
            || std::isnan(mdcKalTrk->py())
            || std::isnan(mdcKalTrk->pz()) ) {
         Warning("Nan Kalman track");
         continue;
      }

      // Kaon
      mdcKalTrk->setPidType(RecMdcKalTrack::kaon);
      hst[1]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      helix_cor -> calibration( mdcKalTrk );
      hst[2]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );

      // pion must be the last one - it will remain the default PID
      mdcKalTrk->setPidType(RecMdcKalTrack::pion);
      hst[5]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
      helix_cor -> calibration( mdcKalTrk );
      hst[6]->Fill( mdcKalTrk->charge()*mdcKalTrk->p() );
   } // end of charged track cycle
   return false;
}

// {{{1 EndJob
//--------------------------------------------------------------------
void PsipHelixParCorEndJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }
   if ( !isMC ) {
      return;
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

#ifdef __cplusplus
}
#endif
