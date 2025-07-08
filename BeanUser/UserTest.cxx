//////////////////////////////////////////////////////////////////////
//
// UserTest
//
// This is example of user functions (here name is UserTest)
// IMPORTANT: It MUST contain "name"Event() function.
// Two functions "name"StartJob() and "name"EndJob() are optional.
//
// The return value of Event() function used only if
// option -o (define output ROOT tree file name) is defined:
//      true  -- save this event in output ROOT
//      false -- skip event
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>


#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"
#include "RootEventData/THltEvent.h"


#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static std::vector<TH1*> hst;
static TNtuple* fNtp;


//--------------------------------------------------------------------
BeanUserShared_EXPORT
void UserTestStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   // reserve
   hst.resize(100,(TH1D*)0);

   // book histograms
   hst[1] = new TH1D("nmdctrk","Ntracks in MDC", 20,-0.5,19.5);
   hst[2] = new TH1D("p_pos","P positive (2trk)",100,0.,3.);
   hst[3] = new TH1D("p_neg","P negative (2trk)",100,0.,3.);

   hst[4] = new TH2D("phi_pos_neg","Phi_pos vs Phi_neg",
         100,-M_PI,M_PI,100,-M_PI,M_PI);

   fNtp = new TNtuple("ntuple","Demo ntuple","charge:p:phi:theta");

   // register in selector to save in given directory
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,"UserTest");

   VecObj ntuples(1,fNtp);
   selector->RegInDir(ntuples,"UserTest");
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool UserTestEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   const TObjArray* m_mdcTrackCol = m_TDstEvent->getMdcTrackCol();
   int NmdcTracks = m_mdcTrackCol->GetEntries();

   hst[1]->Fill(NmdcTracks);

   if( NmdcTracks != 2 ) return false; // skip event

   int sum_charge = 0;
   vector<double> phi_trk(2,0);
   TIter mdcTrackIter(m_mdcTrackCol);
   TMdcTrack* mdcTrack = 0;
   while ((mdcTrack = (TMdcTrack*)mdcTrackIter.Next())) {
      sum_charge += mdcTrack->charge();

      if( mdcTrack->charge() > 0 ) {
         hst[2]->Fill(mdcTrack->p());
         phi_trk[0] = mdcTrack->phi();
      } else {
         hst[3]->Fill(mdcTrack->p());
         phi_trk[1] = mdcTrack->phi();
      }

      fNtp->Fill( mdcTrack->charge(), mdcTrack->p(), mdcTrack->phi(),
            mdcTrack->theta() );
   }

   if( sum_charge==0 ) hst[4]->Fill(phi_trk[0],phi_trk[1]);

   return (sum_charge==0); // save event if sum_charge is 0
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void UserTestEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }
}

#ifdef __cplusplus
}
#endif
