//////////////////////////////////////////////////////////////////////
//
// EntryList
//
// This simple example illustrates how to use TEntryList in BEAN
// -- set global variable 'action' to ENL_FILL to save list of
//    selected events in file 'FileName'
// -- set action=ENL_READ to read entry list from file 'FileName'
//    and run program EntryListEvent for selected events only
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <cstdlib>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TEntryList.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static std::vector<TH1*> his1;

//
enum ENL { ENL_FILL, ENL_READ, ENL_NULL };
static TEntryList* enlst = 0;
static const char* FileName = "EnLst.root";
static const char* EntryListName = "EnLst";
static enum ENL action = ENL_FILL;
// static enum ENL action = ENL_READ;

inline Double_t DtoR(Double_t ang){ return ang/180*M_PI; }

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void EntryListStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " EntryListStartJob() " << endl;

   // reserve
   his1.resize(10,nullptr);

   // book histograms:
   his1[1] = new TH1D("nmdctrk","Ntracks in MDC", 20,-0.5,19.5);
   his1[2] = new TH1D("p_pos","P positive (2trk)",100,0.,3.);
   his1[3] = new TH1D("p_neg","P negative (2trk)",100,0.,3.);

   his1[4] = new TH2D("phi_pos_neg","Phi_pos vs Phi_neg",
         100,-M_PI,M_PI,100,-M_PI,M_PI);

   // register in selector to save in given directory
   VecObj his1o(his1.begin(),his1.end());
   selector->RegInDir( his1o, "EntryList");

   // ===== Entry List =====
   switch( action ) {
      case ENL_FILL:
         enlst = new TEntryList(EntryListName,"Entry List example");
         break;
      case ENL_READ:
         {
            TFile file(FileName,"READ");
            if( file.IsZombie() ) {
               cout << " can not open " << file.GetName() << endl;
               exit(1);
            }
            // rean Entry List
            enlst = (TEntryList*) file.Get(EntryListName);
            if( !enlst ) {
               cout << " ERROR: Can not fine EnLst! " << endl;
               exit(1);
            }
            enlst->SetDirectory(0);
            // This function is called from SlaveBegin()
            // therefore the fChain is not accessible here.
            // create a trigger to setup TEntryList in Init()
            selector->SetEntryList(enlst);
            break;
         }
      default:
         break;
   }
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool EntryListEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   const TObjArray* m_mdcTrackCol = m_TDstEvent->getMdcTrackCol();
   int NmdcTracks = m_mdcTrackCol->GetEntries();

   his1[1]->Fill(NmdcTracks);

   if( NmdcTracks != 2 ) return false; // skip event

   int sum_charge = 0;
   vector<double> phi_trk(2,0);
   vector<double> p_trk(2,0);
   TIter mdcTrackIter(m_mdcTrackCol);
   TMdcTrack* mdcTrack = 0;
   while ((mdcTrack = (TMdcTrack*)mdcTrackIter.Next())) {
      sum_charge += mdcTrack->charge();

      if( mdcTrack->charge() > 0 ) {
         his1[2]->Fill(mdcTrack->p());
         p_trk[0] = mdcTrack->p();
         phi_trk[0] = mdcTrack->phi();
      } else {
         his1[3]->Fill(mdcTrack->p());
         p_trk[1] = mdcTrack->p();
         phi_trk[1] = mdcTrack->phi();
      }
   }

   // if( sum_charge==0 ) {
   if( sum_charge==0 && p_trk[0] < 1. && p_trk[1] < 1. ) {
      his1[4]->Fill(phi_trk[0],phi_trk[1]);
      // ===== Entry List =====
      if( action == ENL_FILL ) {
         selector->SaveEntryInList(enlst);
      }
   }

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void EntryListEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " EntryListEndJob() " << endl;

   // ===== Entry List =====
   if( action == ENL_FILL ) {
      TFile f(FileName,"RECREATE");
      enlst->Write();
      cout << enlst->GetN() << " entries were selected " << endl;
   }
}

#ifdef __cplusplus
}
#endif
