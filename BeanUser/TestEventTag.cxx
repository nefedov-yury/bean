//////////////////////////////////////////////////////////////////////
//
// TestEventTag
//
// This example shows how to use EventTag library ported from BOSS
//
// NOTE: You MUST specify the paths to pdg and decayCodes files.
//       Use the function AbsPath(file) which returns absolute path
//       in the form: "path_to_base_Bean_dir/" + file
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

#include "EventTag/EventTagSvc.h"

#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static std::vector<TH1D*> his1;
//static TNtuple* fNtp;

static EventTagSvc * m_EventTagSvc;

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestEventTagStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " TestEventTagStartJob() " << endl;
   }

   // reserve
   his1.resize(100,(TH1D*)0);

   // book histograms
   his1[1] = new TH1D("evttagfirst",
         "EventTag second byte (first decay)", 257,-0.5,256.5);

     // register in selector to save in given directory
   VecObj his1o(his1.begin(),his1.end());
   selector->RegInDir(his1o," TestEventTag");

   //initialize EventTag
   m_EventTagSvc = EventTagSvc::instance();
   if ( !m_EventTagSvc->IsInitialized() ) {
      // set paths to pdg & decayCodes files:
      m_EventTagSvc->setPdtFile(
        selector->AbsPath("Analysis/EventTag/share/pdt.table") );
      m_EventTagSvc->setDecayTabsFile(
        selector->AbsPath("Analysis/EventTag/share/decay.codes") );

      // one MUST add this *.codes file to toplevel Makefile in order
      // to work with PROOF (only *.cxx/*.h files are included in
      // PROOF packages by defgault)
      // Or you can use network path here or copy this file to all the
      // PROOF worker nodes
      m_EventTagSvc->setUserDecayTabsFile(
        selector->AbsPath("BeanUser/mydecay.codes") );

      if( selector->Verbose() ) m_EventTagSvc->setVerbose(1);

      m_EventTagSvc->initialize();
   }

}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool TestEventTagEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestEventTagEvent() " << endl;

   m_EventTagSvc->setMcEvent(m_TMcEvent);
   if ( selector->Verbose() ) {
      cout  << "TAG:"  << hex << m_EventTagSvc->getEventTag() << endl;
   }
   int decCode = ((m_EventTagSvc->getEventTag())&0xFF00) >> 8;
   his1[1]->Fill( (double) decCode);

   return (0);
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestEventTagEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestEventTagEndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
