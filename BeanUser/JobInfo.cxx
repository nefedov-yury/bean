//////////////////////////////////////////////////////////////////////
//
// JobInfo
//
// Print information stored in JobInfoTree of dst.
// It does not need event loop.
// Use "bean.exe -N1 -u JobInfo -h /dev/null file_name.dst"
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"
#include "RootEventData/THltEvent.h"

#include "RootEventData/TJobInfo.h"

#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void JobInfoStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if ( selector->Verbose() ) cout << " JobInfoStartJob() " << endl;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool JobInfoEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if ( selector->Verbose() ) cout << " JobInfoEvent() " << endl;

   static TJobInfo*     m_TJobInfo = 0;
   static TBranch*      b_TJobInfo = 0;

   if ( m_TJobInfo ) return false; // print only once
   TFile* tf = selector->GetCurrentFile();
   if ( !tf ) {
      cout << " JobInf::ERROR: no open dst file" << endl;
      return false;
   }

   cout << " ====================== JOB INFO ====================== "
      << endl;

   tf->cd();
   TTree* tree = (TTree*)gDirectory->Get("JobInfoTree");
   if ( !tree ) {
      cout << " JobInf::ERROR: can not find JobInfoTree" << endl;
      return false;
   }
   tree->SetMakeClass(0);
   tree->SetBranchAddress("JobInfo", &m_TJobInfo, &b_TJobInfo);

   Long64_t ientry = tree->LoadTree(0);
   if ( ientry < 0 ) {
      cout << " JobInf::ERROR: can not load tree, ientry= "
         << ientry << endl;
      return false;
   }
   tree->GetEntry(0);

   if ( m_TJobInfo ) {
      cout << " BOSS version: " << m_TJobInfo->getBossVer() << endl;
      vector<string> JobOptions = m_TJobInfo->getJobOptions();
      cout << " Job Options: " << endl;
      for(unsigned int j = 0; j < JobOptions.size(); j++ ) {
         cout << JobOptions[j] << endl;
      }
      cout << " Decay options: " << m_TJobInfo->getDecayOptions()
         << endl;
      vector<int> TotEvents = m_TJobInfo->getTotEvtNo();
      cout << " Events number information:" << endl;
      int sum = 0;
      size_t nruns2 = TotEvents.size();
      for(size_t i = 0; i < nruns2-1; i+=2 ) {
         cout << " Run# " << TotEvents[i];
         int nevents =  TotEvents[i+1];
         cout << " N(events)= " << nevents << endl;
         sum += nevents;
      }
      cout << " ---------------------------- " << endl;
      if( sum > 0 )
         cout << " Total " << sum << " events " << endl;
   }
   cout << " ====================================================== "
      << endl;

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void JobInfoEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " JobInfoEndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
