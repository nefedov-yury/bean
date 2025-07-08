//////////////////////////////////////////////////////////////////////
//
// TestAbsCor
//
// This is a test example of using functions of the
// Absorption Correction algorithm `boss/Analysis/PhotonCor/AbsCor/`
//
// The main differences from BOSS version:
// 1) You MUST point on directory with data-files.
//    In most cases it is:
//              selector->AbsPath("Analysis/AbsCor")
//    Here AbsPath(dir_name) function returns full path to this
//    directory: "path_to_base_Bean_dir/" + dir_name
//
// 2) The AbsCor algorithm in BEAN does not use the database, unlike
//    the BOSS version. However, you can read calibration files using
//    the functions:
//
//    void ReadDatac3p(std::string DataPathc3p, bool info=true);
//    void ReadParMcCor(std::string paraPath, bool info=true);
//    void ReadCorFunpara(std::string CorFunparaPath, bool info=true);
//
//    The path to the calibration files in CernVM `/cvmfs/` file
//    system can be obtained using the python3 program
//    `GetAbsCorFiles.py` located in `bean/scripts` directory.
//
// 3) I recommend copying calibration files to a local filesystem to
//    avoid dependency on the network file system (cvmfs).
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>

#include <TH1.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "AbsCor/AbsCor.h"

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

// Pointer to the class for Absorption Corrections.
// The 'static' limits the scope of the pointer to this file only.
// The object is created in the TestAbsCorStartJob() function
// and used in the TestAbsCorEvent() function.
static AbsCor* m_abscor = nullptr;

// histograms
static vector<TH1*> hst;

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestAbsCorStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   // init Absorption Corrections
   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

   // read files for pi0 and 'evsetTofCor' from 'local' files.
   // you can use either an absolute path to the file:
   // string data_c3p = "/home/nefedov/study/bes3/bean_67_new/"
   // "Analysis/AbsCor/dat/00-00-41/c3ptof2021psip.txt";
   // or a path relative to the BEAN folder.

   // Psi(2S) 2021: run66257-run69292
   string data_c3p = selector->AbsPath(
         "Analysis/AbsCor/dat/00-00-41/c3ptof2021psip.txt");
   m_abscor->ReadDatac3p( data_c3p );

   string cor_evsetTof = selector->AbsPath(
         "Analysis/AbsCor/dat/00-00-41/"
         "evsetTofCorFunctionPar2021psip.txt");
   m_abscor->ReadCorFunpara(cor_evsetTof);

   // Book histograms
   hst.resize(50,nullptr);
   hst[1] = new TH1D("G_ang","angle with closest chg.trk",
         180,0.,180.);
   hst[2] = new TH1D("G_n","N_{#gamma} in event", 11,-0.5,10.5);

   hst[3] = new TH1D("E_M2pi0","M^{2}(2#gamma)", 200,0.,0.04);
   hst[4] = new TH1D("E_M2all","M^{2}(2#gamma)", 500,0.,2.);
   hst[5] = new TH1D("E_Npi0","N(#pi^{0})", 10,-0.5,9.5);

   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,"TestAbsCor");
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool TestAbsCorEvent( ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   // Get event information
   int runNo   = m_TEvtHeader->getRunId();
   // int eventNo = m_TEvtHeader->getEventId();

   // check that it is Psi(2S) 2021: run66257-run69292
   int run = abs(runNo);
   if ( !(66257 <= run && run <= 69292) ) {
      cerr << " ERROR: this is not Psi(2S) 2021 data "
         "Absorption Corrections are probably wrong"
         << endl;
      exit(EXIT_FAILURE);
   }

   // Apply absorption corrections to all neutral tracks in the event
   m_abscor->AbsorptionCorrection(selector);

   // Select gammas candidates
   // parameters of reconstruction
   static const double min_angle = 10 * M_PI/180; // 10 grad

   // const TEvtRecObject* m_TEvtRecObject = selector->GetEvtRecObject();
   const TEvtRecEvent* evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   Hep3Vector xorig {0.,0.,0.};
   vector<RecEmcShower*> gtrk;      // gamma tracks
   vector<HepLorentzVector> Pg;     // 4-momentum of gammas

   for ( int i = evtRecEvent->totalCharged();
         i < evtRecEvent->totalTracks(); i++ ) {
      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if ( !itTrk->isEmcShowerValid() ) {
         continue;
      }
      RecEmcShower* emcTrk = itTrk->emcShower();

      // *) good EMC time:
      if ( emcTrk->time() < 0 || emcTrk->time() > 14 ) {
         continue;
      }

      // *) good EMC energy deposited in the barrel (endcap) part of EMC
      double eraw = emcTrk->energy();
      double absCosTheta = fabs(  cos(emcTrk->theta()) );

      bool GoodCluster=false;
      if ( absCosTheta < 0.8 ) {  //barrel
         GoodCluster = eraw > 25E-3;
      } else if ( absCosTheta > 0.85 && absCosTheta < 0.92 ) { //endcap
         GoodCluster = eraw > 50E-3;
      }
      if ( !GoodCluster ) {
         continue;
      }

      // *) the nearest charged track is far from cluster
      Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

      double tang = 200.; // min angle between cluster and track
      for ( int j = 0; j < evtRecEvent->totalCharged(); j++ ) {
         DstEvtRecTracks* jtTrk =
            static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(j));
         if ( !jtTrk->isExtTrackValid() ) {
            continue;
         }
         RecExtTrack* extTrk = jtTrk->extTrack();
         if ( extTrk->emcVolumeNumber() == -1 ) {
            continue;   //track does not hit EMC
         }
         Hep3Vector extpos = extTrk->emcPosition();
         double angd = fabs(extpos.angle(emcpos)); // [0,pi]
         if ( angd < tang ) {
            tang = angd;
         }
      } //--------------------------------------------End for(j)

      hst[1]->Fill( tang*180/M_PI );
      if ( tang < min_angle ) {
         continue;
      }

      gtrk.push_back(emcTrk);

      // (E,Px,Py,Pz) for good gammas
      Hep3Vector p3 = emcpos - xorig;
      p3 *= eraw / p3.mag();
      Pg.push_back( HepLorentzVector(p3,eraw) );
   } //-----------------------------------------------End for (i)

   size_t Ng = gtrk.size();
   hst[2]->Fill(Ng);

   if ( Ng < 2 ) {
      return false;
   }

   // analyse gg pairs: search for the decay pi0->2gamma
   size_t Npi0 = 0;
   for ( size_t i = 0; i < Ng-1; ++i ) {
      const auto& LVgi = Pg[i];
      for ( size_t j = i+1; j < Ng; ++j ) {
         const auto& LVgj = Pg[j];

         double Mgg2 = (LVgi+LVgj).m2();
         hst[3]->Fill(Mgg2);
         hst[4]->Fill(Mgg2);

         if ( 0.013 < Mgg2 && Mgg2 < 0.022 ) { // CUT pi0
            Npi0 += 1;
         }
      } // end of for(j)
   } // end of for(i)
   hst[5]->Fill(Npi0);

   return false;
}

#ifdef __cplusplus
}
#endif
