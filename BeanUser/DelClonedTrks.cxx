//////////////////////////////////////////////////////////////////////
//
// DelClonedTrks
//
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <TH1.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "MagneticField/MagneticFieldSvc.h"

#include "DstEvtRecTracks.h"
#include "ReadDst.h"

using namespace std;

// {{{1 Global variables
//--------------------------------------------------------------------
// static bool isMC = false;

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
void DelClonedTrksStartJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
   }

   // initialize DatabaseSvc (for VertexDbSvc) -----------------------
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if ( (dbs->GetDBFilePath()).empty() ) {
      // set path to directory with databases:
      dbs->SetDBFilePath(
            selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   // initialize Magnetic field (for VertexFitBField) ----------------
   MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   if ( (mf->GetPath()).empty() ) {
      // set path to directory with magnetic fields tables
      mf->SetPath(selector->AbsPath("Analysis/MagneticField"));
      mf->UseDBFlag(false); // like in the boss program
      mf->RunMode(3); // like in the boss program
   }

   // Book histograms ------------------------------------------------
   hst.resize(20,nullptr);

   hst[1] = new TH1D("ct_Ntrk","Ntrk initial", 21,-10.5,+10.5);
   hst[2] = new TH1D("ct_pp","cos(ang) pairs of +", 200,-1.,1.);
   hst[3] = new TH1D("ct_mm","cos(ang) pairs of -", 200,-1.,1.);
   hst[4] = new TH1D("ct_ang","cos(ang) + and -", 100,0.998,1.);
   hst[5] = new TH1D("ct_dpp","dP for + and -", 200,-0.02,0.02);
   hst[7] = new TH1D("ct_nrm","Ntrk to delete", 21,-10.5,+10.5);
   hst[8] = new TH1D("ct_ncl","N cloned", 10,0.5,+10.5);
   hst[9] = new TH1D("ct_status","status of tracks", 1001,-0.5,+1000.5);

   // register in selector to save in given directory
   const char* SaveDir = "DelClonedTrks";
   VecObj hsto(hst.begin(),hst.end());
   selector->RegInDir(hsto,SaveDir);
}

// {{{1 getVertexOrigin()
//--------------------------------------------------------------------
static Hep3Vector getVertexOrigin(int runNo, bool verbose = false) {
//--------------------------------------------------------------------
   static int save_runNo = 0;
   static Hep3Vector xorigin;

   if ( runNo == save_runNo ) {
      return xorigin;
   }

   // update vertex for new run
   xorigin.set(0.,0.,0.);
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();

#if BOSS_VER == 711
   string BossVer("7.1.1");     // just for test
#elif BOSS_VER > 700
   string BossVer("7.0.9");     // 2009,2012,2021 Psi(2S)
#else
   string BossVer("6.6.4");     // 2009 Psi(2S)
   int run = abs(runNo);
   if ( (run >= 25338 && run <= 27090)) {  // 2012 Psi(2S)
      BossVer = "6.6.4.p03";
   }
#endif
   vtxsvc->SetBossVer(BossVer);
   vtxsvc->handle(runNo);

   if ( vtxsvc->isVertexValid() ) {
      double* dbv = vtxsvc->PrimaryVertex();
      xorigin.set(dbv[0],dbv[1],dbv[2]);
      if( verbose ) {
         cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
      }
   } else {
      cerr << " FATAL ERROR:"
         " Cannot obtain vertex information for run#" << runNo
         << endl;
      Warning("Cannot obtain vertex information");
      exit(EXIT_FAILURE);
   }

   save_runNo = runNo;
   return xorigin;
}

// {{{1 MAIN: Event()
//--------------------------------------------------------------------
bool DelClonedTrksEvent( ReadDst* selector,
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
   // ++ Get event information
   //-----------------------------------------------------------------
   int runNo   = m_TEvtHeader->getRunId();
   // bool isMC = (runNo < 0);

   //-----------------------------------------------------------------
   // ++ Parameters for selecting a track and deleting cloned tracks
   //-----------------------------------------------------------------
   static const double Rvxy0_max = 1.0;
   static const double Rvz0_max = 10.0;
   // static const double cosTheta_max = 0.80;  // barrel only
   static const double cosTheta_max = 0.93;

   double maxCos = 0.9998; // ~1. degrees
   double maxDp  = 0.005;  // 5 MeV

   // ++ get interaction point
   Hep3Vector xorigin = getVertexOrigin(runNo);

   //-----------------------------------------------------------------
   // ++ loop over all charged tracks
   //-----------------------------------------------------------------
   const TEvtRecEvent* evtRecEvent =m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   vector<DstEvtRecTracks*> trk_p; // plus
   vector<DstEvtRecTracks*> trk_m; // minus
   trk_p.reserve(8);
   trk_m.reserve(8);
   for ( int i = 0; i < evtRecEvent->totalCharged(); i++ ) {
      DstEvtRecTracks* itTrk =
         static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
      if( !itTrk->isMdcTrackValid() ) {
         continue;
      }
      if( !itTrk->isMdcKalTrackValid() ) {
         continue;
      }

      RecMdcTrack* mdcTrk = itTrk->mdcTrack();
      hst[9]->Fill( double(mdcTrk->stat()) );

      double theta = mdcTrk->theta();
      double cosTheta = cos(theta);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // initial point for MDC rec.
      HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
      VFHelix helixip(point0,a,Ea);
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      // the nearest distance to IP
      double Rvxy0 = vecipa[0]; // in xy plane
      double Rvz0  = vecipa[3]; // in z direction

      if( fabs(Rvxy0) >= Rvxy0_max ) {
         continue;
      }
      if( fabs(Rvz0) >= Rvz0_max ) {
         continue;
      }
      if ( fabs(cosTheta) >= cosTheta_max ) {
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
         continue;
      }

      if( mdcKalTrk->charge() > 0 ) {
         trk_p.push_back(itTrk);
      } else {
         trk_m.push_back(itTrk);
      }
   }

   // Search for cloned tracks
   int Ntrkp = trk_p.size();
   int Ntrkm = trk_m.size();
   hst[1]->Fill(Ntrkp);
   hst[1]->Fill(-Ntrkm);
   set<DstEvtRecTracks*> trk_clone;

   for ( int i = 1; i < Ntrkp; ++i ) { // positive trks
      const auto& ti = trk_p[i];
      const auto& t1 = ti->mdcKalTrack();
      Hep3Vector Vp1(t1->px(), t1->py(), t1->pz());
      for ( int j = 0; j < i; ++j ) {
         const auto& tj = trk_p[j];
         const auto& t2 = tj->mdcKalTrack();
         Hep3Vector Vp2(t2->px(), t2->py(), t2->pz());
         double ang = Vp1.angle(Vp2);
         double ca = cos(ang);
         hst[2]->Fill(ca);
         hst[4]->Fill(ca);
         if ( ca > maxCos ) {
            double dp = Vp1.mag()-Vp2.mag();
            hst[5]->Fill(dp);
            if ( fabs(dp) < maxDp ) {
               auto tb = tj;
               if ( abs(t1->ndof() - t2->ndof()) < 4 ) {
                  if ( t1->chi2() > t2->chi2() ) {
                     tb = ti;
                  }
               } else {
                  if ( t1->ndof() < t2->ndof() ) {
                     tb = ti;
                  }
               }
               trk_clone.insert(tb);
               hst[7]->Fill(Ntrkp);
            }
         }
      }
   }
   // remove positive cloned tracks
   // auto check = [&trk_clone](DstEvtRecTracks* tr) {
      // return trk_clone.find(tr) != end(trk_clone);
   // };
   // trk_p.erase(remove_if(begin(trk_p),end(trk_p),check), end(trk_p));
   // trk_clone.clear();

   for ( int i = 1; i < Ntrkm; ++i ) { // negative trks
      const auto& ti = trk_m[i];
      const auto& t1 = ti->mdcKalTrack();
      Hep3Vector Vp1(t1->px(), t1->py(), t1->pz());
      for ( int j = 0; j < i; ++j ) {
         const auto& tj = trk_m[j];
         const auto& t2 = tj->mdcKalTrack();
         Hep3Vector Vp2(t2->px(), t2->py(), t2->pz());
         double ang = Vp1.angle(Vp2);
         double ca = cos(ang);
         hst[3]->Fill(ca);
         hst[4]->Fill(ca);
         if ( ca > maxCos ) {
            double dp = Vp1.mag()-Vp2.mag();
            hst[5]->Fill(dp);
            if ( fabs(dp) < maxDp ) {
               auto tb = tj;
               if ( abs(t1->ndof() - t2->ndof()) < 4 ) {
                  if ( t1->chi2() > t2->chi2() ) {
                     tb = ti;
                  }
               } else {
                  if ( t1->ndof() < t2->ndof() ) {
                     tb = ti;
                  }
               }
               trk_clone.insert(tb);
               hst[7]->Fill(-Ntrkm);
            }
         }
      }
   }
   // remove negative cloned tracks
   // trk_m.erase(remove_if(begin(trk_m),end(trk_m),check), end(trk_m));

   // mark cloned tracks
   if ( !trk_clone.empty() ) {
      hst[8]->Fill( double(trk_clone.size()) );
      for ( int i = 0; i < evtRecEvent->totalCharged(); i++ ) {
         DstEvtRecTracks* itTrk =
            static_cast<DstEvtRecTracks*>(evtRecTrkCol->At(i));
         if ( trk_clone.find(itTrk) != end(trk_clone) ) {
            RecMdcTrack* mdcTrk = itTrk->mdcTrack();
            mdcTrk->setStat(-222); // set label that track clone
            if ( selector->Verbose() ) {
               cout << " +++++ Cloned Track +++++" << endl;
               cout << " trackId= " << mdcTrk->trackId() << endl;
               cout << " FitQuality= " << mdcTrk->stat()
                  << " chi2= " << mdcTrk->chi2()
                  << " ndof = " << mdcTrk->ndof() << endl;
               cout << " charge= " << mdcTrk->charge()
                  << " Pxy= " << mdcTrk->pxy()
                  << " Pz= " << mdcTrk->pz()
                  << " P= " << mdcTrk->p() << endl;
               cout << " Theta= " << mdcTrk->theta()
                  << " Phi= " << mdcTrk->phi() << endl;
            }
         }
      }
   }

   return false;
}

// {{{1 EndJob
//--------------------------------------------------------------------
void DelClonedTrksEndJob(ReadDst* selector) {
//--------------------------------------------------------------------
   if ( selector->Verbose() ) {
      cout << " Start: " << __func__ << "()" << endl;
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
