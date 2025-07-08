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


static THashList  histograms_cache;
static TNtuple* ntp_double;

static AbsCor* m_abscor = 0;

#define hist(name) ( (TH1 *) histograms_cache.FindObject(name) )

static inline double sqr(const double x ) {
   return x *x;
}
#define D_PLUS_MASS 1.86962
#define MASS_KAON 0.493677
#define MASS_PION 0.13957018

// http://stackoverflow.com/questions/1964150/c-test-if-2-sets-are-disjoint
template<class Set1, class Set2>
bool is_disjoint(const Set1& set1, const Set2& set2) {
   if(set1.empty() || set2.empty()) {
      return true;
   }

   typename Set1::const_iterator
   it1 = set1.begin(),
   it1End = set1.end();
   typename Set2::const_iterator
   it2 = set2.begin(),
   it2End = set2.end();

   if(*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) {
      return true;
   }

   while(it1 != it1End && it2 != it2End) {
      if(*it1 == *it2) {
         return false;
      }
      if(*it1 < *it2) {
         it1++;
      } else {
         it2++;
      }
   }

   return true;
}


static EventTagSvc* m_EventTagSvc;


Hep3Vector getShowerVector(RecEmcShower* shower);
bool isGoodShower(RecEmcShower* emcTrk, const TObjArray* evtRecTrkCol,
      const TEvtRecEvent*  m_evtRecEvent, int exclude_track = -1);

bool is_pion_dtag(Int_t index, TEvtRecDTag* evtRecDTag);
bool is_kaon_dtag(int index, TEvtRecDTag* evtRecDTag);
Hep3Vector getShowerVector(RecEmcShower* shower);

HepLorentzVector calc_miss_4momenum(TEvtRecDTag* evtRecDTag,
      HepLorentzVector observed4Momentum ) {
   //~ HepLorentzVector observed4Momentum =
   //~ pion4Momentum +
   //~ kaon4Momentum +
   //~ electron4Momentum;

   HepLorentzVector taggedD4Momentum(evtRecDTag->px(), evtRecDTag->py(),
         evtRecDTag->pz(), evtRecDTag->pe());
   cout << "taggedD4Momentum: " << taggedD4Momentum << endl;
   cout << "observed4Momentum: " << observed4Momentum << endl;

   double beamE = evtRecDTag->beamE();

   Hep3Vector cmsBoost(-11E-3, 0, 0);
   Hep3Vector initialMomentum = -cmsBoost * 2 * beamE;
   HepLorentzVector ecms(initialMomentum, 2 * beamE);

   HepLorentzVector taggedD4Momentum_cms =  taggedD4Momentum;
   HepLorentzVector observed4Momentum_cms = observed4Momentum;


   taggedD4Momentum_cms.boost(cmsBoost);
   observed4Momentum_cms.boost(cmsBoost);

   HepLorentzVector ecms_cms = ecms;
   ecms_cms.boost(cmsBoost);
   double beamE_cms = 0.5 * ecms_cms.e() ;

   HepLorentzVector expectedD4MomentumBC_cms;
   expectedD4MomentumBC_cms.setVectM(-taggedD4Momentum_cms.vect().unit()
         * TMath::Sqrt( sqr(beamE_cms) - sqr(D_PLUS_MASS)),
         D_PLUS_MASS );

   cout << "expectedD4MomentumBC_cms: " << expectedD4MomentumBC_cms <<
         endl;
   cout << "observed4Momentum_cms: " << observed4Momentum_cms << endl;



   HepLorentzVector missing4MomentumBC_cms = expectedD4MomentumBC_cms -
         observed4Momentum_cms;
   return missing4MomentumBC_cms;

}

#ifdef __cplusplus
extern "C" {
#endif





//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void DoubleDTagStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   VecObj his1, his_aux, his1_init, his_dtag;

   //~ his_aux.push_back( new TH1D("", "miss momentum", 200, -1, 1));
   //~ ntp_electrons = new TNtuple("ntp_electrons", "ntp_electrons", "e_emc:p:cost:charge:evt_tot_emc:evt_dphi:pid_prob_e:pid_prob_sum"  );
   //~ his1.push_back( ntp_electrons);
   //~ his1.push_back( new TH1D("h_ncharged", "h_ncharged", 10,-0.5,9.5  ));
   //~ his1.push_back( new TH1D("h_delta_phi", "h_delta_phi", 800,-1.9, 1.9  ));
   //~ his1.push_back( new TH1D("h_delta_phi2", "h_delta_phi2", 800,-TMath::Pi(), TMath::Pi()  ));


   //~ his1.push_back( new TH1D("h_total_emc", "h_total_emc", 400,0, 4  ));
   //~ his1.push_back( new TH1D("h_max_gamma_emc", "h_max_gamma_emc", 400,0, 2  ));
//~
   //~ his1.push_back( new TH1D("h_high_emc_e", "h_high_emc_e", 400,0, 4  ));
   //~ his1.push_back( new TH1D("h_high_tof_e", "h_high_tof_e", 400,0, 4  ));
   //~ his1.push_back( new TH1D("h_high_total_e", "h_high_total_e", 400,0, 4  ));
//~
   //~ his1.push_back( new TH1D("h_high_p_delta", "h_high_p_delta", 400,0, 2  ));
   //~ his1.push_back( new TH1D("h_high_ep", "h_high_ep", 240, 0, 1.2));


   //~ selector->RegInDir(his1,"RadBhabha");
   //~ histograms_cache.AddAll(selector->GetOutputList());


   // We have to initialize DatabaseSvc
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   //~ DatabaseSvc* dbs = DatabaseSvc::instance();
   //~ if( (dbs->GetDBFilePath()).empty() ) {
   //~ dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   //~ }

   // We have to initialize Magnetic field
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   //~ MagneticFieldSvc* mf = MagneticFieldSvc::instance();
   //~ if( !(mf->GetPath()).empty() ) {
   //~ cerr << " RhopiStartJob::WARNING:"
   //~ << " MagneticFieldSvc has already been initialized" << endl
   //~ << "                         path = " << mf->GetPath() << endl;
   //~ }
   //~ // set path to directory with magnetic fields tables
   //~ mf->SetPath( selector->AbsPath("Analysis/MagneticField") );
   //~ mf->UseDBFlag(false); // like in the boss program
   //~ mf->RunMode(3); // like in the boss program

   m_EventTagSvc = EventTagSvc::instance();
   m_EventTagSvc->setPdtFile(
         selector->AbsPath("Analysis/EventTag/share/pdt.table") );
   m_EventTagSvc->setDecayTabsFile(
         selector->AbsPath("Analysis/EventTag/share/decay.codes") );


   his1.push_back( new TH1D("h_ncharged", "h_ncharged", 10,-0.5,9.5  ));
   his1.push_back( new TH1D("h_200_200_dE",
               "delta energy  of selected  tagged Dp;#Delta E (GeV)", 750, -0.15,
               0.15));
   his1.push_back( new TH1D("h_200_200_mBC",
               "Beam constrained mass of selected  tagged Dp;m_{BC} (GeV)", 400,
               1.83, 1.911));
   his1.push_back( new TH1D("h_200_200_maxShowerE",
               "h_200_200_maxShowerE", 100, 0, 0.5));

   ntp_double = new TNtuple("ntp_double", "ntp_double",
         "mode_1:mode_2:dE_1:dE_2:mBC_1:mBC_2:max_shower_e:evt_tag:run:u_miss:p_miss");
   his1.push_back(ntp_double);


   selector->RegInDir(his1,"DoubleDTag");
   histograms_cache.AddAll(selector->GetOutputList());


   m_EventTagSvc -> initialize();

   m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

}

static vector<pair<TEvtRecDTag*, TEvtRecDTag*> >
selectDoubleDTag(const TObjArray* m_evtRecDTagCol,
      const TObjArray* evtRecTrkCol, int decayMode1, int decayMode2 ) {
   vector<pair<TEvtRecDTag*, TEvtRecDTag*> > result ;


   TIter evtRecDTag1Iter(m_evtRecDTagCol);

   while (TEvtRecDTag* curEvtRecDTag1 = (TEvtRecDTag*)
               evtRecDTag1Iter.Next()) {
      if (curEvtRecDTag1->type() != 1) {
         continue;   // require PID
      }

      if (curEvtRecDTag1->decayMode() != decayMode1)  {
         continue;
      }
      // found first

      TIter evtRecDTag2Iter = evtRecDTag1Iter;
      while (TEvtRecDTag* curEvtRecDTag2 = (TEvtRecDTag*)
                  evtRecDTag2Iter.Next()) {
         if (curEvtRecDTag2->type() != 1) {
            continue;   // require PID
         }
         if (curEvtRecDTag2->charm() ==  curEvtRecDTag1->charm()) {
            continue;   // skip same charm
         }

         if (curEvtRecDTag2->decayMode() != decayMode2)  {
            continue;   //skip wrong mode
         }
         // compare tracks

         if (std::set<Int_t>(curEvtRecDTag1->tracks().begin(),
                     curEvtRecDTag1->tracks().end()) !=
               std::set<Int_t>(curEvtRecDTag2->otherTracks().begin(),
                     curEvtRecDTag2->otherTracks().end())) {
            continue;
         }


         // no common photons
         if (!is_disjoint(std::set<Int_t>(curEvtRecDTag1->showers().begin(),
                           curEvtRecDTag1->showers().end()),
                     std::set<Int_t>(curEvtRecDTag2->showers().begin(),
                           curEvtRecDTag2->showers().end()))) {
            cout << "common photon!" << endl;
            continue;
         }


         cout << "Found double dtag:" << endl;
         cout << "  dtag1: mode=" << curEvtRecDTag1->decayMode() << " charm="
               << curEvtRecDTag1->charm() << endl;
         cout << "         mBC=" << curEvtRecDTag1->mBC() << " deltaE=" <<
               curEvtRecDTag1->deltaE() << endl;

         cout << "  dtag2: mode=" << curEvtRecDTag2->decayMode() << " charm="
               << curEvtRecDTag2->charm() << endl;
         cout << "         mBC=" << curEvtRecDTag2->mBC() << " deltaE=" <<
               curEvtRecDTag2->deltaE() << endl;


         result.push_back( make_pair(curEvtRecDTag1, curEvtRecDTag2));





      }
   }

   return result;
};


static void formatTag(char* tag_str, unsigned int eventTag_,
      TMcEvent* m_TMcEvent) {
   unsigned int eventTag = eventTag_;

   unsigned char byte1 =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   unsigned char byte2 =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   unsigned char byte3 =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   if (byte2 > byte3) {
      sprintf(tag_str, "%d-%d-%d", byte1, byte2, byte3);
   } else {
      sprintf(tag_str, "%d-%d-%d", byte1, byte3, byte2);
   }


}


//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
bool DoubleDTagEvent( ReadDst* selector,
      TEvtHeader* m_TEvtHeader,
      TDstEvent* m_TDstEvent,
      TEvtRecObject* m_TEvtRecObject,
      TMcEvent* m_TMcEvent,
      TTrigEvent* m_TTrigEvent,
      TDigiEvent* m_TDigiEvent,
      THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{

   m_abscor->SuppressHotCrystals(selector);

   int event = m_TEvtHeader->getEventId();
   int run = m_TEvtHeader->getRunId();


   char tag_str[12];
   m_EventTagSvc->setMcEvent(m_TMcEvent);
   unsigned int eventTagDecay =  (m_EventTagSvc->getEventTag()  &
               0xFFFFFF00) >> 8;
   formatTag(tag_str, eventTagDecay, m_TMcEvent);
   //~ cout << "tag_str: " << tag_str << endl;



   // evtRecTrkCol
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();
   const TEvtRecEvent*  m_evtRecEvent =
         m_TEvtRecObject->getEvtRecEvent();

   const TObjArray*  m_evtRecDTagCol=
         m_TEvtRecObject->getEvtRecDTagCol();



   auto double_tags_200_201 = selectDoubleDTag(m_evtRecDTagCol,
               evtRecTrkCol, 200, 201);
   auto double_tags_200_200 = selectDoubleDTag(m_evtRecDTagCol,
               evtRecTrkCol, 200, 200);

   auto double_tags = double_tags_200_200;
   double_tags.insert( double_tags.end(), double_tags_200_201.begin(),
         double_tags_200_201.end() );





   for (size_t i = 0; i < double_tags.size(); ++i) {
      auto it = double_tags[i];



      double maxShowerE = 0;
      for (unsigned shower_id = m_evtRecEvent->totalCharged();
            shower_id < m_evtRecEvent->totalTracks(); ++shower_id) {
         //check this shower doesn't belong to any of D
         if (find(it.first->showers().begin(), it.first->showers().end(),
                     shower_id) != it.first->showers().end()) {
            continue;
         }
         if (find(it.second->showers().begin(), it.second->showers().end(),
                     shower_id) != it.second->showers().end()) {
            continue;
         }



         DstEvtRecTracks* track = (DstEvtRecTracks*) evtRecTrkCol->At(
                     shower_id);
         if (track->isEmcShowerValid()) {
            RecEmcShower* emcTrk = track->emcShower();
            double energy = emcTrk->energy();

            if (isGoodShower(emcTrk, evtRecTrkCol, m_evtRecEvent)) {
               double energy = emcTrk->energy();
               cout << "good shower id:" << shower_id << endl;

               if (energy > maxShowerE) {
                  maxShowerE = energy;
               }
            }
         }
      }


      if ((it.first->decayMode()==200) && (it.second->decayMode()==200)) {
         hist("h_200_200_dE")->Fill(it.first->deltaE());
         hist("h_200_200_dE")->Fill(it.second->deltaE());
         hist("h_200_200_mBC")->Fill(it.first->mBC());
         hist("h_200_200_mBC")->Fill(it.second->mBC());
         hist("h_200_200_maxShowerE")->Fill(maxShowerE);
      }


      double u_miss = -100;
      double p_miss = -100;

      if ((it.first->decayMode()==200) && (it.second->decayMode()==201)) {
         TEvtRecDTag* evtRecDTag = it.first;
         TEvtRecDTag* signalD = it.second;
         HepLorentzVector observed4Momentum(0,0,0,0);

         for (unsigned charged_id = 0; charged_id < signalD->tracks().size() ;
               ++charged_id) {
            double charged_mass = -1;

            if (is_pion_dtag( signalD->tracks()[charged_id], signalD)) {
               charged_mass = MASS_PION;
            } else if (is_kaon_dtag( signalD->tracks()[charged_id], signalD)) {
               charged_mass = MASS_KAON;
            } else {
               cerr << "unknown type of charged track! track_id: " <<
                     signalD->tracks()[charged_id] << endl;
            }

            DstEvtRecTracks* charged_track = (DstEvtRecTracks*) evtRecTrkCol->At(
                        signalD->tracks()[charged_id]);
            HepLorentzVector charged4Momentum = charged_track->mdcKalTrack()->p4(
                        charged_mass);
            observed4Momentum += charged4Momentum;
         }

         // add single photon
         //~ for (int photon_num=0; photon_num < 1; ++photon_num)
         int photon_num = event % 2; // kind of deterministic random
         {
            DstEvtRecTracks* photon_track = (DstEvtRecTracks*) evtRecTrkCol->At(
                        signalD->showers()[photon_num]);

            Hep3Vector photonShowerVector = getShowerVector(
                        photon_track->emcShower());
            Hep3Vector photonMomentum = photonShowerVector.unit();
            photonMomentum.setMag(photon_track->emcShower()->energy());

            HepLorentzVector photon4Momentum;
            photon4Momentum.setVectM(photonMomentum, 0);

            observed4Momentum +=  photon4Momentum;
         }

         //~ HepLorentzVector taggedD4Momentum(evtRecDTag->px(), evtRecDTag->py(), evtRecDTag->pz(), evtRecDTag->pe());
         //~ u_miss = calc_umiss_bc(evtRecDTag, taggedD4Momentum);

         HepLorentzVector missing4MomentumBC_cms = calc_miss_4momenum(
                     evtRecDTag, observed4Momentum);
         u_miss = missing4MomentumBC_cms.e() -
               missing4MomentumBC_cms.vect().mag();
         p_miss = missing4MomentumBC_cms.vect().mag();




      }




      //mode_1:mode_2:dE_1:dE_2:mBC_1:mBC_2:max_shower_e
      ntp_double->Fill( it.first->decayMode(), it.second->decayMode(),
            it.first->deltaE(), it.second->deltaE(),
            it.first->mBC(), it.second->mBC(),
            maxShowerE, eventTagDecay, run, u_miss, p_miss);


      if (maxShowerE > 0.55) {
         cout << "huge shower "  << maxShowerE << endl;
      }
   }






   if ( !double_tags.empty()) {
      return true;
   }


   return false;

}


//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void DoubleDTagEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " DoubleDTagEndJob() " << endl;
   }
   delete m_abscor;
}

#ifdef __cplusplus
}
#endif
