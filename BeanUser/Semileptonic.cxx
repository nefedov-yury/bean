#include "DLLDefines.h"         // mandatory!

#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TNtuple.h>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TMap.h>
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

#include "ReadDst.h"
#include <TObjString.h>
#include <assert.h>
// PDG2010
#define MASS_KAON 0.493677
#define MASS_PION 0.13957018
#define MASS_ELECTRON 0.510998910E-3
#define D_PLUS_MASS 1.86962
#define PSI3770_MASS  3.77292
#define MASS_MUON 0.10565836668
using namespace std;


static EventTagSvc* m_EventTagSvc;

//~ typedef vector<TH1*> VecHist1;
static VecObj hists_dtag_decay_good_mBC;
static VecObj hists_dtag_decay_good_dE;
static VecObj hists_dtag_decay_good_Ntags;

static VecObj hists_dtag_decay_bad_mBC;
static VecObj hists_dtag_decay_bad_dE;

static VecObj  hists_dtag_decay_selected_bad_mBC;
static VecObj  hists_dtag_decay_selected_good_mBC;
static VecObj  hists_dtag_decay_selected_mBC;

static VecObj  hists_dtag_decay_good_mBC_vs_dE;
static VecObj  hists_dtag_decay_bad_mBC_vs_dE;
static TNtuple* ntp_features;

static TNtuple* ntp_angular_selected;
static TNtuple*  ntp_angular_selected_wmc;
static TNtuple* ntp_angular_mc;
static TNtuple* ntp_momenta_mc;
static TNtuple* ntp_electrons;

static THashList  histograms_cache;
static TMap* cuts_boundaries;

//~ #define hist(name) ( (TH1 *) selector->GetOutputList()->FindObject(name) )
#define hist(name) ( (TH1 *) histograms_cache.FindObject(name) )

static void smartFillHist(std::string name, double value,
      bool is_good, bool is_excbad = false, bool has_mc = true) {
   hist(name.c_str())->Fill(value);
   if (has_mc) {
      if (is_good) {
         hist((name + "_good").c_str())->Fill(value);
      } else {
         hist((name + "_bad").c_str())->Fill(value);
         if (is_excbad) {
            hist((name + "_exbad").c_str())->Fill(value);
         }
      }
   }
}


static void smartBookHistImpl(VecObj& his1, TString title,
      TString description, int nBins, double low, double high ) {
   //~ cout << title << endl;
   his1.push_back( new TH1D( title, TString::Format(description, "all"),
               nBins, low, high));
   his1.push_back( new TH1D( title + "_good",
               TString::Format(description, "signal"), nBins, low, high));
   his1.push_back( new TH1D( title + "_bad", TString::Format(description,
                     "background"), nBins, low, high));
   his1.push_back( new TH1D( title + "_exbad",
               TString::Format(description, "exclusive background"), nBins, low,
               high));
}

static void smartBookHist(VecObj& his1, VecObj& his_dtag,
      TString title, TString description, int nBins, double low,
      double high ) {
   for (int decay = 200; decay < 206; ++decay) {
      smartBookHistImpl(his_dtag, TString::Format("%s_%d",title.Data(),
                  decay), description, nBins, low, high);
      //~ smartBookHistImpl(his_dtag, TString::Format("%s_-%d",title.Data(), decay), description, nBins, low, high);
   }

   smartBookHistImpl(his_dtag, title + "_+", description, nBins, low,
         high);
   smartBookHistImpl(his_dtag, title + "_-", description, nBins, low,
         high);

   smartBookHistImpl(his1, title, description, nBins, low, high);

}

static void smartBookHist(VecObj& his_init, VecObj& his,
      VecObj& his_dtag, TString title, TString description, int nBins,
      double low, double high ) {
   smartBookHist(his_init, his_init, "i_" + title,
         description + "  [init.]", nBins, low, high);
   smartBookHist(his_init, his_init, "s_" + title,
         description + "  [seq.]", nBins, low, high);
   smartBookHist(his, his_dtag, title, description, nBins, low, high);
}




class TCutsManager {
   private:
      struct TSmartHistValue {
         double value;
         bool is_excbad;
      };
      typedef std::map<std::string, std::vector<TSmartHistValue> >
      THistValues;
      typedef std::vector<std::string> TCutsList;

      TCutsList cuts_applied;
      TCutsList cuts;
      THistValues hist_values;

      std::set<std::string> cuts_set;

      bool is_good;
      bool has_mc;
      TEvtRecDTag* evtRecDTag;
      int charm;

   public:
      TCutsManager(bool mc, bool good, TEvtRecDTag* evtRecDTag_ = 0,
            int charm_ = 0)
         : is_good(good)
         , has_mc(mc)
         , evtRecDTag(evtRecDTag_)
         , charm(charm_)
      {};

      void initCut(std::string name) {
         //save cuts order
         if (cuts_set.count(name) == 0) {
            cuts.push_back(name);
            cuts_set.insert(name);
         }
      }
      void smartFillCutHistImpl(std::string name, double value,
            bool is_excbad = false) {
         initCut(name);

         TSmartHistValue smart_value;
         smart_value.value = value;
         smart_value.is_excbad = is_excbad;

         hist_values[name].push_back(smart_value);
      }

      void smartFillCutHist(std::string name, double value,
            bool is_excbad = false) {
         if (evtRecDTag) {
            smartFillCutHistImpl(TString::Format("%s_%d",name.c_str(),
                        evtRecDTag->decayMode()).Data(), value, is_excbad);
         }

         if (charm != 0) {
            std::string sign_str;
            if (charm == 1) {
               sign_str = "+";
            } else {
               sign_str = "-";
            }

            smartFillCutHistImpl((name + "_" + sign_str).data(), value,
                  is_excbad);
         }


         smartFillCutHistImpl(name, value, is_excbad);
      }


      void rejectByCut(std::string cut_name) {
         cuts_applied.push_back(cut_name);
      }
      void addBoundaries(std::string cut_name, TString boundaries) {
         if (!cuts_boundaries->Contains(cut_name.c_str())) {
            cuts_boundaries->Add(new TObjString(cut_name.c_str()),
                  new TObjString(boundaries));
         }

      }


      bool rejectByCutIfLower(std::string cut_name, double value,
            double left) {
         if (value < left) {
            rejectByCut(cut_name);
            addBoundaries(cut_name, TString::Format("%f", left));
            return true;
         } else {
            return false;
         }
      }

      bool rejectByCutIfHigher(std::string cut_name, double value,
            double right) {
         if (value > right) {
            rejectByCut(cut_name);
            addBoundaries(cut_name, TString::Format("%f", right));
            return true;
         } else {
            return false;
         }
      }

      bool rejectByCutIfOutsideWindow(std::string cut_name, double value,
            double left, double right) {

         if ( (value < left) || (value > right)) {
            rejectByCut(cut_name);

            addBoundaries(cut_name, TString::Format("%f,%f", left, right));


            return true;
         } else {
            return false;
         }

      }

      bool rejectByCutIfInsideWindow(std::string cut_name, double value,
            double left, double right) {
         if ( (value > left) && (value < right)) {
            rejectByCut(cut_name);

            addBoundaries(cut_name, TString::Format("%f,%f", left, right));


            return true;
         } else {
            return false;
         }

      }




      void smartFillHistValues(const std::string& name,
            const std::vector<TSmartHistValue>& values) {
         for (std::vector<TSmartHistValue>::const_iterator value =
                     values.begin();
               value != values.end(); ++value) {
            //~ if (evtRecDTag) {
            //~ smartFillHist(TString::Format("%s_%d",name.c_str(), evtRecDTag->decayMode()).Data(), value->value, is_good, value->is_excbad, has_mc);
            //~ }

            smartFillHist(name.c_str(), value->value, is_good, value->is_excbad,
                  has_mc);
         }
      }
      void fillCutsHisto(const char* name, bool is_good = false ) {
         if (is_good) {
            hist("hCuts_good")->Fill(name, 1.0);
         } else {
            hist("hCuts_bad")->Fill(name,1.0);
         }
      }

      void fillCutsStatistics() {
         fillCutsHisto("all", is_good);
         for (TCutsList::const_iterator cut = cuts.begin();
               cut != cuts.end() && (cuts_applied.empty()
                     || *cut != cuts_applied[0] );
               ++cut) {
            fillCutsHisto(cut->c_str(), is_good);

         }

         if (survived()) {
            fillCutsHisto("survived", is_good);
         }

      }

      void fatalByCut(std::string cut_name) {
         rejectByCut(cut_name);
         processCuts(true);



      }

      void fillCutsSequentalHistos() {
         for (TCutsList::const_iterator cut = cuts.begin();
               cut != cuts.end();
               ++cut) {

            THistValues::const_iterator iter =  hist_values.find(*cut);
            if (iter != hist_values.end()) {
               smartFillHistValues("s_" + *cut, iter->second);
            }

            if  (!cuts_applied.empty() && cuts_applied[0] == *cut ) {
               break;
            }
         }

      }

      bool processCuts(bool fatal = false) {
         for (THistValues::const_iterator iter = hist_values.begin();
               iter != hist_values.end(); ++iter) {
            const std::string& cut_name = iter->first;

            if (!fatal) {
               if (cuts_applied.empty() || ( cuts_applied.size() == 1
                           && cuts_applied[0] == cut_name )) {
                  //fill cut hist
                  smartFillHistValues(cut_name, iter->second);
               }
            }
            smartFillHistValues("i_" + cut_name, iter->second);
         }

         fillCutsSequentalHistos();
         fillCutsStatistics();

         return survived();
      }
      bool survived() {
         return cuts_applied.empty();
      }

};





static inline int get_charge(const TMcParticle* particle) {
   return (particle->getParticleID() < 0) ? -1 : 1;
}

static bool hasDecay(unsigned int eventTag, unsigned char top,
      unsigned char decay ) {
   unsigned char topDecay =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   unsigned char leftDecay =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   unsigned char rightDecay =   (eventTag & 0xFF); //last byte
   eventTag >>= 8;

   if ((top == topDecay) && ( ( leftDecay == decay )
               ||( rightDecay == decay) )) {
      return true;
   } else {
      return false;
   }

}

static bool selectDtagDecay(TEvtRecDTag* evtRecDTag, int decayMode,
      double deltaELow, double deltaEHigh, double mBCLow, double mBCHigh) {
   if ( (evtRecDTag->decayMode() == decayMode) ) {
      if ((evtRecDTag->deltaE() > deltaELow)
            && (evtRecDTag->deltaE() < deltaEHigh)) {
         if  ((evtRecDTag->mBC() > mBCLow) && (evtRecDTag->mBC() < mBCHigh)) {
            return true;
         }
      }
   }
   return false;
}

static bool fillDtagDecayHists(TEvtRecDTag* evtRecDTag, int decayMode,
      unsigned int eventTagDecay, int topEventTag, int decayEventTag,
      double deltaELow, double deltaEHigh, double mBCLow, double mBCHigh) {
   bool selection_passed = selectDtagDecay(evtRecDTag, decayMode,
               deltaELow, deltaEHigh,  mBCLow, mBCHigh);

   if ( (evtRecDTag->decayMode() == decayMode) ) {
      if (hasDecay(eventTagDecay, topEventTag, decayEventTag )) {
         ((TH1*) hists_dtag_decay_good_mBC[decayMode])-> Fill(
               evtRecDTag->mBC() );
         ((TH1*) hists_dtag_decay_good_dE[decayMode]) -> Fill(
               evtRecDTag->deltaE() );

         ((TH2*) hists_dtag_decay_good_mBC_vs_dE[decayMode]) -> Fill(
               evtRecDTag->mBC(), evtRecDTag->deltaE());


         if (selection_passed) {
            ((TH1*) hists_dtag_decay_selected_good_mBC[decayMode]) -> Fill(
                  evtRecDTag->mBC() );
         }
      } else {
         ((TH1*) hists_dtag_decay_bad_mBC[decayMode])-> Fill(
               evtRecDTag->mBC() );
         ((TH1*) hists_dtag_decay_bad_dE[decayMode]) -> Fill(
               evtRecDTag->deltaE() );

         ((TH2*) hists_dtag_decay_bad_mBC_vs_dE[decayMode]) -> Fill(
               evtRecDTag->mBC(), evtRecDTag->deltaE());

         if (selection_passed) {
            ((TH1*) hists_dtag_decay_selected_bad_mBC[decayMode]) -> Fill(
                  evtRecDTag->mBC() );
         }

      }

   }
   return selection_passed;



}
//~ static void initialize_pid(ParticleID* pid, DstEvtRecTracks* dst_track, )
static void init_pid(ParticleID* pid, DstEvtRecTracks* dst_tracks,
      int run) {
   pid->init();
   pid->setMethod(pid->methodProbability());
   pid->setChiMinCut(4); //FIXME: 2013-12-06 which one should we use?

   pid->setRecTrack(dst_tracks);
   //~ pid->usePidSys(pid->useDedx() | pid->useTof()| pid->useTof1() | pid->useTof2() | pid->useEmc() | pid->useMuc());
   //~ pid->usePidSys(0xff);

   pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() |
         pid->useTofE() | pid->useEmc());

}
static bool is_electron(int index, ParticleID* pid,
      const  TObjArray* evtRecTrkCol, int run) {
   DstEvtRecTracks* dst_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
               index);

   init_pid(pid, dst_tracks, run);

   pid->identify(
         pid->onlyPion()    | pid->onlyKaon() |
         pid->onlyElectron()
   );

   pid->calculate(run);


   if ( ( pid->probElectron() > 0.8 * (  pid->probPion() +
                     pid->probKaon() + pid->probElectron())) &&
         ( pid->probElectron() > 0.001)
   ) {
      return true;
   }

   return false;
}

bool is_pion_dtag(Int_t index, TEvtRecDTag* evtRecDTag) {
   const vector<Int_t>& piid = evtRecDTag->pionId();
   return  (find (piid.begin(), piid.end(), index ) != piid.end());
}

bool is_kaon_dtag(int index, TEvtRecDTag* evtRecDTag) {
   const vector<Int_t>& kid = evtRecDTag->kaonId();
   return  (find (kid.begin(), kid.end(), index ) != kid.end());
}

static unsigned int get_num_kaon_pion(const vector<Int_t>& tracks,
      const vector<Int_t>&   kp_tracks) {
   unsigned int counter = 0;

   for (size_t i = 0; i < kp_tracks.size(); ++i)
      for (size_t j = 0; j < tracks.size(); ++j)  {
         if (kp_tracks[i] == tracks[j]) {
            counter += 1;
            break;
         }
      }
   return counter;

}

static bool hasPi0 (TMcEvent* m_TMcEvent) {
   const TObjArray* particles = m_TMcEvent->getMcParticleCol();

   TIter next(particles);
   while (TMcParticle* particle = (TMcParticle*) next()) {
      if (particle->getParticleID() == 111) {
         return true;
      }
   }
   return false;

}


static double getZ0Energy (TMcEvent* m_TMcEvent) {
   const TObjArray* particles = m_TMcEvent->getMcParticleCol();

   TObjArrayIter next(particles);
   while (TMcParticle* particle = (TMcParticle*) next()) {
      if (particle->getParticleID() == 23) {
         return particle->getInitialMomentumE();
      }
   }
   return -1;

};



static inline double sqr(const double x ) {
   return x *x;
}

static  bool isGoodPoca(BMdcKalTrack* mdcKalTrk,
      Hep3Vector& xorigin) {
   HepVector a = mdcKalTrk->getZHelix();
   HepSymMatrix Ea = mdcKalTrk->getZError();
   HepPoint3D point0(0.,0.,0.);
   HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
   VFHelix helixip3(point0,a,Ea);
   helixip3.pivot(IP);
   HepVector  vecipa = helixip3.a();

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


static bool isGoodChargedTrack(DstEvtRecTracks* tracks,
      Hep3Vector& xorigin) {
   if (!tracks->isMdcKalTrackValid()) {
      return false;
   }

   tracks->mdcKalTrack()->setPidType(BMdcKalTrack::pion);

   if (tracks->mdcKalTrack()->charge() == 0) {
      return false;
   }

   if ( fabs(cos(tracks->mdcKalTrack()->theta())) >= 0.93 ) {
      return false;
   }

   if (!isGoodPoca(tracks->mdcKalTrack(), xorigin)) {
      return false;
   }


   return true;
}


static double get_min_track_momentum(const TEvtRecDTag* evtRecDTag,
      const TObjArray* evtRecTrkCol, bool energy = false) {
   double min_p = 1E10;

   const vector<Int_t>& tracks = energy ? evtRecDTag->showers() :
         evtRecDTag->tracks();


   for (unsigned i = 0; i < tracks.size(); ++i) {
      DstEvtRecTracks* dstTracks =  (DstEvtRecTracks*) evtRecTrkCol->At(
                  tracks[i]);
      if (dstTracks->isMdcKalTrackValid()) {
         double p;

         if (energy) {
            p = dstTracks->emcShower()->energy();
         } else {
            p = dstTracks->mdcKalTrack()->p();
         }
         if (p  < min_p) {
            min_p = p;
         }
      }
   }
   return min_p;
}

//returns first particle found with given particle ID
TMcParticle* getParticleById(TMcEvent* m_TMcEvent, int particle_id) {
   TObjArrayIter next (m_TMcEvent->getMcParticleCol());
   while (TMcParticle* particle = (TMcParticle*) next()) {
      if (abs(particle->getParticleID()) == particle_id) {
         return particle;
      }
   }

   return 0;
}

TMcParticle* getParticleByIdAndDecay(TMcEvent* m_TMcEvent,
      int particle_id, int decay) {
   TIter next(m_TMcEvent->getMcParticleCol());
   while (TMcParticle* particle = (TMcParticle*) next()) {
      if (abs(particle->getParticleID()) == particle_id) {
         if (m_EventTagSvc->getDecayCode(particle) == decay) {
            return particle;
         }
      }
   }

   return 0;
}



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

   if ((hasDecay(eventTag_, 80, 8)
               && (m_EventTagSvc->getDecayCode(getParticleById(m_TMcEvent,
                                 313)) == 1)) ||
         (hasDecay(eventTag_, 80, 12)
               && (m_EventTagSvc->getDecayCode(getParticleById(m_TMcEvent,
                                 315)) == 1)) ) {
      sprintf(tag_str + strlen(tag_str), "%s", "-1");
   }
}

inline HepLorentzVector getMcParticleP4(const TMcParticle* particle) {
   return HepLorentzVector(
               particle->getInitialMomentumX(),
               particle->getInitialMomentumY(),
               particle->getInitialMomentumZ(),
               particle->getInitialMomentumE());
}


const vector<const TMcParticle*> getDaughtersById(
      TMcEvent* m_TMcEvent, const TMcParticle* mother, int particle_id) {
   const vector<Int_t>& daughters = mother->getDaughters();
   vector<const TMcParticle*> result;

   for (unsigned index = 0; index < daughters.size(); ++index) {
      const TMcParticle* daughter = m_TMcEvent->getMcParticle(
                  daughters[index]);
      if (abs(daughter->getParticleID()) == particle_id) {
         result.push_back(daughter);
      }
   }

   return result;

}

const TMcParticle* getDaughterById(TMcEvent* m_TMcEvent,
      const TMcParticle* mother, int particle_id) {
   const vector<const TMcParticle*>& daughters = getDaughtersById(
               m_TMcEvent, mother, particle_id);
   if (!daughters.empty()) {
      return daughters[0];
   } else {
      return 0;
   }
}

struct TAngularVars {
   double q2;
   double mkpi;
   double theta_v;
   double theta_l;
   double chi;
   double d_theta;
   double nu_theta;
};




TAngularVars getAngularVars(HepLorentzVector electron4Momentum,
      HepLorentzVector missing4Momentum,
      HepLorentzVector kaon4Momentum, HepLorentzVector pion4Momentum)

{
   TAngularVars vars;
   // calculate kinematic quantities

   vars.nu_theta = missing4Momentum.vect().theta();  //misc

   //masses
   vars.q2 = (electron4Momentum + missing4Momentum).m2();
   vars.mkpi = (kaon4Momentum + pion4Momentum).m();

   //angles
   // 1. go to D rest frame
   HepLorentzVector d4Momentum = electron4Momentum + missing4Momentum +
         kaon4Momentum + pion4Momentum;
   vars.d_theta =d4Momentum.vect().theta(); //misc

   electron4Momentum.boost(-d4Momentum.boostVector());
   missing4Momentum.boost(-d4Momentum.boostVector());
   kaon4Momentum.boost(-d4Momentum.boostVector());
   pion4Momentum.boost(-d4Momentum.boostVector());

   // 2. get k-pi
   HepLorentzVector kpi4Momentum = kaon4Momentum + pion4Momentum;
   //normal vector to the k-pi plane
   Hep3Vector kpiNormal = kaon4Momentum.vect().cross(
               pion4Momentum.vect());

   // 3. K to k-pi rest frame
   kaon4Momentum.boost(-kpi4Momentum.boostVector());
   vars.theta_v = (kpi4Momentum.vect()).angle(kaon4Momentum.vect());

   // 4. W rest frame
   HepLorentzVector lnu4Momentum = electron4Momentum + missing4Momentum;
   Hep3Vector lnuNormal = electron4Momentum.vect().cross(
               missing4Momentum.vect());

   // 5. electron to W rest frame
   electron4Momentum.boost(-lnu4Momentum.boostVector());
   vars.theta_l = (lnu4Momentum.vect()).angle(electron4Momentum.vect());

   // 6. Angle between planes
   vars.chi = lnuNormal.angle(kpiNormal);

   return vars;
}


void fillAngularVars(TNtuple* ntuple, const TAngularVars& vars,
      int tag_passed, bool passed, int decayMode, int charm) {
   //~ cout << "fillAngularVars: vars.mkpi=" << vars.mkpi << ", vars.q2=" << vars.q2 << " ,tag_passed=" << tag_passed << " , passed=" << passed << endl;
   ntuple->Fill(vars.q2, vars.mkpi, vars.theta_v, vars.theta_l, vars.chi,
         tag_passed, passed, vars.d_theta, vars.nu_theta, decayMode, charm);
}

void fillAngularVars(TNtuple* ntuple,
      HepLorentzVector electron4Momentum, HepLorentzVector missing4Momentum,
      HepLorentzVector kaon4Momentum, HepLorentzVector pion4Momentum,
      int tag_passed, bool passed, int decayMode, int charm) {
   TAngularVars vars = getAngularVars(electron4Momentum,
               missing4Momentum, kaon4Momentum, pion4Momentum);
   //   "q2:mkpi:theta_v:theta_l:chi");
   fillAngularVars(ntuple, vars, tag_passed, passed, decayMode, charm);
}


const TMcParticle* get_d_meson(TMcEvent* m_TMcEvent, int charm) {
   TMcParticle* psi_3770 =  getParticleByIdAndDecay(m_TMcEvent, 30443,
               80);
   if (!psi_3770) {
      return 0;
   }

   const vector<const TMcParticle*>& d_mesons =  getDaughtersById(
               m_TMcEvent, psi_3770, 411); // get D mesons
   assert(d_mesons.size() == 2);

   for (unsigned i = 0; i < 2;  ++i) {
      if (get_charge(d_mesons[i]) == charm) {
         return d_mesons[i];
      }
   }

   return 0;
}


bool getMcAngularVars(TMcEvent* m_TMcEvent, int charm,
      TAngularVars& vars) {
   auto d_meson = get_d_meson(m_TMcEvent, charm);
   if (!d_meson) {
      return false;
   }

   const TMcParticle* K = 0;
   const TMcParticle* pi = 0;
   const TMcParticle* nu = 0;
   const TMcParticle* e = 0;

   int decay = m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
               (d_meson));

   if  (decay == 18) { //  K-  pi+  e+ nu_e
      K = getDaughterById(m_TMcEvent, d_meson, 321);
      pi = getDaughterById(m_TMcEvent, d_meson, 211);
   } else if (decay == 8 ) { // K* e nu
      const TMcParticle* k_star = getDaughterById(m_TMcEvent, d_meson, 313);
      K = getDaughterById(m_TMcEvent, k_star, 321);
      pi = getDaughterById(m_TMcEvent, k_star, 211);
   } else if (decay ==
         12) {  // D+ => anti-K_2*0 e+  nu_e   , K_2*0 => K+  pi-
      const TMcParticle* k_star_2 = getDaughterById(m_TMcEvent, d_meson,
                  315);
      K = getDaughterById(m_TMcEvent, k_star_2, 321);
      pi = getDaughterById(m_TMcEvent, k_star_2, 211);

   }

   nu = getDaughterById(m_TMcEvent, d_meson, 12);
   e = getDaughterById(m_TMcEvent, d_meson, 11);


   if (K && pi && nu && e) {
      vars = getAngularVars(getMcParticleP4(e), getMcParticleP4(nu),
                  getMcParticleP4(K), getMcParticleP4(pi));
      return true;
   }

   return false;
}

bool fillMcTracks(TMcEvent* m_TMcEvent, int charm, TNtuple* ntp,
      int tag_passed, bool passed, int decayMode) {
   TAngularVars vars;
   if (getMcAngularVars(m_TMcEvent, charm, vars)) {
      fillAngularVars(ntp, vars, tag_passed, passed, decayMode, charm);
      return true;
   } else {
      return false;
   }
}

bool fillMcMomenta(TMcEvent* m_TMcEvent, int charm, TNtuple* ntp,
      int tag_passed, bool passed) {
   // begin copypaste
   auto d_meson = get_d_meson(m_TMcEvent, charm);
   if (!d_meson) {
      return false;
   }

   const TMcParticle* K = 0;
   const TMcParticle* pi = 0;
   const TMcParticle* nu = 0;
   const TMcParticle* e = 0;

   int decay = m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
               (d_meson));

   if  (decay == 18) { //  K-  pi+  e+ nu_e
      K = getDaughterById(m_TMcEvent, d_meson, 321);
      pi = getDaughterById(m_TMcEvent, d_meson, 211);
   } else if (decay == 8 ) { // K* e nu
      const TMcParticle* k_star = getDaughterById(m_TMcEvent, d_meson, 313);
      K = getDaughterById(m_TMcEvent, k_star, 321);
      pi = getDaughterById(m_TMcEvent, k_star, 211);
   } else if (decay ==
         12) {  // D+ => anti-K_2*0 e+  nu_e   , K_2*0 => K+  pi-
      const TMcParticle* k_star_2 = getDaughterById(m_TMcEvent, d_meson,
                  315);
      K = getDaughterById(m_TMcEvent, k_star_2, 321);
      pi = getDaughterById(m_TMcEvent, k_star_2, 211);

   }

   nu = getDaughterById(m_TMcEvent, d_meson, 12);
   e = getDaughterById(m_TMcEvent, d_meson, 11);

// end of copypaste

   if (K && pi && nu && e) {

      auto e_p4 = getMcParticleP4(e);
      auto nu_p4 = getMcParticleP4(nu);
      auto k_p4 = getMcParticleP4(K);
      auto pi_p4 = getMcParticleP4(pi);


      double vars_d[18] = {e_p4.px(), e_p4.py(), e_p4.pz(), e_p4.e(),
                  nu_p4.px(), nu_p4.py(), nu_p4.pz(), nu_p4.e(),
                  k_p4.px(), k_p4.py(), k_p4.pz(), k_p4.e(),
                  pi_p4.px(), pi_p4.py(), pi_p4.pz(), pi_p4.e(),
                  (double) tag_passed, (double)passed
            };
      float vars[18];
      for (size_t i = 0; i < sizeof(vars)/sizeof(vars[0]); ++i) {
         vars[i] = vars_d[i];
      }

      ntp->Fill(vars);

      return true;
   }
   return false;
}



#define NUM_ELEM(x) (sizeof (x) / sizeof (*(x)))



TEvtRecDTag* selectDTag(const TObjArray* m_evtRecDTagCol,
      const TObjArray* evtRecTrkCol, int charm, int decayMode) {
   double min_distance = 1E10;
   TEvtRecDTag* evtRecDTag = 0;
   TIter evtRecDTagIter(m_evtRecDTagCol);
   while (TEvtRecDTag* curEvtRecDTag = (TEvtRecDTag*)
               evtRecDTagIter.Next()) {
      if (curEvtRecDTag->type() == 1) { // require PID
         if (curEvtRecDTag->charm() == charm) {
            if (curEvtRecDTag->decayMode() == decayMode)  {
               double distance =   sqr( curEvtRecDTag->deltaE() / 4E-3) ;
               if (distance < min_distance) {
                  min_distance = distance;
                  evtRecDTag = curEvtRecDTag;
               }
            }
         }
      }
   }
   return evtRecDTag;
};

Hep3Vector getShowerVector(RecEmcShower* shower) {
   return Hep3Vector(shower->x(), shower->y(), shower->z());
}


bool isGoodShower(RecEmcShower* emcTrk, const TObjArray* evtRecTrkCol,
      const TEvtRecEvent*  m_evtRecEvent, int exclude_track = -1) {
   if (!emcTrk) {
      return false;
   }

   cout << "check emcTrk " << emcTrk->trackId() <<" energy=" <<
         emcTrk->energy() << " status=" << emcTrk->status() << endl;

   if ( (emcTrk->time() < 0) || (emcTrk->time() > 14)) {
      cout << "bad time" << endl;
      return false;
   }

   if (emcTrk->status() == 9999 ) {
      cout << "hot crystal skipped" << endl;
      return false;
   }

   // find minimum angle

   Hep3Vector showerVector = getShowerVector(emcTrk);

   bool tooCloseToCharged = false;

   for (int charged_id = 0; charged_id < m_evtRecEvent->totalCharged();
         ++charged_id) {
      if (charged_id != exclude_track) {
         DstEvtRecTracks* charged = (DstEvtRecTracks*) evtRecTrkCol->At(
                     charged_id);
         if (! charged->isExtTrackValid()) {
            continue;
         }
         if (charged->extTrack()->emcVolumeNumber() == -1) {
            continue;
         }

         const Hep3Vector& chargedMomentum =
               charged->extTrack()->emcPosition();

         // fixme initial momentum
         double angle = showerVector.angle(chargedMomentum);

         if (angle < 15 * TMath::DegToRad()) {
            cout << "charged tr id " << charged_id << "angle: " << angle << endl;
            tooCloseToCharged = true;
            break;
         }
      }
   }

   if (tooCloseToCharged) {
      cout << "too close to charged" << endl;
      return false;
   }

   double energy = emcTrk->energy();
   double absCosTheta = TMath::Abs(  TMath::Cos(emcTrk->theta()) );

   if  (   ( ( absCosTheta < 0.8 ) &&  ( energy > 25E-3) ) ||
         ( ( absCosTheta > 0.86 ) && ( absCosTheta < 0.92 )
               && ( energy > 50E-3)  )   ) {
      return true;
   } else {
      cout << "bad energy" << endl;
      return false;
   }
}



bool get_good_d_to_k_pi_e_nu_mass(TMcEvent* m_TMcEvent,
      const TMcParticle* d_meson, double& mass) {
   int decay = m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
               (d_meson));

   if  (decay == 18) { //  K-  pi+  e+ nu_e
      const TMcParticle* K = getDaughterById(m_TMcEvent, d_meson, 321);
      const TMcParticle* pi = getDaughterById(m_TMcEvent, d_meson, 211);
      mass = (getMcParticleP4(K) + getMcParticleP4(pi)).m();
      return true;
   } else if (decay == 8 ) { // K* e nu
      const TMcParticle* k_star = getDaughterById(m_TMcEvent, d_meson, 313);
      if (m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
                  (k_star)) == 1) { // K+  pi-
         mass = getMcParticleP4(k_star).m();
         return true;
      }
   } else if (decay ==
         12) {  // D+ => anti-K_2*0 e+  nu_e   , K_2*0 => K+  pi-
      const TMcParticle* k_star_2 = getDaughterById(m_TMcEvent, d_meson,
                  315);
      if (m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
                  (k_star_2)) == 1) { // K+  pi-
         mass = getMcParticleP4(k_star_2).m();
         return true;
      }

   }

   mass = -10000;
   return false;
}



bool is_good_k_pi_e_nu(TMcEvent* m_TMcEvent, int charm,
      double mass_low = -1E10, double mass_high = +1E10,
      double* out_mass = 0) {
   auto d_meson = get_d_meson(m_TMcEvent, charm);
   if (!d_meson) {
      return false;
   }

   double mass;
   if (get_good_d_to_k_pi_e_nu_mass(m_TMcEvent, d_meson, mass)) {
      if ( (mass >= mass_low) && (mass <= mass_high)) {
         if (out_mass) {
            *out_mass = mass;
         }

         return true;
      }
   }

   return false;
}

double getTrueKPiMass(TMcEvent* m_TMcEvent, int charm) {
   double mass;
   bool good = is_good_k_pi_e_nu(m_TMcEvent, charm, -1E10, +1E10, &mass);
   assert(good);
   return mass;
}




static multiset<long> to_pdg_set(const vector<string>&
      particle_list) {
   multiset<long> result;
   for (auto it = particle_list.begin(); it != particle_list.end();
         ++it) {
      long int pdg = m_EventTagSvc->name2pdg(*it);
      assert( pdg != 0);
      result.insert(pdg);
   }
   return result;
};

static void  print_pdg_set(const multiset<long>& pdg_set) {
   cout << "{";
   for (auto it = pdg_set.begin(); it != pdg_set.end(); ++it) {
      cout << m_EventTagSvc->pdg2name(*it) <<  " ";
   }
   cout << "}" << endl;
}

static const multiset<long> get_final_state(TMcEvent* m_TMcEvent,
      const TMcParticle* parent, bool norm_sign = true) {
   static const multiset<long> fs_particles = to_pdg_set(vector<string>({"K+", "K-", "K_S0", "pi0", "pi+", "pi-", "e+", "e-", "mu+", "mu-", "nu_e", "anti-nu_e", "nu_mu", "anti-nu_mu"}));
   static const multiset<long> ignore_particles = to_pdg_set(
               vector<string>({"gammaFSR"}));

   multiset<long> final_state;

   if (parent) {
      const vector<Int_t>& daughters = parent->getDaughters();



      if (daughters.empty()
            || fs_particles.count(parent->getParticleID())) {
         // no further decay
         if (ignore_particles.count(parent->getParticleID()) ==
               0) { // not a FSR photon
            long particle_id = parent->getParticleID();
            if (norm_sign) {
               particle_id = abs(particle_id);
            }
            final_state.insert(particle_id);
         }
      } else {
         for (size_t index = 0; index < daughters.size(); ++index) {
            const TMcParticle* daughter = m_TMcEvent->getMcParticle(
                        daughters[index]);
            auto daughter_final_state = get_final_state(m_TMcEvent, daughter,
                        norm_sign);
            final_state.insert(daughter_final_state.begin(),
                  daughter_final_state.end());
         }
      }
   }

   return final_state;
};



int get_mc_tag_mode(TMcEvent* m_TMcEvent,
      int charm) { //charm of tag side
   static multimap<int, multiset<long> > tag_finals {
      { 200, to_pdg_set(vector<string>({"K+","pi+", "pi+"})) },
      { 201, to_pdg_set(vector<string>({"K+", "pi+", "pi+", "pi0"})) },
      { 202, to_pdg_set(vector<string>({"K_S0", "pi+"})) },
      { 203, to_pdg_set(vector<string>({"K_S0", "pi+", "pi0"})) },
      { 203, to_pdg_set(vector<string>({"pi+", "pi+", "pi+", "pi0"})) },
      { 204, to_pdg_set(vector<string>({"K_S0", "pi+", "pi+", "pi+"})) },
      { 204, to_pdg_set(vector<string>({"K_S0", "pi+", "K_S0"})) },
      { 204, to_pdg_set(vector<string>({"pi+", "pi+", "pi+", "pi+", "pi+"})) },
      { 205, to_pdg_set(vector<string>({"K+","K+", "pi+"})) },

   };

   const TMcParticle* tag_d_meson = get_d_meson(m_TMcEvent, charm);
   auto tag_final_state = get_final_state(m_TMcEvent,  tag_d_meson);

   for (auto it = tag_finals.begin(); it != tag_finals.end(); ++it) {
      if (it->second == tag_final_state) {
         return it->first;
      }
   }
   return -1;
}



BMdcKalTrack::PidType get_kal_pid_type(const TMcParticle* particle) {
   int mc_particle_id = abs(particle->getParticleID());

   BMdcKalTrack::PidType kal_pid_type = BMdcKalTrack::null;
   if (mc_particle_id == 11) {
      kal_pid_type = BMdcKalTrack::electron;
   } else if (mc_particle_id == 13) {
      kal_pid_type = BMdcKalTrack::muon;
   } else if (mc_particle_id == 211) {
      kal_pid_type = BMdcKalTrack::pion;
   } else if (mc_particle_id == 321) {
      kal_pid_type = BMdcKalTrack::kaon;
   } else if (mc_particle_id == 2212) {
      kal_pid_type = BMdcKalTrack::proton;
   }

   return kal_pid_type;
}

bool get_track_momentum(const DstEvtRecTracks* tracks,
      BMdcKalTrack::PidType kal_pid_type, Hep3Vector& p3) {
   assert(tracks);

   if (kal_pid_type != BMdcKalTrack::null) {
      // try to use mdc kal track
      if (tracks->isMdcKalTrackValid()) {
         tracks->mdcKalTrack()->setPidType(BMdcKalTrack::pion);
         p3 = tracks->mdcKalTrack()->p3();
         return true;
      }
   }

   // fallback to mdctrack
   if (tracks->isMdcTrackValid()) {
      p3 = tracks->mdcTrack()->p3();
      return true;
   }

   return false;

}

DstEvtRecTracks* get_closest_track(const TObjArray* evtRecTrkCol,
      const vector<Int_t>& track_ids, const TMcParticle* ref_particle,
      double threshold = 1.0) {
   Hep3Vector ref_p3 = getMcParticleP4(ref_particle).vect();
   BMdcKalTrack::PidType kal_pid_type = get_kal_pid_type(ref_particle);

   std::map<double, int> dist_map;

   for (size_t index = 0; index < track_ids.size(); ++index) {
      DstEvtRecTracks* tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                  track_ids[index]);

      Hep3Vector p3;
      if (get_track_momentum(tracks, kal_pid_type, p3)) {
         double distance = (p3 - ref_p3).mag();
         dist_map[distance] = track_ids[index];
      }

   }

   auto first = dist_map.begin();


   if (first != dist_map.end()) {
      // at least one element
      auto second = ++dist_map.begin();

      if (second != dist_map.end()) {
         // second min distance. check it against threshold
         double ratio = (first->first / second->first);
         if (ratio > threshold) {
            cout << "big ratio: " << ratio << endl;
            cout << "track_ids len" << track_ids.size() << endl;
            return 0;
         }
      }

      return (DstEvtRecTracks*) evtRecTrkCol->At(first->second);


   }


   return 0;

}


#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void SemileptonicStartJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{

   if( selector->Verbose() ) {
      cout << " TestPIDStartJob() " << endl;
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


   VecObj his1, his_aux, his1_init, his_dtag;

   his_aux.push_back( new TH1D("hDp_pMiss_x", "miss momentum", 200, -1,
               1));
   his_aux.push_back( new TH1D("hDp_pMiss_y", "miss momentum", 200, -1,
               1));
   his_aux.push_back( new TH1D("hDp_pMiss_z", "miss momentum", 200, -1,
               1));
   his_aux.push_back( new TH1D("hDp_pMiss_e", "miss momentum", 200, -1,
               1));



   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_mBC", "Beam constrained mass of selected  tagged Dp", 100, 1.862, 1.878);
   smartBookHist(his1_init, his1, his_dtag, "hDp_dE",
         "delta energy  of selected  tagged Dp;#Delta E (GeV)", 750, -0.15,
         0.15);
   smartBookHist(his1_init, his1, his_dtag, "hDp_mBC",
         "Beam constrained mass of selected  tagged Dp;m_{BC} (GeV)", 400,
         1.83, 1.911);



   smartBookHist(his1_init, his1, his_dtag, "_hDp_dE_plus",
         "delta energy  of selected  tagged Dp", 750, -0.15, 0.15);
   smartBookHist(his1_init, his1, his_dtag, "_hDp_dE_minus",
         "delta energy  of selected  tagged Dp", 750, -0.15, 0.15);

   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_dE",
         "minimum dE of D0 tag;#Delta E (GeV)", 100, -0.1, 0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_dE_sideband",
         "minimum dE of D0 tag (mass sideband);#Delta E (GeV)", 100, -0.1,
         0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_dE_sideband2",
         "minimum dE of D0 tag (mass sideband);#Delta E (GeV)", 100, -0.1,
         0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_1_dE",
         "minimum dE of D0 tag 1;#Delta E (GeV)", 100, -0.1, 0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_3_dE",
         "minimum dE of D0 tag 3;#Delta E (GeV)", 100, -0.1, 0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_201_3_dE",
         "minimum dE of D0 tag;#Delta E (GeV)", 100, -0.1, 0.1);

   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_mode",
         "Decay mode of tagged D0", 510, 0,509);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_mode_sideband",
         "Decay mode of tagged D0 (mass sideband)", 510, 0,509);
   smartBookHist(his1_init, his1, his_dtag, "hDp_d0_mode_sideband2",
         "Decay mode of tagged D0 (mass sideband)", 510, 0,509);


   //~ for (int decay = 200; decay < 210; ++decay) {
   //~ smartBookHist(his1_init, his1, his_dtag , TString::Format("hDp_dE_%d", decay), "delta energy  of selected  tagged Dp", 200, -0.15, 0.15);
   //~ smartBookHist(his1_init, his1, his_dtag , TString::Format("hDp_mBC_%d", decay), "Beam constrained mass of selected  tagged Dp", 400, 1.83, 1.911);
   //~ }

   smartBookHist(his1_init, his1, his_dtag, "hDp_dE_mcgood",
         "delta energy  of selected  tagged Dp", 750, -0.15, 0.15);

   smartBookHist(his1_init, his_aux, his_dtag, "hDTag_DecayMode",
         "Decay mode of tagged D", 510, 0,509);

   smartBookHist(his1_init, his1, his_dtag, "hDp_mKPi",
         "Missing mass of K-#pi pair of %s events", 300, 0, 2);



   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss",
         "Missing mass squared of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         300, -0.1, 0.2);
   smartBookHist(his1_init, his1, his_dtag, "hDp_mMissBC",
         "BC Missing mass squared of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         300, -0.1, 0.2);

   smartBookHist(his1_init, his1, his_dtag, "hDp_uMiss",
         "Umiss = E - |p| %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         400, -0.2, 0.2);
   smartBookHist(his1_init, his1, his_dtag, "hDp_uMiss_cms",
         "CMS Umiss = E - |p| %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         400, -0.2, 0.2);
   smartBookHist(his1_init, his1, his_dtag, "hDp_uMiss_cms2",
         "CMS Umiss = E - |p| %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         400, -0.2, 0.2);
   smartBookHist(his1_init, his1, his_dtag, "hDp_uMiss_cms_bc",
         "CMS Umiss = E - |p| %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         400, -0.2, 0.2);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_mMissAlt",
         "Missing mass squared(using four-vectors) of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.2, 0.2);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_mMissInitial",
         "Missing mass squared(before any cuts) of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.2, 0.2);


   smartBookHist(his1_init, his1, his_dtag, "hDp_goodOtherSide",
         "Good D+ candidate on the other side of %s ", 2, -.5, 1.5);

   smartBookHist(his1_init, his1, his_dtag, "hDp_nPions",
         "Number of %s pion candidates", 6, 0, 5);
   smartBookHist(his1_init, his1, his_dtag, "hDp_nKaons",
         "Number of %s kaon candidates", 6, 0, 5);

   smartBookHist(his1_init, his1, his_dtag, "hDp_n3tr_comb",
         "Number of %s 3 tracks combinations", 3, -0.5, 2.5);





   smartBookHist(his1_init, his_aux, his_dtag, "hDp_nCharged",
         "Number of %s charged candidates", 16, 0, 15);

   smartBookHist(his1_init, his1, his_dtag, "hDp_TotCharge",
         "Total charge of %s ", 20, -10.5, 9.5);




   smartBookHist(his1_init, his_aux, his_dtag, "hDp_nShowers",
         "Number of %s chower candidates", 16, 0, 15);

   smartBookHist(his1_init, his1, his_dtag, "hDp_pMiss",
         "Missing momentum of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.1, 0.9);
   smartBookHist(his1_init, his1, his_dtag, "hDp_EMiss",
         "Missing energy of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.1, 0.9);


   smartBookHist(his1_init, his1, his_dtag, "hDp_missCosTheta",
         "Cosine of polar angle of missing momentum of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         50, 0, 1);

   //~ smartBookHist(his1, "hDp_Lxy_kpnue", "L_{xy} of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}", 100, 0, 2);
   //~ smartBookHist(his1, "hDp_Lxy_k0nue", "L_{xy} of %s D^{#pm} #Rightarrow Ks^{0} e^{#pm} #nu_{e}", 100, 0, 2);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_chiPrimVertex",
         "#chi^2 of primary vertex fit of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}.",
         201,0,200 );
   smartBookHist(his1_init, his1, his_dtag, "hDp_mPiPi",
         "Invariant mass of %s K^{#pm} #pi^{#mp}  pair with mass of K^{#pm} candidate fixed at pi mass.",
         160,0.2,1.4 );

   smartBookHist(his1_init, his1, his_dtag, "hDp_nPi0",
         "Number of #pi^{0} candidates with #chi^{2} < 20 (%s events)", 7, 0,
         6);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_pi0Chisq",
         "Distribution of #chi^{2} of  #pi^{0} candidates (%s events) ", 1000,
         0, 1000);

   smartBookHist(his1_init, his1, his_dtag, "hDp_nGoodShowers",
         "N of good showers of  %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         7, -0.5, 6.5);
   smartBookHist(his1_init, his1, his_dtag, "hDp_showerE", "%s", 140, 0,
         0.7);
   smartBookHist(his1_init, his1, his_dtag, "hDp_maxShowerE", "%s", 140,
         0, 0.7);
   smartBookHist(his1_init, his1, his_dtag, "hDp_totalShowersE",
         "Total unassosiated shower energy of %s events", 100, 0, 1);


   smartBookHist(his1_init, his1, his_dtag, "hDp_costE",
         "cos(theta electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, -1, 1);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_pxyK",
         "XY Momentum of K^{#pm} candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_pxyE",
         "XY Momentum electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_pxyPi",
         "XY Momentum of pi candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_pE",
         "Momentum of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_pPi",
         "Momentum of pion candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_pK",
         "Momentum of kaon candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         100, 0, 1);

   smartBookHist(his1_init, his1, his_dtag, "hDp_eE",
         "EMC Energy of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ",
         200, 0, 2);




   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_pE_sm_kal", "TMdcKalTrack Momentum of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_pE_sm", "TMdcTrack Momentum of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);

   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_pE_sm_mupi", "TMdcTrack Momentum of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_pE_sm_mue", "TMdcTrack Momentum of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);


   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_sm_pRealMu", "TMdcTrack Momentum of muon misidentified as E or Pi in sm of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}", 100, 0, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_sm_pRealPi", "TMdcTrack Momentum of pion misidentified as Mu or Pi in sm of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}", 100, 0, 1);


   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_thE_sm_mue","Cosine of theta E of semileptonic mue", 100, -1, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_thE_sm_mupi","Cosine of theta E of semileptonic mue", 100, -1, 1);

   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_phiE_sm_mue","phi E of semileptonic mue", 362, -180.5, 180.5);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_phiE_sm_mupi","phi E of semileptonic mue", 362, -180.5, 180.5);



   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_pMu", "MC Momentum of muon in  %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_pMu_mue", " mu=>e MC Momentum of muon in  %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 100, 0, 1);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_eMu", "MC full energy of muon in  %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 150, 0, 1.5);

   smartBookHist(his1_init, his1, his_dtag, "hDp_vi_MC_pMu",
         "initial MC Momentum of muon in  %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)",
         100, 0, 1);

   smartBookHist(his1_init, his1, his_dtag, "hDp_MC_pPi",
         "MC Momentum of pion in semimuonic", 100, 0, 1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_vi_MC_pPi",
         "vi MC Momentum of pion in semimuonic", 100, 0, 1);


   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_pFakeE", "MC Momentum of whatever identified as electron (pion or muon) in semimuonic", 100, 0, 1);



   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pE_withTOF", "Momentum of electron candidates with valid TOF of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, 0, 2);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pE_withTOF_E", "Momentum of electron candidates with valid TOF E of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, 0, 2);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pE_withTOF_Q", "Momentum of electron candidates with valid TOF Q of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, 0, 2);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pE_withTOF_C", "Momentum of electron candidates with valid TOF C of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, 0, 2);

   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_tof", "TOF of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, -20, 20);


   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_EMCeE", "Shower energy of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 150, 0, 1.5);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_EMCeE_sm_mupi", "Shower energy of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 150, 0, 1.5);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_EMCeE_sm_mue", "Shower energy of electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} (semimuonic excbad)", 150, 0, 1.5);





   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_K_as_pi",
         "Missing mass of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} with Kaon mass replaced with Pion mass",
         140, -0.1, 0.6);
   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_e_as_mu",
         "Missing mass of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} with E mass replaced with Mu mass",
         100, -0.08, 0.08);

   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_mMiss_K_as_pi2", "Missing mass of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} with Kaon mass replaced with Pion mass", 400, -2 , 2);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_mMiss_K_as_pi_minus_mMiss", "Missing mass K-as-pi minus missing mass of %s tracks ", 400, -2 , 2);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_eRelProbE",
         "probE normalized to sum PID prob. %s tracks, electron candidate ",
         800, 0, 10000);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_eRelProbMu",
         "probMu normalized to sum PID prob. %s tracks, electron candidate ",
         1000, 0,0.1);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_eLogRelProbMu",
         " log probMu normalized to sum PID prob. %s tracks, electron candidate ",
         1000, -100,100);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_nMucLayers",
         "N of layers in MUC with hits. %s tracks, electron candidate. (inclusive)",
         15, -0.5, 14.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_Ntracks",
         "N of signal side tracks( %s events)", 15, -0.5, 14.5);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTag",
         "Has good D+", 3, -0.5, 2.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTag2",
         "Has good D+", 3, -0.5, 2.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTagFlag",
         "Has good D+", 3, -0.5, 2.5);

   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTagMode",
         "Mode of good D+", 501, -1.5, 499.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTagMode2",
         "Mode of good D+", 501, -1.5, 499.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodTagMode3",
         "Mode of good D+", 501, -1.5, 499.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hDp_goodMcTagMode",
         "Mode of good D+", 501, -1.5, 499.5);




   smartBookHist(his1_init, his1, his_dtag, "hDp_epEMC",
         "EMC energy deposit / momentum of electron candidates", 120, 0, 1.2);

   smartBookHist(his1_init, his1, his_dtag, "hDp_epEMC_sm_mupi",
         "EMC energy deposit / momentum of electron candidates, semileptonic mu as pi",
         120, 0, 1.2);
   smartBookHist(his1_init, his1, his_dtag, "hDp_epEMC_sm_mue",
         "EMC energy deposit / momentum of electron candidates, semileptonic mu as e",
         120, 0, 1.2);

   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_epEMC_sm", "Rec. EMC energy deposit / momentum of muon./, semileptonic, before PID", 120, 0, 1.2);
   //~ smartBookHist(his1_init, his1, his_dtag , "hDp_MC_dedx_sm", "Rec. dedx of muon, semileptonic, before PID", 200, 0, 2);


   //~ smartBookHist(his1, "hDP_dEdx", "dE/dx %s electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 200, 0, 2);

   //~ smartBookHist(his1, "hDP_dEdxChiE", "dE/dx #chi(E) %s electron candidates of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e} ", 120, 0, 6);


   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_mu",
         "Missing mass of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}. (bad is for semimuonic)",
         200, -0.5, 0.5);

   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_vxfit",
         "Missing mass squared (after VertexFit) of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.1, 0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_err",
         "Missing mass squared absolute error of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, 0, 0.04);

   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_alt",
         "Missing mass squared (P4_cms - sum(P4_i) ) of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.1, 0.1);
   smartBookHist(his1_init, his1, his_dtag, "hDp_mMiss_alt2",
         "Missing mass squared (using sum(p3_i) for dtag side ) of %s D^{#pm} #Rightarrow K^{#pm} #pi^{#mp} e^{#pm} #nu_{e}",
         100, -0.1, 0.1);

   // FSR
   smartBookHist(his1_init, his1, his_dtag, "hDp_nFSR",
         "Number of FSR photons in %s ", 6, 0, 5);
   smartBookHist(his1_init, his1, his_dtag, "hDp_FSRE",
         "Total energy of FSR photons in %s", 100, 0, 0.5);

   smartBookHist(his1_init, his_aux, his_dtag, "hViSignalDecays",
         "very initial decay codes", 302, -1.5, 300.5);
   smartBookHist(his1_init, his_aux, his_dtag, "hSignalDecays",
         "very initial decay codes", 302, -1.5, 300.5);
   smartBookHist(his1_init, his_aux, his_dtag,
         "hGoodTagFlagSignalDecays", "very initial decay codes", 302, -1.5,
         300.5);


   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_e_r", "|V_xy| of %s electron tracks ", 200, 0, 10);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_kaon_r", "|V_xy| of %s kaon tracks ", 200, 0, 10);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pion_r", "|V_xy| of %s pion tracks ", 200, 0, 10);
//~
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_e_z", "|V_z| of %s electron tracks ", 200, 0, 20);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_kaon_z", "|V_z| of %s kaon tracks ", 200, 0, 20);
   //~ smartBookHist(his1_init, his_aux, his_dtag , "hDp_pion_z", "|V_z| of %s pion tracks ", 200, 0, 20);


   //~ his_aux.push_back(new TH1D("hDp_chisq", "#chi^{2} of kinematic fitting with missing neutrino mass and reusing D-Tag momentum. D mass is fixed", 600, 0, 200));
   //~ his_aux.push_back(new TH1D("hDp_Z0E", "Energy of Z0", 400, 3.700, 3.800));

   his_aux.push_back(new TH1D("hDp_mKPiTruthInitial",
               "True Missing mass of K-#pi pair of events (intitial histo)", 500, 0,
               2));
   his_aux.push_back(new TH1D("hDp_mKPiTruthInitialGoodTag",
               "True Missing mass of K-#pi pair of events (intitial histo)", 500, 0,
               2));


   his_aux.push_back(new TH1D("hDp_mKPiTruthVeryInitial",
               "True Missing mass of K-#pi pair of events (intitial histo)", 500, 0,
               2));
   his_aux.push_back(new TH1D("hDp_mKPiTruth",
               "True Missing mass of K-#pi pair of events ", 500, 0, 2));
   his_aux.push_back(new TH1D("hDp_mKPiTruthKStar",
               "True Missing mass of K-#pi pair of K star events (intitial histo)",
               500, 0, 2));

   his_aux.push_back(new TH1D("hDp_mBC_mctag_good",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));
   his_aux.push_back(new TH1D("hDp_mBC_mctag_bad",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));
   his_aux.push_back(new TH1D("hDp_mBC_tagmode_good",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));
   his_aux.push_back(new TH1D("hDp_mBC_tagmode_bad",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));
   his_aux.push_back(new TH1D("hDp_mBC_tagflag_good",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));
   his_aux.push_back(new TH1D("hDp_mBC_tagflag_bad",
               "Beam constrained mass of selected  tagged Dp (bad tags)", 400, 1.83,
               1.911));

   his_aux.push_back(new TH1D("hDp_4MDist",
               "4 momentum distance Dtag-MCtag", 600, -0.1, 1));
   his_aux.push_back(new TH1D("hDp_4MDist_mctag_good",
               "4 momentum distance Dtag-MCtag", 600, -0.1, 1));
   his_aux.push_back(new TH1D("hDp_4MDist_mctag_bad",
               "4 momentum distance Dtag-MCtag", 600, -0.1, 1));
   his_aux.push_back(new TH1D("hDp_4MDist_tagmode_good",
               "4 momentum distance Dtag-MCtag", 600, -0.1, 1));
   his_aux.push_back(new TH1D("hDp_4MDist_tagmode_bad",
               "4 momentum distance Dtag-MCtag", 600, -0.1, 1));


   //~ his_aux.push_back( new TH2D("hDp_MC_PvsEMC_sm_-", "mu- EMC deposit vs p of true muons (semileptonic events);p (GeV/c);EMC deposit (GeV)", 100, 0, 1,  100, 0, 1));


   his1.push_back( new TH2D("hDp_e_p_vs_cost",
               "p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)", 11, 0.1,
               1.2, 20, -1.0, 1.0));
   his1.push_back( new TH2D("hDp_e_p_vs_cost_good",
               "(good) p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)",
               11, 0.1, 1.2, 20, -1.0, 1.0));
   his1.push_back( new TH2D("s_hDp_e_p_vs_cost_good",
               "(good) p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)",
               11, 0.1, 1.2, 20, -1.0, 1.0));

   his1.push_back( new TH2D("hDp_e_p_vs_cost_we",
               "p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)", 11, 0.1,
               1.2, 20, -1.0, 1.0));
   his1.push_back( new TH2D("hDp_e_p_vs_cost_we_good",
               "(good) p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)",
               11, 0.1, 1.2, 20, -1.0, 1.0));
   his1.push_back( new TH2D("s_hDp_e_p_vs_cost_we_good",
               "(good) p vs cos(t) for electron candidates;p_e (GeV);cos(#theta)",
               11, 0.1, 1.2, 20, -1.0, 1.0));




   hists_dtag_decay_good_mBC.resize(500);
   hists_dtag_decay_good_dE.resize(500);
   hists_dtag_decay_good_Ntags.resize(500);

   hists_dtag_decay_bad_mBC.resize(500);
   hists_dtag_decay_bad_dE.resize(500);

   hists_dtag_decay_selected_bad_mBC.resize(500);
   hists_dtag_decay_selected_good_mBC.resize(500);
   hists_dtag_decay_selected_mBC.resize(500);
   hists_dtag_decay_good_mBC_vs_dE.resize(500);
   hists_dtag_decay_bad_mBC_vs_dE.resize(500);

   for (int decay = 200; decay < 210; ++decay) {
      hists_dtag_decay_good_mBC[decay]   = new TH1D(
            TString::Format("hDp_decay_%d_good_mBC", decay),
            "Beam constrained mass of good kDptoKPiPi	(200)", 100, 1.83, 1.9);
      hists_dtag_decay_good_dE[decay]    = new TH1D(
            TString::Format("hDp_decay_%d_good_dE", decay),
            "delta energy  of good kDptoKPiPi	(200)", 100, -0.1, 0.1);
      hists_dtag_decay_good_Ntags[decay] = new TH1I(
            TString::Format("hDp_decay_%d_good_Ntags", decay),
            "Count of tags per good event DptoKPiPi	(200)", 11, 0, 10);

      hists_dtag_decay_bad_mBC[decay]   = new TH1D(
            TString::Format("hDp_decay_%d_bad_mBC", decay),
            "Beam constrained mass of bad kDptoKPiPi	", 100, 1.83, 1.9);
      hists_dtag_decay_bad_dE[decay]    = new TH1D(
            TString::Format("hDp_decay_%d_bad_dE", decay),
            "delta energy  of bad kDptoKPiPi", 100, -0.1, 0.1);


      hists_dtag_decay_selected_bad_mBC[decay] = new TH1D(
            TString::Format("hDp_decay_%d_selected_bad_mBC", decay),
            "Beam constrained mass of bad kDptoKPiPi passed the selection of deltaE",
            100, 1.83, 1.9);
      hists_dtag_decay_selected_good_mBC[decay] = new TH1D(
            TString::Format("hDp_decay_%d_selected_good_mBC", decay),
            "Beam constrained mass of good kDptoKPiPi passed the selection of deltaE",
            100, 1.83, 1.9);
      hists_dtag_decay_selected_mBC[decay] = new TH1D(
            TString::Format("hDp_decay_%d_selected_total_mBC", decay),
            "Beam constrained mass of all kDptoKPiPi passed the selection of deltaE",
            100, 1.83, 1.9);

      hists_dtag_decay_good_mBC_vs_dE[decay] = new TH2D(
            TString::Format("hDp_decay_%d_good_mBC_vs_dE", decay),
            "Beam constrained mass of all good tracks versus deltaE", 100, 1.83,
            1.9, 100, -0.1, 0.1);

      hists_dtag_decay_bad_mBC_vs_dE[decay] = new TH2D(
            TString::Format("hDp_decay_%d_bad_mBC_vs_dE", decay),
            "Beam constrained mass of bad tracks versus deltaE", 100, 1.83, 1.9,
            100, -0.1, 0.1);



   }

   selector->RegInDir(hists_dtag_decay_good_mBC,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_good_dE,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_good_Ntags,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_bad_mBC,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_bad_dE,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_selected_bad_mBC,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_selected_good_mBC,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_selected_mBC,
         "Semileptonic_Dtag_decays");

   selector->RegInDir(hists_dtag_decay_good_mBC_vs_dE,
         "Semileptonic_Dtag_decays");
   selector->RegInDir(hists_dtag_decay_bad_mBC_vs_dE,
         "Semileptonic_Dtag_decays");



   his_aux.push_back(new TH1I("hBadTags", "Bad tags", 11, 1,0));
   his_aux.push_back(new TH1I("hAllTags", "All tags", 11, 1,0));
   his_aux.push_back(new TH1I("hInitialAllTags", "All tags", 11, 1,0));
   his_aux.push_back(new TH1I("hInitialBadTags", "hInitialBadTags tags",
               11, 1,0));
   his_aux.push_back(new TH1I("hVeryInitialBadTags",
               "hVeryInitialBadTags tags", 11, 1,0));
   his_aux.push_back(new TH1I("hIntermediateAllTags",
               "All tags after 3tracks", 11, 1,0));


   his_aux.push_back(new TH1I("hRealEventSurvived",
               "Number of \"events\" survived in one real event", 13, -0.5, 12.5));
   his_aux.push_back(new TH1I("hCharmSurvived",
               "Number of \"events\" with same charm survived in one real event", 7,
               -0.5, 6.5));




   his_aux.push_back(new TH1I("hAux", "aux", 11, 1,0));
   his_aux.push_back(new TH1I("hBadTags_e", "Bad tags", 11, 1,0));


   for (int decay = 200; decay < 206; ++decay) {
      his_aux.push_back(new TH1I(TString::Format("hGoodWrongTags_%d",
                        decay), "GoodWrongTags", 11, 1,0));
      his_aux.push_back(new TH1I(TString::Format("hBadTags_%d", decay),
                  "Bad tags for mode", 11, 1,0));
   }


   his_aux.push_back(new TH1I("hCuts_good",
               "Count of good events after each cut", 11, 1,0));
   his_aux.push_back(new TH1I("hCuts_bad",
               "Count of bad events after each cut", 11, 1,0));

   his_aux.push_back(new TH1I("hGoodMass_vi",
               "Is good event with good kpi mass, at the very begining", 2, -0.5,
               1.5));
   his_aux.push_back(new TH1I("hGoodMass_hastag",
               "Is good event with good kpi mass, after some tag selection, before dE and mBC cuts",
               2, -0.5, 1.5));
   his_aux.push_back(new TH1I("hGoodMass_mcgoodtag",
               "Is good event with good kpi mass, after some tag selection, before dE and mBC cuts provided tag was good (MC)",
               2, -0.5, 1.5));
   his_aux.push_back(new TH1I("hGoodMass_mcgoodtag2",
               "Is good event with good kpi mass, after some tag selection, before dE and mBC cuts provided tag was good (MC)",
               2, -0.5, 1.5));
   his_aux.push_back(new TH1I("hGoodMass_good_tag_flag",
               "Is good event with good kpi mass, after some tag selection, before dE and mBC cuts provided tag was good (MC)",
               2, -0.5, 1.5));
   his_aux.push_back(new TH1I("hGoodMass_d",
               "Is good event with good kpi mass, after d selection", 2, -0.5, 1.5));


   cuts_boundaries = new TMap();
   cuts_boundaries->SetName("cuts_boundaries");
   his1.push_back((TNamed*)cuts_boundaries);

   selector->RegInDir(his1,"Semileptonic_Dtag");
   selector->RegInDir(his1_init,"Semileptonic_Dtag_initial");
   selector->RegInDir(his_aux,"Semileptonic_Dtag_Aux");
   selector->RegInDir(his_dtag,"Semileptonic_Dtag_Modes");
   VecObj his2;
   his2.resize(100);



   m_EventTagSvc = EventTagSvc::instance();
   m_EventTagSvc->setPdtFile(
         selector->AbsPath("Analysis/EventTag/share/pdt.table") );
   m_EventTagSvc->setDecayTabsFile(
         selector->AbsPath("Analysis/EventTag/share/decay.codes") );


   m_EventTagSvc->addUserChainTrig("anti-K*0");
   m_EventTagSvc->addUserChainTrig("K*0");
   m_EventTagSvc->addUserChainTrig("anti-K_2*0");
   m_EventTagSvc->addUserChainTrig("K_2*0");

   //~ m_EventTagSvc->setUserDecayTabsFile("BeanUser/semileptonic.codes");

   m_EventTagSvc -> initialize();



   histograms_cache.AddAll(selector->GetOutputList());
   hist("hBadTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hAllTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hInitialAllTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hInitialBadTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hVeryInitialBadTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hIntermediateAllTags")->Fill("blah",
         0.0); //one MUST keep this for correct merging!


   hist("hAux")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hBadTags_e")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hCuts_good")->Fill("blah",
         0.0); //one MUST keep this for correct merging!
   hist("hCuts_bad")->Fill("blah",
         0.0); //one MUST keep this for correct merging!

   //~ hist("hCut")->Fill("blah",0.0); //one MUST keep this for correct merging!
   //~ hist("hCut_good")->Fill("blah",0.0); //one MUST keep this for correct merging!
   //~ hist("hCut_bad")->Fill("blah",0.0); //one MUST keep this for correct merging!

   for (int decay = 200; decay < 206; ++decay) {
      hist(TString::Format("hGoodWrongTags_%d", decay))->Fill("blah", 0.0);
      hist(TString::Format("hBadTags_%d", decay))->Fill("blah", 0.0);

   }




   ntp_features = new TNtuple("ntp_features","Features ntuple",
         "good:ep:pmiss:muclayers:etotal:mmiss_china:mmiss");


   VecObj features_ntuples(1, ntp_features);
   selector->RegInDir(features_ntuples,"Semileptonic_Features");

   char angular_vars[] =
         "q2:mkpi:theta_v:theta_l:chi:tag_passed:passed:d_theta:nu_theta:tag_mode:charm";
   ntp_angular_selected = new TNtuple("ntp_angular_selected",
         "Angular data of the survived itemsntuple",angular_vars);
   ntp_angular_mc = new TNtuple("ntp_angular_mc",
         "Angular data of the MC truth information",angular_vars);


   char momenta_vars[] =
         "l_px:l_py:l_pz:l_pe:nu_px:nu_py:nu_pz:nu_pe:k_px:k_py:k_pz:k_pe:pi_px:pi_py:pi_pz:pi_pe:tag_passed:passed:is_good";

   ntp_momenta_mc = new TNtuple("ntp_momenta_mc",
         "Momenta data of the MC truth information",momenta_vars);

   ntp_angular_selected_wmc = new TNtuple("ntp_angular_selected_wmc", "",
         "q2:mkpi:theta_v:theta_l:chi:mc_q2:mc_mkpi:mc_theta_v:mc_theta_l:mc_chi:tag_passed:good_mc:tag_mode");



   VecObj result_ntuples;
   result_ntuples.push_back(ntp_angular_mc);
   result_ntuples.push_back(ntp_angular_selected);
   result_ntuples.push_back(ntp_angular_selected_wmc);
   result_ntuples.push_back(ntp_momenta_mc);
   selector->RegInDir(result_ntuples, "Semileptonic_Tracks");


   ntp_electrons = new TNtuple("ntp_electrons", "ntp_electrons",
         "e_emc:p:cost:charge:evt_tot_emc:evt_dphi:pid_prob_e:pid_prob_sum:passed"  );
   VecObj e_ntuples;
   e_ntuples.push_back(ntp_electrons);
   selector->RegInDir(e_ntuples, "Semileptonic_Electrons");



}


// <is_mc, mode> => <de_low, de_high>
static map<pair<bool,int>, pair<double, double> > tag_de_cuts {
   // mc
   { {true, 200}, {-0.022, 0.020} },
   { {true, 201}, {-0.060, 0.033} },
   { {true, 202}, {-0.019, 0.020} },
   { {true, 203}, {-0.067, 0.036} },
   { {true, 204}, {-0.023, 0.021} },
   { {true, 205}, {-0.017, 0.015} },
   //data
   { {false, 200}, {-0.027,0.025} },
   { {false, 201}, {-0.063,0.036} },
   { {false, 202}, {-0.024,0.024} },
   { {false, 203}, {-0.071, 0.040} },
   { {false, 204}, {-0.035, 0.032} },
   { {false, 205}, {-0.021, 0.018} },

};


Hep3Vector getVertexOrigin(int runNo) {
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



//~ #include <TObjectTable.h> // debug memeory leak

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
bool SemileptonicEvent( ReadDst* selector,
      TEvtHeader* m_TEvtHeader,
      TDstEvent* m_TDstEvent,
      TEvtRecObject* m_TEvtRecObject,
      TMcEvent* m_TMcEvent,
      TTrigEvent* m_TTrigEvent,
      TDigiEvent* m_TDigiEvent,
      THltEvent* m_THltEvent)
//-----------------------------------------------------------------------------
{
   //~ static int counter = 0;
   //~ ++counter;
   //~ if (counter % 500 == 0) {
   //~ gObjectTable->Print();
   //~ }

   int event = m_TEvtHeader->getEventId();
   int run = m_TEvtHeader->getRunId();

   // evtRecTrkCol
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   // some initializtion
   ParticleID* pid = ParticleID::instance();
#if (BOSS_VER < 700)
   pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
   pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif

   // We have to initialize DatabaseSvc
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   DatabaseSvc* dbs = DatabaseSvc::instance();
   if( (dbs->GetDBFilePath()).empty() ) {
      dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
   }

   Hep3Vector xorigin = getVertexOrigin(run);





   char tag_str[12];
   m_EventTagSvc->setMcEvent(m_TMcEvent);




   unsigned int eventTagDecay =  (m_EventTagSvc->getEventTag()  &
               0xFFFFFF00) >> 8;
   //~ if (!hasDecay(eventTagDecay, 80, 71 )) return 0;

   formatTag(tag_str, eventTagDecay, m_TMcEvent);


   //~ cout <<"ChainCode: " << hex << m_EventTagSvc->getChainCode(getParticleById(m_TMcEvent, 30443)) << dec <<endl;



   bool preselection_passed = 0;


   bool good_oldalgo = false;

   good_oldalgo |= hasDecay(eventTagDecay, 80, 18 ); //  K-  pi+  e+ nu_e
   good_oldalgo |= hasDecay(eventTagDecay, 80, 8)
         && (m_EventTagSvc->getDecayCode(getParticleById(m_TMcEvent,
                           313)) == 1); // K+  pi-
   good_oldalgo |= hasDecay(eventTagDecay, 80, 12)
         && (m_EventTagSvc->getDecayCode(getParticleById(m_TMcEvent,
                           315)) == 1); // D+ => anti-K_2*0 e+  nu_e   , K_2*0 => K+  pi-

   if (!m_TEvtRecObject) {
      // rtraw file, nothing to do...
      for (int charm = -1;  charm <= 1; charm+=2) {
         for (int decayMode = 200; decayMode < 206; ++decayMode) {
            fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, false, false,
                  decayMode);
         }
      }

      //~ fillMcMomenta(m_TMcEvent, -1, ntp_momenta_mc, false, false);
      //~ fillMcMomenta(m_TMcEvent, 1, ntp_momenta_mc, false, false);

      return false;
   }

   const TEvtRecEvent*  m_evtRecEvent =
         m_TEvtRecObject->getEvtRecEvent();
   const TObjArray*  m_evtRecDTagCol=
         m_TEvtRecObject->getEvtRecDTagCol();


   int real_event_survived = 0;
   for (int charm = -1;  charm <= 1; charm+=2) {
      int charm_survived = 0;

      assert((charm == -1) ||(charm == 1));
      bool good = is_good_k_pi_e_nu(m_TMcEvent, charm);

      bool has_good_kpi_mass = is_good_k_pi_e_nu(m_TMcEvent, charm, 0.8,
                  1.0);
      cout << "has_good_kpi_mass: "<< has_good_kpi_mass << endl;
      hist("hGoodMass_vi")->Fill(has_good_kpi_mass);




      // calculate true K-pi pair mass
      double mKPi_true = -1;
      if (good) {
         mKPi_true = getTrueKPiMass(m_TMcEvent, charm);
         hist("hDp_mKPiTruthVeryInitial")->Fill(mKPi_true);
      }

      //~ if (!good) {
      //~ hist("hVeryInitialBadTags")->Fill(tag_str,1.0);
      //~ cout << "is_good: " << good << " for charm " << charm << endl;
      //~ selector->PrintMcDecayTree();
      //~ }

      const TMcParticle* signal_d_meson = get_d_meson(m_TMcEvent, charm);
      bool good_by_final_state = (get_final_state(m_TMcEvent,
                        signal_d_meson) == to_pdg_set(vector<string>({"K+","pi+", "e-", "nu_e"})));

      if (good != good_by_final_state) {
         cout << "good_by_final_state(" << good_by_final_state << ") != good ("
               << good << ") !" <<endl;
         cout << "signal final state";
         print_pdg_set(get_final_state(m_TMcEvent, signal_d_meson));
         selector->PrintMcDecayTree();
      }


      int signal_decay_code = -1;
      if (signal_d_meson) {
         signal_decay_code = m_EventTagSvc->getDecayCode(
                     const_cast<TMcParticle*>(signal_d_meson));
      }





      for (int decayMode = 200; decayMode < 206; ++decayMode) {
         bool event_passed = 0;
         int tag_passed = 0;


         TEvtRecDTag* evtRecDTag = selectDTag(m_evtRecDTagCol, evtRecTrkCol,
                     -1 * charm, decayMode);


         bool is_semimuonic =  (get_final_state(m_TMcEvent,
                           signal_d_meson) == to_pdg_set(vector<string>({"K+","pi+", "mu-", "nu_mu"})));



         TCutsManager cutsManager(!m_TMcEvent->getMcParticleCol()->IsEmpty(),
               good, evtRecDTag, charm);


         cutsManager.initCut("hasDTagCol");
         //~ cout << "m_evtRecDTagCol->GetEntriesFast(): " << m_evtRecDTagCol->GetEntriesFast() << endl;
         //~ cout << "m_evtRecDTagCol->IsEmpty(): " << m_evtRecDTagCol->IsEmpty() << endl;

         const TMcParticle* tag_d_meson = get_d_meson(m_TMcEvent, -1 * charm);
         auto tag_final_state = get_final_state(m_TMcEvent,  tag_d_meson);


         //~ if (evtRecDTag) {
         //~ cout << "Found dtag #" << i
         //~ << " : mode=" << evtRecDTag->decayMode()
         //~ << " , type=" << evtRecDTag->type()
         //~ << " , charm=" << evtRecDTag->charm()
         //~ << " , charge=" << evtRecDTag->charge()
         //~ << endl;
         //~ cout << "final state: ";
         //~ print_pdg_set(tag_final_state);
         //~ } else {
         //~ cout << "null evtRecDTag " << endl;
         //~ }

         int mc_tag_mode = get_mc_tag_mode(m_TMcEvent, -1 * charm);
         bool mc_is_good_tag =  ( (mc_tag_mode >= 200)
                     && (mc_tag_mode <= 205));

         //~ if (mc_is_good_tag && !evtRecDTag) {
         //~ cout << "good tag not found: " << mc_tag_mode << endl;
         //~ cout << "final state: ";
         //~ print_pdg_set(tag_final_state);
         //~ selector->PrintMcDecayTree();
         //~ }



         cutsManager.smartFillCutHist("hViSignalDecays",signal_decay_code);


         cutsManager.initCut("hasTags");
         if (!evtRecDTag) {
            cutsManager.fatalByCut( "hasTags");
            fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, tag_passed, false,
                  decayMode);
            //~ fillMcMomenta(m_TMcEvent, charm, ntp_momenta_mc, tag_passed, false);
            continue;
         }


         //~ fillDtagDecayHists(evtRecDTag, 200, eventTagDecay, 80, 71, -0.02, 0.02, 1.868, 1.872);
         //~ fillDtagDecayHists(evtRecDTag, 201, eventTagDecay, 80, 87, -0.02, 0.02, 1.868, 1.872);



         bool good_tag_mode =  (evtRecDTag->decayMode() == mc_tag_mode);

         // =============================================================
         // пытаемся определить, попадает ли таг в 'сигнальную' часть mBC распределения,
         // сравнивая измеренный и MC 4-х векторы D
         HepLorentzVector taggedD4Momentum(evtRecDTag->px(), evtRecDTag->py(),
               evtRecDTag->pz(), evtRecDTag->pe());

         bool good_tag_flag = false;
         double d4MomentumDist = 1E10;

         if (tag_d_meson) {
            HepLorentzVector MCD4Momentum  = getMcParticleP4(tag_d_meson);

            // считаем разницу как длину разности пространственных составляющих
            //~ auto d4MomentumDiff = taggedD4Momentum. - MCD4Momentum;
            //~ double d4MomentumDist = TMath::Sqrt(d4MomentumDiff.vect().mag2() + sqr(d4MomentumDiff.e()));

            //~ d4MomentumDist = (taggedD4Momentum.vect() - MCD4Momentum.vect()).mag();
            d4MomentumDist = taggedD4Momentum.vect().angle(MCD4Momentum.vect());
            // считаем таг хорошим:
            good_tag_flag = ((d4MomentumDist < 0.20) && good_tag_mode) ;
         }

         cout << "good_tag_flag: " << good_tag_flag << endl;


         // =============================================================
         HepLorentzVector mc_pion_4momentum;
         HepLorentzVector mc_muon_4momentum;
         const TMcParticle* mc_muon = 0;

         if (is_semimuonic && signal_d_meson) {
            //find pion
            const TMcParticle* pi = 0;
            int decay = m_EventTagSvc->getDecayCode(const_cast<TMcParticle*>
                        (signal_d_meson));

            if  (decay == 42) { //     K-  pi+  mu+ nu_e
               pi = getDaughterById(m_TMcEvent, signal_d_meson, 211);
            } else if (decay == 32 ) { //   anti-K*0  mu+  nu_mu
               const TMcParticle* k_star = getDaughterById(m_TMcEvent,
                           signal_d_meson, 313);
               pi = getDaughterById(m_TMcEvent, k_star, 211);
            } else if (decay == 36) {  // 36   anti-K_2*0 mu+  nu_mu
               const TMcParticle* k_star_2 = getDaughterById(m_TMcEvent,
                           signal_d_meson, 315);
               pi = getDaughterById(m_TMcEvent, k_star_2, 211);
            }


            if (!pi) {
               is_semimuonic = false;
            } else {
               mc_pion_4momentum = getMcParticleP4(pi);
               mc_muon = getDaughterById(m_TMcEvent, signal_d_meson, 13); // mu
               mc_muon_4momentum = getMcParticleP4(mc_muon); //mu

               cutsManager.smartFillCutHist( "hDp_vi_MC_pPi",
                     mc_pion_4momentum.vect().mag());
               cutsManager.smartFillCutHist( "hDp_vi_MC_pMu",
                     mc_muon_4momentum.vect().mag());

               //~ preselection_passed = 1;
            }
         }


         //  ===========================================================




         if (mc_is_good_tag && !good_tag_mode) {
            cout << "wrong tag mode: mc " << mc_tag_mode  << "!=" <<
                  evtRecDTag->decayMode() << endl;
            cout << "d4MomentumDist: " << d4MomentumDist << endl;
            cout << "final state: ";
            print_pdg_set(tag_final_state);
            selector->PrintMcDecayTree();
         }




         if (good_tag_mode ) {
            cutsManager.smartFillCutHist("hDp_goodTagMode",
                  evtRecDTag->decayMode());
         } else {
            cutsManager.smartFillCutHist("hDp_goodTagMode", -1);
         }




         cutsManager.smartFillCutHist("hDTag_DecayMode",
               evtRecDTag->decayMode());

         if ( (evtRecDTag->decayMode()  < 200)
               || (evtRecDTag->decayMode()  > 205) ) {
            //~ if ( evtRecDTag->decayMode()  != 200){
            cutsManager.rejectByCut("hDTag_DecayMode");
         }




         cutsManager.smartFillCutHist("hDp_goodTag", mc_is_good_tag ); //201



         if (cutsManager.survived())  {
            hist("hGoodMass_hastag")->Fill(has_good_kpi_mass);

            //~ preselection_passed = 1;

         }

         cutsManager.smartFillCutHist("hDp_dE", evtRecDTag->deltaE());

         bool is_mc = !m_TMcEvent->getMcParticleCol()->IsEmpty();



         auto dE_cuts = tag_de_cuts[make_pair(is_mc, evtRecDTag->decayMode())];
         double dE_low = dE_cuts.first;
         double dE_high = dE_cuts.second;



         cutsManager.rejectByCutIfOutsideWindow(TString::Format("hDp_dE_%d",
                     evtRecDTag->decayMode()).Data(), evtRecDTag->deltaE(), dE_low,
               dE_high);



         // vvv enable for D-tag count survey
         //~ if (cutsManager.survived())  preselection_passed = 1;



         // смотрим на распределение mBC для тегов, которые мы считаем "плохими". Ожидаем видеть только подложку, без пиков
         if (cutsManager.survived()) {
            hist("hDp_4MDist")->Fill( d4MomentumDist);


            if (mc_is_good_tag) {
               hist("hDp_mBC_mctag_good")->Fill( evtRecDTag->mBC());
               hist("hDp_4MDist_mctag_good")->Fill( d4MomentumDist);
            } else {
               hist("hDp_mBC_mctag_bad")->Fill( evtRecDTag->mBC());
               hist("hDp_4MDist_mctag_bad")->Fill( d4MomentumDist);
            }

            if (good_tag_mode) {
               hist("hDp_mBC_tagmode_good")->Fill( evtRecDTag->mBC());
               hist("hDp_4MDist_tagmode_good")->Fill( d4MomentumDist);
            } else {
               hist("hDp_mBC_tagmode_bad")->Fill( evtRecDTag->mBC());
               hist("hDp_4MDist_tagmode_bad")->Fill( d4MomentumDist);
            }

            if (good_tag_flag) {
               hist("hDp_mBC_tagflag_good")->Fill( evtRecDTag->mBC());
            } else {
               hist("hDp_mBC_tagflag_bad")->Fill( evtRecDTag->mBC());
            }


         }


         // ========= mBC cut =============
         cutsManager.smartFillCutHist("hDp_mBC", evtRecDTag->mBC());
         cutsManager.rejectByCutIfOutsideWindow("hDp_mBC", evtRecDTag->mBC(),
               1.863, 1.877);




         cutsManager.smartFillCutHist("hDp_goodTagMode2",
               good_tag_mode ? evtRecDTag->decayMode() : -1);
         cutsManager.smartFillCutHist("hDp_goodTagMode3",
               good_tag_flag ? evtRecDTag->decayMode() : -1);



         cutsManager.smartFillCutHist("hDp_goodMcTagMode",
               mc_is_good_tag ? mc_tag_mode : -1);



         cutsManager.smartFillCutHist("hDp_goodTag2",  mc_is_good_tag);
         cutsManager.smartFillCutHist("hDp_goodTagFlag",  good_tag_flag);


         cutsManager.smartFillCutHist("hSignalDecays",signal_decay_code);
         cutsManager.smartFillCutHist("hGoodTagFlagSignalDecays",
               good_tag_flag ? signal_decay_code : -1);

         if (cutsManager.survived()) {
            if (!mc_is_good_tag) {
               cout << "bad tag survived as " << evtRecDTag->decayMode() << endl;
               cout << "final state: ";
               print_pdg_set(tag_final_state);
               selector->PrintMcDecayTree();
            }


            hist("hGoodMass_d")->Fill(has_good_kpi_mass);


            if (good_tag_mode) {
               hist("hGoodMass_mcgoodtag")->Fill(has_good_kpi_mass);
            }

            if (mc_is_good_tag) {
               hist("hGoodMass_mcgoodtag2")->Fill(has_good_kpi_mass);
            }

            if (good_tag_flag) {
               hist("hGoodMass_good_tag_flag")->Fill(has_good_kpi_mass);
            }


            // vvv enable for deep cuts research
            //~ preselection_passed = 1;


            //~ hist("hInitialAllTags")->Fill(tag_str,1.0);

            if (good) {
               hist("hDp_mKPiTruthInitial")->Fill(mKPi_true);

               if (good_tag_flag) {
                  hist("hDp_mKPiTruthInitialGoodTag")->Fill(mKPi_true);
               }


               if (hasDecay(eventTagDecay, 80, 8)) {
                  hist("hDp_mKPiTruthKStar")->Fill(mKPi_true);
               }

            }
            tag_passed = good_tag_flag ? evtRecDTag->decayMode() : -1;
            cout << "tag_passed: " << tag_passed  << endl;
         } else { //save time
            cutsManager.processCuts(true);
            fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, false, false,
                  decayMode);
            //~ fillMcMomenta(m_TMcEvent, charm, ntp_momenta_mc, false, false);
            continue;
         }






         const vector<Int_t>& otrk = evtRecDTag->otherTracks();
         const vector<Int_t>& oshw = evtRecDTag->otherShowers();

         vector<Int_t> signal_tracks;

         for (unsigned track_i = 0; track_i < otrk.size(); ++track_i) {
            DstEvtRecTracks* tracks =  (DstEvtRecTracks*) evtRecTrkCol->At(
                        otrk[track_i]);
            if (!isGoodChargedTrack(tracks, xorigin)) {
               continue;
            }

            signal_tracks.push_back(otrk[track_i]);
         }




         // simple selection: other tracks == 3
         // one of the other tracks is electron
         // two others should be pion and kaons and so should be contained in Dtag list

         cutsManager.smartFillCutHist("hDp_Ntracks", signal_tracks.size());

         if (signal_tracks.size() != 3) {
            if (good ) {
               cout << "hDp_Ntracks fail" << signal_tracks.size() << endl;
            }
            cutsManager.fatalByCut( "hDp_Ntracks");
            fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, tag_passed, false,
                  decayMode);
            //~ fillMcMomenta(m_TMcEvent, charm, ntp_momenta_mc, tag_passed, false);
            continue;
         }




         int total_charge = evtRecDTag->charge();
         for (int tmpi = 0; tmpi < 3; ++tmpi) {
            DstEvtRecTracks* tracks =  (DstEvtRecTracks*) evtRecTrkCol->At(
                        signal_tracks[tmpi]);
            total_charge  += tracks->mdcTrack()->charge();
         }

         cutsManager.smartFillCutHist("hDp_TotCharge", total_charge);

         if (cutsManager.survived()) {
            cout << "total charge: " << total_charge << endl;
         }



         DstEvtRecTracks* kaon_tracks = 0;
         DstEvtRecTracks* pion_tracks = 0;
         DstEvtRecTracks* e_tracks = 0;
         bool found_combination = false;
         int n_3tracks_combinations = 0;


         for (unsigned pion_id = 0; pion_id < signal_tracks.size() ;
               ++pion_id) {
            for (unsigned kaon_id = 0; kaon_id < signal_tracks.size() ;
                  ++kaon_id) {
               if (pion_id != kaon_id)
                  for (unsigned electron_id = 0; electron_id < signal_tracks.size() ;
                        ++electron_id)
                     if ( (electron_id != kaon_id) && (electron_id != pion_id)) {

                        DstEvtRecTracks* _kaon_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                                    signal_tracks[kaon_id]);
                        DstEvtRecTracks* _pion_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                                    signal_tracks[pion_id]);
                        DstEvtRecTracks* _e_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                                    signal_tracks[electron_id]);

                        //~ if ( ! _kaon_tracks->isMdcTrackValid())
                        //~ continue;
                        //~ if ( ! _pion_tracks->isMdcTrackValid())
                        //~ continue;
                        //~ if ( ! _e_tracks->isMdcTrackValid())
                        //~ continue;



                        //FIXME: проверить что в pions нормальные критерии на пионы

                        if (   is_pion_dtag( signal_tracks[pion_id], evtRecDTag) &&
                              is_kaon_dtag( signal_tracks[kaon_id], evtRecDTag) &&
                              is_electron(  signal_tracks[electron_id], pid, evtRecTrkCol, run) &&
                              ( _e_tracks   ->mdcTrack()->charge() == -evtRecDTag->charge() ) &&
                              ( _kaon_tracks->mdcTrack()->charge() ==
                                    -_pion_tracks->mdcTrack()->charge()) &&
                              ( _e_tracks->mdcTrack()->charge() ==
                                    _pion_tracks->mdcTrack()->charge())

                        )

                        {
                           //~ found_combination = true;
                           n_3tracks_combinations += 1;

                           e_tracks = _e_tracks;
                           pion_tracks = _pion_tracks;
                           kaon_tracks = _kaon_tracks;


                           //~ break;
                        }
                     }
            }
         }

         if (n_3tracks_combinations > 0) {
            found_combination = true;
         }


         cutsManager.initCut("3tracks");

         if (!found_combination) {
            cutsManager.fatalByCut( "3tracks");
            fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, tag_passed, false,
                  decayMode);
            //~ fillMcMomenta(m_TMcEvent, charm, ntp_momenta_mc, tag_passed, false);
            cout << "preselection_passed: " << preselection_passed << endl;
            continue;
         }


         // for cut study 20121210
         if (cutsManager.survived()) {
            preselection_passed = 1;
         }
         //~ preselection_passed

         cout << "========================================================================================"
               << endl;
         cout  << "Found combination" << endl;
         //~ selector->PrintEvtRec();
         selector->PrintMcDecayTree();
         cout << "is_good: " << good << endl;
         cout << "decayMode: "<< evtRecDTag->decayMode();



         cutsManager.smartFillCutHist( "hDp_n3tr_comb",
               n_3tracks_combinations);
         cutsManager.rejectByCutIfHigher( "hDp_n3tr_comb",
               n_3tracks_combinations, 1.01);






         e_tracks->mdcKalTrack()->setPidType(BMdcKalTrack::electron);
         kaon_tracks->mdcKalTrack()->setPidType(BMdcKalTrack::kaon);
         pion_tracks->mdcKalTrack()->setPidType(BMdcKalTrack::pion);


         HepLorentzVector kaon4Momentum = kaon_tracks->mdcKalTrack()->p4(
                     MASS_KAON);
         HepLorentzVector pion4Momentum = pion_tracks->mdcKalTrack()->p4(
                     MASS_PION);
         double mKPi = (pion4Momentum + kaon4Momentum).m();
         if (good) {
            double diff  = (mKPi_true - mKPi) * 100.0 / mKPi;
            cout << tag_str << " mKPi_true: " << mKPi_true << ", mKPi: " << mKPi
                  << ", diff: " <<  diff <<  "%" << endl;
            if (diff > 10) {
               cout << "big diff" << endl;
            }


         }






         // check whether this tracks found as another tagged D
         bool goodOtherTag = false;

         TIter otherEvtRecDTagIter(m_evtRecDTagCol);
         while (TEvtRecDTag* otherEvtRecDTag = (TEvtRecDTag*)
                     otherEvtRecDTagIter.Next()) {
            // all tracks of otherEvtRecDTag should present in evtRecDTag's otherTracks
            // all showers of otherEvtRecDTag should present in evtRecDTag's otherShowers
            const vector<Int_t>& otherSideTrk = otherEvtRecDTag->tracks();
            const vector<Int_t>& otherSideShw = otherEvtRecDTag->showers();

            bool goodIntersection = 1;

            for (unsigned i = 0; i < otherSideTrk.size(); ++i) {
               if (find(otrk.begin(), otrk.end(), otherSideTrk[i]) == otrk.end()) {
                  goodIntersection = 0;
                  break;
               }
            }
            for (unsigned i = 0; i < otherSideShw.size(); ++i) {
               if (find(oshw.begin(), oshw.end(), otherSideShw[i]) == oshw.end()) {
                  goodIntersection = 0;
                  break;
               }
            }
            if (goodIntersection) {
               // check whether this tag has good deltaE and mBC
               goodOtherTag |= selectDtagDecay(otherEvtRecDTag, 200, -0.015, 0.015,
                           1.866, 1.875);
            }
         }

         cutsManager.smartFillCutHist("hDp_goodOtherSide", (int)goodOtherTag);
         if (goodOtherTag) {
            //~ cout << "goodOtherTag! skipping" << endl;
            //~ cutsManager.rejectByCut("hDp_goodOtherSide");
         }


         // FSR recovery tests
         cout << "e_tracks->mdcKalTrack()->p3(): " <<
               e_tracks->mdcKalTrack()->p3() << endl;
         cout << "e kal th: " << e_tracks->mdcKalTrack()->theta() << " phi: "
               << e_tracks->mdcKalTrack()->phi() << endl;
         cout << "e kal p3 th: " << e_tracks->mdcKalTrack()->p3().theta() <<
               " phi: " << e_tracks->mdcKalTrack()->p3().phi() << endl;
         cout << "e mdc th: " << e_tracks->mdcTrack()->theta() << " phi: " <<
               e_tracks->mdcTrack()->phi() << endl;

         if (e_tracks->isEmcShowerValid()) {
            cout << "e shower th: " << e_tracks->emcShower()->theta() << " phi: "
                  << e_tracks->emcShower()->phi() << endl;
            if (e_tracks->isExtTrackValid()) {
               cout << "e emc ext th: " <<
                     e_tracks->extTrack()->emcPosition().theta() << " phi: " <<
                     e_tracks->extTrack()->emcPosition().phi() << endl;
            }
         }

         cout << "e recdsttracks track_id " <<  e_tracks->trackId() << endl;


         set<Int_t> fsr_showers;
         HepLorentzVector totalFsr4Momentum;

         // loop over EMC showers
         for (unsigned shower_id = 0; shower_id < oshw.size(); ++shower_id) {
            DstEvtRecTracks* shower_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                        oshw[shower_id]);
            RecEmcShower* emcTrk = shower_tracks->emcShower();

            if (!isGoodShower(emcTrk, evtRecTrkCol,
                        m_evtRecEvent, /*exclude_track=*/e_tracks->trackId())) {
               continue;
            }

            Hep3Vector showerVector = getShowerVector(emcTrk);

            double initial_angle =  showerVector.angle(
                        e_tracks->mdcKalTrack()->p3());
            cout << "shower id: " << oshw[shower_id] <<
                  " angle to initial electron: " <<  initial_angle / TMath::DegToRad()
                  << " degrees" << endl;
            if (initial_angle < TMath::DegToRad() * 5.0 ) {
               // found fsr photon
               cout << "found FSR photon: " << oshw[shower_id] << ", energy: " <<
                     emcTrk->energy() * 1E3<< " MeV"<< endl;
               cout << "FSR shower th: " << emcTrk->theta() << " phi: " <<
                     emcTrk->phi() << endl;

               Hep3Vector fsrMomentum = showerVector.unit();
               fsrMomentum.setMag(emcTrk->energy());

               HepLorentzVector fsr4Momentum;
               fsr4Momentum.setVectM(fsrMomentum, 0);

               cout << "fsr4momentum: " << fsr4Momentum << endl;
               totalFsr4Momentum += fsr4Momentum;

               fsr_showers.insert(oshw[shower_id]);
            }

         }

         cutsManager.smartFillCutHist( "hDp_nFSR", fsr_showers.size());
         cutsManager.smartFillCutHist( "hDp_FSRE", totalFsr4Momentum.e());
         cout << "totalFsr4Momentum: " << totalFsr4Momentum << endl;

         HepLorentzVector electron4Momentum;
         electron4Momentum.setVectM(
               e_tracks->mdcKalTrack()->p3() + totalFsr4Momentum.vect(),
               MASS_ELECTRON);

         // vvv outdated
         //~ electron4Momentum.setVectM(
         //~ e_tracks->mdcKalTrack()->p3(),
         //~ MASS_ELECTRON);
         //~ electron4Momentum += totalFsr4Momentum;




         // end of FSR tests

         // epEMC

         double e_shower_energy = -1;
         double ep_ratio = -1;
         if (e_tracks->isEmcShowerValid()) {
            e_shower_energy = e_tracks->emcShower()->energy();
            ep_ratio = (e_shower_energy / e_tracks->mdcKalTrack()->p());
         }



         bool with_electron = good ||  hasDecay(eventTagDecay, 80, 9 )
               ||  hasDecay(eventTagDecay, 80, 10 )
               ||  hasDecay(eventTagDecay, 80, 11)
               ||  hasDecay(eventTagDecay, 80, 12)
               ||  hasDecay(eventTagDecay, 80, 13)
               ||  hasDecay(eventTagDecay, 80, 14)
               ||  hasDecay(eventTagDecay, 80, 15)
               ||  hasDecay(eventTagDecay, 80, 16)
               ||  hasDecay(eventTagDecay, 80, 17)
               ||  hasDecay(eventTagDecay, 80, 18)
               ||  hasDecay(eventTagDecay, 80, 19)
               ||  hasDecay(eventTagDecay, 80, 20) ;


         double e_cost = TMath::Cos(e_tracks->mdcKalTrack()->theta());
         double e_p = e_tracks->mdcKalTrack()-> p();
         cutsManager.smartFillCutHist( "hDp_pxyE",
               e_tracks->mdcKalTrack()-> pxy());
         cutsManager.smartFillCutHist( "hDp_pE", e_p, !with_electron);
         cutsManager.smartFillCutHist( "hDp_eE", e_shower_energy);
         cutsManager.smartFillCutHist( "hDp_costE", e_cost);

         if (good) {
            ((TH2D*) hist("s_hDp_e_p_vs_cost_good"))->Fill(e_p, e_cost);

            if (e_shower_energy > 0) {
               ((TH2D*) hist("s_hDp_e_p_vs_cost_we_good"))->Fill(e_p, e_cost);
            }



         };




         cutsManager.smartFillCutHist( "hDp_epEMC", ep_ratio, is_semimuonic);

         //~ cutsManager.rejectByCutIfOutsideWindow("hDp_epEMC", ep_ratio, 0.8, 1.1);
         cutsManager.rejectByCutIfLower("hDp_epEMC", ep_ratio, 0.8);


         // study of different epEMC distributions for different signs

         //~ if (is_semimuonic && signal_d_meson) {
         //~ DstEvtRecTracks * rec_mu_tracks = get_closest_track(evtRecTrkCol, otrk, mc_muon, 0.1);
         //~ cout << "rec_mu_tracks: " << rec_mu_tracks << endl;
         //~ cout << "e_tracks: " << e_tracks << endl;
         //~ cout << "pion_tracks: " << pion_tracks << endl;
         //~ cout << "kaon_tracks: " << kaon_tracks << endl;



         //~ cutsManager.smartFillCutHist( "hDp_EMCeE", e_shower_energy , is_semimuonic);
         //~ cutsManager.smartFillCutHist( "hDp_pE_sm_kal", e_tracks->mdcKalTrack()->p() , is_semimuonic);
         //~ cutsManager.smartFillCutHist( "hDp_pE_sm", e_tracks->mdcTrack()->p() , is_semimuonic);

         // pi-pi, mu-e
         //~ double pipi_mue_dist = sqr(mc_pion_4momentum.angle(pion4Momentum)) + sqr(mc_muon_4momentum.angle(electron4Momentum));

         // pi-e, mu-pi
         //~ double pie_mupi_dist = sqr(mc_pion_4momentum.angle(electron4Momentum)) + sqr(mc_muon_4momentum.angle(pion4Momentum));

         //~ double fake_e_momentum;
         //~ cout << "pipi_mue_dist: " << pipi_mue_dist << " pie_mupi_dist: " << pie_mupi_dist << endl;
         //~ if  (pipi_mue_dist < pie_mupi_dist) { //mu as e, pi as pi
         //~ fake_e_momentum = mc_muon_4momentum.vect().mag();
         //~
         //~ cutsManager.smartFillCutHist( "hDp_epEMC_sm_mue", ep_ratio);
         //~ cutsManager.smartFillCutHist( "hDp_EMCeE_sm_mue", e_shower_energy);
         //~ cutsManager.smartFillCutHist( "hDp_pE_sm_mue", e_tracks->mdcKalTrack()->p());
         //~
         //~ cutsManager.smartFillCutHist( "hDp_sm_pRealMu", e_tracks->mdcKalTrack()->p());
         //~ cutsManager.smartFillCutHist( "hDp_sm_pRealPi", pion_tracks->mdcKalTrack()->p());
         //~
         //~ if (e_tracks->isEmcShowerValid()) {
         //~ cutsManager.smartFillCutHist( "hDp_thE_sm_mue", TMath::Cos(e_tracks->emcShower()->theta()));
         //~ cutsManager.smartFillCutHist( "hDp_phiE_sm_mue", TMath::RadToDeg()*e_tracks->emcShower()->phi());
         //~ }
         //~
         //~ cutsManager.smartFillCutHist( "hDp_MC_pMu_mue", mc_muon_4momentum.vect().mag());
         //~
         //~
         //~ } else { // mu as pi, pi as e
         //~ fake_e_momentum = mc_pion_4momentum.vect().mag();
         //~ cutsManager.smartFillCutHist( "hDp_epEMC_sm_mupi", ep_ratio);
         //~ cutsManager.smartFillCutHist( "hDp_EMCeE_sm_mupi", e_shower_energy);
         //~ cutsManager.smartFillCutHist( "hDp_pE_sm_mupi", e_tracks->mdcKalTrack()->p());
         //~
         //~ cutsManager.smartFillCutHist( "hDp_sm_pRealMu", pion_tracks->mdcKalTrack()->p());
         //~ cutsManager.smartFillCutHist( "hDp_sm_pRealPi", e_tracks->mdcKalTrack()->p());
         //~
         //~ if (e_tracks->isEmcShowerValid()) {
         //~ cutsManager.smartFillCutHist( "hDp_thE_sm_mupi", TMath::Cos(e_tracks->emcShower()->theta()));
         //~ cutsManager.smartFillCutHist( "hDp_phiE_sm_mupi", TMath::RadToDeg()*e_tracks->emcShower()->phi());
         //~ }
         //~ }

         //~ cutsManager.smartFillCutHist( "hDp_MC_pFakeE", fake_e_momentum);


         //~ cutsManager.smartFillCutHist( "hDp_MC_pPi", mc_pion_4momentum.vect().mag());
         //~ cutsManager.smartFillCutHist( "hDp_MC_pMu", mc_muon_4momentum.vect().mag());

         //~ cutsManager.smartFillCutHist( "hDp_MC_eMu", mc_muon_4momentum.e());

         //~ if ((mc_muon_4momentum.vect().mag() -e_tracks->mdcKalTrack()->p()) / e_tracks->mdcKalTrack()->p() > 0.3) {
         //~ cout << "mc_muon_4momentum.vect().mag(): " << mc_muon_4momentum.vect().mag() << endl;
         //~ cout << "e_tracks->mdcKalTrack()->p(): " << e_tracks->mdcKalTrack()->p() << endl;
         //~ cout << "electron4Momentum: " << electron4Momentum << endl;
         //~ cout << "electron4Momentum (w/o fsr): " << electron4Momentum -totalFsr4Momentum << endl;
         //~ cout << "pion4Momentum: " << pion4Momentum << endl;
         //~
         //~ cout << "mc_muon_4momentum: " <<  mc_muon_4momentum << endl;
         //~
         //~ selector->PrintMcDecayTree();
         //~ }

         //~ }


         //~ cutsManager.smartFillCutHist(  "hDp_pxyPi", pion_tracks->mdcKalTrack()-> pxy());
         //~ cutsManager.rejectByCutIfLower("hDp_pxyPi", pion_tracks->mdcKalTrack()-> pxy(), 0.012);
         //~
         //~
         //~ cutsManager.smartFillCutHist( "hDp_pxyK", kaon_tracks->mdcKalTrack()-> pxy());
         //~ cutsManager.rejectByCutIfLower("hDp_pxyK", kaon_tracks->mdcKalTrack()-> pxy(), 0.012);
         //~
         //~
         //~ cutsManager.smartFillCutHist( "hDp_pxyE", e_tracks->mdcKalTrack()-> pxy());
         //~ cutsManager.rejectByCutIfLower("hDp_pxyE", e_tracks->mdcKalTrack()-> pxy(), 0.012);



         cutsManager.smartFillCutHist(  "hDp_pPi",
               pion_tracks->mdcKalTrack()-> p());
         cutsManager.smartFillCutHist(  "hDp_pK",
               kaon_tracks->mdcKalTrack()-> p());






         KinematicFit* kmfit = KinematicFit::instance();
         kmfit->init();
         Hep3Vector taggedDMomentum(evtRecDTag->px(), evtRecDTag->py(),
               evtRecDTag->pz());

         //~ HepLorentzVector taggedD4Momentum(evtRecDTag->px(), evtRecDTag->py(), evtRecDTag->pz(), evtRecDTag->pe());

         cout << "evtRecDTag->pe(): " << evtRecDTag->pe() << endl;
         cout << "taggedD4Momentum.mass" << (-taggedD4Momentum).m2() << endl;
         cout << "taggedD4Momentum" <<   taggedD4Momentum << endl;
         cout << "eventTagDecay: " << eventTagDecay << " text:   " << tag_str
               << endl;





         HepLorentzVector observed4Momentum =
               pion4Momentum +
               kaon4Momentum +
               electron4Momentum;


         Hep3Vector observedMomentum = observed4Momentum.vect();

         double beamE = evtRecDTag->beamE();

         Hep3Vector cmsBoost(-11E-3, 0, 0);

         Hep3Vector initialMomentum = -cmsBoost * 2 * beamE;
         HepLorentzVector ecms(initialMomentum, 2 * beamE);

         HepLorentzVector expected4Momentum;

         expected4Momentum.setVectM(initialMomentum -taggedDMomentum,
               D_PLUS_MASS);

         //~ HepLorentzVector expected4MomentumBC(initialMomentum -taggedDMomentumBC,
         //~ evtRecDTag->beamE())
         //~ ;



         cout << "missing mass: " <<  (ecms -  taggedD4Momentum -
                     observed4Momentum).m2() << endl;
         cout << "missing mass2: " <<  (expected4Momentum -
                     observed4Momentum).m2() << endl;
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





         HepLorentzVector expectedD4Momentum_cms;
         expectedD4Momentum_cms.setVectM(-taggedD4Momentum_cms.vect(),
               D_PLUS_MASS );


         HepLorentzVector expectedD4Momentum_cms2;
         HepLorentzVector taggedD4Momentum_cms2;
         taggedD4Momentum_cms2.setVectM(taggedD4Momentum.vect(), D_PLUS_MASS);
         taggedD4Momentum_cms2.boost(cmsBoost);
         expectedD4Momentum_cms2.setVectM(-taggedD4Momentum_cms2.vect(),
               D_PLUS_MASS );




         HepLorentzVector missing4MomentumBC_cms = expectedD4MomentumBC_cms -
               observed4Momentum_cms;
         HepLorentzVector missing4Momentum_cms = expectedD4Momentum_cms -
               observed4Momentum_cms;
         HepLorentzVector missing4Momentum_cms2 = expectedD4Momentum_cms2 -
               observed4Momentum_cms;


         cout << "expectedD4MomentumBC_cms: " << expectedD4MomentumBC_cms <<
               endl;
         cout << "observed4Momentum_cms: " << observed4Momentum_cms << endl;
         cout << "missing4MomentumBC_cms: " << missing4MomentumBC_cms << endl;
         cout << "missing4Momentum_cms: " << missing4Momentum_cms << endl;
         cout << "missing4Momentum_cms2: " << missing4Momentum_cms2 << endl;

         cout << "missing mass3: " <<  (missing4MomentumBC_cms).m2() << endl;

         double missingEnergy = evtRecDTag->beamE() - observed4Momentum.e();
         double missingEnergy_cms = evtRecDTag->beamE() -
               observed4Momentum_cms.e();


         //~ hist("hDp_mMiss_kpnue")->Fill( (expected4Momentum - observed4Momentum).m());


         cout << "missingEnergy: " << missingEnergy << endl;
         cout << "missingEnergy_cms: " << missingEnergy_cms << endl;
         cout << "missingEnergy_cms rel delta: " << (missingEnergy_cms
                     -missingEnergy)/missingEnergy << endl;

         cout << "observedMomentum: " << observedMomentum << endl;
         cout << "observed4Momentum:  " << observed4Momentum << endl;
         cout << "expected4Momentum:  " << expected4Momentum << endl;
         //~ Hep3Vector missingMomentum = (0 - taggedDMomentum - observedMomentum);
         Hep3Vector missingMomentum = (expected4Momentum.vect() -
                     observedMomentum);

         HepLorentzVector missing4Momentum = (expected4Momentum -
                     observed4Momentum);

         cout << "missing momentum: " <<missingMomentum << endl;
         cout << "missing momentum(2): " <<missing4Momentum.vect() << endl;
         assert(missing4Momentum.vect() == missingMomentum);

         cout << "kaon 4-momentum: " << kaon4Momentum << endl;
         cout << "pion 4-momentum: " << pion4Momentum << endl;
         cout << "electron 4-momentum: " <<electron4Momentum << endl;



         HepLorentzVector observed4Momentum_KaonAsPion =
               kaon_tracks->mdcKalTrack()->p4(MASS_PION) +
               pion4Momentum +
               electron4Momentum;

         smartFillHist("hDp_mMissInitial",
               (expected4Momentum - observed4Momentum).m2(), good);


         // ========== missing momentum  =================
         cutsManager.smartFillCutHist( "hDp_pMiss", missingMomentum.mag() );

         // mostly kinematic cut. consider deletion
         // to suppress background from K pi pi (with no nu)

         //~ cutsManager.rejectByCutIfLower("hDp_pMiss", missingMomentum.mag() , 0.1); // this suppresses sigificant part of  K*_20 decay channel!

         // ========== missing energy    =================
         cutsManager.smartFillCutHist( "hDp_EMiss", missingEnergy_cms);
         cutsManager.rejectByCutIfLower("hDp_EMiss", missingEnergy_cms,
               0.04);  // this suppresses sigificant part of  K*_20 decay channel!

         //~ cutsManager.rejectByCutIfLower("hDp_pMiss", missingMomentum.mag() , 0.00); // to preserve K*_20



         // ========== electron momentum cut =================




         //~ const std::vector<RecTofTrack* >& tof_tracks = e_tracks->tofTrack();
         //~ for (size_t tof_index = 0; tof_index < tof_tracks.size(); ++tof_index) {
         //~ RecTofTrack* tof_track = tof_tracks[tof_index];
         //~ if ((tof_track->status() & 0x00000004) >> 2) {
         //~ smartFillHist("hDp_tof", tof_track->tof(), good);
         //~ }
         //~ }





         //~ if (pid->IsTofEInfoValid()) {
         //~ cutsManager.smartFillCutHist( "hDp_pE_withTOF_E", e_tracks->mdcKalTrack()-> p(), !with_electron);
         //~ }
         //~
         //~ if (pid->IsTofCInfoValid()) {
         //~ cutsManager.smartFillCutHist( "hDp_pE_withTOF_C", e_tracks->mdcKalTrack()-> p(), !with_electron);
         //~ }




         // ========== MUC hits
         int nMucLayers =  e_tracks->isMucTrackValid() ?
               e_tracks->mucTrack()->numLayers() : 0;

         cutsManager.smartFillCutHist( "hDp_nMucLayers", nMucLayers,
               is_semimuonic);
         //~ cutsManager.rejectByCutIfHigher("hDp_nMucLayers", nMucLayers , 1.01);



         //~ smartFillHist("hDp_mMiss_mu",  (expected4Momentum - observed4Momentum).m2(), good, is_semimuonic);
         // ========== EMC deposit

         //~ smartFillHist("hDp_eEMC_inc", e_tracks->emcShower()->energy(), good);
         //~ smartFillHist("hDp_eEMC_exc", e_tracks->emcShower()->energy(), good, is_semimuonic);


         // === showers

         double massKPiAsPiPi = (kaon_tracks->mdcTrack()->p4(
                           MASS_PION) + pion_tracks->mdcTrack()->p4(MASS_PION)).m();

         cutsManager.smartFillCutHist( "hDp_mPiPi", massKPiAsPiPi );

         int nGoodShowers = 0;
         double totalShowersE =0;
         double maxShowerE = 0;

         for (unsigned shower_id = 0; shower_id < oshw.size(); ++shower_id) {
            if (fsr_showers.count(oshw[shower_id])) {
               continue;
            }

            DstEvtRecTracks* shower_tracks = (DstEvtRecTracks*) evtRecTrkCol->At(
                        oshw[shower_id]);
            RecEmcShower* emcTrk = shower_tracks->emcShower();

            if (isGoodShower(emcTrk, evtRecTrkCol, m_evtRecEvent)) {
               double energy = emcTrk->energy();
               nGoodShowers += 1;
               cout << "good shower id:" << oshw[shower_id] << endl;

               totalShowersE += energy;
               if (energy > maxShowerE) {
                  maxShowerE = energy;
               }
            }
         }


         cutsManager.smartFillCutHist( "hDp_nPions",
               get_num_kaon_pion(evtRecDTag->otherTracks(),evtRecDTag->pionId()));

         // pointless cut
         //~ if (get_num_kaon_pion(evtRecDTag->otherTracks(),evtRecDTag->pionId()) > 2) {
         //~ cutsManager.rejectByCut( "hDp_nPions");
         //~ }


         cutsManager.smartFillCutHist( "hDp_nKaons",
               get_num_kaon_pion(evtRecDTag->otherTracks(),evtRecDTag->kaonId()));

         // pointless cut (it is overlaped by extremely low pMiss cut)
         //~ if (get_num_kaon_pion(evtRecDTag->otherTracks(),evtRecDTag->kaonId()) > 1) {
         //~ cutsManager.rejectByCut( "hDp_nKaons");
         //~ }

         cutsManager.smartFillCutHist( "hDp_nCharged",
               evtRecDTag->otherTracks().size());
         cutsManager.smartFillCutHist( "hDp_nShowers",
               evtRecDTag->otherShowers().size(), hasPi0(m_TMcEvent));


         int nPi0 = 0;

         const TObjArray*  m_evtRecPi0Col = m_TEvtRecObject->getEvtRecPi0Col();
         if( m_evtRecPi0Col->GetEntries() ) {
            TIter evtRecPi0Iter(m_evtRecPi0Col);
            TEvtRecPi0* evtRecPi0 = NULL;
            while ((evtRecPi0 = (TEvtRecPi0*) evtRecPi0Iter.Next())) {
               if ( evtRecPi0->chisq() < 20 ) {
                  nPi0 += 1;
               };
               //~ smartFillHist("hDp_pi0Chisq", evtRecPi0->chisq(), good);
               //FIXME: check gammas against otherShowers
               //~ cout << "      hiEnGamma= " << evtRecPi0->hiEnGamma()
               //~ << " loEnGamma= " << evtRecPi0->loEnGamma() << endl;
            }
         }

         cutsManager.smartFillCutHist( "hDp_nPi0", nPi0, hasPi0(m_TMcEvent));
         cutsManager.smartFillCutHist( "hDp_nGoodShowers", nGoodShowers,
               hasPi0(m_TMcEvent));


         cutsManager.smartFillCutHist( "hDp_maxShowerE", maxShowerE,
               hasPi0(m_TMcEvent));
         cutsManager.rejectByCutIfHigher("hDp_maxShowerE", maxShowerE, 0.250);



         cutsManager.smartFillCutHist( "hDp_totalShowersE", totalShowersE,
               hasPi0(m_TMcEvent));
         //~ cutsManager.rejectByCutIfHigher("hDp_totalShowersE", totalShowersE , 0.25);



         //============ Chineese mass cut =============

         // M_miss with K as Pi. See http://www.sciencedirect.com/science?_ob=ArticleURL&_udi=B6TVN-4CXMXN8-2&_user=1562599&_coverDate=09/02/2004&_rdoc=6&_fmt=high&_orig=browse&_origin=browse&_zone=rslt_list_item&_srch=doc-info(%23toc%235539%232004%23994029998%23513688%23FLA%23display%23Volume)&_cdi=5539&_sort=d&_docanchor=&_ct=15&_acct=C000053742&_version=1&_urlVersion=0&_userid=1562599&md5=fa68574d6ccfc939c42796441dad3c3b&searchtype=a
         // doi:10.1016/j.physletb.2004.07.004 |

         double mMiss_K_as_pi = (expected4Momentum -
                     observed4Momentum_KaonAsPion).m2();


         cutsManager.smartFillCutHist( "hDp_mMiss_K_as_pi", mMiss_K_as_pi);
         //~ cutsManager.rejectByCutIfLower("hDp_mMiss_K_as_pi", mMiss_K_as_pi , 0.025);


         // ============ 2013.12.02 test muon bkg ========

         HepLorentzVector observed4Momentum_e_as_mu =
               kaon4Momentum +
               pion4Momentum +
               electron4Momentum
               - e_tracks->mdcKalTrack()->p4(MASS_ELECTRON)
               + e_tracks->mdcKalTrack()->p4(MASS_MUON);

         double mMiss_e_as_mu = (expected4Momentum -
                     observed4Momentum_e_as_mu).m2();
         cutsManager.smartFillCutHist( "hDp_mMiss_e_as_mu", mMiss_e_as_mu);


         // ===== Try to suppress D0 bkg ==== //

         // 1. there is Pi0 in our Dtag mode

         bool d0_201_3_found = 0;
         double d0_201_3_dE = 1;
         double d0_dE = 1;
         double d0_dE_sideband = 1;
         double d0_dE_sideband2 = 1;
         int d0_mode = -1;
         int d0_mode_sideband = -1;
         int d0_mode_sideband2 = -1;

         double d0_1_dE = 1;
         double d0_3_dE = 1;
         double d0_mbc = 1E10;

         if ( (evtRecDTag->decayMode() == 201) ||
               (evtRecDTag->decayMode() == 203) ||
               (evtRecDTag->decayMode() == 204) ) {

            TIter d0EvtRecDTagIter(m_evtRecDTagCol);
            while (TEvtRecDTag* d0EvtRecDTag = (TEvtRecDTag*)
                        d0EvtRecDTagIter.Next()) {
               if (d0EvtRecDTag->type() != 1) {
                  continue;
               }
               //~ if (d0EvtRecDTag->charm() != evtRecDTag->charm()) continue;
               if (d0EvtRecDTag->decayMode() > 126) {
                  continue;
               }

               if ((d0EvtRecDTag->mBC() > 1.860) && (d0EvtRecDTag->mBC() < 1.875) ) {
                  // all
                  if (fabs(d0EvtRecDTag->deltaE()) < fabs(d0_dE)) {
                     d0_dE = d0EvtRecDTag->deltaE();
                     d0_mode = d0EvtRecDTag->decayMode();
                  }
               }

               if ((d0EvtRecDTag->mBC() > 1.845) && (d0EvtRecDTag->mBC() < 1.860) ) {
                  // mass sideband
                  if (fabs(d0EvtRecDTag->deltaE()) < fabs(d0_dE_sideband)) {
                     d0_dE_sideband = d0EvtRecDTag->deltaE();
                     d0_mode_sideband = d0EvtRecDTag->decayMode();
                  }
               }

               if ((d0EvtRecDTag->mBC() > 1.83) && (d0EvtRecDTag->mBC() < 1.845) ) {
                  // mass sideband 2
                  if (fabs(d0EvtRecDTag->deltaE()) < fabs(d0_dE_sideband2)) {
                     d0_dE_sideband2 = d0EvtRecDTag->deltaE();
                     d0_mode_sideband2 = d0EvtRecDTag->decayMode();
                  }
               }




            }


         }


         cutsManager.smartFillCutHist( "hDp_d0_mode_sideband",
               d0_mode_sideband, hasDecay(eventTagDecay, 81, 8 ) );
         cutsManager.smartFillCutHist( "hDp_d0_dE_sideband", d0_dE_sideband,
               hasDecay(eventTagDecay, 81, 8 ) );
         cutsManager.smartFillCutHist( "hDp_d0_mode_sideband2",
               d0_mode_sideband2, hasDecay(eventTagDecay, 81, 8 ) );
         cutsManager.smartFillCutHist( "hDp_d0_dE_sideband2", d0_dE_sideband2,
               hasDecay(eventTagDecay, 81, 8 ) );

         cutsManager.smartFillCutHist( "hDp_d0_mode", d0_mode,
               hasDecay(eventTagDecay, 81, 8 ) );
         cutsManager.smartFillCutHist( "hDp_d0_dE", d0_dE,
               hasDecay(eventTagDecay, 81, 8 ) );


         cutsManager.rejectByCutIfInsideWindow("hDp_d0_dE", d0_dE, -0.01,
               0.01);


         cout << "pi trackid: " << pion_tracks->trackId() << endl;
         cout << "K trackid: " << kaon_tracks->trackId() << endl;
         cout << "e trackid: " << e_tracks->trackId() << endl;



         //============ Missing mass cut =============
         double mMiss =  (missing4Momentum).m2();

         double uMiss = (expected4Momentum - observed4Momentum).e() -
               (expected4Momentum - observed4Momentum).vect().mag();
         double uMiss_cms_bc = missing4MomentumBC_cms.e() -
               missing4MomentumBC_cms.vect().mag();
         double uMiss_cms = missing4Momentum_cms.e() -
               missing4Momentum_cms.vect().mag();
         double uMiss_cms2 = missing4Momentum_cms2.e() -
               missing4Momentum_cms2.vect().mag();


         cutsManager.smartFillCutHist( "hDp_mMissBC",
               missing4MomentumBC_cms.m2());

         cutsManager.smartFillCutHist( "hDp_mMiss", mMiss,
               hasDecay(eventTagDecay, 80, 11 ) );
         //~ cutsManager.rejectByCutIfOutsideWindow("hDp_mMiss", mMiss , -0.04, 0.04);


         cout << "U_miss_cms_BC= " << uMiss_cms << endl;
         if (cutsManager.survived()) {
            cout << "run= " << run << " event= " << event << " U_miss=" << uMiss
                  << endl;
         }



         cutsManager.smartFillCutHist( "hDp_uMiss_cms", uMiss_cms );
         cutsManager.smartFillCutHist( "hDp_uMiss_cms2", uMiss_cms2 );

         //~ cutsManager.rejectByCutIfOutsideWindow("hDp_uMiss_cms", uMiss_cms , -0.04, 0.04);

         cutsManager.smartFillCutHist( "hDp_uMiss", uMiss );
         //~ cutsManager.rejectByCutIfOutsideWindow("hDp_uMiss", uMiss , -0.04, 0.04);

         cutsManager.smartFillCutHist( "hDp_uMiss_cms_bc", uMiss_cms_bc );
         cutsManager.rejectByCutIfOutsideWindow("hDp_uMiss_cms_bc",
               uMiss_cms_bc, -0.04, 0.04);


         //~ cout << "ep_ratio: " << ep_ratio << endl;
         ntp_features->Fill(good, ep_ratio, missingMomentum.mag(), nMucLayers,
               totalShowersE,mMiss_K_as_pi,mMiss);


         if (cutsManager.survived() && good) {
            hist("hDp_mKPiTruth")->Fill(mKPi_true);
         }

         cutsManager.smartFillCutHist( "hDp_mKPi",   mKPi);

         if (cutsManager.survived()) {
            cout << "event has passed to mkpi cut" << endl;
            fillAngularVars(ntp_angular_selected, electron4Momentum,
                  missing4Momentum, kaon4Momentum, pion4Momentum, tag_passed, true,
                  decayMode, charm);



            electron4Momentum.boost(cmsBoost);
            kaon4Momentum.boost(cmsBoost);
            pion4Momentum.boost(cmsBoost);
            missing4Momentum.boost(cmsBoost);
            //~ TAngularVars rec_vars = getAngularVars(electron4Momentum, missing4MomentumBC_cms, kaon4Momentum, pion4Momentum);
            TAngularVars rec_vars = getAngularVars(electron4Momentum,
                        missing4Momentum, kaon4Momentum, pion4Momentum);

            cout << "electron4Momentum: " << electron4Momentum << endl;
            cout << "missing4momentumBC_cms: " << missing4MomentumBC_cms << endl;
            cout << "kaon4Momentum: " << kaon4Momentum << endl;
            cout << "pion4Momentum: " << pion4Momentum << endl;


            cout << "rec_vars: q2= " << rec_vars.q2 << " , mkpi= " <<
                  rec_vars.mkpi << " , theta_v= " << rec_vars.theta_v << " ,theta_l=" <<
                  rec_vars.theta_l << " ,chi= " << rec_vars.chi << endl;

            TAngularVars mc_vars;
            bool good_mc = getMcAngularVars(m_TMcEvent, charm, mc_vars);

            ntp_angular_selected_wmc->Fill(rec_vars.q2, rec_vars.mkpi,
                  rec_vars.theta_v, rec_vars.theta_l, rec_vars.chi,
                  mc_vars.q2, mc_vars.mkpi, mc_vars.theta_v, mc_vars.theta_l,
                  mc_vars.chi,
                  tag_passed, good_mc, decayMode);


            if (hasDecay(eventTagDecay, 81, 8 )) {
               preselection_passed = true;
            }





            event_passed = true;
         }



         //~ cutsManager.rejectByCutIfOutsideWindow("hDp_mKPi", mKPi , 0.8, 1.0);

         if (cutsManager.survived()) {
            real_event_survived += 1;
            charm_survived += 1;

            cout << "event survived!" << endl;


            smartFillHist("hDp_missCosTheta",
                  TMath::Abs(missingMomentum.cosTheta()), good);

            ((TH2D*) hist("hDp_e_p_vs_cost"))->Fill(e_p, e_cost);
            if (good) {
               ((TH2D*) hist("hDp_e_p_vs_cost_good"))->Fill(e_p, e_cost);
            };

            if (e_shower_energy > 0) {
               ((TH2D*) hist("hDp_e_p_vs_cost_we"))->Fill(e_p, e_cost);
               if (good) {
                  ((TH2D*) hist("hDp_e_p_vs_cost_we_good"))->Fill(e_p, e_cost);
               };
            }





            hist("hAllTags")->Fill(tag_str,1.0);
            if (! good) {
               cout << "not good, hDp_mKPi: " << mKPi  << endl;
               hist("hBadTags")->Fill(tag_str,1.0);
               hist(TString::Format("hBadTags_%d",
                           evtRecDTag->decayMode()))->Fill(tag_str, 1.0);
            }

            if (good && !good_tag_mode) {
               cout << "good with bad tag" << endl;
               hist(TString::Format("hGoodWrongTags_%d",
                           evtRecDTag->decayMode()))->Fill(tag_str, 1.0);
            }

            //TODO: ошибки из vertexfit


            HepPoint3D vx(0., 0., 0.);
            HepSymMatrix Evx(3, 0);
            double bx = 1E+6;
            double by = 1E+6;
            double bz = 1E+6;
            Evx[0][0] = bx*bx;
            Evx[1][1] = by*by;
            Evx[2][2] = bz*bz;



            VertexParameter vxpar;
            vxpar.setVx(vx);
            vxpar.setEvx(Evx);

            VertexFit* vtxfit = VertexFit::instance();
            vtxfit->init();

            HepLorentzVector tag4Momentum_alt;
            Hep3Vector tagMomentum_alt;

            const vector<Int_t>& dtag_tracks = evtRecDTag->tracks();
            //~
            for (unsigned id = 0; id < dtag_tracks.size(); ++id) {
               RecMdcKalTrack* kalTrk = ((DstEvtRecTracks*) evtRecTrkCol->At(
                                 dtag_tracks[id]))->mdcKalTrack();
               tagMomentum_alt += kalTrk->p3();

               double mass = -1;

               if ( is_pion_dtag( dtag_tracks[id], evtRecDTag) ) {
                  mass = MASS_PION;
               } else if ( is_kaon_dtag( dtag_tracks[id], evtRecDTag) ) {
                  mass = MASS_KAON;
               }

               WTrackParameter wvTrk = WTrackParameter(mass, kalTrk->getZHelix(),
                           kalTrk->getZError());
               tag4Momentum_alt += wvTrk.p();

               vtxfit->AddTrack(id,  wvTrk);
            }
            //~ cout << "tag4Momentum_alt:" << tag4Momentum_alt << endl;

            vtxfit->AddTrack(dtag_tracks.size(), WTrackParameter(MASS_ELECTRON,
                        e_tracks->mdcKalTrack()->getZHelix(),
                        e_tracks->mdcKalTrack()->getZError()));
            vtxfit->AddTrack(dtag_tracks.size() + 1, WTrackParameter(MASS_PION,
                        pion_tracks->mdcKalTrack()->getZHelix(),
                        pion_tracks->mdcKalTrack()->getZError()));
            vtxfit->AddTrack(dtag_tracks.size() + 2, WTrackParameter(MASS_KAON,
                        kaon_tracks->mdcKalTrack()->getZHelix(),
                        kaon_tracks->mdcKalTrack()->getZError()));

            vtxfit->AddVertex(0, vxpar, 0, 1, 2, 3, 4, 5);


            if(!vtxfit->Fit(0)) {
               cout << "bad fit!!" << endl;
            }

            vtxfit->BuildVirtualParticle(0);

            //~ cout << "vtxfit0->chisq(0): " << vtxfit->chisq(0) << endl;

            WTrackParameter virt_vertex = vtxfit->wVirtualTrack(0);

            //~ vtxfit->Swim(0);

            HepLorentzVector missing4Momentum_vxfit = ecms - virt_vertex.p();
            cout << missing4Momentum_vxfit << endl;
            cout << (expected4Momentum - observed4Momentum) << endl;


            HepVector D(4, 0); //derivative vector d(m^2)/d(p_i)
            D[0] = 2 * missing4Momentum_vxfit.px();
            D[1] = 2 * missing4Momentum_vxfit.py();
            D[2] = 2 * missing4Momentum_vxfit.pz();
            D[3] = -2 * missing4Momentum_vxfit.e();


            cout << "missing mass (vxfit): ("
                  << missing4Momentum_vxfit.m2()
                  << " +/- "
                  << sqrt((D.T() * virt_vertex.Ep() * D)[0])
                  << " ) GeV^2/c^2" << endl;

            //~ smartFillHist("hDp_mMiss_vxfit",  (missing4Momentum_vxfit).m2(), good);
            //~ smartFillHist("hDp_mMiss_err",  sqrt((D.T() * virt_vertex.Ep() * D)[0]), good);

            //~ HepLorentzVector ecms2(0, 0,0 , getZ0Energy(m_TMcEvent));

            //~ smartFillHist("hDp_mMiss_alt",  (ecms - tag4Momentum_alt - observed4Momentum ).m2(), good);


            HepLorentzVector expected4Momentum_alt;

            expected4Momentum_alt.setVectM(ecms.vect() -tagMomentum_alt,
                  D_PLUS_MASS);

            //~ cout << "taggedDMomentum: " << taggedDMomentum << endl;
            //~ cout << "tagMomentum_alt: " << tagMomentum_alt << endl;

            //~ smartFillHist("hDp_mMiss_alt2",  (expected4Momentum_alt - observed4Momentum).m2(), good);


            //~ hist("hDp_Z0E")->Fill(getZ0Energy(m_TMcEvent));
            //we allow only one combination per event
         }

         cutsManager.processCuts();

         // electron PID systematics
         //~ if (good) {
//~
         //~ is_electron( e_tracks->trackId(), pid, evtRecTrkCol, run);
//~
         //~ double pid_prob_e = pid->probElectron() ;
         //~ double pid_prob_sum = pid->probPion() + pid->probKaon() + pid->probElectron();
//~
//~
//~
         //~ ntp_electrons->Fill(e_shower_energy,
         //~ e_tracks->mdcKalTrack()->p(),
         //~ TMath::Cos(e_tracks->mdcKalTrack()->theta()),
         //~ e_tracks->mdcKalTrack()->charge(),
         //~ 0,0,pid_prob_e,pid_prob_sum,event_passed);
         //~ }


         fillMcTracks(m_TMcEvent, charm, ntp_angular_mc, tag_passed,
               event_passed, decayMode);
         //~ fillMcMomenta(m_TMcEvent, charm, ntp_momenta_mc, tag_passed, event_passed);

      }
      hist("hCharmSurvived")->Fill(charm_survived);
   }

   hist("hRealEventSurvived")->Fill(real_event_survived);

   return preselection_passed;

}

//-----------------------------------------------------------------------------
BeanUserShared_EXPORT
void SemileptonicEndJob(ReadDst* selector)
//-----------------------------------------------------------------------------
{
   if( selector->Verbose() ) {
      cout << " TestPIDEndJob() " << endl;
   }
}

#ifdef __cplusplus
}
#endif
//~ root [71] for (int i=0; i< hDp_pE_eff->GetNbinsX()  ; ++i) { hDp_pE_eff->SetBinContent(i, hDp_pE_bad->Integral(i,hDp_pE_good->GetNbinsX() ) ? hDp_pE_good->Integral(i,hDp_pE_good->GetNbinsX())*1.0 / sqrt(hDp_pE_bad->Integral(i,hDp_pE_good->GetNbinsX() )) : 0 );}

//TODO: оценить эффективность оценки трёх треков
//+абсолютные числа
//+другие таги
//test2






// ======= вопросы 19.09.12
// 1. в определении uMiss/mMiss и т.д., точнее в определении missing4momentum используются измеренные импульсы+массы.
// Внимание, вопрос: чем чунлею мешают FSR фотоны?



// K-star shape factor =0.893084
