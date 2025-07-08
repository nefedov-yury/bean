#ifndef DstEvtRecTracks_H
#define DstEvtRecTracks_H
//////////////////////////////////////////////////////////////////////
//                                                                  //
// DstEvtRecTracks                                                  //
// * This is a transitional class intended to be used in place of   //
//   "EvtRecEvent/EvtRecTrack.h"                                    //
// * SmartRefVector<RecTofTrack> -> std::vector<RecTofTrack*>&      //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <vector>

#include <TObject.h>

#include "RootEventData/TDstEvent.h"

#include "BMdcTrack.h"
#include "BMdcKalTrack.h"
#include "BTofTrack.h"
#include "BExtTrack.h"
#include "BEmcTrack.h"

class TEvtRecTrack;

typedef  BMdcTrack        RecMdcTrack;
typedef  TMdcDedx         RecMdcDedx;
typedef  BMdcKalTrack     RecMdcKalTrack;
typedef  BTofTrack        RecTofTrack;
typedef  BExtTrack        RecExtTrack;
typedef  BEmcTrack        RecEmcShower;
typedef  TMucTrack        RecMucTrack;

class DstEvtRecTracks : public TObject
{
   public:
      DstEvtRecTracks(TEvtRecTrack* rec_trk,TDstEvent* dst_event);
      ~DstEvtRecTracks();

      int trackId() const;

      // TMdcTrack
      bool isMdcTrackValid() const {return (mdc_trk != nullptr);}
      RecMdcTrack* mdcTrack() const {return mdc_trk;}

      // TMdcDedx
      bool isMdcDedxValid() const {return (mdc_dedx != nullptr);}
      RecMdcDedx* mdcDedx() const {return mdc_dedx;}

      // TMdcKalTrack
      bool isMdcKalTrackValid() const {return (mdc_kal_trk!=nullptr);}
      RecMdcKalTrack* mdcKalTrack() const {return mdc_kal_trk;}

      // TofTrack
      bool isTofTrackValid() const {return (!tof_trk.empty());}
      const std::vector<RecTofTrack* >& tofTrack() const {
         return tof_trk;
      }

      // TExtTrack
      bool isExtTrackValid() const {return (ext_trk != nullptr);}
      RecExtTrack* extTrack() const {return ext_trk;}

      // TEmcTrack
      bool isEmcShowerValid() const {return (emc_trk != nullptr);}
      RecEmcShower* emcShower() const {return emc_trk;}

      // TMucTrack
      bool isMucTrackValid() const {return (muc_trk != nullptr);}
      RecMucTrack* mucTrack() const {return muc_trk;}

   private:
      TEvtRecTrack*              m_rec_trk;
      RecMdcTrack*               mdc_trk;
      RecMdcDedx*                mdc_dedx;
      RecMdcKalTrack*            mdc_kal_trk;
      RecExtTrack*               ext_trk;
      std::vector<RecTofTrack* > tof_trk;
      RecEmcShower*              emc_trk;
      RecMucTrack*               muc_trk;
};

typedef  DstEvtRecTracks  EvtRecTrack;
#endif
