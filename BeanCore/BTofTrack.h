#ifndef _BTofTrack_h
#define _BTofTrack_h 1
//////////////////////////////////////////////////////////////////////
//                                                                  //
// BTofTrack                                                        //
// * This is adapter between RootEventData/TTofTrack.h and          //
//   DstEvent/DstTofTrack.h                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "RootEventData/TTofTrack.h"
#include "TofHitStatus.h"

class BTofTrack : public TTofTrack
{
   public:
      BTofTrack(const TTofTrack* trk) : TTofTrack(*trk) {}

      // two new function in BOSS-7 (DstEvent-00-02-51)
      int tofID() {
         auto t_tofID = TTofTrack::tofID();
         if( t_tofID < 0 ) {
            return t_tofID;
         }
         if( TofHitStatus::is_mrpc( status() ) ) {
            return static_cast<int>( t_tofID/12 );
         } else {
            return t_tofID;
         }
         return -1;
      }

      int strip() {
         auto t_tofID = TTofTrack::tofID();
         if( t_tofID<0 ) {
            return -1;
         }
         if( TofHitStatus::is_mrpc( status() ) ) {
            return static_cast<int>( t_tofID%12 );
         }
         return -1;
      }

};
#endif
