// #include <iostream>

#include "RootEventData/TEvtRecTrack.h"
#include "RootEventData/TDstEvent.h"

#include "DstEvtRecTracks.h"

//--------------------------------------------------------------------
DstEvtRecTracks::DstEvtRecTracks(TEvtRecTrack* rec_trk,
      TDstEvent* dst_event)
//--------------------------------------------------------------------
{
   mdc_trk = nullptr;
   mdc_dedx = nullptr;
   mdc_kal_trk = nullptr;
   ext_trk = nullptr;
   emc_trk = nullptr;
   muc_trk = nullptr;

   m_rec_trk = rec_trk;
   if( !rec_trk ) {
      return;
   }

   int mdc_trk_id = rec_trk->mdcTrackId();
   if( mdc_trk_id >= 0 ) {
      const TObjArray* m_mdcTrackCol = dst_event->getMdcTrackCol();
      TMdcTrack* tmp = (TMdcTrack*) m_mdcTrackCol->At(mdc_trk_id);
      if( tmp ) {
         mdc_trk = new RecMdcTrack(tmp);
      }
   }

   int mdc_dedx_trk_id = rec_trk->mdcDedxId();
   if( mdc_dedx_trk_id >= 0 ) {
      const TObjArray* m_mdcDedxCol = dst_event->getMdcDedxCol();
      TMdcDedx* tmp = (TMdcDedx*) m_mdcDedxCol->At(mdc_dedx_trk_id);
      if( tmp ) {
         mdc_dedx = new RecMdcDedx(*tmp); // just copy
      }
   }

   int mdc_kal_trk_id = rec_trk->mdcKalTrackId();
   if( mdc_kal_trk_id >= 0 ) {
      const TObjArray* m_mdcKalTrackCol =
         dst_event->getMdcKalTrackCol();
      TMdcKalTrack* tmp =
         (TMdcKalTrack*) m_mdcKalTrackCol->At(mdc_kal_trk_id);
      if( tmp ) {
         mdc_kal_trk = new RecMdcKalTrack(tmp);
      }
   }

   const std::vector<int>& tof_trk_id = rec_trk->tofTrackIds();
   size_t id_size = tof_trk_id.size();
   if( id_size > 0 ) {
      const TObjArray* m_tofTrackCol = dst_event->getTofTrackCol();
      for( size_t i = 0; i < id_size; i++ ) {
         TTofTrack* tmp =
            (TTofTrack*) m_tofTrackCol->At(tof_trk_id[i]);
         if( tmp ) {
            tof_trk.push_back( new RecTofTrack(tmp) );
         }
      }
   }

   int ext_trk_id = rec_trk->extTrackId();
   if( ext_trk_id >= 0 ) {
      const TObjArray* m_extTrackCol = dst_event->getExtTrackCol();
      TExtTrack* tmp = (TExtTrack*) m_extTrackCol->At(ext_trk_id);
      if( tmp ) {
         ext_trk = new RecExtTrack(tmp);
      }
   }

   int emc_trk_id = rec_trk->emcShowerId();
   if( emc_trk_id >= 0 ) {
      const TObjArray* m_emcTrackCol = dst_event->getEmcTrackCol();
      TEmcTrack* tmp = (TEmcTrack*) m_emcTrackCol->At(emc_trk_id);
      if( tmp ) {
         emc_trk = new RecEmcShower(tmp);
      }
   }

   int muc_trk_id = rec_trk->mucTrackId();
   if( muc_trk_id >= 0 ) {
      const TObjArray* m_mucTrackCol = dst_event->getMucTrackCol();
      TMucTrack* tmp = (TMucTrack*) m_mucTrackCol->At(muc_trk_id);
      if( tmp ) {
         muc_trk = new RecMucTrack(*tmp); // just copy
      }
   }
}

//--------------------------------------------------------------------
DstEvtRecTracks::~DstEvtRecTracks()
//--------------------------------------------------------------------
{
   delete(mdc_trk);
   delete(mdc_dedx);
   delete(mdc_kal_trk);
   for( auto& tt : tof_trk ) {
      delete(tt);
   }
   delete(ext_trk);
   delete(emc_trk);
   delete(muc_trk);
}

//--------------------------------------------------------------------
int DstEvtRecTracks::trackId() const
//--------------------------------------------------------------------
{
   return m_rec_trk->trackId();
}
