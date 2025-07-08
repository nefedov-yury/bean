// #include <iostream>
#include "BExtTrack.h"

// None: see ExtTrackCnv.cxx

//--------------------------------------------------------------------
CLHEP::HepSymMatrix BExtTrack::tof1ErrorMatrix() const
//--------------------------------------------------------------------
{
   CLHEP::HepSymMatrix tmp(6);
   for(int i = 0; i < 6; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = m_trk->GetTof1ErrorMatrix(i,j);
      }
   }
   return tmp;
}

//--------------------------------------------------------------------
CLHEP::HepSymMatrix BExtTrack::tof2ErrorMatrix() const
//--------------------------------------------------------------------
{
   CLHEP::HepSymMatrix tmp(6);
   for(int i = 0; i < 6; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = m_trk->GetTof2ErrorMatrix(i,j);
      }
   }
   return tmp;
}

//--------------------------------------------------------------------
CLHEP::HepSymMatrix BExtTrack::emcErrorMatrix() const
//--------------------------------------------------------------------
{
   CLHEP::HepSymMatrix tmp(6);
   for(int i = 0; i < 6; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = m_trk->GetEmcErrorMatrix(i,j);
      }
   }
   return tmp;
}

//--------------------------------------------------------------------
CLHEP::HepSymMatrix BExtTrack::mucErrorMatrix() const
//--------------------------------------------------------------------
{
   CLHEP::HepSymMatrix tmp(6);
   for(int i = 0; i < 6; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = m_trk->GetMucErrorMatrix(i,j);
      }
   }
   return tmp;
}
