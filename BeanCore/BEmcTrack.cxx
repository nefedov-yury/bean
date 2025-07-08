// #include <iostream>
#include "BEmcTrack.h"

//--------------------------------------------------------------------
HepPoint3D BEmcTrack::position() const
//--------------------------------------------------------------------
{
   return HepPoint3D(x(),y(),z());
}

//--------------------------------------------------------------------
CLHEP::HepSymMatrix BEmcTrack::errorMatrix() const
//--------------------------------------------------------------------
{
   CLHEP::HepSymMatrix tmp(3);
   // Error Matrix: 0:dxx, 1:dyy, 2:dzz
   //               3:dxy, 4:dxz, 5:dyz
   tmp[0][0] = err(0);
   tmp[1][1] = err(1);
   tmp[2][2] = err(2);
   tmp[0][1] = err(3);
   tmp[0][2] = err(4);
   tmp[1][2] = err(5);

   return tmp;
}
