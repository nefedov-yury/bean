#ifndef _BEmcTrack_h
#define _bEmcTrack_h 1
//////////////////////////////////////////////////////////////////////
//                                                                  //
// BEmcTrack                                                        //
// * This is adapter between RootEventData/TEmcTrack.h and          //
//   DstEvent/DstEmcTrack.h                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Geometry/Point3D.h>
typedef HepGeom::Point3D<double> HepPoint3D;

#include "RootEventData/TEmcTrack.h"

class BEmcTrack : public TEmcTrack
{
   public:
      BEmcTrack(const TEmcTrack* trk) : TEmcTrack(*trk) {}
      ~BEmcTrack() = default;

      // "CLHEP-functions" absent in TEmcTrack
      HepPoint3D position() const;
      CLHEP::HepSymMatrix errorMatrix() const;
};
#endif
