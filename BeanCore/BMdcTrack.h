#ifndef _BMdcTrack_h
#define _BMdcTrack_h 1
//////////////////////////////////////////////////////////////////////
//                                                                  //
// BMdcTrack                                                        //
// * This is adapter between RootEventData/TMdcTrack.h and          //
//   DstEvent/DstMdcTrack.h                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
typedef HepGeom::Point3D<double> HepPoint3D;

#include "RootEventData/TMdcTrack.h"

class BMdcTrack : public TMdcTrack
{
   public:
      BMdcTrack(const TMdcTrack* trk) : TMdcTrack(*trk) {}
      ~BMdcTrack() = default;

      // this is to avoid name mismatch:
      using TMdcTrack::helix;
      using TMdcTrack::err;

      // "CLHEP-functions" absent in TMdcTrack
      CLHEP::HepVector helix() const;
      CLHEP::HepSymMatrix err() const;
      CLHEP::Hep3Vector p3() const;
      CLHEP::HepLorentzVector p4(const double mass) const;
      // poca - position of closest approach to origin
      HepPoint3D poca() const;
      // The Coordinate of Track origin(m)
      HepPoint3D x3() const;
};
#endif
