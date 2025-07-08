#ifndef _BExtTrack_h
#define _BExtTrack_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// BExtTrack                                                        //
// * This is adapter between RootEventData/TExtTrack.h and          //
//   DstEvent/DstExtTrack.h                                         //
// * Root file contains information only for pion extrapolation,    //
//   therefore all functions with 'int parID' can not be usable     //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <string>

#include <CLHEP/Matrix/SymMatrix.h>
// #include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include "RootEventData/TExtTrack.h"

using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

class BExtTrack
{
   public:
      BExtTrack(TExtTrack* trk) {m_trk = trk;}

      const TExtTrack* getTExtTrack() const {return m_trk;}

      int GetTrackId() const {return m_trk->GetTrackId();}
      int trackId() const    {return m_trk->GetTrackId();}

      //Get track extrapolation data @ Tof layer1.
      CLHEP::Hep3Vector tof1Position() const {
         return CLHEP::Hep3Vector(
               m_trk->GetTof1PositionX(),
               m_trk->GetTof1PositionY(),
               m_trk->GetTof1PositionZ());
      }
      CLHEP::Hep3Vector tof1Momentum() const {
         return CLHEP::Hep3Vector(
               m_trk->GetTof1MomentumX(),
               m_trk->GetTof1MomentumY(),
               m_trk->GetTof1MomentumZ());
      }
      std::string tof1VolumeName() const {
         return std::string((m_trk->GetTof1VolumeName()).Data());}
      int tof1VolumeNumber() const {
         return m_trk->GetTof1VolumeNumber();
      }
      double tof1() const {return m_trk->GetTof1();}
      double tof1Path() const {return m_trk->GetTof1Path();}
      double tof1PosSigmaAlongZ() const {
         return m_trk->GetTof1PosSigmaAlongZ();
      }
      double tof1PosSigmaAlongT() const {
         return m_trk->GetTof1PosSigmaAlongT();
      }
      double tof1PosSigmaAlongX() const {
         return m_trk->GetTof1PosSigmaAlongX();
      }
      double tof1PosSigmaAlongY() const {
         return m_trk->GetTof1PosSigmaAlongY();
      }
      CLHEP::HepSymMatrix tof1ErrorMatrix() const;

      //Get track extrapolation data @ Tof layer2.
      CLHEP::Hep3Vector tof2Position() const {
         return CLHEP::Hep3Vector(
               m_trk->GetTof2PositionX(),
               m_trk->GetTof2PositionY(),
               m_trk->GetTof2PositionZ());
      }
      CLHEP::Hep3Vector tof2Momentum() const {
         return CLHEP::Hep3Vector(
               m_trk->GetTof2MomentumX(),
               m_trk->GetTof2MomentumY(),
               m_trk->GetTof2MomentumZ());
      }
      std::string tof2VolumeName() const {
         return std::string((m_trk->GetTof2VolumeName()).Data());
      }
      int tof2VolumeNumber() const {
         return m_trk->GetTof2VolumeNumber();
      }
      double tof2() const {return m_trk->GetTof2();}
      double tof2Path() const {return m_trk->GetTof2Path();}
      double tof2PosSigmaAlongZ() const {
         return m_trk->GetTof2PosSigmaAlongZ();
      }
      double tof2PosSigmaAlongT() const {
         return m_trk->GetTof2PosSigmaAlongT();
      }
      double tof2PosSigmaAlongX() const {
         return m_trk->GetTof2PosSigmaAlongX();
      }
      double tof2PosSigmaAlongY() const {
         return m_trk->GetTof2PosSigmaAlongY();
      }
      CLHEP::HepSymMatrix tof2ErrorMatrix() const;

      //Get track extrapolation data @ EMC.
      CLHEP::Hep3Vector emcPosition() const {
         return CLHEP::Hep3Vector(
               m_trk->GetEmcPositionX(),
               m_trk->GetEmcPositionY(),
               m_trk->GetEmcPositionZ());
      }
      CLHEP::Hep3Vector emcMomentum() const {
         return CLHEP::Hep3Vector(
               m_trk->GetEmcMomentumX(),
               m_trk->GetEmcMomentumY(),
               m_trk->GetEmcMomentumZ());
      }
      std::string emcVolumeName() const {
         return std::string((m_trk->GetEmcVolumeName()).Data());
      }
      int emcVolumeNumber() const {
         return m_trk->GetEmcVolumeNumber();
      }
      double emcPosSigmaAlongTheta() const {
         return m_trk->GetEmcPosSigmaAlongTheta();
      }
      double emcPosSigmaAlongPhi() const {
         return m_trk->GetEmcPosSigmaAlongPhi();
      }
      CLHEP::HepSymMatrix emcErrorMatrix() const;
      double emcPath() const {return m_trk->emcPath();}

      //Get track extrapolation data @ MUC.
      CLHEP::Hep3Vector mucPosition() const {
         return CLHEP::Hep3Vector(
               m_trk->GetMucPositionX(),
               m_trk->GetMucPositionY(),
               m_trk->GetMucPositionZ());
      }
      CLHEP::Hep3Vector mucMomentum() const {
         return CLHEP::Hep3Vector(
               m_trk->GetMucMomentumX(),
               m_trk->GetMucMomentumY(),
               m_trk->GetMucMomentumZ());
      }
      std::string mucVolumeName() const {
         return std::string((m_trk->GetMucVolumeName()).Data());
      }
      int mucVolumeNumber() const {
         return m_trk->GetMucVolumeNumber();
      }
      double mucPosSigmaAlongZ() const {
         return m_trk->GetMucPosSigmaAlongZ();
      }
      double mucPosSigmaAlongT() const {
         return m_trk->GetMucPosSigmaAlongT();
      }
      double mucPosSigmaAlongX() const {
         return m_trk->GetMucPosSigmaAlongX();
      }
      double mucPosSigmaAlongY() const {
         return m_trk->GetMucPosSigmaAlongY();
      }
      CLHEP::HepSymMatrix mucErrorMatrix() const;

   private:
      TExtTrack* m_trk;
      // int myParticleType;// it is always = 2
};
#endif
