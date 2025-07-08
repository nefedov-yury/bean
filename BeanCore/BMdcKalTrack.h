#ifndef _BMdcKalTrack_h
#define _BMdcKalTrack_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// BMdcKalTrack                                                     //
// * This is adapter between RootEventData/TMdcKalTrack.h           //
//   and DstEvent/DstMdcKalTrack.h                                  //
//   and partialy MdcRecEvent/RecMdcKalTrack.h                      //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
typedef HepGeom::Point3D<double> HepPoint3D;

using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "RootEventData/TMdcKalTrack.h"

class BMdcKalTrack
{
   public:
      enum PidType
      {
         null = -1,
         electron = 0,
         muon = 1,
         pion = 2,
         kaon = 3,
         proton = 4
      };

   public:
      BMdcKalTrack(TMdcKalTrack* root_trk) {
         m_trk = root_trk; m_pid = pion;
      }

      const TMdcKalTrack* getTMdcKalTrack() const {return m_trk;}

      void setPidType(PidType pidType = pion) {m_pid = pidType;}
      PidType getPidType() const {return m_pid;}

      int trackId() const {return m_trk->getTrackId();}
      // const double mass() const { return m_mass[m_pid];  }
      int charge() const {return getCharge(m_pid);}
      double pxy() const;
      double px() const;
      double py() const;
      double pz() const;
      double p() const;
      double theta() const;
      double phi() const;

      double x() const {return x(m_pid);}
      double y() const {return y(m_pid);}
      double z() const {return z(m_pid);}

      double x(int pid) const;
      double y(int pid) const;
      double z(int pid) const;
      double r() const;

      int stat() const {return m_trk->getStat(m_pid);}
      double chi2() const {return m_trk->getChisq(m_pid);}
      int ndof() const {return m_trk->getNdf(m_pid);}
#if BOSS_VER < 663
      int nster() const {return m_trk->getNster(m_pid);}
      int firstLayer() const {return m_trk->getFirstLayer(m_pid);}
      int lastLayer() const {return m_trk->getLastLayer(m_pid);}
#endif

      double dr()    const {return getZHelix_pid(0,m_pid);}
      double fi0()   const {return getZHelix_pid(1,m_pid);}
      double kappa() const {return getZHelix_pid(2,m_pid);}
      double dz()    const {return getZHelix_pid(3,m_pid);}
      double tanl()  const {return getZHelix_pid(4,m_pid);}

      HepVector helix()   const {return getZHelix(m_pid);}
      HepSymMatrix err()  const {return getZError(m_pid);}
      HepVector fhelix()  const {return getFHelix(m_pid);}
      HepSymMatrix ferr() const {return getFError(m_pid);}
      // HepPoint3D poca()   const {return getPoca(m_pid);}
      Hep3Vector p3() const;
      HepPoint3D x3() const;

      // HepLorentzVector p4() const;
      HepLorentzVector p4(double mass) const;

      int getTrackId() const {return m_trk->getTrackId();}
      int getCharge(const int pid) const;
      int getStat(const int pid) const {return m_trk->getStat(pid);}
      double getChisq(const int pid) const {
         return m_trk->getChisq(pid);
      }
      int getNdf(const int pid) const {return m_trk->getNdf(pid);}
#if BOSS_VER < 663
      int getNster(const int pid) const {return m_trk->getNster(pid);}
      int getFirstLayer(const int pid) const {
         return m_trk->getFirstLayer(pid);
      }
      int getLastLayer(const int pid) const {
         return m_trk->getLastLayer(pid);
      }
#endif

      HepVector getZHelix(const int pid) const;
      HepSymMatrix getZError(const int pid) const;
      HepVector getFHelix(const int pid) const;
      HepSymMatrix getFError(const int pid) const;
      // HepPoint3D getPoca(const int pid) const;

      // this is functions from MdcRecEvent/RecMdcKalTrack.h
      HepVector    getZHelix() const {return getZHelix(pion);}
      HepSymMatrix getZError() const {return getZError(pion);}
      HepVector    getFHelix() const {return getFHelix(pion);}
      HepSymMatrix getFError() const {return getFError(pion);}

      HepVector    getZHelixE() const {return getZHelix(electron);}
      HepSymMatrix getZErrorE() const {return getZError(electron);}
      HepVector    getFHelixE() const {return getFHelix(electron);}
      HepSymMatrix getFErrorE() const {return getFError(electron);}

      HepVector    getZHelixMu() const {return getZHelix(muon);}
      HepSymMatrix getZErrorMu() const {return getZError(muon);}
      HepVector    getFHelixMu() const {return getFHelix(muon);}
      HepSymMatrix getFErrorMu() const {return getFError(muon);}

      HepVector    getZHelixK() const {return getZHelix(kaon);}
      HepSymMatrix getZErrorK() const {return getZError(kaon);}
      HepVector    getFHelixK() const {return getFHelix(kaon);}
      HepSymMatrix getFErrorK() const {return getFError(kaon);}

      HepVector    getZHelixP() const {return getZHelix(proton);}
      HepSymMatrix getZErrorP() const {return getZError(proton);}
      HepVector    getFHelixP() const {return getFHelix(proton);}
      HepSymMatrix getFErrorP() const {return getFError(proton);}

      // HepPoint3D getPocaE()  const { return getPoca(electron); }
      // HepPoint3D getPocaMu() const { return getPoca(muon); }
      // HepPoint3D getPoca()   const { return getPoca(pion); }
      // HepPoint3D getPocaK()  const { return getPoca(kaon); }
      // HepPoint3D getPocaP()  const { return getPoca(proton); }

      // set functions: see DstEvent/DstMdcKalTrack.h
      void setZHelix(const HepVector& helix, const int pid);
      void setZHelix(double* helix, const int pid);
      void setZError(const HepSymMatrix& error, const int pid);
      void setZError(double* error, const int pid);

      void setFHelix(const HepVector& fhelix, const int pid);
      void setFHelix(double* fhelix, const int pid);
      void setFError(const HepSymMatrix& ferror, const int pid);
      void setFError(double* ferror, const int pid);

   private:
      TMdcKalTrack* m_trk;
      PidType m_pid;

      double getZHelix_pid(int i,int pid = 2) const {
         switch( pid ) {
            case 0: return m_trk->getZHelixE(i);
            case 1: return m_trk->getZHelixMu(i);
            case 2: return m_trk->getZHelix(i);
            case 3: return m_trk->getZHelixK(i);
            case 4: return m_trk->getZHelixP(i);
         }
         return m_trk->getZHelix(i);
      }

      double getZError_pid(int i,int j,int pid = 2) const {
         switch( pid ) {
            case 0: return m_trk->getZErrorE(i,j);
            case 1: return m_trk->getZErrorMu(i,j);
            case 2: return m_trk->getZError(i,j);
            case 3: return m_trk->getZErrorK(i,j);
            case 4: return m_trk->getZErrorP(i,j);
         }
         return m_trk->getZError(i,j);
      }

      // ther are no getPoca.. functions for RootEventData >= 6.5.5
      // double getPoca_pid(int i,int pid = 2) const {
         // switch( pid ) {
            // case 0: return m_trk->getPocaE(i);
            // case 1: return m_trk->getPocaMu(i);
            // case 2: return m_trk->getPoca(i);
            // case 3: return m_trk->getPocaK(i);
            // case 4: return m_trk->getPocaP(i);
         // }
         // return m_trk->getPoca(i);
      // }

      double getFHelix_pid(int i,int pid = 2) const {
         switch( pid ) {
            case 0: return m_trk->getFHelixE(i);
            case 1: return m_trk->getFHelixMu(i);
            case 2: return m_trk->getFHelix(i);
            case 3: return m_trk->getFHelixK(i);
            case 4: return m_trk->getFHelixP(i);
         }
         return m_trk->getFHelix(i);
      }

      double getFError_pid(int i,int j,int pid = 2) const {
         switch( pid ) {
            case 0: return m_trk->getFErrorE(i,j);
            case 1: return m_trk->getFErrorMu(i,j);
            case 2: return m_trk->getFError(i,j);
            case 3: return m_trk->getFErrorK(i,j);
            case 4: return m_trk->getFErrorP(i,j);
         }
         return m_trk->getFError(i,j);
      }
};
#endif
