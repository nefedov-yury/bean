// #include <iostream>
#include <cmath>

#include "BMdcKalTrack.h"

//--------------------------------------------------------------------
int BMdcKalTrack::getCharge(const int pid) const
//--------------------------------------------------------------------
{
   int charge = 0;
   double kappa = getZHelix_pid(2,pid);
   if( kappa > 0.0000000001 ) {
      charge = 1 ;
   } else if( kappa < -0.0000000001 ) {
      charge = -1;
   } else {
      charge = 0;
   }
   return charge;
}

//--------------------------------------------------------------------
double BMdcKalTrack::pxy() const
//--------------------------------------------------------------------
{
   double kappa = getZHelix_pid(2,m_pid);
   if( kappa != 0 ) {
      return 1./fabs(kappa);
   }
   return 0.;
}

//--------------------------------------------------------------------
double BMdcKalTrack::px() const
//--------------------------------------------------------------------
{
   double phi0 = getZHelix_pid(1,m_pid);
   return pxy()*(-sin(phi0));
}

//--------------------------------------------------------------------
double BMdcKalTrack::py() const
//--------------------------------------------------------------------
{
   double phi0 = getZHelix_pid(1,m_pid);
   return pxy()*cos(phi0);
}

//--------------------------------------------------------------------
double BMdcKalTrack::pz() const
//--------------------------------------------------------------------
{
   double tanl = getZHelix_pid(4,m_pid);
   return pxy()*tanl;
}

//--------------------------------------------------------------------
double BMdcKalTrack::p() const
//--------------------------------------------------------------------
{
   double tanl = getZHelix_pid(4,m_pid);
   return pxy()*sqrt(1. + tanl*tanl);
}

//--------------------------------------------------------------------
double BMdcKalTrack::theta() const
//--------------------------------------------------------------------
{
   double tanl = getZHelix_pid(4,m_pid);
   return acos(tanl/sqrt(1. + tanl*tanl));
}

//--------------------------------------------------------------------
double BMdcKalTrack::phi() const
//--------------------------------------------------------------------
{
   // return atan2(py(),px());
   double phi0 = getZHelix_pid(1,m_pid);
   return atan2( cos(phi0),-sin(phi0) );
}

//--------------------------------------------------------------------
double BMdcKalTrack::x(int pid) const
//--------------------------------------------------------------------
{
   double dr   = getZHelix_pid(0,pid);
   double phi0 = getZHelix_pid(1,pid);
   return dr*cos(phi0);
}

//--------------------------------------------------------------------
double BMdcKalTrack::y(int pid) const
//--------------------------------------------------------------------
{
   double dr   = getZHelix_pid(0,pid);
   double phi0 = getZHelix_pid(1,pid);
   return dr*sin(phi0);
}

//--------------------------------------------------------------------
double BMdcKalTrack::z(int pid) const
//--------------------------------------------------------------------
{
   double dz = getZHelix_pid(3,pid);
   return dz;
}

//--------------------------------------------------------------------
double BMdcKalTrack::r() const
//--------------------------------------------------------------------
{
   double dr = getZHelix_pid(0,m_pid);
   return fabs(dr);
}

//--------------------------------------------------------------------
HepVector BMdcKalTrack::getZHelix(const int pid) const
//--------------------------------------------------------------------
{
   HepVector tmp(5);
   for(int i = 0; i < 5; i++) {
      tmp[i] = getZHelix_pid(i,pid);
   }
   return tmp;
}

//--------------------------------------------------------------------
HepSymMatrix BMdcKalTrack::getZError(const int pid) const
//--------------------------------------------------------------------
{
   HepSymMatrix tmp(5);
   for(int i = 0; i < 5; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = getZError_pid(i,j,pid);
      }
   }
   return tmp;
}

//--------------------------------------------------------------------
HepVector BMdcKalTrack::getFHelix(const int pid) const
//--------------------------------------------------------------------
{
   HepVector tmp(5);
   for(int i = 0; i < 5; i++) {
      tmp[i] = getFHelix_pid(i,pid);
   }
   return tmp;
}

//--------------------------------------------------------------------
HepSymMatrix BMdcKalTrack::getFError(const int pid) const
//--------------------------------------------------------------------
{
   HepSymMatrix tmp(5);
   for(int i = 0; i < 5; i++) {
      for(int j = 0; j <= i; j++) {
         tmp[i][j] = getFError_pid(i,j,pid);
      }
   }
   return tmp;
}

// ther are no getPoca.. functions for RootEventData >= 6.5.5
//-----------------------------------------------------------------
// HepPoint3D BMdcKalTrack::getPoca(const int pid) const
//-----------------------------------------------------------------
// {
   // return HepPoint3D(
         // getPoca_pid(0,pid),
         // getPoca_pid(1,pid),
         // getPoca_pid(2,pid));
// }

//--------------------------------------------------------------------
Hep3Vector BMdcKalTrack::p3() const
//--------------------------------------------------------------------
{
   double pt = pxy();
   double phi0 = getZHelix_pid(1,m_pid);
   double tanl = getZHelix_pid(4,m_pid);
   return Hep3Vector(-pt*sin(phi0),pt*cos(phi0),pt*tanl);
}

//--------------------------------------------------------------------
HepPoint3D BMdcKalTrack::x3() const
//--------------------------------------------------------------------
{
   return HepPoint3D(x(),y(),z());
}

//--------------------------------------------------------------------
HepLorentzVector BMdcKalTrack::p4(double mass) const
//--------------------------------------------------------------------
{
   double pt = pxy();
   double phi0 = getZHelix_pid(1,m_pid);
   double tanl = getZHelix_pid(4,m_pid);
   double E = sqrt(pt*pt*(1. + tanl*tanl) + mass*mass);
   return HepLorentzVector(-pt*sin(phi0),pt*cos(phi0),pt*tanl,E);
}

//--------------------------------------------------------------------
void BMdcKalTrack::setZHelix(const HepVector& helix, const int pid)
//--------------------------------------------------------------------
{
   double zhelix[5] {helix[0],helix[1],helix[2],helix[3],helix[4]};
   setZHelix(zhelix,pid);
}

//--------------------------------------------------------------------
void BMdcKalTrack::setZHelix(double* helix, const int pid)
//--------------------------------------------------------------------
{
   switch( pid ) {
      case 0: m_trk->setZHelixE(helix);
              break;
      case 1: m_trk->setZHelixMu(helix);
              break;
      case 2: m_trk->setZHelix(helix);
              break;
      case 3: m_trk->setZHelixK(helix);
              break;
      case 4: m_trk->setZHelixP(helix);
              break;
   }
}

//--------------------------------------------------------------------
void BMdcKalTrack::setZError(const HepSymMatrix& error, const int pid)
//--------------------------------------------------------------------
{
   double zerror[5][5];
   for ( int i = 0; i < 5; i++ ) {
      for ( int j = 0; j < 5; j++ ) {
         zerror[i][j] = error[i][j];
      }
   }
   switch( pid ) {
      case 0: m_trk->setZErrorE(zerror);
              break;
      case 1: m_trk->setZErrorMu(zerror);
              break;
      case 2: m_trk->setZError(zerror);
              break;
      case 3: m_trk->setZErrorK(zerror);
              break;
      case 4: m_trk->setZErrorP(zerror);
              break;
   }
}

//--------------------------------------------------------------------
void BMdcKalTrack::setZError(double* error, const int pid)
//--------------------------------------------------------------------
{
   double zerror[5][5];
   int k=0;
   for ( int i = 0; i < 5 ; i++) {
      for( int j = 0; j <= i; j++,k++) {
         zerror[i][j] = error[k];
         zerror[j][i] = error[k];
      }
   }
   switch( pid ) {
      case 0: m_trk->setZErrorE(zerror);
              break;
      case 1: m_trk->setZErrorMu(zerror);
              break;
      case 2: m_trk->setZError(zerror);
              break;
      case 3: m_trk->setZErrorK(zerror);
              break;
      case 4: m_trk->setZErrorP(zerror);
              break;
   }
}

//--------------------------------------------------------------------
void BMdcKalTrack::setFHelix(const HepVector& helix, const int pid)
//--------------------------------------------------------------------
{
   double fhelix[5] {helix[0],helix[1],helix[2],helix[3],helix[4]};
   setFHelix(fhelix,pid);
}

//--------------------------------------------------------------------
void BMdcKalTrack::setFHelix(double* helix, const int pid)
//--------------------------------------------------------------------
{
   switch( pid ) {
      case 0: m_trk->setFHelixE(helix);
              break;
      case 1: m_trk->setFHelixMu(helix);
              break;
      case 2: m_trk->setFHelix(helix);
              break;
      case 3: m_trk->setFHelixK(helix);
              break;
      case 4: m_trk->setFHelixP(helix);
              break;
   }
}

//--------------------------------------------------------------------
void BMdcKalTrack::setFError(const HepSymMatrix& error, const int pid)
//--------------------------------------------------------------------
{
   double ferror[5][5];
   for ( int i = 0; i < 5; i++ ) {
      for ( int j = 0; j < 5; j++ ) {
         ferror[i][j] = error[i][j];
      }
   }
   switch( pid ) {
      case 0: m_trk->setFErrorE(ferror);
              break;
      case 1: m_trk->setFErrorMu(ferror);
              break;
      case 2: m_trk->setFError(ferror);
              break;
      case 3: m_trk->setFErrorK(ferror);
              break;
      case 4: m_trk->setFErrorP(ferror);
              break;
   }
}

//--------------------------------------------------------------------
void BMdcKalTrack::setFError(double* error, const int pid)
//--------------------------------------------------------------------
{
   double ferror[5][5];
   int k=0;
   for ( int i = 0; i < 5 ; i++) {
      for( int j = 0; j <= i; j++,k++) {
         ferror[i][j] = error[k];
         ferror[j][i] = error[k];
      }
   }
   switch( pid ) {
      case 0: m_trk->setFErrorE(ferror);
              break;
      case 1: m_trk->setFErrorMu(ferror);
              break;
      case 2: m_trk->setFError(ferror);
              break;
      case 3: m_trk->setFErrorK(ferror);
              break;
      case 4: m_trk->setFErrorP(ferror);
              break;
   }
}
