//////////////////////////////////////////////////////////////////////
//
// MagField
//
// This is a test example of using the MagneticField library functions
// It does not need event loop. Use "bean.exe -N 1 -u MagField ..."
// to test it.
//
// NOTE: You MUST specify the paths:
//       1) to the directory with magnetic fields tables
//       2) to the database directory (see TestDb)
//       Use the function AbsPath(file) which returns absolute path
//       in the form: "path_to_base_Bean_dir/" + file
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cstdlib>

#include <TH1.h>
#include <TH2.h>

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Geometry/Point3D.h"

typedef HepGeom::Point3D<double> HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "ReadDst.h"
#include "MagneticField/MagneticFieldSvc.h"

using namespace std;
using namespace CLHEP;

#ifdef __cplusplus
extern "C" {
#endif

static std::vector<TH1D*> his1;
static std::vector<TH2D*> his2;

static MagneticFieldSvc* mag_field;
static double max_x = 1.8*m, min_x = -1.8*m;
static double max_y = 1.8*m, min_y = -1.8*m;
static double max_z = 2.1*m, min_z = -2.1*m;


//--------------------------------------------------------------------
BeanUserShared_EXPORT
void MagFieldStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " MagFieldStartJob() " << endl;

   // create MagneticFieldSvc
   mag_field = MagneticFieldSvc::instance();
   if( !(mag_field->GetPath()).empty() ) {
      cerr << " MagFieldStartJob::WARNING:"
         << " MagneticFieldSvc has already been initialized" << endl
         << "                         path = "
         << mag_field->GetPath() <<endl;
   }

   // set paths to directory with magnetic fields tables
   // and to database directory
   mag_field->SetPaths( selector->AbsPath("Analysis/MagneticField"),
         selector->AbsPath("Analysis/DatabaseSvc/dat") );

   if( !mag_field->initialize() ) {
      cout << " Can not initialize MagneticField. Stop." << endl;
      exit(1);
   }

   // book histograms
   his1.resize(100,(TH1D*)0);
   his2.resize(100,(TH2D*)0);

   his1[1] = new TH1D("Bx_z","Bx(z) x=0 y=0", 420, min_z, max_z);
   his1[2] = new TH1D("By_z","By(z) x=0 y=0", 420, min_z, max_z);
   his1[3] = new TH1D("Bz_z","Bz(z) x=0 y=0", 420, min_z, max_z);
   his1[4] = new TH1D("Bx_z1","Bx(z) x=0.25 y=0", 420, min_z, max_z);
   his1[5] = new TH1D("By_z1","By(z) x=0.25 y=0", 420, min_z, max_z);
   his1[6] = new TH1D("Bz_z1","Bz(z) x=0.25 y=0", 420, min_z, max_z);
   his1[7] = new TH1D("Bx_z2","Bx(z) x=0.25 y=0.25", 420, min_z, max_z);
   his1[8] = new TH1D("By_z2","By(z) x=0.25 y=0.25", 420, min_z, max_z);
   his1[9] = new TH1D("Bz_z2","Bz(z) x=0.25 y=0.25", 420, min_z, max_z);

   his2[1] = new TH2D("Bx_xy","Bx(x,y) z=0",72,min_x,max_x,72,min_y,max_y);
   his2[2] = new TH2D("By_xy","By(x,y) z=0",72,min_x,max_x,72,min_y,max_y);
   his2[3] = new TH2D("Bz_xy","Bz(x,y) z=0",72,min_x,max_x,72,min_y,max_y);

   // register in selector to save in given directory
   VecObj his1o(his1.begin(),his1.end());
   selector->RegInDir( his1o     ,"MagField");

   VecObj his2o(his2.begin(),his2.end());
   selector->RegInDir( his2o     ,"MagField2");
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool MagFieldEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " MagFieldEvent() " << endl;

   int run = m_TEvtHeader->getRunId();
   mag_field->handle(run); // MUST BE !

   double bx = 0,by = 0,bz = 0;

   // center of detector:
   double pz = 0;
   double dz = 0.01*m;
   for(int i = 0; i < 420; i++) {
      pz = min_z + dz*(i+0.5);
      HepPoint3D r(0.,0.,pz);
      HepVector3D b;
      mag_field->fieldVector(r,b);
      bx = b.x()/tesla;
      by = b.y()/tesla;
      bz = b.z()/tesla;
      his1[1]->Fill(pz,bx);
      his1[2]->Fill(pz,by);
      his1[3]->Fill(pz,bz);

      HepPoint3D r2(0.25,0.,pz);
      mag_field->fieldVector(r2,b);
      bx = b.x()/tesla;
      by = b.y()/tesla;
      bz = b.z()/tesla;
      his1[4]->Fill(pz,bx);
      his1[5]->Fill(pz,by);
      his1[6]->Fill(pz,bz);

      HepPoint3D r3(0.25,0.25,pz);
      mag_field->fieldVector(r3,b);
      bx = b.x()/tesla;
      by = b.y()/tesla;
      bz = b.z()/tesla;
      his1[7]->Fill(pz,bx);
      his1[8]->Fill(pz,by);
      his1[9]->Fill(pz,bz);
   }

   pz = 0;
   double px = 0, py = 0;
   double dx = 0.05*m;
   double dy = 0.05*m;
   for(int i = 0; i < 72; i++) {
      px = min_x + dx*(i+0.5);
      for(int j = 0; j < 72; j++) {
         py = min_y + dy*(j+0.5);
         HepPoint3D r(px,py,pz);
         HepVector3D b;
         mag_field->fieldVector(r,b);
         bx = b.x()/tesla;
         by = b.y()/tesla;
         bz = b.z()/tesla;
         his2[1]->Fill(px,py,bx);
         his2[2]->Fill(px,py,by);
         his2[3]->Fill(px,py,bz);
//          if( abs(i-36) < 3 && abs(j-36) < 3 ) {
//            cout << i << ":" << j
//                 << " (x,y,z)= " << r << " vecB= " << b << endl;
//          }
      }
   }

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void MagFieldEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " MagFieldEndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
