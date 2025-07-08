//--------------------------------------------------------------------
// AbsCor: <=> boss/Analysis/PhotonCor/AbsCor
//
//    Nefedov: The original program from the BOSS can not be used
//    directly in the BEAN:
//    1) Corrections 'edgecor' use functions from the
//       detector description package (see EmcRecGeoSvc)
//       This is definitely not for BEAN.
//    2) BEAN does not use the database, but read files whose
//       names can be specified in functions: SetCorFunparaPath()
//       and SetDataPathc3ptof()
//    3) The number of modifications to be applied exceeds
//       reasonable limits. The program becomes unreadable.
//
//--------------------------------------------------------------------

#include "AbsCor/AbsCor.h"   // must be first!

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"

#include "DstEvtRecTracks.h"

#if (BOSS_VER > 700)
#include "Identifier/EmcID.h"
#endif

#if   (BOSS_VER >= 661 && BOSS_VER <= 665)
#define DATAPATH "00-00-28"
#elif (BOSS_VER >= 702 && BOSS_VER <= 707)
#define DATAPATH "00-00-36"
#else
// TODO: it seems like DATAPATH will never change again
#define DATAPATH "00-00-41"
#endif

using namespace std;

//--------------------------------------------------------------------
AbsCor::AbsCor(const string& path, bool _usetof, bool _dodatacor) {
//--------------------------------------------------------------------
   usetof = _usetof;
   dodatacor = _dodatacor;

   // default values (taken from BOSS);
   hotcellmask = false;
   dopi0Cor = true;
   MCuseTof = true;

   MCCorUseFunction = true;

   dt = nullptr;
   dtErr = nullptr;

   //-----------------------------------------------------------------
   // read the parameter of pi0 calibration for data
   //-----------------------------------------------------------------
   string DataPathc3p = path;
   if( !usetof ) {
      DataPathc3p += "/dat/" DATAPATH "/c3p.txt";
   } else {
      DataPathc3p += "/dat/" DATAPATH "/c3ptof.txt";
   }
   ReadDatac3p(DataPathc3p);

   // Single gamma corrections:
#if (BOSS_VER >= 700)
   //-----------------------------------------------------------------
   // Read energy correction parameters from PhotonCor/McCor
   //-----------------------------------------------------------------
   string paraPath = path;
   if ( !MCuseTof ) {
      paraPath += "/dat/" DATAPATH "/evset.txt";
   } else {
      paraPath += "/dat/" DATAPATH "/evsetTof.txt";
   }
   ReadParMcCor(paraPath);
#endif

#if (BOSS_VER >= 708)
   //-----------------------------------------------------------------
   // Read energy correction Function parameters
   //-----------------------------------------------------------------
   string CorFunparaPath = path;
   if ( !MCuseTof ) {
      CorFunparaPath += "/dat/" DATAPATH "/evsetCorFunctionPar.txt";
   } else {
      CorFunparaPath += "/dat/" DATAPATH "/evsetTofCorFunctionPar.txt";
   }
   ReadCorFunpara(CorFunparaPath);

   string info709 = R"(
AbsCor-Info: if DST was made by BOSS >= 7.0.8:
AbsCor in Bean does not use the database, but you can read files with
calibrations in functions:
      void ReadDatac3p(std::string DataPathc3p, bool info=true);
      void ReadCorFunpara(std::string CorFunparaPath, bool info=true);
An example of using these functions is contained in the file:
BeanUser/TestAbsCor.cxx

The path to the calibration files (on /cvmfs) can be obtained using
the python3 program 'GetAbsCorFiles.py', see folder bean/scripts/
)";
   cout << info709 << endl;
#endif

#if BOSS_VER > 713
   // message for untested version of BOSS
   cout << "AbsCor-WARNING: untested BOSS= " << BOSS_VER << endl
      << "  check AbsCor correction files in\n"
      << "  /cvmfs/bes3.ihep.ac.cn/CalibConst/emc/ShEnCalib/\n"
      << endl;
#endif

   //-----------------------------------------------------------------
   // Suppression of hot crystals
   // Reading the map from hotcry.txt (Hajime, Jul 2013)
   for(int ih=0; ih<10; ih++) {
      hrunstart[ih]=-1;
      hrunend[ih]=-1;
      hcell[ih]=-1;
   }
   int numhots=4; // numbers of hot crystals
   int dumring,dumphi,dummod,dumid;
   string HotList = path + "/dat/hotcry.txt";

   ifstream hotcrys;
   hotcrys.exceptions( ifstream::failbit | ifstream::badbit );
   hotcrys.open(HotList.c_str(),ios::in);
   for(int il=0; il<numhots; il++) {
      hotcrys>>hrunstart[il];
      hotcrys>>hrunend[il];
      hotcrys>>hcell[il];
      hotcrys>>dumring;
      hotcrys>>dumphi;
      hotcrys>>dummod;
      hotcrys>>dumid;
   }
   hotcrys.close();
}

//--------------------------------------------------------------------
void AbsCor::ReadDatac3p(std::string DataPathc3p, bool info) {
//--------------------------------------------------------------------
   ifstream inc3p(DataPathc3p.c_str(),ios::in);
   if ( !inc3p.is_open() ) {
      cout << " can not open: " << DataPathc3p << endl;
      exit(EXIT_FAILURE);
   }

   for(int i = 0; i < 4; i++) {
      double am,ame;
      inc3p >> am;
      inc3p >> ame;
      ai[i] = am;
      // cerr << " ai[" << i << "]= " << am << endl;
   }

   if( inc3p.fail() ) {
      cout << " error while reading file " << DataPathc3p << endl;
      exit(EXIT_FAILURE);
   }
   inc3p.close();

   if ( info ) {
      cout << __func__ << "(): successfully read file:" << endl
         << DataPathc3p << endl;
   }
}

//--------------------------------------------------------------------
void AbsCor::ReadParMcCor(std::string paraPath, bool info) {
//--------------------------------------------------------------------
   ifstream in2(paraPath.c_str(),ios::in);
   if ( !in2.is_open() ) {
      cout << " can not open: " << paraPath << endl;
      exit(EXIT_FAILURE);
   }

   double energy,thetaid,peak1,peakerr1,res,reserr;
   dt = new TGraph2DErrors();
   dtErr = new TGraph2DErrors();
   //for(int i=0;i<560;i++){
   for(int i=0; i<1484; i++) { //53*28
      in2>>energy;
      in2>>thetaid;
      in2>>peak1;
      in2>>peakerr1;
      in2>>res;
      in2>>reserr;
      dt->SetPoint(i,energy,thetaid,peak1);
      dt->SetPointError(i,0,0,peakerr1);
      dtErr->SetPoint(i,energy,thetaid,res);
      dtErr->SetPointError(i,0,0,reserr);
      if(i<28) {
         e25min[int(thetaid)]=energy;
      }
      if(i>=1484-28) {
         e25max[int(thetaid)]=energy;
      }
      // if(i>=560-28) e25max[int(thetaid)]=energy;
   }

   if( in2.fail() ) {
      cout << " error while reading file " << paraPath << endl;
      exit(EXIT_FAILURE);
   }
   in2.close();

   if ( info ) {
      cout << __func__ << "(): successfully read file:" << endl
         << paraPath << endl;
   }
}

//--------------------------------------------------------------------
void AbsCor::ReadCorFunpara(std::string CorFunparaPath, bool info) {
//--------------------------------------------------------------------
   ifstream in2corfun(CorFunparaPath.c_str(),ios::in);
   if ( !in2corfun.is_open() ) {
      cout << " can not open: " << CorFunparaPath << endl;
      exit(EXIT_FAILURE);
   }

   for(int i=0; i<28; i++) {
      in2corfun>>m_corFunPar[i][0];
      in2corfun>>m_corFunPar[i][1];
      in2corfun>>m_corFunPar[i][2];
      in2corfun>>m_corFunPar[i][3];
      in2corfun>>m_corFunPar[i][4];
      in2corfun>>m_corFunPar[i][5];
   }

   if( in2corfun.fail() ) {
      cout << " error while reading file " << CorFunparaPath << endl;
      exit(EXIT_FAILURE);
   }
   in2corfun.close();

   if ( info ) {
      cout << __func__ << "(): successfully read file:" << endl
         << CorFunparaPath << endl;
   }
}

// Nefedov: I keep this function for back compatibility
//          see BeanUser/ DaubleDTag & RadBhabha
//--------------------------------------------------------------------
void AbsCor::SuppressHotCrystals(ReadDst* selector) {
//--------------------------------------------------------------------
   if( selector->Verbose() ) {
      cout << " SuppressHotCrystals() " << endl;
   }

   int runNo = const_cast<TEvtHeader* >
         (selector->GetEvtHeader())->getRunId();

   const TEvtRecEvent* evtRecEvent =
         selector->GetEvtRecObject()->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   if( evtRecEvent->totalTracks() > evtRecTrkCol->GetSize() ) {
      return;
   }
   if( evtRecEvent->totalTracks() > 50 ) {
      return;
   }

   for(int i = 0; i < evtRecEvent->totalTracks(); i++) {
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
      if( !itTrk->isEmcShowerValid() ) {
         continue;
      }

      RecEmcShower* emcTrk = itTrk->emcShower();

      // If it is "hot", return "9999" (Hajime, Jul 2013)

      for ( int ih = 0; ih < 10; ih++ ) {
         if ( (hrunstart[ih] == -1) ||
               (hrunend[ih] == -1) ||
               (hcell[ih] == -1) ) {
            continue;
         }
         if ( (abs(runNo) < hrunstart[ih]) ||
               (abs(runNo) > hrunend[ih]) ) {
            continue;
         }

         if ( emcTrk->cellId() == hcell[ih] ) {
            emcTrk->setStatus(9999);

            if (selector->Verbose()) {
               cout << " AbsCor::SuppressHotCrystals():"
                  " suppressed cell=" << emcTrk->cellId()
                     << " energy=" << emcTrk->energy() << endl;
            }
         }
      }
   }

}


//--------------------------------------------------------------------
void AbsCor::AbsorptionCorrection(ReadDst* selector) {
//--------------------------------------------------------------------
   if( selector->Verbose() ) {
      cout << " AbsorptionCorrection() " << endl;
   }

   int runNo = const_cast<TEvtHeader* >
         (selector->GetEvtHeader())->getRunId();

   const TEvtRecEvent* evtRecEvent =
         selector->GetEvtRecObject()->getEvtRecEvent();
   const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

   if( evtRecEvent->totalTracks() > evtRecTrkCol->GetSize() ) {
      return;
   }
   if( evtRecEvent->totalTracks() > 50 ) {
      return;
   }

   for(int i = 0; i < evtRecEvent->totalTracks(); i++) {
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
      if( !itTrk->isEmcShowerValid() ) {
         continue;
      }

      RecEmcShower* emcTrk = itTrk->emcShower();

#if (BOSS_VER < 700)
      // this piece of code prevents that AbsCor corrections
      // from being applied twice
      int st = emcTrk->status();
      if( selector->Verbose() ) {
         cout << " AbsCor:: EMC status= " << emcTrk->status();
      }
      if( st > 10000 ) {
         if( selector->Verbose() ) {
            cout << " skip ... " << endl;
         }
         continue;
      }
      emcTrk->setStatus(st+20000);
      if( selector->Verbose() ) {
         cout << " set to " << emcTrk->status() << endl;
      }

      if( emcTrk->e5x5() < 0.015 ) {
         continue;
      }
#endif

      double etof=0;
      if( usetof && itTrk->isTofTrackValid() ) {
         const std::vector<RecTofTrack* >&
            recTofTrackVec = itTrk->tofTrack();
         if( !recTofTrackVec.empty() ) {
            etof = recTofTrackVec[0]->energy();
         }
         if( etof>100. ) {
            etof = 0;
         }
      }

      if( selector->Verbose() ) {
         cout << " AbsCor:: etof = " << etof << endl;
      }

      double energyC;
#if (BOSS_VER < 700)
      energyC = emcTrk->energy() + etof;

#else
      double e5x5=emcTrk->e5x5();

      Identifier id(emcTrk->cellId());

      unsigned int module = EmcID::barrel_ec(id);
//     unsigned int ntheta = EmcID::theta_module(id);
//     unsigned int nphi = EmcID::phi_module(id);

      // id=EmcID::crystal_id(module,ntheta,nphi);

      unsigned int thetaModule = EmcID::theta_module(id);
//     unsigned int phiModule = EmcID::phi_module(id);

      int thetaId = thetaModule;
      if (module==0||module==2) {
         thetaId = thetaModule;
      }
      if (module==1 && thetaModule<=21) {
         thetaId = thetaModule + 6;
      }
      if (module==1 && thetaModule>21) {
         thetaId = 43 - thetaModule + 6;
      }

#if (BOSS_VER <= 707)
      if ( MCuseTof ) {
         energyC=ECorrMC(e5x5+etof,thetaId);
      } else {
         energyC=ECorrMC(e5x5,thetaId);
      }
#else
      double DthetaId = double(thetaId);
      if( MCuseTof ) {
         if ( thetaId < 6 ) {
            etof = 0.0;   // in EMC endcap
         }
         if ( MCCorUseFunction ) {
            energyC = ECorrFunctionMC(e5x5+etof,DthetaId);
         } else {
            energyC = ECorrMC(e5x5+etof,DthetaId);
         }
      } else {
         if ( MCCorUseFunction ) {
            energyC = ECorrFunctionMC(e5x5,DthetaId);
         } else {
            energyC = ECorrMC(e5x5,DthetaId);
         }
      }
#endif
#endif

      if( selector->Verbose() ) {
         cout << " AbsCor:: energyC= " << energyC << endl;
      }

      double lnEcor=1.0;
      if ( dopi0Cor ) {
         if ( runNo > 0 && dodatacor ) {
            double lnE = std::log(energyC);
#if (BOSS_VER > 700)
            if ( energyC>1.0 ) {
               lnE=std::log(1.0);
            }
#if (BOSS_VER <= 700)
            if ( energyC<0.07 ) {
               lnE=std::log(0.07);
            }
#else
            if ( energyC<0.05 ) {
               lnE=std::log(0.05);
            }
#endif
#endif
            lnEcor = ai[0] + lnE*(ai[1] + lnE*(ai[2]+lnE*ai[3]));
         }
      }

      if( lnEcor < 0.5 ) {
         continue;
      }

      // Nefedov: function SuppressHotCrystals() is available in BEAN
      //          regardless of version
      //
      // If it is "hot", return "9999" (Hajime, Jul 2013)
      if ( hotcellmask ) {
         for ( int ih = 0; ih < 10; ih++ ) {
            if ( (hrunstart[ih] == -1) ||
                  (hrunend[ih] == -1) ||
                  (hcell[ih] == -1) ) {
               continue;
            }
            if ( (abs(runNo) < hrunstart[ih]) ||
                  (abs(runNo) > hrunend[ih]) ) {
               continue;
            }

            if ( emcTrk->cellId() == hcell[ih] ) {
               emcTrk->setStatus(9999);
            }
         }
      }

      double enecor=1.;
      // remove part for "edgecor"

      double energyCC = energyC/(lnEcor*enecor);
      emcTrk->setEnergy(energyCC);

      if( selector->Verbose() ) {
         cout << " AbsCor:: energyCC = " << energyCC << endl;
      }
   }
}


//--------------------------------------------------------------------
//The following function is copied from PhotonCor/McCor
double AbsCor::ECorrMC(double eg, double theid) const
//--------------------------------------------------------------------
{
   double Energy5x5=eg;
   if(eg<E25min(int(theid))) {
      eg=E25min(int(theid));
   }
   if(eg>E25max(int(theid))) {
      eg=E25max(int(theid))-0.001;
   }

   if(theid<=0) {
      theid=0.001;
   }
   if(theid>=27) {
      theid=26.999;
   }
   Float_t einter = eg + 0.00001;
   Float_t tinter = theid+0.0001;
   //cout<<"inter="<< einter<<"   "<<tinter<<endl;
   double ecor=dt->Interpolate(einter,tinter);
   // cout<<"ecor="<<ecor<<endl;
   if(!(ecor)) {
      return Energy5x5;
   }
   if(ecor<0.5) {
      return Energy5x5;
   }
   double EnergyCor=Energy5x5/ecor;
   return EnergyCor;
}

//--------------------------------------------------------------------
// Get energy error
double AbsCor::ErrMC(double eg, double theid) const
//--------------------------------------------------------------------
{
   if(eg<E25min(int(theid))) {
      eg=E25min(int(theid));
   }
   if(eg>E25max(int(theid))) {
      eg=E25max(int(theid))-0.001;
   }
   if(theid<=0) {
      theid=0.001;
   }
   if(theid>=27) {
      theid=26.999;
   }
   Float_t einter = eg + 0.00001;
   Float_t tinter = theid+0.0001;
   double err=dtErr->Interpolate(einter,tinter);
   return err;
}

//--------------------------------------------------------------------
double AbsCor::ECorrFunctionMC(double eg, double theid) const
//--------------------------------------------------------------------
{
   if(theid<0||theid>27) {
      cout << "in AbsCor EcorrFunctionMC error::"
            "thetaId is out of the region [0,27]" << endl;
   }
   double Energy5x5=eg;
   double x=Energy5x5;
   //if(Energy5x5>E25max(int(theid))) x=E25max(int(theid));

   int ith=int(theid);
   double ecor;
   ecor = m_corFunPar[ith][0]/(m_corFunPar[ith][1]+exp(-x))
         + m_corFunPar[ith][2]
         + m_corFunPar[ith][3]*x
         + m_corFunPar[ith][4]*x*x
         + m_corFunPar[ith][5]*x*x*x;

   double EnergyCor=Energy5x5/ecor;
   return EnergyCor;
}
