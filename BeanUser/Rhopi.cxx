//////////////////////////////////////////////////////////////////////
//
// Rhopi
//
// This is Rhopi programm from Analysis/Physics/RhopiAlg
// modified (accommodated) for BEAN.
//
// Main points:
//      1) we use root-ntuples here
//      2) we create and give parameters of magnetic fild in StartJob
//      3) we initialize Absorption Correction alg. in StartJob
//         and use it in Event functions, see TestAbsCor.cxx
//      4) we use ParticleID, see TestPID.cxx for details
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <cmath>
#include <vector>

#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TMath.h>

#include "CLHEP/Units/PhysicalConstants.h"
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/TwoVector.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Geometry/Point3D.h>
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using CLHEP::pi;

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "AbsCor/AbsCor.h"
#include "VertexFit/VertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "MagneticField/MagneticFieldSvc.h"

#include "DstEvtRecTracks.h"
#include "TofHitStatus.h"
#include "ReadDst.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif


const double mpi = 0.13957;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double velc = 299.792458;   // tof path unit in mm

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

static int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6;

static std::vector<TNtuple* > m_tuple;

//ReadBeamParFromDb m_reader;
// Declare r0, z0 cut for charged tracks
static double m_vr0cut;
static double m_vz0cut;

//Declare energy, dphi, dthe cuts for fake gamma's
static double m_energyThreshold;
static double m_gammaPhiCut;
static double m_gammaThetaCut;
static double m_gammaAngleCut;

static int m_test4C;
static int m_test5C;

static int m_checkDedx;
static int m_checkTof;

static AbsCor* m_abscor = 0;

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void RhopiStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " RhopiStartJob() " << endl;

  // init Absorption Correction
  m_abscor = new AbsCor(selector->AbsPath("Analysis/AbsCor"));

  // Initial values
  m_vr0cut=1.0;
  m_vz0cut=5.0;

  m_energyThreshold=0.04;
  m_gammaPhiCut=20.0;
  m_gammaThetaCut=20.0;
  m_gammaAngleCut=20.0;

  m_test4C = 1;
  m_test5C = 1;

  m_checkDedx = 1;
  m_checkTof = 1;

  if( selector->Verbose() ) {
    cout << endl;
    cout << " Rhopi internal parameters: " << endl;
    cout << " Rhopi.Vr0cut =            " << m_vr0cut << endl;
    cout << " Rhopi.Vz0cut =            " << m_vz0cut << endl;
    cout << " Rhopi.EnergyThreshold =   " << m_energyThreshold << endl;
    cout << " Rhopi.GammaPhiCut =       " << m_gammaPhiCut << endl;
    cout << " Rhopi.GammaThetaCut =     " << m_gammaThetaCut << endl;
    cout << " Rhopi.GammaAngleCut =     " << m_gammaAngleCut << endl;
    cout << " Rhopi.Test4C =            " << m_test4C << endl;
    cout << " Rhopi.Test5C =            " << m_test5C << endl;
    cout << " Rhopi.CheckDedx =         " << m_checkDedx << endl;
    cout << " Rhopi.CheckTof =          " << m_checkTof << endl;
    cout << endl;
  }

  // We have to initialize DatabaseSvc
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  DatabaseSvc* dbs = DatabaseSvc::instance();
  if( (dbs->GetDBFilePath()).empty() ) {
     // set path to directory with databases:
     dbs->SetDBFilePath(selector->AbsPath("Analysis/DatabaseSvc/dat"));
  }

  // We have to initialize Magnetic field
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  MagneticFieldSvc* mf = MagneticFieldSvc::instance();
  if( !(mf->GetPath()).empty() ) {
    cerr << " RhopiStartJob::WARNING:"
         << " MagneticFieldSvc has already been initialized" << endl
         << "                         path = " << mf->GetPath() << endl;
  }
  // set path to directory with magnetic fields tables
  mf->SetPath( selector->AbsPath("Analysis/MagneticField") );
  mf->UseDBFlag(false); // like in the boss program
  mf->RunMode(3); // like in the boss program


  // Book n-tuples
  m_tuple.resize(15,0);
  m_tuple[1] = new TNtuple("vxyz","vxyz", "vx0:vy0:vz0:vr0:rvxy0:rvz0:rvphi0");

  m_tuple[2] = new TNtuple("photon","photon", "dthe:dphi:dang:eraw");

  m_tuple[3] = new TNtuple("etot","etot", "m2gg:etot");

  if(m_test4C==1) {
    m_tuple[4] = new TNtuple("fit4c","fit4c", "chi2:mpi0");
  } // test 4C


  if(m_test5C==1) {
    m_tuple[5] = new TNtuple("fit5c","fit5c", "chi2:mrh0:mrhp:mrhm");

    m_tuple[6] = new TNtuple("geff","geff", "fcos:elow");
  } // test 5c

  if(m_checkDedx == 1) {
    m_tuple[7] = new TNtuple("dedx","dedx",
                   "ptrk:chie:chimu:chipi:chik:chip:probPH:normPH:ghit:thit");
  } // check dE/dx

  if(m_checkTof == 1) {
    m_tuple[8] = new TNtuple("tofe","tofe",
                           "ptrk:cntr:ph:rhit:qual:te:tmu:tpi:tk:tp");
  } // check Tof:endcap


  if(m_checkTof == 1) {
    m_tuple[9] = new TNtuple("tof1","tof1",
                           "ptrk:cntr:ph:zhit:qual:te:tmu:tpi:tk:tp");
  } // check Tof:barrel inner Tof


  if(m_checkTof == 1) {
    m_tuple[10] = new TNtuple("tof2","tof2",
                            "ptrk:cntr:ph:zhit:qual:te:tmu:tpi:tk:tp");
  } // check Tof:barrel outter Tof


  m_tuple[11] = new TNtuple("pid","pid", "ptrk:cost:dedx:tof1:tof2:prob");

  // register in selector to save in given directory
  VecObj ntuples(m_tuple.begin(),m_tuple.end());
  selector->RegInDir(ntuples,"Rhopi");
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
bool RhopiEvent(ReadDst* selector,
                   TEvtHeader* m_TEvtHeader,
                   TDstEvent* m_TDstEvent,
                   TEvtRecObject* m_TEvtRecObject,
                   TMcEvent* m_TMcEvent,
                   TTrigEvent* m_TTrigEvent,
                   TDigiEvent* m_TDigiEvent,
                   THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " RhopiEvent() " << endl;

  m_abscor->AbsorptionCorrection(selector);

  // eventHeader <=> m_TEvtHeader;
  int event = m_TEvtHeader->getEventId();
  int runNo = m_TEvtHeader->getRunId();
  if( selector->Verbose() )
    cout << "run, evtnum = " << runNo << " , " << event <<endl;

  //~ cout<<"event "<<event<<endl;
  Ncut0++;

  //evtRecEvent
  const TEvtRecEvent*  evtRecEvent = m_TEvtRecObject->getEvtRecEvent();
  if( selector->Verbose() ) {
    cout <<"ncharg, nneu, tottks = "
      << evtRecEvent->totalCharged() << " , "
      << evtRecEvent->totalNeutral() << " , "
      << evtRecEvent->totalTracks() <<endl;
  }

  // evtRecTrkCol
  const TObjArray* evtRecTrkCol = selector->GetEvtRecTrkCol();

  //
  // check x0, y0, z0, r0
  // suggest cut: |z0|<5 && r0<1
  //
  Vint iGood, ipip, ipim;
  iGood.clear();
  ipip.clear();
  ipim.clear();
  Vp4 ppip, ppim;
  ppip.clear();
  ppim.clear();

  int nCharge = 0;

  // this is output of boss-6.5.2 Rhopi for sim/rec jobOpt.
  Hep3Vector xorigin(0.192687,-0.306578,-0.096103);

  VertexDbSvc* vtxsvc = VertexDbSvc::instance();
  // vtxsvc->SetBossVer("6.5.5");
  vtxsvc->SetBossVer("6.6.3");
  vtxsvc->handle(runNo); // MUST BE !
  if(vtxsvc->isVertexValid()){
    double* dbv = vtxsvc->PrimaryVertex();
    // double*  vv = vtxsvc->SigmaPrimaryVertex();
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
    // cout << " sqlite-db vertex: (x,y,z)= " << xorigin << endl;
  }

  for(int i = 0; i < evtRecEvent->totalCharged(); i++){
//    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if( !itTrk->isMdcTrackValid() ) continue;
    RecMdcTrack *mdcTrk = itTrk->mdcTrack();
//     double pch=mdcTrk->p();
    double x0=mdcTrk->x();
    double y0=mdcTrk->y();
    double z0=mdcTrk->z();
    double phi0=mdcTrk->helix(1);
    double xv=xorigin.x();
    double yv=xorigin.y();
    double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);

    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
    VFHelix helixip(point0,a,Ea);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
    double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
    double  Rvphi0=vecipa[1];

    m_tuple[1]->Fill(x0,y0,z0,Rxy,Rvxy0,Rvz0,Rvphi0);
//    if(fabs(z0) >= m_vz0cut) continue;
//    if(fabs(Rxy) >= m_vr0cut) continue;

    if(fabs(Rvz0) >= 10.0) continue;
    if(fabs(Rvxy0) >= 1.0) continue;

    iGood.push_back(i);
    nCharge += mdcTrk->charge();
  }

  //
  // Finish Good Charged Track Selection
  //
  int nGood = iGood.size();
  if( selector->Verbose() )
    cout << "ngood, totcharge = " << nGood << " , " << nCharge << endl;
  if((nGood != 2)||(nCharge!=0)){
    return false; // skip event
  }
  Ncut1++;

  Vint iGam;
  iGam.clear();
  for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
//     EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(i);
    if(!itTrk->isEmcShowerValid()) continue;
    RecEmcShower *emcTrk = itTrk->emcShower();
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.;
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
//       EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      DstEvtRecTracks* jtTrk = (DstEvtRecTracks*) evtRecTrkCol->At(j);
      if(!jtTrk->isExtTrackValid()) continue;
      RecExtTrack *extTrk = jtTrk->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      //      double ctht = extpos.cosTheta(emcpos);
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if(angd < dang){
        dang = angd;
        dthe = thed;
        dphi = phid;
      }
    }
    if(dang>=200) continue;
    double eraw = emcTrk->energy();
    dthe = dthe * 180 / (CLHEP::pi);
    dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);

    m_tuple[2]->Fill(dthe,dphi,dang,eraw);

    if(eraw < m_energyThreshold) continue;
//    if((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
    if(fabs(dang) < m_gammaAngleCut) continue;
    //
    // good photon cut will be set here
    //
    iGam.push_back(i);
  }

  //
  // Finish Good Photon Selection
  //
  int nGam = iGam.size();

  if( selector->Verbose() ) {
    cout << "num Good Photon " << nGam  << " , "
         <<evtRecEvent->totalNeutral()<<endl;
  }
  if(nGam<2){
    return false; // skip event
  }
  Ncut2++;



  //
  //
  // check dedx infomation
  //
  //

  if(m_checkDedx == 1) {
    for(int i = 0; i < nGood; i++) {
//       EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(iGood[i]);
      if(!itTrk->isMdcTrackValid()) continue;
      if(!itTrk->isMdcDedxValid())continue;
      RecMdcTrack* mdcTrk = itTrk->mdcTrack();
      RecMdcDedx* dedxTrk = itTrk->mdcDedx();
      double m_ptrk = mdcTrk->p();

      double m_chie = dedxTrk->chiE();
      double m_chimu = dedxTrk->chiMu();
      double m_chipi = dedxTrk->chiPi();
      double m_chik = dedxTrk->chiK();
      double m_chip = dedxTrk->chiP();
      double m_ghit = dedxTrk->numGoodHits();
      double m_thit = dedxTrk->numTotalHits();
      double m_probPH = dedxTrk->probPH();
      double m_normPH = dedxTrk->normPH();

      m_tuple[7]->Fill(m_ptrk,m_chie,m_chimu,m_chipi,m_chik,m_chip,
                     m_probPH,m_normPH,m_ghit,m_thit);
    }
  }

  //
  // check TOF infomation
  //


  if(m_checkTof == 1) {
    for(int i = 0; i < nGood; i++) {
//       EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + iGood[i];
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(iGood[i]);
      if(!itTrk->isMdcTrackValid()) continue;
      if(!itTrk->isTofTrackValid()) continue;

      RecMdcTrack * mdcTrk = itTrk->mdcTrack();
//       SmartRefVector<RecTofTrack> tofTrkCol = itTrk->tofTrack();
      const std::vector<RecTofTrack* >& tofTrkCol = itTrk->tofTrack();

      double ptrk = mdcTrk->p();

//       SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
      std::vector<RecTofTrack* >::const_iterator iter_tof = tofTrkCol.begin();
      for(;iter_tof != tofTrkCol.end(); iter_tof++ ) {
        TofHitStatus *status = new TofHitStatus;
        status->setStatus((*iter_tof)->status());
        if(!(status->is_barrel())){//endcap
          if( !(status->is_counter()) ) continue; // ?
          if( status->layer()!=0 ) continue;//layer1
          double path=(*iter_tof)->path(); // ?
          double tof  = (*iter_tof)->tof();
          double ph   = (*iter_tof)->ph();
          double rhit = (*iter_tof)->zrhit();
          double qual = 0.0 + (*iter_tof)->quality();
          double cntr = 0.0 + (*iter_tof)->tofID();
          double texp[5];
          for(int j = 0; j < 5; j++) {
            double gb = ptrk/xmass[j];
            double beta = gb/sqrt(1+gb*gb);
            texp[j] = 10 * path /beta/velc;
          }
          double m_te_etof    = tof - texp[0];
          double m_tmu_etof   = tof - texp[1];
          double m_tpi_etof   = tof - texp[2];
          double m_tk_etof    = tof - texp[3];
          double m_tp_etof    = tof - texp[4];

          m_tuple[8]->Fill(ptrk,cntr,ph,rhit,qual,
                         m_te_etof,m_tmu_etof,m_tpi_etof,m_tk_etof,m_tp_etof);
        }
        else {//barrel
          if( !(status->is_counter()) ) continue; // ?
          if(status->layer()==1){ //layer1
            double path=(*iter_tof)->path(); // ?
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
              double gb = ptrk/xmass[j];
              double beta = gb/sqrt(1+gb*gb);
              texp[j] = 10 * path /beta/velc;
            }

            double m_te_btof1    = tof - texp[0];
            double m_tmu_btof1   = tof - texp[1];
            double m_tpi_btof1   = tof - texp[2];
            double m_tk_btof1    = tof - texp[3];
            double m_tp_btof1    = tof - texp[4];

            m_tuple[9]->Fill(ptrk,cntr,ph,rhit,qual,
                  m_te_btof1,m_tmu_btof1,m_tpi_btof1,m_tk_btof1,m_tp_btof1);
          }

          if(status->layer()==2){//layer2
            double path=(*iter_tof)->path(); // ?
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            for(int j = 0; j < 5; j++) {
              double gb = ptrk/xmass[j];
              double beta = gb/sqrt(1+gb*gb);
              texp[j] = 10 * path /beta/velc;
            }

            double m_te_btof2    = tof - texp[0];
            double m_tmu_btof2   = tof - texp[1];
            double m_tpi_btof2   = tof - texp[2];
            double m_tk_btof2    = tof - texp[3];
            double m_tp_btof2    = tof - texp[4];

            m_tuple[10]->Fill(ptrk,cntr,ph,rhit,qual,
                  m_te_btof2,m_tmu_btof2,m_tpi_btof2,m_tk_btof2,m_tp_btof2);
          }
        }

        delete status;
      }
    } // loop all charged track
  }  // check tof


  //
  // Assign 4-momentum to each photon
  //

  Vp4 pGam;
  pGam.clear();
  for(int i = 0; i < nGam; i++) {
//     EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
      DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(iGam[i]);
    RecEmcShower* emcTrk = itTrk->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector ptrk;
    ptrk.setPx(eraw*sin(the)*cos(phi));
    ptrk.setPy(eraw*sin(the)*sin(phi));
    ptrk.setPz(eraw*cos(the));
    ptrk.setE(eraw);

//    ptrk = ptrk.boost(-0.011,0,0);// boost to cms

    pGam.push_back(ptrk);
  }
  cout<<"before pid"<<endl;
  //
  // Assign 4-momentum to each charged track
  //
  ParticleID *pid = ParticleID::instance();
#if (BOSS_VER < 700)
  pid->set_path(selector->AbsPath("Analysis/ParticleID_boss6"));
#else
  pid->set_path(selector->AbsPath("Analysis/ParticleID"));
#endif
  for(int i = 0; i < nGood; i++) {
//     EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
    DstEvtRecTracks* itTrk = (DstEvtRecTracks*) evtRecTrkCol->At(iGood[i]);
    //    if(pid) delete pid;
    pid->init();
    pid->setMethod(pid->methodProbability());
//    pid->setMethod(pid->methodLikelihood());  //for Likelihood Method

    pid->setChiMinCut(4);
    pid->setRecTrack(itTrk);
    pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()
                                 | pid->useTofE()); // use PID sub-system
    pid->identify(pid->onlyPion() | pid->onlyKaon());    // seperater Pion/Kaon
    //    pid->identify(pid->onlyPion());
    //  pid->identify(pid->onlyKaon());
    pid->calculate(runNo);
    if(!(pid->IsPidInfoValid())) continue;
    RecMdcTrack* mdcTrk = itTrk->mdcTrack();

    double m_ptrk_pid = mdcTrk->p();
    double m_cost_pid = cos(mdcTrk->theta());
    double m_dedx_pid = pid->chiDedx(2);
    double m_tof1_pid = pid->chiTof1(2);
    double m_tof2_pid = pid->chiTof2(2);
    double m_prob_pid = pid->probPion();

    m_tuple[11]->Fill( m_ptrk_pid,m_cost_pid,m_dedx_pid,
                     m_tof1_pid,m_tof2_pid,m_prob_pid );

//  if(pid->probPion() < 0.001 || (pid->probPion() < pid->probKaon())) continue;
    if(pid->probPion() < 0.001) continue;
//    if(pid->pdf(2)<pid->pdf(3)) continue; //  for Likelihood Method(0=electron 1=muon 2=pion 3=kaon 4=proton)

    RecMdcKalTrack* mdcKalTrk = itTrk->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
//     RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion
    if( !mdcKalTrk ) {
      cout << " No Kalman track ! " << endl;
      return false;
    }
    mdcKalTrk->setPidType(RecMdcKalTrack::pion);

    if(mdcKalTrk->charge() >0) {
      ipip.push_back(iGood[i]);
      HepLorentzVector ptrk;
      ptrk.setPx(mdcKalTrk->px());
      ptrk.setPy(mdcKalTrk->py());
      ptrk.setPz(mdcKalTrk->pz());
      double p3 = ptrk.mag();
      ptrk.setE(sqrt(p3*p3+mpi*mpi));

//      ptrk = ptrk.boost(-0.011,0,0);//boost to cms

      ppip.push_back(ptrk);
    } else {
      ipim.push_back(iGood[i]);
      HepLorentzVector ptrk;
      ptrk.setPx(mdcKalTrk->px());
      ptrk.setPy(mdcKalTrk->py());
      ptrk.setPz(mdcKalTrk->pz());
      double p3 = ptrk.mag();
      ptrk.setE(sqrt(p3*p3+mpi*mpi));

//      ptrk = ptrk.boost(-0.011,0,0);//boost to cms

      ppim.push_back(ptrk);
    }
  }

/*
  for(int i = 0; i < nGood; i++) {//for rhopi without PID
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
    RecMdcTrack* mdcTrk = itTrk->mdcTrack();
    if(mdcTrk->charge() >0) {
      ipip.push_back(iGood[i]);
      HepLorentzVector ptrk;
      ptrk.setPx(mdcTrk->px());
      ptrk.setPy(mdcTrk->py());
      ptrk.setPz(mdcTrk->pz());
      double p3 = ptrk.mag();
      ptrk.setE(sqrt(p3*p3+mpi*mpi));
      ppip.push_back(ptrk);
    } else {
      ipim.push_back(iGood[i]);
      HepLorentzVector ptrk;
      ptrk.setPx(mdcTrk->px());
      ptrk.setPy(mdcTrk->py());
      ptrk.setPz(mdcTrk->pz());
      double p3 = ptrk.mag();
      ptrk.setE(sqrt(p3*p3+mpi*mpi));
      ppim.push_back(ptrk);
    }
  }// without PID
*/

  int npip = ipip.size();
  int npim = ipim.size();
  if(npip*npim != 1) return false; // skip event

  Ncut3++;


  //
  // Loop each gamma pair, check ppi0 and pTot
  //

  HepLorentzVector pTot;
  for(int i = 0; i < nGam - 1; i++){
    for(int j = i+1; j < nGam; j++) {
      HepLorentzVector p2g = pGam[i] + pGam[j];
      pTot = ppip[0] + ppim[0];
      pTot += p2g;
      double m_m2gg = p2g.m();
      double m_etot = pTot.e();
      m_tuple[3] -> Fill(m_m2gg,m_etot);

    }
  }


//   RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ipip[0]))->mdcKalTrack();
//   RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+ipim[0]))->mdcKalTrack();

  RecMdcKalTrack *pipTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(ipip[0]))
                                                             ->mdcKalTrack();
  RecMdcKalTrack *pimTrk = ((DstEvtRecTracks*)evtRecTrkCol->At(ipim[0]))
                                                             ->mdcKalTrack();

  WTrackParameter wvpipTrk, wvpimTrk;
  wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
  wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());

/* Default is pion, for other particles:
  wvppTrk = WTrackParameter(mp, pipTrk->getZHelixP(), pipTrk->getZErrorP());//proton
  wvmupTrk = WTrackParameter(mmu, pipTrk->getZHelixMu(), pipTrk->getZErrorMu());//muon
  wvepTrk = WTrackParameter(me, pipTrk->getZHelixE(), pipTrk->getZErrorE());//electron
  wvkpTrk = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK());//kaon
*/
  //
  //    Test vertex fit
  //

  HepPoint3D vx(0., 0., 0.);
  HepSymMatrix Evx(3, 0);
  double bx = 1E+6;
  double by = 1E+6;
  double bz = 1E+6;
  Evx[0][0] = bx*bx;
  Evx[1][1] = by*by;
  Evx[2][2] = bz*bz;

  VertexParameter vxpar;
  vxpar.setVx(vx);
  vxpar.setEvx(Evx);

  VertexFit* vtxfit = VertexFit::instance();
  vtxfit->init();
  vtxfit->AddTrack(0,  wvpipTrk);
  vtxfit->AddTrack(1,  wvpimTrk);
  vtxfit->AddVertex(0, vxpar,0, 1);
  if(!vtxfit->Fit(0)) return false; //skip event

   vtxfit->BuildVirtualParticle(0);

   cout << "vtxfit0->chisq(0): " << vtxfit->chisq(0) << endl;



  vtxfit->Swim(0);


  WTrackParameter wpip = vtxfit->wtrk(0);
  WTrackParameter wpim = vtxfit->wtrk(1);

  KinematicFit * kmfit = KinematicFit::instance();

  //
  //  Apply Kinematic 4C fit
  //
  cout<<"before 4c"<<endl;
  if(m_test4C==1) {
//    double ecms = 3.097;
    HepLorentzVector ecms(0.034,0,0,3.097);

    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    for(int i = 0; i < nGam-1; i++) {
//       RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
      RecEmcShower *g1Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(iGam[i]))
                                                             ->emcShower();

      for(int j = i+1; j < nGam; j++) {
//         RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
        RecEmcShower *g2Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(iGam[j]))
                                                               ->emcShower();
        kmfit->init();
        kmfit->AddTrack(0, wpip);
        kmfit->AddTrack(1, wpim);
        kmfit->AddTrack(2, 0.0, g1Trk);
        kmfit->AddTrack(3, 0.0, g2Trk);
        kmfit->AddFourMomentum(0, ecms);
        bool oksq = kmfit->Fit();
        if(oksq) {
          double chi2 = kmfit->chisq();
          if(chi2 < chisq) {
            chisq = chi2;
            ig1 = iGam[i];
            ig2 = iGam[j];
          }
        }
      }
    }

    if(chisq < 200) {

//       RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
//       RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
      RecEmcShower *g1Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(ig1))
                                                             ->emcShower();
      RecEmcShower *g2Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(ig2))
                                                             ->emcShower();

      kmfit->init();
      kmfit->AddTrack(0, wpip);
      kmfit->AddTrack(1, wpim);
      kmfit->AddTrack(2, 0.0, g1Trk);
      kmfit->AddTrack(3, 0.0, g2Trk);
      kmfit->AddFourMomentum(0, ecms);
      bool oksq = kmfit->Fit();
      if(oksq) {
        HepLorentzVector ppi0 = kmfit->pfit(2) + kmfit->pfit(3);
        double m_mpi0 = ppi0.m();
        double m_chi1 = kmfit->chisq();
        m_tuple[4]->Fill(m_chi1,m_mpi0);
        Ncut4++;
      }
    }
  }

  //
  //  Apply Kinematic 5C Fit
  //

  // find the best combination over all possible pi+ pi- gamma gamma pair
  if(m_test5C==1) {
//    double ecms = 3.097;
    HepLorentzVector ecms(0.034,0,0,3.097);
    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    for(int i = 0; i < nGam-1; i++) {
//       RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+iGam[i]))->emcShower();
      RecEmcShower *g1Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(iGam[i]))
                                                             ->emcShower();
      for(int j = i+1; j < nGam; j++) {
//         RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+iGam[j]))->emcShower();
        RecEmcShower *g2Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(iGam[j]))
                                                               ->emcShower();
        kmfit->init();
        kmfit->AddTrack(0, wpip);
        kmfit->AddTrack(1, wpim);
        kmfit->AddTrack(2, 0.0, g1Trk);
        kmfit->AddTrack(3, 0.0, g2Trk);
        kmfit->AddResonance(0, 0.135, 2, 3);
        kmfit->AddFourMomentum(1, ecms);
        if(!kmfit->Fit(0)) continue;
        if(!kmfit->Fit(1)) continue;
        bool oksq = kmfit->Fit();
        if(oksq) {
          double chi2 = kmfit->chisq();
          if(chi2 < chisq) {
            chisq = chi2;
            ig1 = iGam[i];
            ig2 = iGam[j];
          }
        }
      }
    }


    cout << " chisq = " << chisq <<endl;

    if(chisq < 200) {
//       RecEmcShower* g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
//       RecEmcShower* g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
      RecEmcShower *g1Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(ig1))
                                                             ->emcShower();
      RecEmcShower *g2Trk = ((DstEvtRecTracks*)evtRecTrkCol->At(ig2))
                                                             ->emcShower();

      kmfit->init();
      kmfit->AddTrack(0, wpip);
      kmfit->AddTrack(1, wpim);
      kmfit->AddTrack(2, 0.0, g1Trk);
      kmfit->AddTrack(3, 0.0, g2Trk);
      kmfit->AddResonance(0, 0.135, 2, 3);
      kmfit->AddFourMomentum(1, ecms);
      bool oksq = kmfit->Fit();
      if(oksq){
        HepLorentzVector ppi0 = kmfit->pfit(2) + kmfit->pfit(3);
        HepLorentzVector prho0 = kmfit->pfit(0) + kmfit->pfit(1);
        HepLorentzVector prhop = ppi0 + kmfit->pfit(0);
        HepLorentzVector prhom = ppi0 + kmfit->pfit(1);

        double m_chi2  = kmfit->chisq();
        double m_mrh0 = prho0.m();
        double m_mrhp = prhop.m();
        double m_mrhm = prhom.m();
        double eg1 = (kmfit->pfit(2)).e();
        double eg2 = (kmfit->pfit(3)).e();
        double fcos = abs(eg1-eg2)/ppi0.rho();
        m_tuple[5]->Fill(m_chi2,m_mrh0,m_mrhp,m_mrhm);
        Ncut5++;
        //
        //  Measure the photon detection efficiences via
        //          J/psi -> rho0 pi0
        //
        if(fabs(prho0.m()-0.770)<0.150) {
          if(fabs(fcos)<0.99) {
            double m_fcos = (eg1-eg2)/ppi0.rho();
            double m_elow =  (eg1 < eg2) ? eg1 : eg2;
            m_tuple[6]->Fill(m_fcos,m_elow);
            Ncut6++;
          }
        } // rho0 cut
      }  //oksq
    }
  }

  return false; // skip event
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void RhopiEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
  if( selector->Verbose() ) cout << " RhopiEndJob() " << endl;

  cout<<"total number:         "<<Ncut0<<endl;
  cout<<"nGood==2, nCharge==0: "<<Ncut1<<endl;
  cout<<"nGam>=2:              "<<Ncut2<<endl;
  cout<<"Pass Pid:             "<<Ncut3<<endl;
  cout<<"Pass 4C:              "<<Ncut4<<endl;
  cout<<"Pass 5C:              "<<Ncut5<<endl;
  cout<<"J/psi->rho0 pi0:      "<<Ncut6<<endl;

  cout << "in finalize()" << endl;
}

#ifdef __cplusplus
}
#endif
