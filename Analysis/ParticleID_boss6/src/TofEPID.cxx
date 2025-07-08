#include <cmath>

#include "ParticleID/TofEPID.h"

#ifndef BEAN
#include "MdcRecEvent/RecMdcTrack.h"
#include "TofRecEvent/RecTofTrack.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#else
#include "TofHitStatus.h"
#endif

TofEPID * TofEPID::m_pointer = 0;

TofEPID * TofEPID::instance() {
   if(!m_pointer) m_pointer = new TofEPID();
   return m_pointer;
}

TofEPID::TofEPID():ParticleIDBase() {
   //readSigPar();
}

void TofEPID::init() {
   for(int i = 0; i < 5; i++) {
      m_chi[i] = 99.0;
      m_prob[i] = -1.0;
      m_offset[i] = 99.0;
      m_sigma[i] = 1.0;
   }
   m_chimin = 99.;
   m_pdfmin =99.;
   m_ndof = 0;
   m_mass2 = -999;
   m_rhit = -99;
}

void TofEPID::calculate() {
   if(particleIDCalculation()==0) m_ndof=1;
}

int TofEPID::particleIDCalculation() {
   int irc = -1;
   EvtRecTrack* recTrk = PidTrk();
   if(!(recTrk->isMdcTrackValid())) return irc;
   RecMdcTrack* mdcTrk = recTrk->mdcTrack();

   double ptrk = mdcTrk->p();
   //   double cost = cos(mdcTrk->theta());
   //   double charge = mdcTrk->charge();

   if(!(recTrk->isTofTrackValid())) return irc;

#ifndef BEAN
   SmartRefVector<RecTofTrack> tofTrk = recTrk->tofTrack();
   SmartRefVector<RecTofTrack>::iterator it;//=tofTrk.begin();
#else
   const std::vector<RecTofTrack* >& tofTrk = recTrk->tofTrack();
   std::vector<RecTofTrack* >::const_iterator it;//=tofTrk.begin();
#endif

   TofHitStatus *hitst = new TofHitStatus;
   std::vector<int> tofecount;
   int goodtofetrk=0;
   for(it = tofTrk.begin(); it!=tofTrk.end(); it++,goodtofetrk++) {
      unsigned int st = (*it)->status();
      hitst->setStatus(st);
      if(  (hitst->is_barrel()) ) continue;
      if( !(hitst->is_counter()) ) continue;
      if( hitst->layer()==1 )  tofecount.push_back(goodtofetrk);
   }
   delete hitst;
   if(tofecount.size()!=1) return irc;//not tof2 track or more than 1 tracks
   it = tofTrk.begin()+tofecount[0];


   //   int qual = (*it)->quality();
   //  if(qual != 1) return irc;
   //   int cntr = (*it)->tofID();
   double tof  = (*it)->tof();
   if(tof <=0 ) return irc;
   double path = (*it)->path();
   //   double m_ph   = (*it)->ph();
   m_rhit = (*it)->zrhit();
   // m_part = tofTrk->getPart();

   double beta2 = path*path/velc()/velc()/tof/tof;
   m_mass2 = ptrk * ptrk * (1/beta2 -1);
   // if ((m_mass2>20)||(m_mass2<-1)) return irc;


   double chitemp = 99.;
   double pdftemp = 0;
   for(int i = 0; i < 5; i++) {
      //     double gb = ptrk/xmass(i);
      //     double beta = gb/sqrt(1+gb*gb);
      double texp = (*it)->texp(i);
      //      path /beta/velc();
      m_offset[i] = tof - texp-(*it)->toffset(i);// - offsetTofE(i, cntr, ptrk, m_rhit, m_ph,charge);
      double sigma_tmp= (*it)->sigma(0);
      /* if(fabs(sigma_tmp/1000.)>0.3) {
           int tofid=(*it)->tofID();
       //       sigma_tmp=getSigbyRun(tofid)*0.1;
        }*/

      if(sigma_tmp!=0) {
         m_sigma[i] = 1.2*sigma_tmp;
         if(i<2) m_sigma[i]=sigma_tmp;
      }
      else m_sigma[i]=0.12;
      //    m_chi[i] = m_offset[i]/m_sigma[i];

      //    m_sigma[i]=sigma_tmp;
      //   if(i!=0) m_sigma[i] = sigma_tmp*1.1;
      //    m_sigma[i] = sigmaTofE(i, cntr,ptrk,m_rhit, m_ph,charge);
      m_chi[i] = m_offset[i]/m_sigma[i];
      if(fabs(m_chi[i]) < chitemp) chitemp = fabs(m_chi[i]);
      double ppp = pdfCalculate(m_chi[i],1);
      if(fabs(ppp) > pdftemp) pdftemp = fabs(ppp);
   }
   m_chimin = chitemp;
   // if(m_chimin > chiMinCut() ) return irc;
   // if(pdftemp < pdfCalculate(pdfMinSigmaCut(),1.0)) return irc;

   // calculate prob

   for(int i = 0; i < 5; i++)
      m_prob[i] = probCalculate(m_chi[i]*m_chi[i], 1);

   m_ndof = 1;
   irc = 0;
   return irc;
}


//
//  TOF endcap: Correction routines
//

double TofEPID::offsetTofE(int n, int cntr, double ptrk, double rtof, double ph,double charge) {
   double offset;
   //   double gb = ptrk/xmass(n);
   switch(n) {
   case 0: {  // Electron
      double ptemp = ptrk;
      if(ptrk<0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      offset = 0.001*(-28.8481+138.159*plog-249.334*plog*plog);
      break;
   }

   case 1: { // Muon
      double ptemp = ptrk;
      if(ptrk<0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      offset = 0.001*(-33.6966+1.91915*plog-0.592320*plog*plog);
      break;
   }
   case 2: { // Pion
      double ptemp = ptrk;
      if(ptrk<0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      offset = 0.001*(-27.9965 + 1.213 * plog - 2.02740 * plog * plog);
      break;
   }
   case 3: { // Kaon
      double ptemp = ptrk;
      if(ptrk<0.3) ptemp = 0.3;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      offset = 0.001*(-23.4842 -28.7569 * plog + 78.21* plog *plog);
      break;
   }

   case 4: { // Proton
      double ptemp = ptrk;
      if(ptrk<0.4) ptemp = 0.4;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      if(charge>0)
         offset = 0.001*(-4.854-110.540*plog+99.8732*plog*plog);
      if(charge<0)
         offset = 0.001*(27.047-145.120*plog+167.014*plog*plog);
      break;
   }

   default:
      offset = 0.0;
      break;
   }
   //  offset = 0.0;
   return offset;
}

double TofEPID::sigmaTofE(int n,  int cntr, double ptrk, double rtof, double ph,double charge) {

   double sigma;
   //   double gb = ptrk/xmass(n);
   switch(n) {

   case 0: { // Electron
      double ptemp = ptrk;
      if(ptrk < 0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      sigma = 0.001 * (109.974 +15.2457 * plog + 36.8139 * plog * plog);

      break;
   }

   case 1: { // Muon
      double ptemp = ptrk;
      if(ptrk < 0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      sigma = 0.001 * (96.5077 -2.96232 * plog + 3.12910 * plog * plog);
      break;
   }

   case 2: { // pion
      double ptemp = ptrk;
      if(ptrk < 0.2) ptemp = 0.2;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      sigma = 0.001 * (105.447 - 2.08044 * plog + 3.44846 * plog * plog);
      break;
   }

   case 3: { // Kaon
      double ptemp = ptrk;
      if(ptrk < 0.3) ptemp = 0.3;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      sigma = 0.001*(88.8806 - 26.8464 * plog + 113.672 * plog * plog);
      break;
   }
   case 4: { // Proton
      double ptemp = ptrk;
      if(ptrk < 0.5) ptemp = 0.5;
      if(ptrk > 2.1) ptemp = 2.1;
      double plog = log(ptemp);
      if(charge>0)
         sigma = 0.001 * (96.3534 -44.1139 * plog + 53.9805 * plog * plog);
      if(charge<0)
         sigma = 0.001 * (157.345 -98.7357 * plog + 55.1145 * plog * plog);
      break;
   }

   default:
      sigma = 0.100;

      break;
   }
   // sigma =1;
   return sigma;
}

