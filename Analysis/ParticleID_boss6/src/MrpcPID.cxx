#include <cmath>

#include "ParticleID/MrpcPID.h"

#ifndef BEAN
#include "MdcRecEvent/RecMdcTrack.h"
#include "TofRecEvent/RecTofTrack.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#else
#include "TofHitStatus.h"
#endif

MrpcPID * MrpcPID::m_pointer = 0;

MrpcPID * MrpcPID::instance() {
   if(!m_pointer) m_pointer = new MrpcPID();
   return m_pointer;
}

MrpcPID::MrpcPID():ParticleIDBase() {
   //readSigPar();
}

void MrpcPID::init() {
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

void MrpcPID::calculate() {
   if(particleIDCalculation()==0) m_ndof=1;
}

int MrpcPID::particleIDCalculation() {
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


   double tof  = (*it)->tof();
   if(tof <=0 ) return irc;
   double path = (*it)->path();
   m_rhit = (*it)->zrhit();


   double beta2 = path*path/velc()/velc()/tof/tof;
   m_mass2 = ptrk * ptrk * (1/beta2 -1);



   double chitemp = 99.;
   double pdftemp = 0;


   

   double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};

   for(int i = 0; i < 5; i++) {


      double texp = (*it)->texp(i);
      m_offset[i] = tof - texp;


   


      double sigma_tmp= (*it)->sigma(i);


      if (sigma_tmp!=0)  m_sigma[i]=sigma_tmp;
      else m_sigma[i]=0.08;


      m_chi[i] = m_offset[i]/m_sigma[i];
      
      
      if(fabs(m_chi[i]) < chitemp) chitemp = fabs(m_chi[i]);
      double ppp = pdfCalculate(m_chi[i],1);
      if(fabs(ppp) > pdftemp) pdftemp = fabs(ppp);
   }
   m_chimin = chitemp;

   // calculate prob

   for(int i = 0; i < 5; i++)
      m_prob[i] = probCalculate(m_chi[i]*m_chi[i], 1);



   m_ndof = 1;
   irc = 0;
   return irc;
}


