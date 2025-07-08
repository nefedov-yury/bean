#include <fstream>
#include <cmath>
#include <cstdlib>
#include "Math/ChebyshevPol.h"
#include "ParticleID/DedxPID.h"
#ifndef BEAN
#include "EvtRecEvent/EvtRecTrack.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcDedx.h"
#else
#include "DstEvtRecTracks.h"
#endif

DedxPID * DedxPID::m_pointer = 0;

DedxPID * DedxPID::instance() {
   if(!m_pointer) m_pointer = new DedxPID();
   return m_pointer;
}

DedxPID::DedxPID():ParticleIDBase() {
   m_readstate=0;
}

void DedxPID::init() {
   for(int i = 0; i < 5; i++) {
      m_chi[i] = 99.0;
      m_prob[i] = -1.0;
      m_offset[i] = 99.0;
      m_sigma[i] = 1.0;
   }
   m_chimin = 99.;
   m_pdfmin =99;
   m_ndof = 0;
   m_goodHits = -99;
   m_normPH = -99;
   m_probPH = -99;
   m_nhitcutdx=5;
}

void DedxPID::calculate() {
   //    int rundedx = getRunNo();
   if(!m_readstate) {
      inputpar();
      m_readstate=1;
   }
   if(particleIDCalculation() == 0) m_ndof=1;
}

int DedxPID::particleIDCalculation() {
   //    int rundedx2 = getRunNo();
   int runNo = getRunNo();
   int nhitcutdedx=getNhitCutDx();
   int irc = -1;
   EvtRecTrack* recTrk = PidTrk();
   if(!(recTrk->isMdcTrackValid())) return irc;
   RecMdcTrack* mdcTrk = recTrk->mdcTrack();

   double ptrk = mdcTrk->p();
   int charge = mdcTrk->charge();
   if(ptrk>5) return irc;
   double cost = cos(mdcTrk->theta());
   //   double sig_the= sin(mdcTrk->theta());

   if(!(recTrk->isMdcDedxValid())) return irc;
   RecMdcDedx* dedxTrk = recTrk->mdcDedx();

#ifndef BEAN
   //if((dedxTrk->normPH()>30)||(dedxTrk->normPH()<0)) return irc;
   m_goodHits = dedxTrk->numGoodHits();
   //if(dedxTrk->numGoodHits()<nhitcutdedx) return irc;
#else
// This is BEAN part:
#if ( BOSS_VER <= 705 )
   if((dedxTrk->normPH()>30)||(dedxTrk->normPH()<0)) return irc;
#endif
   m_goodHits = dedxTrk->numGoodHits();
#if ( BOSS_VER <= 705 )
   if(dedxTrk->numGoodHits()<nhitcutdedx) return irc;
#endif
#endif
   m_normPH = dedxTrk->normPH();
   m_probPH = dedxTrk->probPH();
   // calculate chi and min chi
   double chitemp = 99.;
   double pdftemp = 0;
   //  double testchi[5];
   //  double testptrk[5];
   //  double testcost[5];
   for(int i = 0; i < 5; i++) {
#if !defined(BEAN) || (BOSS_VER >= 708)
        //add by CHEN Tong
        if (recTrk->isMdcKalTrackValid())
        {
            RecMdcKalTrack *kalTrk = recTrk->mdcKalTrack();
            if (i == 0)
            {
                kalTrk->setPidType(RecMdcKalTrack::electron);
            }
            else if (i == 1)
            {
                kalTrk->setPidType(RecMdcKalTrack::muon);
            }
            else if (i == 2)
            {
                kalTrk->setPidType(RecMdcKalTrack::pion);
            }
            else if (i == 3)
            {
                kalTrk->setPidType(RecMdcKalTrack::kaon);
            }
            else if (i == 4)
            {
                kalTrk->setPidType(RecMdcKalTrack::proton);
            }
            ptrk = kalTrk->p();
            cost = cos(kalTrk->theta());
        }
        //CT
#endif
      double sep = dedxTrk->chi(i);

#ifndef BEAN
      string sftver = getenv("BES_RELEASE");
      string sft;
      sft.assign(sftver,0,5);
      if(sft=="6.6.0"||sft=="6.5.5") {
         m_chi[i] = CorrDedx(i,ptrk,cost,sep,charge);
      }
      else
         m_chi[i]=sep;
#else
// This is BEAN part:
#if ( BOSS_VER == 660 || BOSS_VER == 655 )
      m_chi[i] = CorrDedx(i,ptrk,cost,sep,charge);
#else
      m_chi[i]=sep;
#endif
#endif

      m_offset[i] = offsetDedx(i, ptrk, cost);
      m_sigma[i] = sigmaDedx(i, ptrk, cost);
#if !defined(BEAN) || (BOSS_VER >= 708)
      //add by CHEN Tong
      // nefedov: I added parentheses to avoid compiler warnings
        if ( (runNo >= 9947 && runNo <=10878) ||
             (runNo>= 27255 && runNo <= 28236) ||
             (runNo>= 52940 && runNo <= 54976) ||
             (runNo>= 55861 && runNo <= 56546) ||
             (runNo>= 56788 && runNo <= 59015)    )
        {
            m_offsetCorr[i] = offsetCorr(i, charge, ptrk, cost);
            m_sigmaCorr[i] = sigmaCorr(i, charge, ptrk, cost);
            m_chi[i] = (sep - m_offsetCorr[i]) / m_sigmaCorr[i];
        }
        //CT
#endif
      if(fabs(m_chi[i]) < chitemp) chitemp = fabs(m_chi[i]);
      double ppp = pdfCalculate(m_chi[i],1);
      if(fabs(ppp) > pdftemp) pdftemp = fabs(ppp);

   }
   m_chimin = chitemp;
   m_pdfmin = pdftemp;
   if(m_chimin > chiMinCut()) return irc;
   if(pdftemp < pdfCalculate(pdfMinSigmaCut(),1)) return irc;


   // calculate prob

   for(int i = 0; i < 5; i++)
      m_prob[i] = probCalculate(m_chi[i]*m_chi[i], 1);

   m_ndof = 1;
   irc = 0;
   return irc;
}

//
//  dE/dx Correction routines
//



double DedxPID::offsetDedx(int n, double ptrk, double cost) {
   return 0;
}

double DedxPID::CorrDedx(int n, double ptrk, double cost,double chi,int charge) {
   int rundedx2 = getRunNo();
   double offset = 0.0;
   double offsetp = 0.0;
   double offsetc = 0.0;
   double sigcos = 1;
   double sigp = 1;
   double chicor=chi;
   //   double gb = ptrk/xmass(n);

   switch(n) {
   case 0: { // Electron
      break;
   }

   case 1: {// Muon
      break;
   }

   case 2: {// Pion
      //     double  ptemp = ptrk;
      double  costm = cost;
      if(ptrk<0.1||ptrk>1) break;
      int index = int((ptrk-0.1)/0.05);
      if(index<=0) index=1;
      if(index>=17) index=16;

      if(fabs(costm)>=0.8) break;
      int index1 = int((costm+0.8)/0.1);
      if(index1<=0) index1=1;
      if(index1>=15) index1=14;

      //psipp data
      if(rundedx2>=11414&&rundedx2<=14604) {
         offsetp = cal_par(index,m_psipp_pi_ptrk_offset,ptrk,0.125,0.05);
         sigp = cal_par(index,m_psipp_pi_ptrk_sigma,ptrk,0.125,0.05);
         offsetc = cal_par(index1,m_psipp_pi_theta_offset,costm,-0.75,0.1);
         sigcos = cal_par(index1,m_psipp_pi_theta_sigma,costm,-0.75,0.1);
      }
      //psipp mc
      if(rundedx2<=-11414&&rundedx2>=-14604) {
         offsetp = cal_par(index,m_psipp_mc_pi_ptrk_offset,ptrk,0.125,0.05);
         sigp = cal_par(index,m_psipp_mc_pi_ptrk_sigma,ptrk,0.125,0.05);
         offsetc = cal_par(index1,m_psipp_mc_pi_theta_offset,costm,-0.75,0.1);
         sigcos = cal_par(index1,m_psipp_mc_pi_theta_sigma,costm,-0.75,0.1);
      }

      offset=offsetp+sigp*offsetc;
      chicor=(chicor-offset)/(sigcos*sigp);
      break;
   }

   case 3: {// Kaon
      //     double  ptemp = ptrk;
      double  costm = cost;
      if(ptrk<0.3||ptrk>0.8) break;
      offset=0;
      int index = int((ptrk-0.3)/0.1);
      if(index<=0) index=1;
      if(index>=4) index=3;

      int index1 = int((costm+0.9)/0.1);
      if(index1<=0) index1=1;
      if(index1>=17) index1=16;
      //data Jpsi
      if(rundedx2>=9947&&rundedx2<=10878) {
         if(charge>0) {
            offsetp = cal_par(index,m_jpsi_kap_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_kap_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_kap_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_kap_theta_sigma,costm,-0.85,0.1);
            }
         }
         if(charge<0) {
            offsetp = cal_par(index,m_jpsi_kam_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_kam_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_kam_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_kam_theta_sigma,costm,-0.85,0.1);
            }
         }
      }

      //mc Jpsi
      if(rundedx2<=-9947&&rundedx2>=-10878) {
         if(charge>0) {
            offsetp = cal_par(index,m_jpsi_mc_kap_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_mc_kap_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_mc_kap_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_mc_kap_theta_sigma,costm,-0.85,0.1);
            }
         }
         if(charge<0) {
            offsetp = cal_par(index,m_jpsi_mc_kam_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_mc_kam_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_mc_kam_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_mc_kam_theta_sigma,costm,-0.85,0.1);
            }
         }
      }

      //data Psip
      if(rundedx2>=8093&&rundedx2<=9025) {
         if(ptrk<0.3||ptrk>1.2) break;
         index = int((ptrk-0.3)/0.1);
         if(index<=0) index=1;
         if(index>=8) index=7;
         if(charge>0) {
            offsetp = cal_par(index,m_psip_kap_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_kap_ptrk_sigma,ptrk,0.35,0.1);
         }
         if(charge<0) {
            offsetp = cal_par(index,m_psip_kam_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_kam_ptrk_sigma,ptrk,0.35,0.1);
         }
      }

      //mc Psip
      if(rundedx2<=-8093&&rundedx2>=-9025) {
         //    if(ptrk < 0.4) ptrk = 0.4;
         if(ptrk<0.3||ptrk>1.2) break;
         index = int((ptrk-0.3)/0.1);
         if(index<=0) index=1;
         if(index>=8) index=7;
         if(charge>0) {
            offsetp = cal_par(index,m_psip_mc_kap_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_mc_kap_ptrk_sigma,ptrk,0.35,0.1);
         }
         if(charge<0) {
            offsetp = cal_par(index,m_psip_mc_kam_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_mc_kam_ptrk_sigma,ptrk,0.35,0.1);
         }
      }


      //psipp kaon data
      if(rundedx2>=11414&&rundedx2<=14604) {
         if(ptrk<0.15||ptrk>1) break;
         index = int((ptrk-0.15)/0.05);
         if(index<=0) index=1;
         if(index>=16) index=15;
         if(fabs(costm)>=0.8) break;
         index1 = int((costm+0.8)/0.1);
         if(index1<=0) index1=1;
         if(index1>=15) index1=14;

         offsetp = cal_par(index,m_psipp_ka_ptrk_offset,ptrk,0.175,0.05);
         sigp = cal_par(index,m_psipp_ka_ptrk_sigma,ptrk,0.175,0.05);
         offsetc = cal_par(index1,m_psipp_ka_theta_offset,costm,-0.75,0.1);
         sigcos = cal_par(index1,m_psipp_ka_theta_sigma,costm,-0.75,0.1);
      }
      //psipp kaon mc
      if(rundedx2<=-11414&&rundedx2>=-14604) {
         if(ptrk<0.15||ptrk>1) break;
         index = int((ptrk-0.15)/0.05);
         if(index<=0) index=1;
         if(index>=16) index=15;
         if(fabs(costm)>=0.8) break;
         index1 = int((costm+0.8)/0.1);
         if(index1<=0) index1=1;
         if(index1>=15) index1=14;
         offsetp = cal_par(index,m_psipp_mc_ka_ptrk_offset,ptrk,0.175,0.05);
         sigp = cal_par(index,m_psipp_mc_ka_ptrk_sigma,ptrk,0.175,0.05);
         offsetc = cal_par(index1,m_psipp_mc_ka_theta_offset,costm,-0.75,0.1);
         sigcos = cal_par(index1,m_psipp_mc_ka_theta_sigma,costm,-0.75,0.1);
      }

      offset=offsetp+sigp*offsetc;
      chicor=(chicor-offset)/(sigcos*sigp);
      break;
   }

   case 4 : { // Proton
      //     double  ptemp = ptrk;
      double  costm = cost;
      if(ptrk<0.3||ptrk>1.1) break;
      int index = int((ptrk-0.3)/0.1);
      if(index<=0) index=1;
      if(index>=7) index=6;

      int index1 = int((costm+0.9)/0.1);
      if(index1<=0) index1=1;
      if(index1>=17) index1=16;

      //    double plog = log(ptemp);
      offset=0;
      if(rundedx2>=9947&&rundedx2<=10878) {
         if(charge>0) {
            offsetp = cal_par(index,m_jpsi_protonp_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_protonp_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_protonp_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_protonp_theta_sigma,costm,-0.85,0.1);
            }
         }
         if(charge<0) {
            offsetp = cal_par(index,m_jpsi_protonm_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_protonm_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_protonm_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_protonm_theta_sigma,costm,-0.85,0.1);
            }
         }
      }

      //mc JPsi
      if(rundedx2<=-9947&&rundedx2>=-10878) {
         if(charge>0) {
            offsetp = cal_par(index,m_jpsi_mc_protonp_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_mc_protonp_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_mc_protonp_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_mc_protonp_theta_sigma,costm,-0.85,0.1);
            }
         }
         if(charge<0) {
            offsetp = cal_par(index,m_jpsi_mc_protonm_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_jpsi_mc_protonm_ptrk_sigma,ptrk,0.35,0.1);
            if(fabs(costm)<=0.83) {
               offsetc = cal_par(index1,m_jpsi_mc_protonm_theta_offset,costm,-0.85,0.1);
               sigcos = cal_par(index1,m_jpsi_mc_protonm_theta_sigma,costm,-0.85,0.1);
            }
         }
      }

      //data Psip
      if(rundedx2>=8093&&rundedx2<=9025) {
         if(charge>0) {
            offsetp = cal_par(index,m_psip_protonp_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_protonp_ptrk_sigma,ptrk,0.35,0.1);
         }
         if(charge<0) {
            offsetp = cal_par(index,m_psip_protonm_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_protonm_ptrk_sigma,ptrk,0.35,0.1);
         }
      }

      //mc Psip
      if(rundedx2<=-8093&&rundedx2>=-9025) {
         if(charge>0) {
            offsetp = cal_par(index,m_psip_mc_protonp_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_mc_protonp_ptrk_sigma,ptrk,0.35,0.1);
         }
         if(charge<0) {
            offsetp = cal_par(index,m_psip_mc_protonm_ptrk_offset,ptrk,0.35,0.1);
            sigp = cal_par(index,m_psip_mc_protonm_ptrk_sigma,ptrk,0.35,0.1);
         }
      }

      //psipp proton data
      if(rundedx2>=11414&&rundedx2<=14604) {
         if(ptrk<0.2||ptrk>1.1) break;
         index = int((ptrk-0.2)/0.05);
         if(index<=0) index=1;
         if(index>=17) index=16;
         if(fabs(costm)>=0.83) break;
         index1 = int((costm+0.9)/0.1);
         if(index1<=0) index1=1;
         if(index1>=17) index1=16;

         offsetp = cal_par(index,m_psipp_proton_ptrk_offset,ptrk,0.225,0.05);
         sigp = cal_par(index,m_psipp_proton_ptrk_sigma,ptrk,0.225,0.05);
         offsetc = cal_par(index1,m_psipp_proton_theta_offset,costm,-0.85,0.1);
         sigcos = cal_par(index1,m_psipp_proton_theta_sigma,costm,-0.85,0.1);
      }
      //psipp proton mc
      if(rundedx2<=-11414&&rundedx2>=-14604) {
         if(ptrk<0.2||ptrk>1.1) break;
         index = int((ptrk-0.2)/0.1);
         if(index<=0) index=1;
         if(index>=8) index=7;
         if(fabs(costm)>=0.83) break;
         index1 = int((costm+0.9)/0.1);
         if(index1<=0) index1=1;
         if(index1>=17) index1=16;
         offsetp = cal_par(index,m_psipp_mc_proton_ptrk_offset,ptrk,0.25,0.1);
         sigp = cal_par(index,m_psipp_mc_proton_ptrk_sigma,ptrk,0.25,0.1);
         offsetc = cal_par(index1,m_psipp_mc_proton_theta_offset,costm,-0.85,0.1);
         sigcos = cal_par(index1,m_psipp_mc_proton_theta_sigma,costm,-0.85,0.1);
      }
      offset=offsetp+sigp*offsetc;
      chicor=(chicor-offset)/(sigcos*sigp);
      break;
   }

   default:
      offset = 0.0;
      break;
   }
   //  offset = 0.0;
   return chicor;
}

double DedxPID::sigmaDedx(int n, double ptrk, double cost) {

   /*  int rundedx3 = getRunNo();
     double sigma = 1.0;
    double sigp = 1.0;
    double sigmac = 1.0;
    double gb = ptrk/xmass(n);
    switch(n) {

    case 0: {// Electron
      double  ptemp = ptrk;
      double  costm = cost;
   break;
    }

    case 1: {// Muon
      double  ptemp = ptrk;
      double  costm = cost;
   break;
    }

    case 2: {// Pion
      double  ptemp = ptrk;
      double  costm = cost;
   break;
    }

    case 3: { // Kaon
      double  ptemp = ptrk;
      double  costm = cost;
   break;
    }


    case 4: {// Proton
      double  ptemp = ptrk;
      double  costm = cost;
   break;
    }

    default:
      sigma = 1.0;
      break;
    }
   */
   //  sigma = 1.2;
   //  sigma =1.0;
   return 1;
   //  return sigma;
}

double DedxPID::mypol3(double x, double par0, double par1, double par2, double par3)
{
   double y = x;
   return par0 + (par1 * y) +(par2 * y * y) + (par3 * y * y * y);

}

double DedxPID::mypol5(double x, double par0, double par1, double par2, double par3, double par4, double par5)
{
   double y = x;
   return par0 + (par1 * y) +(par2 * y * y) + (par3 * y * y * y) + (par4 * y * y * y *y)+ (par5 * y * y * y * y * y);

}

void DedxPID::inputpar() {

   //Jpsi ka+ momentum correction
   std::string jpsi_kap_mom = path + "/share/JPsi/kaon/dedx_kap.txt";
   std::string jpsi_kap_mom_mc = path + "/share/JPsi/kaon/dedx_kap_mc.txt";
   ifstream inputmomdata6(jpsi_kap_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata6 ) {
      cout << " can not open: " << jpsi_kap_mom << endl;
      exit(1);
   }
   ifstream inputmomdata6mc(jpsi_kap_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata6mc ) {
      cout << " can not open: " << jpsi_kap_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<12; i++) {
      inputmomdata6>>m_jpsi_kap_ptrk_offset[i];
      inputmomdata6>>m_jpsi_kap_ptrk_sigma[i];
      inputmomdata6mc>>m_jpsi_mc_kap_ptrk_offset[i];
      inputmomdata6mc>>m_jpsi_mc_kap_ptrk_sigma[i];
   }

   //Jpsi ka- momentum correction
   std::string jpsi_kam_mom = path + "/share/JPsi/kaon/dedx_kam.txt";
   std::string jpsi_kam_mom_mc =  path + "/share/JPsi/kaon/dedx_kam_mc.txt";
   ifstream inputmomdata7(jpsi_kam_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata7 ) {
      cout << " can not open: " << jpsi_kam_mom << endl;
      exit(1);
   }
   ifstream inputmomdata7mc(jpsi_kam_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata7mc ) {
      cout << " can not open: " << jpsi_kam_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<12; i++) {
      inputmomdata7>>m_jpsi_kam_ptrk_offset[i];
      inputmomdata7>>m_jpsi_kam_ptrk_sigma[i];
      inputmomdata7mc>>m_jpsi_mc_kam_ptrk_offset[i];
      inputmomdata7mc>>m_jpsi_mc_kam_ptrk_sigma[i];

   }


   //Jpsi ka+ theta correction
   std::string jpsi_kap_the =  path + "/share/JPsi/kaon/dedx_kap_theta.txt";
   std::string jpsi_kap_the_mc =  path + "/share/JPsi/kaon/dedx_kap_theta_mc.txt";
   ifstream inputmomdata8(jpsi_kap_the.c_str(),std::ios_base::in);
   if ( !inputmomdata8 ) {
      cout << " can not open: " << jpsi_kap_the << endl;
      exit(1);
   }
   ifstream inputmomdata8mc(jpsi_kap_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata8mc ) {
      cout << " can not open: " << jpsi_kap_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata8>>m_jpsi_kap_theta_offset[i];
      inputmomdata8>>m_jpsi_kap_theta_sigma[i];
      inputmomdata8mc>>m_jpsi_mc_kap_theta_offset[i];
      inputmomdata8mc>>m_jpsi_mc_kap_theta_sigma[i];
   }

   //Jpsi ka- theta correction
   std::string jpsi_kam_the =  path + "/share/JPsi/kaon/dedx_kam_theta.txt";
   std::string jpsi_kam_the_mc =  path + "/share/JPsi/kaon/dedx_kam_theta_mc.txt";
   ifstream inputmomdata9(jpsi_kam_the.c_str(),std::ios_base::in);
   if ( !inputmomdata9 ) {
      cout << " can not open: " << jpsi_kam_the << endl;
      exit(1);
   }
   ifstream inputmomdata9mc(jpsi_kam_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata9mc ) {
      cout << " can not open: " << jpsi_kam_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata9>>m_jpsi_kam_theta_offset[i];
      inputmomdata9>>m_jpsi_kam_theta_sigma[i];
      inputmomdata9mc>>m_jpsi_mc_kam_theta_offset[i];
      inputmomdata9mc>>m_jpsi_mc_kam_theta_sigma[i];
   }

   //Jpsi proton+ momentum correction
   std::string jpsi_protonp_mom =  path + "/share/JPsi/proton/dedx_protonp.txt";
   std::string jpsi_protonp_mom_mc =  path + "/share/JPsi/proton/dedx_protonp_mc.txt";
   ifstream inputmomdata12(jpsi_protonp_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata12 ) {
      cout << " can not open: " << jpsi_protonp_mom << endl;
      exit(1);
   }
   ifstream inputmomdata12mc(jpsi_protonp_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata12mc ) {
      cout << " can not open: " << jpsi_protonp_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<8; i++) {
      inputmomdata12>>m_jpsi_protonp_ptrk_offset[i];
      inputmomdata12>>m_jpsi_protonp_ptrk_sigma[i];
      inputmomdata12mc>>m_jpsi_mc_protonp_ptrk_offset[i];
      inputmomdata12mc>>m_jpsi_mc_protonp_ptrk_sigma[i];
   }

   //Jpsi proton- momentum correction
   std::string jpsi_protonm_mom =  path + "/share/JPsi/proton/dedx_protonm.txt";
   std::string jpsi_protonm_mom_mc =  path + "/share/JPsi/proton/dedx_protonm_mc.txt";
   ifstream inputmomdata13(jpsi_protonm_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata13 ) {
      cout << " can not open: " << jpsi_protonm_mom << endl;
      exit(1);
   }
   ifstream inputmomdata13mc(jpsi_protonm_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata13mc ) {
      cout << " can not open: " << jpsi_protonm_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<8; i++) {
      inputmomdata13>>m_jpsi_protonm_ptrk_offset[i];
      inputmomdata13>>m_jpsi_protonm_ptrk_sigma[i];
      inputmomdata13mc>>m_jpsi_mc_protonm_ptrk_offset[i];
      inputmomdata13mc>>m_jpsi_mc_protonm_ptrk_sigma[i];
   }

   //Jpsi proton+ theta correction
   std::string jpsi_protonp_the = path + "/share/JPsi/proton/dedx_protonp_theta.txt";
   std::string jpsi_protonp_the_mc = path + "/share/JPsi/proton/dedx_protonp_theta_mc.txt";

   ifstream inputmomdata14(jpsi_protonp_the.c_str(),std::ios_base::in);
   if ( !inputmomdata14 ) {
      cout << " can not open: " << jpsi_protonp_the << endl;
      exit(1);
   }
   ifstream inputmomdata14mc(jpsi_protonp_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata14mc ) {
      cout << " can not open: " << jpsi_protonp_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata14>>m_jpsi_protonp_theta_offset[i];
      inputmomdata14>>m_jpsi_protonp_theta_sigma[i];
      inputmomdata14mc>>m_jpsi_mc_protonp_theta_offset[i];
      inputmomdata14mc>>m_jpsi_mc_protonp_theta_sigma[i];
   }

   //Jpsi proton- theta correction
   std::string jpsi_protonm_the = path + "/share/JPsi/proton/dedx_protonm_theta.txt";
   std::string jpsi_protonm_the_mc = path + "/share/JPsi/proton/dedx_protonm_theta_mc.txt";
   ifstream inputmomdata15(jpsi_protonm_the.c_str(),std::ios_base::in);
   if ( !inputmomdata15 ) {
      cout << " can not open: " << jpsi_protonm_the << endl;
      exit(1);
   }
   ifstream inputmomdata15mc(jpsi_protonm_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata15mc ) {
      cout << " can not open: " << jpsi_protonm_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata15>>m_jpsi_protonm_theta_offset[i];
      inputmomdata15>>m_jpsi_protonm_theta_sigma[i];
      inputmomdata15mc>>m_jpsi_mc_protonm_theta_offset[i];
      inputmomdata15mc>>m_jpsi_mc_protonm_theta_sigma[i];
   }




   // Psip ka+ momentum correction
   std::string psip_kap_mom = path + "/share/Psip/kaon/dedx_kap.txt";
   std::string psip_kap_mom_mc = path + "/share/Psip/kaon/dedx_kap_mc.txt";
   ifstream inputmomdata24(psip_kap_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata24 ) {
      cout << " can not open: " << psip_kap_mom << endl;
      exit(1);
   }
   ifstream inputmomdata24mc(psip_kap_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata24mc ) {
      cout << " can not open: " << psip_kap_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<9; i++) {
      inputmomdata24>>m_psip_kap_ptrk_offset[i];
      inputmomdata24>>m_psip_kap_ptrk_sigma[i];
      inputmomdata24mc>>m_psip_mc_kap_ptrk_offset[i];
      inputmomdata24mc>>m_psip_mc_kap_ptrk_sigma[i];
   }

   //Psip ka- momentum correction
   std::string psip_kam_mom = path + "/share/Psip/kaon/dedx_kam.txt";
   std::string psip_kam_mom_mc = path + "/share/Psip/kaon/dedx_kam_mc.txt";
   ifstream inputmomdata25(psip_kam_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata25 ) {
      cout << " can not open: " << psip_kam_mom << endl;
      exit(1);
   }
   ifstream inputmomdata25mc(psip_kam_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata25mc ) {
      cout << " can not open: " << psip_kam_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<9; i++) {
      inputmomdata25>>m_psip_kam_ptrk_offset[i];
      inputmomdata25>>m_psip_kam_ptrk_sigma[i];
      inputmomdata25mc>>m_psip_mc_kam_ptrk_offset[i];
      inputmomdata25mc>>m_psip_mc_kam_ptrk_sigma[i];
   }


   // Psip proton+ momentum correction
   std::string psip_protonp_mom = path + "/share/Psip/proton/dedx_protonp.txt";
   std::string psip_protonp_mom_mc = path + "/share/Psip/proton/dedx_protonp_mc.txt";
   ifstream inputmomdata26(psip_protonp_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata26 ) {
      cout << " can not open: " << psip_protonp_mom << endl;
      exit(1);
   }
   ifstream inputmomdata26mc(psip_protonp_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata26mc ) {
      cout << " can not open: " << psip_protonp_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<9; i++) {
      inputmomdata26>>m_psip_protonp_ptrk_offset[i];
      inputmomdata26>>m_psip_protonp_ptrk_sigma[i];
      inputmomdata26mc>>m_psip_mc_protonp_ptrk_offset[i];
      inputmomdata26mc>>m_psip_mc_protonp_ptrk_sigma[i];
   }

   //Psip proton- momentum correction
   std::string psip_protonm_mom = path + "/share/Psip/proton/dedx_protonm.txt";
   std::string psip_protonm_mom_mc = path + "/share/Psip/proton/dedx_protonm_mc.txt";
   ifstream inputmomdata27(psip_protonm_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata27 ) {
      cout << " can not open: " << psip_protonm_mom << endl;
      exit(1);
   }
   ifstream inputmomdata27mc(psip_protonm_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata27mc ) {
      cout << " can not open: " << psip_protonm_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<9; i++) {
      inputmomdata27>>m_psip_protonm_ptrk_offset[i];
      inputmomdata27>>m_psip_protonm_ptrk_sigma[i];
      inputmomdata27mc>>m_psip_mc_protonm_ptrk_offset[i];
      inputmomdata27mc>>m_psip_mc_protonm_ptrk_sigma[i];
   }

   //Psipp pi momentum correction
   std::string psipp_pi_mom = path + "/share/Psipp/pion/dedx_pi.txt";
   std::string psipp_pi_mom_mc = path + "/share/Psipp/pion/dedx_pi_mc.txt";
   ifstream inputmomdata28(psipp_pi_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata28 ) {
      cout << " can not open: " << psipp_pi_mom << endl;
      exit(1);
   }
   ifstream inputmomdata28mc(psipp_pi_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata28mc ) {
      cout << " can not open: " << psipp_pi_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata28>>m_psipp_pi_ptrk_offset[i];
      inputmomdata28>>m_psipp_pi_ptrk_sigma[i];
      inputmomdata28mc>>m_psipp_mc_pi_ptrk_offset[i];
      inputmomdata28mc>>m_psipp_mc_pi_ptrk_sigma[i];
   }

   //Psipp pi theta correction
   std::string psipp_pi_the = path + "/share/Psipp/pion/dedx_pi_theta.txt";
   std::string psipp_pi_the_mc = path + "/share/Psipp/pion/dedx_pi_theta_mc.txt";
   ifstream inputmomdata29(psipp_pi_the.c_str(),std::ios_base::in);
   if ( !inputmomdata29 ) {
      cout << " can not open: " << psipp_pi_the << endl;
      exit(1);
   }
   ifstream inputmomdata29mc(psipp_pi_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata29mc ) {
      cout << " can not open: " << psipp_pi_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<16; i++) {
      inputmomdata29>>m_psipp_pi_theta_offset[i];
      inputmomdata29>>m_psipp_pi_theta_sigma[i];
      inputmomdata29mc>>m_psipp_mc_pi_theta_offset[i];
      inputmomdata29mc>>m_psipp_mc_pi_theta_sigma[i];
   }

   //Psipp ka momentum correction
   std::string psipp_ka_mom = path + "/share/Psipp/kaon/dedx_ka.txt";
   std::string psipp_ka_mom_mc = path + "/share/Psipp/kaon/dedx_ka_mc.txt";
   ifstream inputmomdata30(psipp_ka_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata30 ) {
      cout << " can not open: " << psipp_ka_mom << endl;
      exit(1);
   }
   ifstream inputmomdata30mc(psipp_ka_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata30mc ) {
      cout << " can not open: " << psipp_ka_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<17; i++) {
      inputmomdata30>>m_psipp_ka_ptrk_offset[i];
      inputmomdata30>>m_psipp_ka_ptrk_sigma[i];
      inputmomdata30mc>>m_psipp_mc_ka_ptrk_offset[i];
      inputmomdata30mc>>m_psipp_mc_ka_ptrk_sigma[i];
   }

   //Psipp ka theta correction
   std::string psipp_ka_the = path + "/share/Psipp/kaon/dedx_ka_theta.txt";
   std::string psipp_ka_the_mc = path + "/share/Psipp/kaon/dedx_ka_theta_mc.txt";
   ifstream inputmomdata31(psipp_ka_the.c_str(),std::ios_base::in);
   if ( !inputmomdata31 ) {
      cout << " can not open: " << psipp_ka_the << endl;
      exit(1);
   }
   ifstream inputmomdata31mc(psipp_ka_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata31mc ) {
      cout << " can not open: " << psipp_ka_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<16; i++) {
      inputmomdata31>>m_psipp_ka_theta_offset[i];
      inputmomdata31>>m_psipp_ka_theta_sigma[i];
      inputmomdata31mc>>m_psipp_mc_ka_theta_offset[i];
      inputmomdata31mc>>m_psipp_mc_ka_theta_sigma[i];
   }


   //Psipp proton momentum correction
   std::string psipp_proton_mom = path + "/share/Psipp/proton/dedx_proton.txt";
   std::string psipp_proton_mom_mc = path + "/share/Psipp/proton/dedx_proton_mc.txt";
   ifstream inputmomdata32(psipp_proton_mom.c_str(),std::ios_base::in);
   if ( !inputmomdata32 ) {
      cout << " can not open: " << psipp_proton_mom << endl;
      exit(1);
   }
   ifstream inputmomdata32mc(psipp_proton_mom_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata32mc ) {
      cout << " can not open: " << psipp_proton_mom_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata32>>m_psipp_proton_ptrk_offset[i];
      inputmomdata32>>m_psipp_proton_ptrk_sigma[i];
   }
   for(int i=0; i<9; i++) {
      inputmomdata32mc>>m_psipp_mc_proton_ptrk_offset[i];
      inputmomdata32mc>>m_psipp_mc_proton_ptrk_sigma[i];
   }

   //Psipp proton theta correction
   std::string psipp_proton_the = path + "/share/Psipp/proton/dedx_proton_theta.txt";
   std::string psipp_proton_the_mc = path + "/share/Psipp/proton/dedx_proton_theta_mc.txt";
   ifstream inputmomdata33(psipp_proton_the.c_str(),std::ios_base::in);
   if ( !inputmomdata33 ) {
      cout << " can not open: " << psipp_proton_the << endl;
      exit(1);
   }
   ifstream inputmomdata33mc(psipp_proton_the_mc.c_str(),std::ios_base::in);
   if ( !inputmomdata33mc ) {
      cout << " can not open: " << psipp_proton_the_mc << endl;
      exit(1);
   }
   for(int i=0; i<18; i++) {
      inputmomdata33>>m_psipp_proton_theta_offset[i];
      inputmomdata33>>m_psipp_proton_theta_sigma[i];
      inputmomdata33mc>>m_psipp_mc_proton_theta_offset[i];
      inputmomdata33mc>>m_psipp_mc_proton_theta_sigma[i];
   }

}

double DedxPID::iterate(double ptrk,double *mean,double *p) {
   double p1,p2,p3;
   p2=((mean[0]-mean[1])*(p[1]*p[1]-p[2]*p[2])-(mean[1]-mean[2])*(p[0]*p[0]-p[1]*p[1]))/((p[0]-p[1])*(p[1]*p[1]-p[2]*p[2])-(p[1]-p[2])*(p[0]*p[0]-p[1]*p[1]));
   p3=((p[0]-p[1])*(mean[1]-mean[2])-(p[1]-p[2])*(mean[0]-mean[1]))/((p[0]-p[1])*(p[1]*p[1]-p[2]*p[2])-(p[1]-p[2])*(p[0]*p[0]-p[1]*p[1]));
   p1=mean[0]-p2*p[0]-p3*p[0]*p[0];
   double mean1 = p1+p2*ptrk+p3*ptrk*ptrk;
   return mean1;
}

double DedxPID::cal_par(int index1,double *m_jpsi_pip_ptrk_offset,double ptrk,double begin,double bin) {
   double mean1[3],p[3];
   p[0]=begin+(index1-1)*bin;
   p[1]=begin+index1*bin;
   p[2]=begin+(index1+1)*bin;
   mean1[0]=m_jpsi_pip_ptrk_offset[index1-1];
   mean1[1]=m_jpsi_pip_ptrk_offset[index1];
   mean1[2]=m_jpsi_pip_ptrk_offset[index1+1];
   double res=iterate(ptrk,mean1,p);
   return res;
}

//add by CHEN Tong
double DedxPID::offsetCorr(int n, int charge, double ptrk, double cost)
{
    double cosbin[30] = {-0.865, -0.7, -0.5, -0.325, -0.225, -0.19, -0.17, -0.15, -0.13, -0.11, -0.09, -0.07, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.225, 0.325, 0.5, 0.7, 0.865};
    double corr_offset = 0;
    if (n == 0 || n == 1 || n == 2 && ptrk > 0.20 || n == 3 && (ptrk > 0.43 || ptrk < 0.2 || fabs(cost) > 0.2) || n == 4 && (ptrk > 0.43 || ptrk < 0.2))
        return 0;
    else if (n == 2)
    {
        int tmp = 0;
        double par_pip_offset[30][7] = {
            {-0.273711, 0.418852,   -0.111059,  0.011915,   0.0595808,  -0.064644,  -0.00930414},
            {-0.269974, 0.295548,   -0.0805116, 0.112013,   -0.0844919, -0.0101669, 0.0472545},
            {-0.248076, 0.230898,   -0.0954416, 0.00461465, -0.0513009, 0.0841044,  -0.0320536},
            {-0.20576,  0.15908,    -0.0827095, -0.0534723, 0.011227,   0.0325395,  -0.0518127},
            {-0.168736, 0.129392,   -0.05107,   -0.126828,  0.126172,   -0.0489699, -0.0373492},
            {-0.176773, 0.163439,   -0.100698,  -0.111413,  0.115099,   -0.0726195, -0.00490596},
            {-0.0385498,    -0.0807258, 0.0897135,  -0.237934,  0.243725,   -0.108172,  -0.0188745},
            {-0.109386, 0.119386,   -0.154018,  -0.119937,  0.110189,   -0.0662777, -0.0217286},
            {-0.0942173,    0.171144,   -0.0806635, -0.182887,  0.22318,    -0.137083,  0.0174856},
            {-0.0209689,    0.0959425,  -0.00407716,    -0.186367,  0.309657,   -0.158974,  0.0899474},
            {-0.0742976,    0.194453,   -0.245653,  -0.119011,  0.169814,   -0.0985074, -0.0145252},
            {-0.00874818,   0.381593,   -0.398284,  0.133183,   0.0303175,  -0.024684,  -0.0535397},
            {0.264664,  0.128337,   -0.322216,  0.00757683, 0.072456,   -0.0374562, -0.0347229},
            {0.756099,  -0.426524,  0.0622856,  -0.194491,  0.148391,   -0.027508,  -0.0244935},
            {1.01453,   -0.255592,  -0.28767,   0.285283,   -0.249035,  0.257705,   -0.188833},
            {0.992926,  -0.470497,  -0.193231,  0.07265,    -0.14322,   0.124669,   -0.147956},
            {0.461868,  0.0132674,  -0.211087,  0.057901,   0.0263037,  0.0836409,  -0.108405},
            {0.1606,    0.0849392,  -0.189776,  -0.0816648, 0.0954246,  -0.00679478,    -0.0341808},
            {0.00423913,    0.0757551,  0.0444685,  -0.281896,  0.379112,   -0.212121,  0.0594902},
            {-0.253426, 0.326755,   -0.225738,  -0.0691578, 0.109788,   -0.0812613, -0.0482903},
            {-0.344417, 0.464815,   -0.268857,  0.0137453,  0.146423,   -0.0908125, 0.0131566},
            {0.180595,  -0.471044,  0.600536,   -0.680218,  0.666876,   -0.345305,  0.215572},
            {-0.202407, 0.290958,   -0.0843309, -0.0730275, 0.157111,   -0.0836953, 0.0374609},
            {-0.222727, 0.258844,   -0.0676259, -0.123459,  0.18615,    -0.0793805, 0.0245286},
            {-0.320709, 0.351355,   -0.22756,   0.019593,   0.0261016,  0.0191958,  -0.0449903},
            {-0.275517, 0.305623,   -0.147753,  -0.0221624, 0.0685188,  -0.00113907,    -0.0110694},
            {-0.237632, 0.222056,   -0.0967139, -0.0467247, 0.0196372,  0.013321,   -0.0537844},
            {-0.254455, 0.247849,   -0.0727208, -0.0359449, -0.0293587, 0.0650545,  -0.0368777},
            {-0.260665, 0.296247,   -0.0244744, 0.0393139,  -0.0729511, 0.0271136,  0.00654422},
            {-0.312108, 0.443109,   -0.0993222, 0.0326532,  4.76077e-05,    -0.026312,  0.00211329}
        };
        while (cost >= cosbin[tmp] && tmp != 28)
        {
            tmp++;
        }
        if (tmp == 0)
            tmp += 1;
        double par_cos[7];
        for (int j = 0; j < 7; j++)
        {
            double cosbin_tmp[3] = {cosbin[tmp - 1], cosbin[tmp], cosbin[tmp + 1]};
            double par_pip_offset_tmp[3] = {par_pip_offset[tmp - 1][j], par_pip_offset[tmp][j], par_pip_offset[tmp + 1][j]};
            par_cos[j] = interpolation(cost, cosbin_tmp, par_pip_offset_tmp);
        }
        double ptrk_tmp = (ptrk - 0.17) / 0.1;
        corr_offset = ROOT::Math::ChebyshevN(6, ptrk_tmp, par_cos);
        return corr_offset;
    }
    else if (n == 3)
    {
        double par_kp_offset[3][5] = {
            {0.00,  0.00,   0.00,   0.00,   0.00},
            {0.90,  0.90,   0.75,   0.15,   0.00},
            {0.00,  0.00,   0.00,   0.00,   0.00}
        };
        double p_bin[6] = {0.175, 0.225, 0.275, 0.325, 0.375, 0.425};
        int bin_p = (ptrk - 0.175) / 0.05;
        double int_p1 = (par_kp_offset[0][bin_p] - par_kp_offset[1][bin_p]) * fabs(cost) / 0.1 + par_kp_offset[1][bin_p];
        double int_p2 = (par_kp_offset[0][bin_p + 1] - par_kp_offset[1][bin_p + 1]) * fabs(cost) / 0.1 + par_kp_offset[1][bin_p + 1];
        corr_offset = (int_p2 - int_p1) * (ptrk - p_bin[bin_p]) + int_p1;
        return corr_offset;
    }
    else if (n == 4)
    {
        int tmp = 0;
        double par_p_offset[30][6] = {
            {-0.826976, 1.26319,    0.0168621,  -0.350471,  0.208162,   -0.0422268},
            {-0.655279, 1.17488,    -0.624155,  0.0140827,  0.16105,    -0.119258},
            {-0.316389, 0.678636,   -0.746556,  0.476035,   -0.189224,  0.0435277},
            {0.151886,  -0.0224519, -0.287477,  0.274142,   -0.153749,  0.0734597},
            {0.448262,  -0.362655,  -0.0415588, 0.110434,   -0.0833664, 0.0432422},
            {0.622144,  -0.564482,  0.13891,    -0.0124711, -0.00501423,    -0.00480475},
            {0.70422,   -0.680696,  0.15964,    -0.0148526, -0.0430492, -0.0110849},
            {0.910677,  -0.964731,  0.367566,   -0.107582,  0.0527411,  0},
            {1.0095,    -1.10729,   0.394171,   -0.119731,  0.0157713,  0},
            {1.22613,   -1.46063,   0.636896,   -0.195013,  0.0624278,  0},
            {1.54499,   -1.90653,   0.818578,   -0.312619,  0.0763852,  0},
            {1.98669,   -2.52962,   1.10066,    -0.335902,  0.0635311,  0},
            {2.58488,   -3.35617,   1.53409,    -0.555987,  0.160924,   0},
            {3.2245,    -4.24687,   1.85847,    -0.643344,  0.1679, 0},
            {3.88496,   -4.6193,    1.77549,    -0.482142,  0.0735315,  0},
            {3.72702,   -4.64294,   1.89108,    -0.510733,  0.0914982,  0},
            {2.44759,   -3.709, 1.62889,    -0.530121,  0.122125,   0},
            {1.48439,   -2.52461,   1.10631,    -0.38799,   0.0801175,  0},
            {0.90893,   -1.66597,   0.723416,   -0.209917,  0.0630468,  0},
            {0.560033,  -1.07386,   0.398224,   -0.103849,  0.00455556, 0},
            {0.419805,  -0.809094,  0.318519,   -0.0969289, 0.0118289,  0},
            {0.287319,  -0.436034,  0.0837204,  0.0581752,  -0.0808111, 0.0359385},
            {0.258884,  -0.357818,  0.0475784,  0.0384534,  -0.0402014, -0.000481745},
            {0.21925,   -0.246146,  0.00484353, 0.0531687,  -0.031646,  0.0224279},
            {0.16258,   -0.193137,  -0.0671331, 0.0707033,  -0.0645643, 0.00244481},
            {0.145212,  -0.095416,  -0.152871,  0.155015,   -0.0852023, 0.0387718},
            {0.00824374,    0.111452,   -0.341958,  0.296775,   -0.163263,  0.0740253},
            {-0.348253, 0.699121,   -0.755615,  0.479183,   -0.188644,  0.043546},
            {-0.662414, 1.16992,    -0.664637,  0.0413632,  0.145736,   -0.110265},
            {-0.852959, 1.28192,    -0.0518724, -0.297676,  0.164896,   -0.0247413}
        };
        while (cost >= cosbin[tmp] && tmp != 28)
        {
            tmp++;
        }
        if (tmp == 0)
            tmp += 1;
        double par_cos[6];
        for (int j = 0; j < 6; j++)
        {
            double cosbin_tmp[3] = {cosbin[tmp - 1], cosbin[tmp], cosbin[tmp + 1]};
            double par_p_offset_tmp[3] = {par_p_offset[tmp - 1][j], par_p_offset[tmp][j], par_p_offset[tmp + 1][j]};
            par_cos[j] = interpolation(cost, cosbin_tmp, par_p_offset_tmp);
        }
        double ptrk_tmp = (ptrk - 0.33) / 0.1;
        corr_offset = ROOT::Math::Chebyshev5(ptrk_tmp, par_cos[0], par_cos[1], par_cos[2], par_cos[3], par_cos[4], par_cos[5]);
        if (cost > 0.83 && ptrk < 0.3)
            return 2 * corr_offset;
        return corr_offset;
    }
}

double DedxPID::sigmaCorr(int n, int charge, double ptrk, double cost)
{
    double cosbin[30] = {-0.865, -0.7, -0.5, -0.325, -0.225, -0.19, -0.17, -0.15, -0.13, -0.11, -0.09, -0.07, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.225, 0.325, 0.5, 0.7, 0.865};
    double corr_sigma = 1;
    if (n == 0 || n == 1 || n == 2 && ptrk > 0.20 || n == 3 && (ptrk > 0.43 || ptrk < 0.2 || fabs(cost) > 0.2) || n == 4 && (ptrk > 0.43 || ptrk < 0.2))
        return 1;
    else if (n == 2)
    {
        int tmp = 0;
        double par_pip_sigma[30][8] = {
            {0.980745,  -0.0320824, 0.148076,   -0.0185231, -0.0287245, -0.0146609, 0.0291458,  -0.0213545},
            {0.909279,  0.171988,   0.0389074,  -0.0877263, 0.0104833,  0.0392911,  -0.0306531, 0.00307295},
            {0.952309,  0.168354,   0.0100906,  -0.0781473, 0.0874947,  -0.0399515, -0.0040687, 0.0273449},
            {1.01305,   0.110082,   0.0356032,  -0.0320006, 0.0398265,  0.00283998, -0.00353269,    0.0323046},
            {1.10214,   -0.0464266, 0.0844326,  -0.0133869, -0.0233443, 0.0470427,  -0.0233875, 0.00796884},
            {1.13722,   -0.0694168, 0.0810207,  -0.00461111,    -0.0657331, 0.0435981,  -0.0305237, 0.000891762},
            {1.12494,   -0.0402269, 0.039213,   0.0410197,  -0.0849904, 0.0588127,  -0.0110273, -0.0102289},
            {1.19926,   -0.124029,  0.0915294,  0.0051942,  -0.0546915, 0.0774849,  -0.0186289, 0.00478784},
            {1.29678,   -0.285694,  0.230571,   -0.127148,  0.0295359,  -0.0319002, 0.0336125,  -0.0387693},
            {1.34107,   -0.299083,  0.21921,    -0.115439,  0.0098274,  -0.0225279, 0.0177997,  -0.0240218},
            {1.37443,   -0.292003,  0.170567,   -0.0571655, -0.0168999, 0.0336943,  0.00351935, -0.0213522},
            {1.53528,   -0.469584,  0.286359,   -0.136798,  -0.0158706, -0.00297735,    0.00902674, -0.0365105},
            {1.76254,   -0.697432,  0.422368,   -0.120523,  -0.070091,  0.094326,   -0.030649,  -0.0066399},
            {1.94569,   -0.911353,  0.454219,   0.0156298,  -0.319211,  0.379627,   -0.238878,  0.0963599},
            {2.43193,   -1.66075,   1.11878,    -0.488709,  0.0340719,  0.160446,   -0.0920272, 0.0514607},
            {2.0932,    -1.05913,   0.580453,   -0.0946158, -0.236143,  0.29226,    -0.176364,  0.0320375},
            {1.92749,   -0.929276,  0.553641,   -0.170987,  -0.098254,  0.135239,   -0.075694,  -0.0032753},
            {1.875, -0.988737,  0.679303,   -0.365717,  0.0721803,  0.0117188,  -0.00415287,    -0.0189484},
            {1.58964,   -0.578412,  0.43666,    -0.221821,  0.051573,   -0.00910371,    0.0148332,  -0.036281},
            {1.37654,   -0.303804,  0.192889,   -0.0902226, -0.00981288,    -0.0260383, 0.0428074,  -0.0312006},
            {1.33936,   -0.328366,  0.265274,   -0.150611,  0.0536288,  -0.0520476, 0.0561347,  -0.0397721},
            {1.21193,   -0.144586,  0.0934788,  -0.0178466, -0.0531063, 0.0124997,  0.0153325,  -0.0307769},
            {1.17235,   -0.0870064, 0.0758076,  0.00109802, -0.0561468, 0.030234,   0.011407,   -0.0287855},
            {1.16807,   -0.127506,  0.118003,   -0.0407013, -0.00810034,    -0.00531038,    0.0271216,  -0.0245021},
            {1.10302,   -0.0128453, 0.0366757,  0.0430047,  -0.0677478, 0.0658191,  -0.0213085, 0.0231517},
            {1.09825,   -0.0253014, 0.0807453,  -0.0315216, -0.00550889,    0.0126133,  0.00452082, 0.00183966},
            {1.0199,    0.0937959,  0.0449512,  -0.0539498, 0.0544941,  -0.0364448, 0.0205853,  0.0103639},
            {0.964461,  0.159884,   0.0289211,  -0.0837419, 0.0984221,  -0.0551503, 0.00340403, 0.0213683},
            {0.92457,   0.159631,   0.0423011,  -0.07931,   0.0114676,  0.0488719,  -0.0348463, 0.00461647},
            {1.00597,   -0.0291862, 0.131163,   -0.03148,   -0.024414,  -0.039103,  0.0488831,  -0.0346602}
        };
        while (cost >= cosbin[tmp] && tmp != 28)
        {
            tmp++;
        }
        if (tmp == 0)
            tmp += 1;
        double par_cos[8];
        for (int j = 0; j < 8; j++)
        {
            double cosbin_tmp[3] = {cosbin[tmp - 1], cosbin[tmp], cosbin[tmp + 1]};
            double par_pip_sigma_tmp[3] = {par_pip_sigma[tmp - 1][j], par_pip_sigma[tmp][j], par_pip_sigma[tmp + 1][j]};
            par_cos[j] = interpolation(cost, cosbin_tmp, par_pip_sigma_tmp);
        }
        double ptrk_tmp = (ptrk - 0.17) / 0.1;
        double corr_sigma = ROOT::Math::ChebyshevN(7, ptrk_tmp, par_cos);
        if (corr_sigma < 1)
            return 1;
        else
            return corr_sigma;
    }
    else if (n == 3)
    {
        double par_kp_sigma[3][5] = {
            {1.00,  1.00,   1.00,   1.00,   1.00},
            {1.80,  1.80,   1.51,   1.41,   1.00},
            {1.00,  1.00,   1.00,   1.00,   1.00}
        };
        double p_bin[6] = {0.175, 0.225, 0.275, 0.325, 0.375, 0.425};
        int bin_p = (ptrk - 0.175) / 0.05;
        double int_p1 = (par_kp_sigma[0][bin_p] - par_kp_sigma[1][bin_p]) * fabs(cost) / 0.1 + par_kp_sigma[1][bin_p];
        double int_p2 = (par_kp_sigma[0][bin_p + 1] - par_kp_sigma[1][bin_p + 1]) * fabs(cost) / 0.1 + par_kp_sigma[1][bin_p + 1];
        corr_sigma = (int_p2 - int_p1) * (ptrk - p_bin[bin_p]) + int_p1;
        return corr_sigma;
    }
    else if (n == 4)
    {
        int tmp = 0;
        double par_p_sigma[30][8] = {
            {0.794024,  0.0425693,  0.0236678,  -0.0382406, 0.0695961,  -0.0580967, 0.035697,   -0.0135807},
            {0.832773,  -0.00113245,    -0.031817,  0.0606602,  -0.0447306, -0.00903627,    0.025789,   -0.0195913},
            {0.908858,  -0.087108,  0.0549567,  0.00174534, -0.0270899, 0.0429156,  -0.0280865, 0.0188789},
            {1.04046,   -0.246353,  0.133491,   -0.049544,  0.0180147,  0,  0,  0},
            {1.25697,   -0.492783,  0.244496,   -0.0930121, 0.0267921,  0,  0,  0},
            {1.40495,   -0.656157,  0.341844,   -0.13557,   0.0444445,  0,  0,  0},
            {1.48819,   -0.722884,  0.375376,   -0.133594,  0.0550627,  0,  0,  0},
            {1.73349,   -1.02811,   0.545484,   -0.22501,   0.0867905,  0,  0,  0},
            {1.86727,   -1.11375,   0.566508,   -0.209777,  0.0683113,  0,  0,  0},
            {2.17391,   -1.50475,   0.78278,    -0.317744,  0.0926452,  0,  0,  0},
            {2.4923,    -1.78499,   0.944323,   -0.412239,  0.144967,   0,  0,  0},
            {2.96861,   -2.37577,   1.24553,    -0.50482,   0.135875,   0,  0,  0},
            {3.31789,   -2.67592,   1.3589, -0.552132,  0.210676,   0,  0,  0},
            {3.7896,    -3.26956,   1.68685,    -0.702016,  0.155152,   0,  0,  0},
            {3.86579,   -3.22667,   1.54792,    -0.607399,  0.174962,   0,  0,  0},
            {3.91034,   -3.35332,   1.66152,    -0.618642,  0.152329,   0,  0,  0},
            {3.33904,   -2.66857,   1.27733,    -0.460009,  0.0800364,  0,  0,  0},
            {2.97639,   -2.27119,   1.1101, -0.410801,  0.120688,   0,  0,  0},
            {2.55881,   -1.83093,   0.905367,   -0.339415,  0.107886,   0,  0,  0},
            {2.34426,   -1.713, 0.866848,   -0.370726,  0.087723,   0,  0,  0},
            {1.98031,   -1.24099,   0.640849,   -0.232078,  0.0726222,  0,  0,  0},
            {1.74302,   -0.984163,  0.490949,   -0.172, 0.0443975,  0,  0,  0},
            {1.56317,   -0.802742,  0.388115,   -0.14842,   0.0359668,  0,  0,  0},
            {1.44037,   -0.668254,  0.352312,   -0.120142,  0.0549672,  0,  0,  0},
            {1.34493,   -0.583195,  0.310501,   -0.130395,  0.0447765,  0,  0,  0},
            {1.22836,   -0.433327,  0.229097,   -0.0728195, 0.022962,   0,  0,  0},
            {1.05117,   -0.246895,  0.142671,   -0.0529643, 0.016318,   0,  0,  0},
            {0.909469,  -0.0691198, 0.0377954,  0.019234,   -0.0322931, 0.0460066,  -0.0270032, 0.02252},
            {0.843402,  -0.0106399, -0.0217012, 0.0502854,  -0.0341327, -0.0117776, 0.0292822,  -0.0190088},
            {0.826268,  -0.00178627,    0.0679738,  -0.065918,  0.0696007,  -0.0648257, 0.0328222,  -0.00459817}
        };
        while (cost >= cosbin[tmp] && tmp != 28)
        {
            tmp++;
        }
        if (tmp == 0)
            tmp += 1;
        double par_cos[8];
        for (int j = 0; j < 8; j++)
        {
            double cosbin_tmp[3] = {cosbin[tmp - 1], cosbin[tmp], cosbin[tmp + 1]};
            double par_p_sigma_tmp[3] = {par_p_sigma[tmp - 1][j], par_p_sigma[tmp][j], par_p_sigma[tmp + 1][j]};
            par_cos[j] = interpolation(cost, cosbin_tmp, par_p_sigma_tmp);
        }
        double ptrk_tmp = (ptrk - 0.33) / 0.1;
        corr_sigma = ROOT::Math::ChebyshevN(7, ptrk_tmp, par_cos);
        if (corr_sigma < 1)
            return 1;
        else
            return corr_sigma;
    }
}

double DedxPID::interpolation(double cost, double *costheta, double *par)
{
    double ux = (costheta[0] - costheta[1]) * (costheta[0] - costheta[2]);
    double uy = (costheta[1] - costheta[0]) * (costheta[1] - costheta[2]);
    double uz = (costheta[2] - costheta[0]) * (costheta[2] - costheta[1]);
    double bx = par[0] / ux;
    double by = par[1] / uy;
    double bz = par[2] / uz;
    double c1 = bx + by + bz;
    double c2 = bx * (costheta[1] + costheta[2]) + by * (costheta[0] + costheta[2]) + bz * (costheta[0] + costheta[1]);
    double c3 = bx * costheta[1] * costheta[2] + by * costheta[0] * costheta[2] + bz * costheta[0] * costheta[1];
    return c1 * cost * cost - c2 * cost + c3;
}
//CT
