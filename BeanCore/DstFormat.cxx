#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>

#include <TChain.h>
#include <TFile.h>
#include <TFileCacheRead.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TVector3.h>

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"
#if BOSS_VER >= 661
#include "RootEventData/TEvtNavigator.h"
#endif

#include "DstFormat.h"

ClassImp(DstFormat);

using namespace std;

//--------------------------------------------------------------------
DstFormat::DstFormat()
//--------------------------------------------------------------------
{
   fChain = nullptr;
   m_TEvtHeader    = nullptr;
   m_TDstEvent     = nullptr;
   m_TEvtRecObject = nullptr;
   m_TMcEvent      = nullptr;
   m_TTrigEvent    = nullptr;
   m_TDigiEvent    = nullptr;
   m_THltEvent     = nullptr;
#if BOSS_VER >= 661
   m_TEvtNavigator = nullptr;
#endif
}

//--------------------------------------------------------------------
DstFormat::~DstFormat()
//--------------------------------------------------------------------
{
   delete m_TEvtHeader;
   delete m_TDstEvent;
   delete m_TEvtRecObject;
   delete m_TMcEvent;
   delete m_TTrigEvent;
   delete m_TDigiEvent;
   delete m_THltEvent;
#if BOSS_VER >= 661
   delete m_TEvtNavigator;
#endif
}

//--------------------------------------------------------------------
void DstFormat::Init(TTree *tree)
//--------------------------------------------------------------------
{
   // The Init() function is called when the selector needs to
   // initialize a new tree or chain. Typically here the branch
   // addresses and branch pointers of the tree will be set.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   if( Verbose() ) {
      cout << "DstFormat::Init(tree= " << tree << ")" << endl;
   }

   if (!tree) return;

   fChain = tree;
   fChain->SetMakeClass(0); //!!! this is not generated code !!!
   SetBranches(fChain);

}

//--------------------------------------------------------------------
Bool_t DstFormat::Notify()
//--------------------------------------------------------------------
{
   // The Notify() function is called when a new file is opened.
   // This can be either for a new TTree in a TChain or when
   // a new TTree is started when using PROOF.
   // The return value is currently not used.

   // disable prefetch when reading TTree
   bool disable_pref = true;
   if( disable_pref ) {
      TFile* tf = fChain->GetCurrentFile();
      if( tf ) {
         TFileCacheRead* fc = tf->GetCacheRead();
         if( fc ) {
            fc->SetEnablePrefetching(false);
         }
      }
   }

   if( Verbose() ) {
      cout << "DstFormat::Notify()" << endl;
      TFile* tf = fChain->GetCurrentFile();
      if( tf ) {
         cout << "   -> open file: '" << tf->GetName() << "'" << endl
            << "   -> option= " << tf->GetOption()
            << "; Tree number= " << fChain->GetTreeNumber() << endl;
         // list and contents of StreamerInfo for all objects in file
         // tf->ShowStreamerInfo();
      }
   }

   return kTRUE;
}

//--------------------------------------------------------------------
Int_t DstFormat::GetEntry(Long64_t entry, Int_t getall)
//--------------------------------------------------------------------
{
   // Read contents of entry.
   return (fChain) ? fChain->GetTree()->GetEntry(entry, getall) : 0;
}

//--------------------------------------------------------------------
void DstFormat::SetBranches(TTree *tree)
//--------------------------------------------------------------------
{
   // auxiliary function to set branch addresses and branch pointers
   if( Verbose() ) {
      cout << "DstFormat::SetBranches(tree= " << tree << ")" << endl;
   }

   tree->SetBranchAddress("TEvtHeader", &m_TEvtHeader, &b_TEvtHeader);
   tree->SetBranchAddress("TDstEvent", &m_TDstEvent, &b_TDstEvent);
   tree->SetBranchAddress("TEvtRecObject",
         &m_TEvtRecObject, &b_TEvtRecObject);

   tree->SetBranchAddress("TMcEvent", &m_TMcEvent, &b_TMcEvent);

   tree->SetBranchAddress("TTrigEvent", &m_TTrigEvent, &b_TTrigEvent);

   tree->SetBranchAddress("TDigiEvent", &m_TDigiEvent, &b_TDigiEvent);
   tree->SetBranchAddress("THltEvent", &m_THltEvent, &b_THltEvent);

#if BOSS_VER >= 661
   tree->SetBranchAddress("EventNavigator",
         &m_TEvtNavigator, &b_TEvtNavigator);
#endif
}

//--------------------------------------------------------------------
void DstFormat::Show(Long64_t entry)
//--------------------------------------------------------------------
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if( fChain ) {
      fChain->Show(entry);
   }
}

//--------------------------------------------------------------------
TFile* DstFormat::GetCurrentFile() const
//--------------------------------------------------------------------
{
   return (fChain != 0) ? fChain->GetCurrentFile() : 0;
}

//--------------------------------------------------------------------
void DstFormat::PrintHeader() const
//--------------------------------------------------------------------
{
   cout << " ==================== EVENT HEADER ==================== "
      << endl;
   if( m_TEvtHeader ) {
      int run = m_TEvtHeader->getRunId();
      bool isMC = run < 0;
      cout << " Event# " << m_TEvtHeader->getEventId()
         << " Run# " << run
         << " Tag: " << m_TEvtHeader->getEventTag();
      if( !isMC ) {
         time_t t_event = m_TEvtHeader->time();
         cout << " Time: " << ctime(&t_event);
      } else {
         cout << endl;
      }
#if BOSS_VER >= 655
      if( Verbose() ) {
         cout << " Flag1: " << m_TEvtHeader->getFlag1()
            << " Flag2: " << m_TEvtHeader->getFlag2() << endl;
      }
#endif
#if BOSS_VER >= 705
      if( Verbose() ) {
         cout << " ETS T1: " << m_TEvtHeader->getEtsT1()
            << " ETS T2: " << m_TEvtHeader->getEtsT2() << endl;
      }
#endif
   }
   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintDst() const
//--------------------------------------------------------------------
{
   if( !m_TDstEvent ) return;

   const TObjArray* m_mdcTrackCol = m_TDstEvent->getMdcTrackCol();
   const TObjArray* m_mdcDedxCol = m_TDstEvent->getMdcDedxCol();
   const TObjArray* m_mdcKalTrackCol=m_TDstEvent->getMdcKalTrackCol();
   const TObjArray* m_tofTrackCol = m_TDstEvent->getTofTrackCol();
   const TObjArray* m_extTrackCol = m_TDstEvent->getExtTrackCol();
   const TObjArray* m_emcTrackCol = m_TDstEvent->getEmcTrackCol();
   const TObjArray* m_mucTrackCol = m_TDstEvent->getMucTrackCol();

   int NmdcTracks    = m_mdcTrackCol->GetEntries();
   int NmdcDedx      = m_mdcDedxCol->GetEntries();
   int NmdcKalTracks = m_mdcKalTrackCol->GetEntries();
   int NtofTracks    = m_tofTrackCol->GetEntries();
   int NextTracks    = m_extTrackCol->GetEntries();
   int NemcTracks    = m_emcTrackCol->GetEntries();
   int NmucTracks    = m_mucTrackCol->GetEntries();

   bool isNULL = (NemcTracks==0 && NmdcTracks==0 && NtofTracks==0 &&
         NmucTracks==0 && NmdcDedx==0 && NextTracks==0 &&
         NmdcKalTracks==0);

   if( !isNULL || Verbose() ) {
      cout
         << " ===================== Dst Event ====================== "
         << endl;
      cout << " mdcTrackColSize= "    << NmdcTracks    << endl;
      cout << " mdcDedxColSize= "     << NmdcDedx      << endl;
      cout << " mdcKalTrackColSize= " << NmdcKalTracks << endl;
      cout << " tofTrackColSize= "    << NtofTracks    << endl;
      cout << " extTrackColSize= "    << NextTracks    << endl;
      cout << " emcTrackColSize= "    << NemcTracks    << endl;
      cout << " mucTrackColSize= "    << NmucTracks    << endl;
   }

   if( isNULL )         return;
   if( NmdcTracks )     PrintMdcTracks(m_mdcTrackCol);
   if( NmdcDedx )       PrintMdcDedx(m_mdcDedxCol);
   if( NmdcKalTracks )  PrintMdcKalTracks(m_mdcKalTrackCol);
   if( NtofTracks )     PrintTofTrack(m_tofTrackCol);
   if( NextTracks )     PrintExtTrack(m_extTrackCol);
   if( NemcTracks )     PrintEmcTrack(m_emcTrackCol);
   if( NmucTracks )     PrintMucTrack(m_mucTrackCol);

   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintMdcTracks(const TObjArray* m_mdcTrackCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_mdcTrackCol: " << endl;
   TIter mdcTrackIter(m_mdcTrackCol);
   while( TMdcTrack* mdcTrack = (TMdcTrack*)mdcTrackIter.Next() ) {
      cout << " +++++ trackId= " << mdcTrack->trackId() << endl;
      cout << " FitQuality= " << mdcTrack->stat()
         << " chi2= " << mdcTrack->chi2()
         << " ndof = " << mdcTrack->ndof() << endl;
      cout << " charge= " << mdcTrack->charge()
         << " Pxy= " << mdcTrack->pxy()
         << " Pz= " << mdcTrack->pz()
         << " P= " << mdcTrack->p() << endl;
      cout << " Theta= " << mdcTrack->theta()
         << " Phi= " << mdcTrack->phi() << endl;

      if( Verbose() ) {
         cout << " Helix parameters: d0, phi0, 1/Pt, dz, tan(lambda)"
            << endl;
         for(int i = 0; i < 5; i++) {
            cout << setw(12) << mdcTrack->helix(i) << " : ";
         }
         cout << endl;
         cout << " Helix error matrix: " << endl;
         for(int k = 0, i = 0; i < 5; i++) {
            if( i != 0 ) cout << setw(15*i) << " : ";
            for(int j = i; j < 5; j++,k++) {
               cout << setw(12) << mdcTrack->err(k) << " : ";
            }
            cout << endl;
         }
         cout << " Nster= " << mdcTrack->nster()
            << " FirstLayer= " << mdcTrack->firstLayer()
            << " LastLayer= " << mdcTrack->lastLayer() << endl;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMdcDedx(const TObjArray* m_mdcDedxCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_mdcDedxCol: " << endl;
   TIter mdcDedxIter(m_mdcDedxCol);
   while( TMdcDedx* mdcDedx = (TMdcDedx*)mdcDedxIter.Next() ) {
      cout << " +++++ trackId= " << mdcDedx->trackId() << endl;
      cout << " Status= " << mdcDedx->status()
         << " numGoodHits= " << mdcDedx->numGoodHits()
         << " numTotalHits= " << mdcDedx->numTotalHits() << endl;
      cout << " probPH= " << mdcDedx->probPH()
         << " normPH= " << mdcDedx->normPH()
         << " particleId= " << mdcDedx->particleId()
         << " truncAlg= " << mdcDedx->truncAlg() << endl;

      if( Verbose() ) {
         cout << " errorPH= " << mdcDedx->errorPH()
            << " twentyPH= " << mdcDedx->twentyPH() << endl;
         cout << " chi2:  e, mu, pi, K, prot " << endl;
         for(int i = 0; i < 5; i++) {
            cout << setw(12) << mdcDedx->chi(i) << " : ";
         }
         cout << endl;
         // cout << " chi2= "
            // << mdcDedx->chiE() << " : "
            // << mdcDedx->chiMu() << " : "
            // << mdcDedx->chiPi() << " : "
            // << mdcDedx->chiK() << " : "
            // << mdcDedx->chiP() << endl;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMdcKalTracks(const TObjArray*
      m_mdcKalTrackCol) const
//--------------------------------------------------------------------
{
   // pointers to finctions:
   typedef Double_t (TMdcKalTrack::*Helix)(Int_t i) const;
   typedef Double_t (TMdcKalTrack::*Error)(Int_t i, Int_t j) const;

   // arrays of function for easy access
   const static string pid_name[5] =
   {"Electron","Muon","Pion","Kaon","Proton"};
   const static Helix Zhel[5] = {       &TMdcKalTrack::getZHelixE,
      &TMdcKalTrack::getZHelixMu,     &TMdcKalTrack::getZHelix,
      &TMdcKalTrack::getZHelixK,      &TMdcKalTrack::getZHelixP };
   const static Error Zerr[5] = {       &TMdcKalTrack::getZErrorE,
      &TMdcKalTrack::getZErrorMu,     &TMdcKalTrack::getZError,
      &TMdcKalTrack::getZErrorK,      &TMdcKalTrack::getZErrorP };

   const static Helix Fhel[5] = {       &TMdcKalTrack::getFHelixE,
      &TMdcKalTrack::getFHelixMu,     &TMdcKalTrack::getFHelix,
      &TMdcKalTrack::getFHelixK,      &TMdcKalTrack::getFHelixP };
   const static Error Ferr[5] = {       &TMdcKalTrack::getFErrorE,
      &TMdcKalTrack::getFErrorMu,     &TMdcKalTrack::getFError,
      &TMdcKalTrack::getFErrorK,      &TMdcKalTrack::getFErrorP };


   cout << endl << " L________> m_mdcKalTrackCol: " << endl;
   TIter mdcKalTrackIter(m_mdcKalTrackCol);
   TMdcKalTrack* mdcKalTrack = nullptr;
   while( (mdcKalTrack = (TMdcKalTrack*)mdcKalTrackIter.Next()) ) {
      cout << " +++++ trackId= " << mdcKalTrack->getTrackId() << endl;

#if BOSS_VER < 663
      if( Verbose() ) {
         // use Pion (ip=2) for this variables
         cout << " Nster= " << mdcKalTrack->getNster(2)
            << " FirstLayer= " << mdcKalTrack->getFirstLayer(2)
            << " LastLayer= " << mdcKalTrack->getLastLayer(2) << endl;
      }
#endif

      for(int ip = 0; ip < 5; ip++) { // e/mu/pi/K/p
         if( ip == 2 || Verbose() ) { // always print for pi
            cout << "  ----------------- " << endl;
            cout << " | " << setw(8) << pid_name[ip] << " fit    |\n";
            cout << "  ----------------- " << endl;
            cout
#if BOSS_VER >= 655
               << " FitQuality(1:2) = " << mdcKalTrack->getStat(ip)
               << " : " << mdcKalTrack->getStat2(ip)
#endif
               << "   chi2= " << mdcKalTrack->getChisq(ip)
               << " ndof = " << mdcKalTrack->getNdf(ip) << endl;

            cout << " Helix parameters at zero point:" << endl;
            Helix fp = Zhel[ip];// pointer to TMdcKalTrack::getZHelix
            for(int i = 0; i < 5; i++) {
               cout << setw(12) << (mdcKalTrack->*fp)(i) << " : ";
            }
            cout << endl;
            Error fpe = Zerr[ip];// pointer to TMdcKalTrack::getZError
            cout << " Helix error matrix at zero point: " << endl;
            for(int i = 0; i < 5; i++) {
               if( i != 0 ) cout << setw(15*i) << " : ";
               for(int j = i; j < 5; j++) {
                  cout << setw(12) << (mdcKalTrack->*fpe)(i,j) << " : ";
               }
               cout << endl;
            }
         }

         if( Verbose() ) {
            cout << " Helix parameters at first Mdc hit:" << endl;
            Helix fp = Fhel[ip];// pointer to TMdcKalTrack::getFHelix
            for(int i = 0; i < 5; i++) {
               cout << setw(12) << (mdcKalTrack->*fp)(i) << " : ";
            }
            cout << endl;
            Error fpe = Ferr[ip];// pointer to TMdcKalTrack::getFError
            cout << " Helix error matrix at zero point: " << endl;
            for(int i = 0; i < 5; i++) {
               if( i != 0 ) cout << setw(15*i) << " : ";
               for(int j = i; j < 5; j++) {
                  cout << setw(12) << (mdcKalTrack->*fpe)(i,j)
                     << " : ";
               }
               cout << endl;
            }
         }
      } // end of for(e/mu/pi/K/p)
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintTofTrack(const TObjArray* m_tofTrackCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_tofTrackCol: " << endl;
   TIter tofTrackIter(m_tofTrackCol);
   while( TTofTrack* tofTrack = (TTofTrack*)tofTrackIter.Next() ) {
      cout << " +++++ trackId= " << tofTrack->trackID()
         << " tofTrackID= " << tofTrack->tofTrackID()
         << " tofID= " << tofTrack->tofID() << endl;
      cout << " Status= " << hex << (tofTrack->status()) << endl
         << "\t\t" << " layers = "
         << ( (tofTrack->status() & 0x000000C0) >> 6) << endl
         << "\t\t" << " is east = "
         << ( (tofTrack->status() & 0x00000020) >> 5) << endl
         << "\t\t" << " is barrel = "
         << ( (tofTrack->status() & 0x00000010) >> 4) << endl
         << "\t\t" << " is counter (0 - bad) = "
         <<  ( (tofTrack->status() & 0x00000004) >> 2) << endl
         << dec << " quality= " << tofTrack->quality() << endl;

      cout << " tof= " << tofTrack->tof()
         // << " +/- " << sqrt(tofTrack->errtof())
         << " path= " << tofTrack->path()
         << " beta= " << tofTrack->beta() << endl;

      cout << " texp:  e, mu, pi, K, prot " << endl;
      for(int i = 0; i < 5; i++) {
         cout << setw(12) << tofTrack->texp(i) << " : ";
      }
      cout << endl;

      if( Verbose() ) {
         streamsize oldprc = cout.precision(4);
         cout << " Time offset: e,mu,pi,K,prot,pbar" << endl;
         for(int i = 0; i < 6; i++) {
            cout << setw(10) << tofTrack->toffset(i) << " : ";
         }
         cout << endl;
         cout << " Time resolution: e,mu,pi,K,prot,pbar" << endl;
         for(int i = 0; i < 6; i++) {
            cout << setw(10) << tofTrack->sigma(i) << " : ";
         }
         cout << endl;
         cout.precision(oldprc);

         cout << " t0= " << tofTrack->t0()
            << " errt0= " << tofTrack->errt0() << endl;

         cout << " ph= " << tofTrack->ph() << endl;

         cout << " zrhit= " << tofTrack->zrhit()
            << " errzrhit= " << tofTrack->errz() << endl;
         cout << " phi= " << tofTrack->phi()
            << " errphi= " << tofTrack->errphi() << endl;
         cout << " energy= " << tofTrack->energy()
            << " errenergy= " << tofTrack->errenergy() << endl;
      }
      // cout << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintExtTrack(const TObjArray* m_extTrackCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_extTrackCol: " << endl;
   TIter extTrackIter(m_extTrackCol);
   while( TExtTrack* extTrack = (TExtTrack*)extTrackIter.Next() ) {
      cout << " +++++ trackId= " << extTrack->GetTrackId() << endl;
      cout << " Tof1: x,y,z= " << extTrack->GetTof1PositionX()
         << ", " << extTrack->GetTof1PositionY()
         << ", " << extTrack->GetTof1PositionZ() << endl;
      cout << "       Px,Py,Pz= " << extTrack->GetTof1MomentumX()
         << ", " << extTrack->GetTof1MomentumY()
         << ", " << extTrack->GetTof1MomentumZ() << endl;
      cout << "       Tof= " << extTrack->GetTof1()
         << "  Path= " << extTrack->GetTof1Path() << endl;
      if( Verbose() ) {
         cout << "       VolumeName= " <<extTrack->GetTof1VolumeName()
            << " VolumeNumber= " << extTrack->GetTof1VolumeNumber()
            << endl;
         cout << "       SigmaAlong: X,Y,Z,T= "
            << extTrack->GetTof1PosSigmaAlongX() << ", "
            << extTrack->GetTof1PosSigmaAlongY() << ", "
            << extTrack->GetTof1PosSigmaAlongZ() << ", "
            << extTrack->GetTof1PosSigmaAlongT() << endl;

         streamsize oldprc = cout.precision(4);
         cout << " Error matrix 6x6 (x,p): " << endl;
         for(int i = 0; i < 6; i++) {
            if( i != 0 ) cout << setw(13*i) << " : ";
            for(int j = i; j < 6; j++) {
               cout << setw(10) << extTrack->GetTof1ErrorMatrix(i,j)
                  << " : ";
            }
            cout << endl;
         }
         cout.precision(oldprc);
      }

      cout << " Tof2: x,y,z= " << extTrack->GetTof2PositionX()
         << ", " << extTrack->GetTof2PositionY()
         << ", " << extTrack->GetTof2PositionZ() << endl;
      cout << "       Px,Py,Pz= " << extTrack->GetTof2MomentumX()
         << ", " << extTrack->GetTof2MomentumY()
         << ", " << extTrack->GetTof2MomentumZ() << endl;
      cout << "       Tof= " << extTrack->GetTof2()
         << "  Path= " << extTrack->GetTof2Path() << endl;
      if( Verbose() ) {
         cout << "       VolumeName= " <<extTrack->GetTof2VolumeName()
            << " VolumeNumber= " << extTrack->GetTof2VolumeNumber()
            << endl;
         cout << "       SigmaAlong: X,Y,Z,T= "
            << extTrack->GetTof2PosSigmaAlongX() << ", "
            << extTrack->GetTof2PosSigmaAlongY() << ", "
            << extTrack->GetTof2PosSigmaAlongZ() << ", "
            << extTrack->GetTof2PosSigmaAlongT() << endl;

         streamsize oldprc = cout.precision(4);
         cout << " Error matrix 6x6 (x,p): " << endl;
         for(int i = 0; i < 6; i++) {
            if( i != 0 ) cout << setw(13*i) << " : ";
            for(int j = i; j < 6; j++) {
               cout << setw(10) << extTrack->GetTof2ErrorMatrix(i,j)
                  << " : ";
            }
            cout << endl;
         }
         cout.precision(oldprc);
      }

      cout << "  EMC: x,y,z= " << extTrack->GetEmcPositionX()
         << ", " << extTrack->GetEmcPositionY()
         << ", " << extTrack->GetEmcPositionZ() << endl;
      cout << "       Px,Py,Pz= " << extTrack->GetEmcMomentumX()
         << ", " << extTrack->GetEmcMomentumY()
         << ", " << extTrack->GetEmcMomentumZ() << endl;
      cout << "       emcPath= " << extTrack->emcPath() << endl;
      if( Verbose() ) {
         cout << "       VolumeName= " <<extTrack->GetEmcVolumeName()
            << " VolumeNumber= " << extTrack->GetEmcVolumeNumber()
            << endl;
         cout << "       SigmaAlong: Theta,Phi= "
            << extTrack->GetEmcPosSigmaAlongTheta() << ", "
            << extTrack->GetEmcPosSigmaAlongPhi() << endl;

         streamsize oldprc = cout.precision(4);
         cout << " Error matrix 6x6 (x,p): " << endl;
         for(int i = 0; i < 6; i++) {
            if( i != 0 ) cout << setw(13*i) << " : ";
            for(int j = i; j < 6; j++) {
               cout << setw(10) << extTrack->GetEmcErrorMatrix(i,j)
                  << " : ";
            }
            cout << endl;
         }
         cout.precision(oldprc);
      }

      cout << "  Muc: x,y,z= " << extTrack->GetMucPositionX()
         << ", " << extTrack->GetMucPositionY()
         << ", " << extTrack->GetMucPositionZ() << endl;
      cout << "       Px,Py,Pz= " << extTrack->GetMucMomentumX()
         << ", " << extTrack->GetMucMomentumY()
         << ", " << extTrack->GetMucMomentumZ() << endl;
      if( Verbose() ) {
         cout << "       VolumeName= " <<extTrack->GetMucVolumeName()
            << " VolumeNumber= " << extTrack->GetMucVolumeNumber()
            << endl;
         cout << "       SigmaAlong: X,Y,Z,T= "
            << extTrack->GetMucPosSigmaAlongX() << ", "
            << extTrack->GetMucPosSigmaAlongY() << ", "
            << extTrack->GetMucPosSigmaAlongZ() << ", "
            << extTrack->GetMucPosSigmaAlongT() << endl;

         streamsize oldprc = cout.precision(4);
         cout << " Error matrix 6x6 (x,p): " << endl;
         for(int i = 0; i < 6; i++) {
            if( i != 0 ) cout << setw(13*i) << " : ";
            for(int j = i; j < 6; j++) {
               cout << setw(10) << extTrack->GetMucErrorMatrix(i,j)
                  << " : ";
            }
            cout << endl;
         }
         cout.precision(oldprc);
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintEmcTrack(const TObjArray* m_emcTrackCol) const
//--------------------------------------------------------------------
{
   int MaxEmcTracks = m_emcTrackCol->GetEntries();
   int i_emc_trk = 0;
   if( !Verbose() ) {
      int NmdcTracks = m_TDstEvent->getMdcTrackCol()->GetEntries();
      // max number of emc clusters to print
      MaxEmcTracks = NmdcTracks+3;
   }

   cout << endl << " L________> m_emcTrackCol: " << endl;

   // print first tracks with trackId() >= 0
   for( int condition = 0; condition < 2; condition++) {
      TIter emcTrackIter(m_emcTrackCol);
      while( TEmcTrack* emcTrack = (TEmcTrack*)emcTrackIter.Next() ) {
         if( condition == 0 ) {
            if( emcTrack->trackId() < 0 ) continue;
         } else {
            if( emcTrack->trackId() >= 0 ) continue;
            if( i_emc_trk >= MaxEmcTracks ) {
               cout << " ... skip "
                  << m_emcTrackCol->GetEntries() - i_emc_trk
                  << " EMC neutral tracks ..." << endl;
               break;
            }
         }
         i_emc_trk++;

         cout << " +++++ trackId= " << emcTrack->trackId()
            << " Status: ";
         if( emcTrack->status() == 1 ) {
            cout << "single seed cluster; ";
         } else if( emcTrack->status() == 2 ) {
            cout << "splitted from multi-seeds cluster; ";
         } else {
            cout << "unknown= " << emcTrack->status();
         }
         cout << endl;

         cout << " numHits= " << emcTrack->numHits()
            << " cell= " << emcTrack->cellId()
            << " module= " << emcTrack->module() << endl;
         // Module: 0:east endcap; 1:barrel; 2:west endcap
         cout << " x,y,z= " << emcTrack->x()
            << ", " << emcTrack->y()
            << ", " << emcTrack->z() << endl;
         cout << " energy= " << emcTrack->energy()
            << " +/- " << emcTrack->dE() << endl;

         if( Verbose() ) {
            cout << " eSeed= " << emcTrack->eSeed()
               << " e3x3= " << emcTrack->e3x3()
               << " e5x5= " << emcTrack->e5x5() << endl;
            cout << " dTheta= " << emcTrack->dtheta()
               << " dPhi= " << emcTrack->dphi() << endl;
            cout << " time= " << emcTrack->time() << endl; // It is not set !!
            cout << " secondMoment= " << emcTrack->secondMoment()
               << " latMom= " << emcTrack->latMoment()
               << " a20Mom= " << emcTrack->a20Moment()
               << " a42Mom= " << emcTrack->a42Moment()  << endl;
            cout << " Error matrix 3x3 (x,y,z): " << endl;
            double errmat[3][3] = {
               {emcTrack->err(0),emcTrack->err(3),emcTrack->err(4)},
               {              0 ,emcTrack->err(1),emcTrack->err(5)},
               {              0 ,              0 ,emcTrack->err(2)}
            };
            for(int i = 0; i < 3; i++) {
               if( i != 0 ) cout << setw(15*i) << " : ";
               for(int j = i; j < 3; j++) {
                  cout << setw(12) << errmat[i][j] << " : ";
               }
               cout << endl;
            }
         }
      } // end of while
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMucTrack(const TObjArray* m_mucTrackCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_mucTrackCol: " << endl;
   TIter mucTrackIter(m_mucTrackCol);
   while( TMucTrack* mucTrack = (TMucTrack*)mucTrackIter.Next() ) {
      cout << " +++++ trackId= " << mucTrack->trackId() << endl;
      cout << " Status= " << mucTrack->status()
         << " id= " << mucTrack->id()
         << " type= " << mucTrack->type()
         << " numHits= " << mucTrack->numHits() << endl;

      if( Verbose() ) {
         cout << " StartPart= " << mucTrack->startPart()
            << " EndPart= " << mucTrack->endPart() << endl;
         cout << " BarrelLastLayer " << mucTrack->brLastLayer()
            << " EndcapLastLayer= " << mucTrack->ecLastLayer()
            << endl;
         cout << " numLayers= " << mucTrack->numLayers()
            << " maxHitsInLayer= " << mucTrack->maxHitsInLayer()
            << endl;
      }

      cout << " depth= " << mucTrack->depth()
         << " chi2= " << mucTrack->chi2()
         << " dof= " << mucTrack->dof()
         << " rms= " << mucTrack->rms() << endl;
      cout << " x,y,z= " << mucTrack->xPos()
         << ", " << mucTrack->yPos()
         << ", " << mucTrack->zPos() << endl;
      cout << " Px,Py,Pz= " << mucTrack->px()
         << ", " << mucTrack->py()
         << ", " << mucTrack->pz() << endl;

      if( Verbose() ) {
         cout << " Sigma x,y,z= " << mucTrack->xPosSigma()
            << ", " << mucTrack->yPosSigma()
            << ", " << mucTrack->zPosSigma() << endl;
         cout << " distance= " << mucTrack->distance()
            << " deltaPhi= " << mucTrack->deltaPhi() << endl;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintEvtRec() const
//--------------------------------------------------------------------
{
   if( !m_TEvtRecObject ) return;

   const TEvtRecEvent* m_evtRecEvent =
      m_TEvtRecObject->getEvtRecEvent();
   const TObjArray* m_evtRecTrackCol =
      m_TEvtRecObject->getEvtRecTrackCol();

   const TEvtRecPrimaryVertex* m_evtRecPrimaryVertex =
      m_TEvtRecObject->getEvtRecPrimaryVertex();
   const TObjArray* m_evtRecVeeVertexCol =
      m_TEvtRecObject->getEvtRecVeeVertexCol();

   const TObjArray* m_evtRecPi0Col =
      m_TEvtRecObject->getEvtRecPi0Col();
   const TObjArray* m_evtRecEtaToGGCol =
      m_TEvtRecObject->getEvtRecEtaToGGCol();
   const TObjArray* m_evtRecDTagCol =
      m_TEvtRecObject->getEvtRecDTagCol();

   int NevtRecTracks = m_evtRecTrackCol->GetEntries();
   int NevtRecVeeVertexs = m_evtRecVeeVertexCol->GetEntries();
   int NevtRecPi0 = m_evtRecPi0Col->GetEntries();
   int NevtRecEtaToGG = m_evtRecEtaToGGCol->GetEntries();
   int NevtRecDTag = m_evtRecDTagCol->GetEntries();

   bool isNULL = (NevtRecTracks==0 && NevtRecVeeVertexs==0 &&
         NevtRecPi0==0 && NevtRecEtaToGG==0 && NevtRecDTag==0);
   bool isVtxValid = m_evtRecPrimaryVertex->isValid();

   if( !isNULL || Verbose() ) {
      cout << " =================== EvtRec Object ==================== "
         << endl;
      cout << " evtRecTrackColSize= "     << NevtRecTracks     << endl;
      cout << " valid primary vertex= "   << isVtxValid        << endl;
      cout << " evtRecVeeVertexColSize= " << NevtRecVeeVertexs << endl;
      cout << " evtRecPi0ColSize= "       << NevtRecPi0        << endl;
      cout << " evtRecEtaToGGColSize= "   << NevtRecEtaToGG    << endl;
      cout << " evtRecDTagColSize= "      << NevtRecDTag       << endl;
   }

   if( isNULL ) return; // nothing to print

   // TEvtRecEvent
   cout << endl << " L________> m_evtRecEvent: " << endl;
   cout << " Tracks: total= " << m_evtRecEvent->totalTracks()
      << " charged= " << m_evtRecEvent->totalCharged()
      << " neutral= " << m_evtRecEvent->totalNeutral() << endl;
   if( Verbose() ) { // uninitialized variables ??
      cout << " Numbers of: Vee= " << m_evtRecEvent->numberOfVee()
         << " Pi0= " << m_evtRecEvent->numberOfPi0()
         << " EtaToGG= " << m_evtRecEvent->numberOfEtaToGG()
         << " DTag= " << m_evtRecEvent->numberOfDTag()
         << " ?? probably wrong ?? " << endl;
   }

   if( isVtxValid ) {
      cout << endl << " L________> m_evtRecPrimaryVertex: " << endl;
      cout << " +++++ nTracks= " << m_evtRecPrimaryVertex->nTracks()
         << " trkID's= ";
      for( const auto& ti : m_evtRecPrimaryVertex->trackIdList() ) {
         cout << ti << " ";
      }
      cout << endl;
      cout << "       fitMethod= "
         << m_evtRecPrimaryVertex->fitMethod()
         << " chi2= " << m_evtRecPrimaryVertex->chi2()
         << " ndof= " << m_evtRecPrimaryVertex->ndof() << endl;
      cout << "       (X,Y,Z)= "
         << m_evtRecPrimaryVertex->vertex(0) << ", "
         << m_evtRecPrimaryVertex->vertex(1) << ", "
         << m_evtRecPrimaryVertex->vertex(2) << endl;

      if( Verbose() ) {
         cout << "       Error matrix 3x3 (X,Y,Z): " << endl;
         for(int i = 0, k = 0; i < 3; i++) {
            cout << setw(6) << " ";
            if( i != 0 ) cout << setw(15*i) << " : ";
            for(int j = i; j < 3; j++,k++) {
               cout << setw(12)
                  << m_evtRecPrimaryVertex->errorVertex(k) << " : ";
            }
            cout << endl;
         }
      }
   }

   if( NevtRecTracks )     PrintRecTrack(m_evtRecTrackCol);
   if( NevtRecVeeVertexs ) PrintRecVeeVertex(m_evtRecVeeVertexCol);
   if( NevtRecPi0 )        PrintRecPi0(m_evtRecPi0Col);
   if( NevtRecEtaToGG )    PrintRecEtaToGG(m_evtRecEtaToGGCol);
   if( NevtRecDTag )       PrintRecDTag(m_evtRecDTagCol);
   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintRecTrack(const TObjArray* m_evtRecTrackCol) const
//--------------------------------------------------------------------
{
   // max number of rec. tracks to print
   const static int max_rectrk = 20;

   cout << endl << " L________> m_evtRecTrackCol: " << endl;
   TIter evtRecTrackIter(m_evtRecTrackCol);
   int i_rec_trk = 0;
   TEvtRecTrack* evtRecTrack = nullptr;
   while( (evtRecTrack = (TEvtRecTrack*)evtRecTrackIter.Next()) ) {
      cout << " +++++ trackId= " << evtRecTrack->trackId()
         << " partId= " << evtRecTrack->partId()  // -1 ???
         << " quality= " << evtRecTrack->quality()
         << endl;
      cout << "  Id:"
         << " mdcTrack="       << evtRecTrack->mdcTrackId()
         << " mdcDedx="        << evtRecTrack->mdcDedxId()
         << " mdcKalTrack="    << evtRecTrack->mdcKalTrackId()
         << " extTrack="       << evtRecTrack->extTrackId()
         << " emcShower="      << evtRecTrack->emcShowerId()
         << " mucTrack="       << evtRecTrack->mucTrackId()
         << endl;
      const vector<Int_t>& tof_id = evtRecTrack->tofTrackIds();
      if( !tof_id.empty() ) {
         cout << "      tofTracks= ";
         copy(tof_id.begin(),tof_id.end(),
               ostream_iterator<int>(cout, " "));
         cout << endl;
      }
      if( ++i_rec_trk >= max_rectrk ) {
         int n_skip = m_evtRecTrackCol->GetEntries() - i_rec_trk;
         if( n_skip > 0 ) {
            cout << " ... skip " << n_skip << " tracks" << endl;
         }
         break;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintRecVeeVertex(const TObjArray*
      m_evtRecVeeVertexCol) const
//--------------------------------------------------------------------
{
   const static string pid_name[6] = {
      "Electron", "Muon", "Pion", "Kaon", "Proton", "---" };

   cout << endl << " L________> m_evtRecVeeVertexCol: " << endl;
   TIter evtRecVeeVertexIter(m_evtRecVeeVertexCol);
   TEvtRecVeeVertex* evtRecVeeVertex = nullptr;
   while( (evtRecVeeVertex =
            (TEvtRecVeeVertex*)evtRecVeeVertexIter.Next()) ) {
      cout << " +++++ vertexId= " << evtRecVeeVertex->vertexId()
         << endl;
      cout << " V-type= " << setw(8);
      switch( evtRecVeeVertex->vertexType() ) {
         case  0:   cout << "Ks";      break;
         case  1:   cout << "Lambda";  break;
         case  2:   cout << "gamma";   break;
         default:   cout << "unknown"; break;
      }
      cout << " : mass= " << evtRecVeeVertex->mass()
         << " chi2= " << evtRecVeeVertex->chi2()
         << " ndof= " << evtRecVeeVertex->ndof() << endl;

      streamsize oldprc = cout.precision(3);
      cout << "  parameters: (px, py, pz, E, x, y, z)" << endl;
      for(int i = 0; i < 7; i++) {
         cout << setw(9) << evtRecVeeVertex->w(i) << ": ";
      }
      cout << endl;

      if( Verbose() ) {
         cout << " Error matrix (7x7): " << endl;
         for(int i = 0, k = 0; i < 7; i++) {
            if( i != 0 ) cout << setw(11*i) << ": ";
            for(int j = i; j < 7; j++, k++) {
               cout << setw(9) << evtRecVeeVertex->Ew(k) << ": ";
            }
            cout << endl;
         }
      }
      cout.precision(oldprc);

      cout << " Daughters: " << pid_name[evtRecVeeVertex->pair(0)]
         << " + " << pid_name[evtRecVeeVertex->pair(1)]
         << "   with EvtRecTrkId: " << evtRecVeeVertex->daughter(0)
         << " + " << evtRecVeeVertex->daughter(1) << endl;

      if( Verbose() ) {
         cout << "            nTracks= "
            << evtRecVeeVertex->nTracks()
            << " nCharge= " << evtRecVeeVertex->nCharge()
            << " nNeutral= " << evtRecVeeVertex->nTracks() -
            evtRecVeeVertex->nCharge() << endl;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintRecPi0(const TObjArray* m_evtRecPi0Col) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_evtRecPi0Col: " << endl;
   TIter evtRecPi0Iter(m_evtRecPi0Col);
   while( TEvtRecPi0* evtRecPi0=(TEvtRecPi0*)evtRecPi0Iter.Next() ) {
      cout << " +++++ unconMass= " << evtRecPi0->unconMass()
         << " chi2= " << evtRecPi0->chisq() << endl;
      cout << "      hiP(x,y,z,e)_fit= "
         << evtRecPi0-> hiPxfit() << ", "
         << evtRecPi0-> hiPyfit() << ", "
         << evtRecPi0-> hiPzfit() << ", "
         << evtRecPi0-> hiPefit() << endl;
      cout << "      loP(x,y,z,e)_fit= "
         << evtRecPi0-> loPxfit() << ", "
         << evtRecPi0-> loPyfit() << ", "
         << evtRecPi0-> loPzfit() << ", "
         << evtRecPi0-> loPefit() << endl;
      cout << "      hiEnGamma= " << evtRecPi0->hiEnGamma()
         << " loEnGamma= " << evtRecPi0->loEnGamma() << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintRecEtaToGG(const TObjArray*
      m_evtRecEtaToGGCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_evtRecEtaToGGCol: " << endl;
   TIter evtRecEtaToGGIter(m_evtRecEtaToGGCol);
   TEvtRecEtaToGG* evtRecEtaToGG = nullptr;
   while((evtRecEtaToGG=(TEvtRecEtaToGG*)evtRecEtaToGGIter.Next())) {
      cout << " +++++ unconMass= " << evtRecEtaToGG->unconMass()
         << " chi2= " << evtRecEtaToGG->chisq() << endl;
      cout << "      hiP(x,y,z,e)_fit= "
         << evtRecEtaToGG-> hiPxfit() << ", "
         << evtRecEtaToGG-> hiPyfit() << ", "
         << evtRecEtaToGG-> hiPzfit() << ", "
         << evtRecEtaToGG-> hiPefit() << endl;
      cout << "      loP(x,y,z,e)_fit= "
         << evtRecEtaToGG-> loPxfit() << ", "
         << evtRecEtaToGG-> loPyfit() << ", "
         << evtRecEtaToGG-> loPzfit() << ", "
         << evtRecEtaToGG-> loPefit() << endl;
      cout << "      hiEnGamma= " << evtRecEtaToGG->hiEnGamma()
         << " loEnGamma= " << evtRecEtaToGG->loEnGamma() << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintRecDTag(const TObjArray* m_evtRecDTagCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_evtRecDTagCol: " << endl;
   TIter evtRecDTagIter(m_evtRecDTagCol);
   TEvtRecDTag* evtRecDTag = 0;
   while( (evtRecDTag = (TEvtRecDTag*) evtRecDTagIter.Next()) ) {
      if (evtRecDTag->type() == 0) continue;
      cout << " +++++ decayMode= " << evtRecDTag->decayMode()
         << " type= " << evtRecDTag->type() << endl;
      cout << "      beamE= " << evtRecDTag->beamE()
         << " mass= " << evtRecDTag->mass()
         << " mBC= " << evtRecDTag->mBC()
         << " deltaE= " << evtRecDTag->deltaE() << endl;
      cout << "      charge= " << evtRecDTag->charge()
         << " charm= " << evtRecDTag->charm()
         << " numOfChildren= " << evtRecDTag->numOfChildren() <<endl;
      cout << "      P(x,y,z,e)= "
         << evtRecDTag->px() << ", "
         << evtRecDTag->py() << ", "
         << evtRecDTag->pz() << ", "
         << evtRecDTag->pe() << endl;

      const vector<Int_t>& trk = evtRecDTag->tracks();
      if( !trk.empty() ) {
         cout << "      vector<tracks>= ";
         copy(trk.begin(),trk.end(),ostream_iterator<int>(cout," "));
         cout << endl;
      }

      const vector<Int_t>& shw = evtRecDTag->showers();
      if( !shw.empty() ) {
         cout << "      vector<showers>= ";
         copy(shw.begin(),shw.end(),ostream_iterator<int>(cout," "));
         cout << endl;
      }

      const vector<Int_t>& otrk = evtRecDTag->otherTracks();
      if( !otrk.empty() ) {
         cout << "      vector<otherTracks>= ";
         copy( otrk.begin(), otrk.end(),
               ostream_iterator<int>(cout, " ") );
         cout << endl;
      }

      const vector<Int_t>& oshw = evtRecDTag->otherShowers();
      if( !oshw.empty() ) {
         cout << "      vector<otherShowers>= ";
         copy( oshw.begin(), oshw.end(),
               ostream_iterator<int>(cout, " ") );
         cout << endl;
      }

      const vector<Int_t>& piid = evtRecDTag->pionId();
      if( !piid.empty() ) {
         cout << "      vector<pionId>= ";
         copy( piid.begin(), piid.end(),
               ostream_iterator<int>(cout, " ") );
         cout << endl;
      }

      const vector<Int_t>& kid = evtRecDTag->kaonId();
      if( !kid.empty() ) {
         cout << "      vector<kaonId>= ";
         copy(kid.begin(),kid.end(),ostream_iterator<int>(cout," "));
         cout << endl;
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMcEvent() const
//--------------------------------------------------------------------
{
   if( !m_TMcEvent ) return;

   const TObjArray* m_mdcMcHitCol   = m_TMcEvent->getMdcMcHitCol();
   const TObjArray* m_emcMcHitCol   = m_TMcEvent->getEmcMcHitCol();
   const TObjArray* m_tofMcHitCol   = m_TMcEvent->getTofMcHitCol();
   const TObjArray* m_mucMcHitCol   = m_TMcEvent->getMucMcHitCol();
   const TObjArray* m_mcParticleCol = m_TMcEvent->getMcParticleCol();

   int NmdcMcHits   = m_mdcMcHitCol->GetEntries();
   int NemcMcHits   = m_emcMcHitCol->GetEntries();
   int NtofMcHits   = m_tofMcHitCol->GetEntries();
   int NmucMcHits   = m_mucMcHitCol->GetEntries();
   int NmcParticles = m_mcParticleCol->GetEntries();

   bool isNULL = (NmdcMcHits==0 && NemcMcHits==0 && NtofMcHits==0 &&
         NmucMcHits==0 && NmcParticles==0);

   if( !isNULL || Verbose() ) {
      cout
         << " ====================== Mc event ====================== "
         << endl;
      cout << " mdcMcHitColSize= "   << NmdcMcHits   << endl;
      cout << " emcMcHitColSize= "   << NemcMcHits   << endl;
      cout << " tofMcHitColSize= "   << NtofMcHits   << endl;
      cout << " mucMcHitColSize= "   << NmucMcHits   << endl;
      cout << " mcParticleColSize= " << NmcParticles << endl;
   }

   if( isNULL ) return; // nothing to print

   if( NmdcMcHits )   PrintMdcMC(m_mdcMcHitCol);
   if( NemcMcHits )   PrintEmcMC(m_emcMcHitCol);
   if( NtofMcHits )   PrintTofMC(m_tofMcHitCol);
   if( NmucMcHits )   PrintMucMC(m_mucMcHitCol);
   if( NmcParticles ) PrintMCParticle(m_mcParticleCol);

   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintMdcMC(const TObjArray* m_mdcMcHitCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_mdcMcHitCol: " << endl;
   TIter mdcMcHitIter(m_mdcMcHitCol);
   while( TMdcMc* mdcMc = (TMdcMc*)mdcMcHitIter.Next() ) {
      cout << " +++++ mdcMcHitId= " << mdcMc->getId()
         << " trackID= " << mdcMc->getTrackIndex() << endl;
      cout << "       (x,y,z)= " << mdcMc->getPositionX()
         << ", " << mdcMc->getPositionY()
         << ", " << mdcMc->getPositionZ()  << endl;
      cout << "       driftDistance= " << mdcMc->getDriftDistance()
         << " depositEnergy= " << mdcMc->getDepositEnergy()
         << " positionFlag= " << mdcMc->getPositionFlag() << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintEmcMC(const TObjArray* m_emcMcHitCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_emcMcHitCol: " << endl;
   TIter emcMcHitIter(m_emcMcHitCol);
   while( TEmcMc* emcMc = (TEmcMc*)emcMcHitIter.Next() ) {
      cout << " +++++ emcMcHitId= " << emcMc->getId()
         << " trackID= " << emcMc->getTrackIndex() << endl;
      cout << "       isHitEmc= " << emcMc->getHitEmc()
         << " PDGCode= " << emcMc->getPDGCode()
         << " PDGCharge= " << emcMc->getPDGCharge()  << endl;
      cout << "       (x,y,z)= " << emcMc->getPositionX()
         << ", " << emcMc->getPositionY()
         << ", " << emcMc->getPositionZ()  << endl;
      cout << "       Px,Py,Pz= " << emcMc->getPx()
         << ", " << emcMc->getPy()
         << ", " << emcMc->getPz()  << endl;
      cout << "       crystal or dead time= " << emcMc->getTime()
         << " depositEnergy= " << emcMc->getDepositEnergy() << endl;

      if( Verbose() ) {
         cout << " EmcMcHitMap: " << endl;
         typedef std::map<Int_t, Double_t> HitMap;
         const HitMap emcMcHitMap = emcMc->getHitMap();
         for( const auto& hm : emcMcHitMap ) {
            cout << hm.first << " -> " << hm.second << endl;
         }
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintTofMC(const TObjArray* m_tofMcHitCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_tofMcHitCol: " << endl;
   TIter tofMcHitIter(m_tofMcHitCol);
   while( TTofMc* tofMc = (TTofMc*)tofMcHitIter.Next() ) {
      cout << " +++++ tofMcHitId= " << tofMc->getId()
         << " trackID= " << tofMc->getTrackIndex() << endl;
      cout << "       (x,y,z)= " << tofMc->getPositionX()
         << ", " << tofMc->getPositionY()
         << ", " << tofMc->getPositionZ()  << endl;
      cout << "       Px,Py,Pz= " << tofMc->getPx()
         << ", " << tofMc->getPy()
         << ", " << tofMc->getPz()  << endl;
      cout << "       TrackLength= " << tofMc->getFlightTime()
         << " FlightTime= " << tofMc->getFlightTime() << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMucMC(const TObjArray* m_mucMcHitCol) const
//--------------------------------------------------------------------
{
   cout << endl << " L________> m_mucMcHitCol: " << endl;
   TIter mucMcHitIter(m_mucMcHitCol);
   while( TMucMc* mucMc = (TMucMc*)mucMcHitIter.Next() ) {
      cout << " +++++ mucMcHitId= " << mucMc->getId()
         << " trackID= " << mucMc->getTrackIndex() << endl;
      cout << "       (x,y,z)= " << mucMc->getPositionX()
         << ", " << mucMc->getPositionY()
         << ", " << mucMc->getPositionZ()  << endl;
      cout << "       Px,Py,Pz= " << mucMc->getPx()
         << ", " << mucMc->getPy()
         << ", " << mucMc->getPz()  << endl;
      // cout << "       DepositEnergy= "
         // << mucMc->getDepositEnergy() << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMCParticle(const TObjArray* m_mcParticleCol)const
//--------------------------------------------------------------------
{
   if( !Verbose() ) { // Short print
      cout
         << " =================== Mc Decay Tree ==================== "
         << endl;
      PrintMcDecayTree(-99,0);
      return;
   }

   // Long print
   cout << endl << " L________> m_mcParticleCol: " << endl;
   TIter mcParticleIter(m_mcParticleCol);
   TMcParticle* mcParticle = nullptr;
   while( (mcParticle = (TMcParticle*)mcParticleIter.Next()) ) {
      int pId = mcParticle->getParticleID();
      cout << " +++++ mcParticle: " ;
      TParticlePDG* pdgParticle =
         TDatabasePDG::Instance()->GetParticle(pId);
      if( pdgParticle ) {
         cout << '\"' << pdgParticle->GetTitle() << '\"';
      } else {
         cout << "??[" <<  pId  << "]";
      }
      cout << "  trkID= " << mcParticle->getTrackIndex()
         << " vtxID0= " << mcParticle->getVertexIndex0()
         << " vtxID1= " << mcParticle->getVertexIndex1();
#if BOSS_VER >= 653
      if( mcParticle->primaryParticle() ) {
         cout << " primary particle";
      } else if( mcParticle->leafParticle() ) {
         cout << " leaf in the tree";
      } else if( mcParticle->decayFromGenerator() ) {
         cout << " decayed by generator";
      } else if( mcParticle->decayInFlight() ) {
         cout << " decayed in flight";
      } else {
         cout << " unknown= " << mcParticle->getStatusFlags();
      }
      cout << endl;
#endif
      cout << "       ini (x,y,z,t)= "
         << mcParticle->getInitialPositionX()
         << ", " << mcParticle->getInitialPositionY()
         << ", " << mcParticle->getInitialPositionZ()
         << ", " << mcParticle->getInitialPositionT()
         << endl;
      cout << "       ini (Px,Py,Pz,E)= "
         << mcParticle->getInitialMomentumX()
         << ", " << mcParticle->getInitialMomentumY()
         << ", " << mcParticle->getInitialMomentumZ()
         << ", " << mcParticle->getInitialMomentumE()
         << endl;
      cout << "       fin (x,y,z,t)= "
         << mcParticle->getFinalPositionX()
         << ", " << mcParticle->getFinalPositionY()
         << ", " << mcParticle->getFinalPositionZ()
         << ", " << mcParticle->getFinalPositionT()
         << endl;
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintMcDecayTree(int root, int shift) const
//--------------------------------------------------------------------
{
   const TObjArray* particles = m_TMcEvent->getMcParticleCol();
   TIter next(particles);
   while( TMcParticle* particle = (TMcParticle *) next() ) {
      if( particle->getMother() == root ) { // daughter
         cout << string(shift+1,' ');

         TParticlePDG * pdgParticle = TDatabasePDG::Instance()
            ->GetParticle(particle->getParticleID());
         if( pdgParticle ) {
            cout << pdgParticle->GetTitle() << " ["
               << particle->getParticleID() << "]";
         } else {
            cout << "?? [" <<  particle->getParticleID() << "]";
         }

         cout << " "
            << particle->getInitialMomentumE() * 1000 << " MeV"
            << " p: ("
            << particle->getInitialMomentumX() * 1000 << ", "
            << particle->getInitialMomentumY() * 1000 << ", "
            << particle->getInitialMomentumZ() * 1000
            << ") MeV/c"
            << " th: "
            << TVector3(particle->getInitialMomentumX(),
                  particle->getInitialMomentumY(),
                  particle->getInitialMomentumZ()).Theta()
            * (180/M_PI) << " deg"
            << endl;

         PrintMcDecayTree(particle->getTrackIndex(), shift+1);
      }
   }
}

//--------------------------------------------------------------------
void DstFormat::PrintTrigEvent() const
//--------------------------------------------------------------------
{
   if( !m_TTrigEvent ) return;

   cout << " ==================== Triger Event ==================== "
      << endl;
   const TTrigData*  m_trigData = m_TTrigEvent->getTrigData();
   // m_trigData->Dump();
   if( m_trigData && m_trigData->IsA() == TTrigData::Class() ) {
      cout << " L________> m_trigData: " << endl;
      cout << " TimeWindow= " <<  m_trigData->getTimeWindow()
         << " TimingType= " <<  m_trigData->getTimingType()
         << " PreScale= " <<  m_trigData->getPreScale() << endl;

      if( Verbose() ) {
         const int* trg_cnd = m_trigData->getTrigCondition();
         if( trg_cnd ) {
            cout << " TrigCondition= ";
            copy( trg_cnd, trg_cnd+48,
                  ostream_iterator<int>(cout," ") );
            cout << endl;
         }
         const int* trg_chan = m_trigData->getTrigChannel();
         if( trg_cnd ) {
            cout << " TrigChannel= ";
            copy( trg_chan, trg_chan+16,
                  ostream_iterator<int>(cout," ") );
            cout << endl;
         }
      }
   }
   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintDigiEvent() const
//--------------------------------------------------------------------
{
   if( !m_TDigiEvent ) return;

   Bool_t  m_fromMc = m_TDigiEvent->getFromMc();
   const TObjArray*  m_mdcDigiCol = m_TDigiEvent->getMdcDigiCol();
   const TObjArray*  m_emcDigiCol = m_TDigiEvent->getEmcDigiCol();
   const TObjArray*  m_tofDigiCol = m_TDigiEvent->getTofDigiCol();
   const TObjArray*  m_mucDigiCol = m_TDigiEvent->getMucDigiCol();

   int NmdcDigi = m_mdcDigiCol->GetEntries();
   int NemcDigi = m_emcDigiCol->GetEntries();
   int NtofDigi = m_tofDigiCol->GetEntries();
   int NmucDigi = m_mucDigiCol->GetEntries();

   if( NmdcDigi==0 && NemcDigi==0 &&
         NtofDigi==0 && NmucDigi==0 ) return;

   cout << " ===================== Digi Event ===================== "
      << endl;
   cout << " m_fromMc= " << m_fromMc << endl;
   cout << " mdcDigiColSize= " << NmdcDigi << endl;
   cout << " emcDigiColSize= " << NemcDigi << endl;
   cout << " tofDigiColSize= " << NtofDigi << endl;
   cout << " mucDigiColSize= " << NmucDigi << endl;
   cout << " ====================================================== "
      << endl;
}

//--------------------------------------------------------------------
void DstFormat::PrintHltEvent() const
//--------------------------------------------------------------------
{
   if( !m_THltEvent ) return;

   const TObjArray*  m_hltRawCol = m_THltEvent->getHltRawCol();
   const THltInf*  m_hltInf = m_THltEvent->getHltInf();
   const TDstHltInf*  m_dstHltInf = m_THltEvent->getDstHltInf();

   int NhltRaw = m_hltRawCol->GetEntries();

   cout << " ===================== Hlt Event ====================== "
      << endl;
   cout << " hltRawColSize= " << NhltRaw << endl;

   if( NhltRaw ) {
      cout << endl << " L________> m_hltRawCol: " << endl;
      TIter hltRawIter(m_hltRawCol);
      while( THltRaw* hltRaw = (THltRaw*)hltRawIter.Next() ) {
         cout << " +++++ hltRawId= " << hltRaw->getIntId()
            << " trkID= " << hltRaw->getTrackIndex()
            << " TimeChannel= " << hltRaw->getTimeChannel()
            << " ChargeChannel= " << hltRaw->getChargeChannel()
            << endl;
      }
   }

   if( m_hltInf && m_hltInf->IsA() == THltInf::Class() ) {
      cout << " L________> m_hltInf: " << endl;
      cout << " EventType= " << m_hltInf->getEventType()
         << " AlgProcess= " << m_hltInf->getAlgProcess()
         << " CriteriaTable= " << m_hltInf->getCriteriaTable()
         << " Version= " << m_hltInf->getVersion() << endl;
      cout << " TotalEnergy= " << m_hltInf->getTotalEnergy()
         << " Number= " << m_hltInf->getNumber()
         << " NCON= " << m_hltInf->getNCON() << endl;

      if( Verbose() ) {
         std::vector<Int_t> mdcData = m_hltInf->getMdcData();
         if( !mdcData.empty() ) {
            cout << " MdcData= ";
            copy(mdcData.begin(), mdcData.end(),
                  ostream_iterator<int>(cout, " "));
            cout << endl;
         }

         std::vector<Int_t> tofData = m_hltInf->getTofData();
         if( !tofData.empty() ) {
            cout << " TofData= ";
            copy(tofData.begin(), tofData.end(),
                  ostream_iterator<int>(cout, " "));
            cout << endl;
         }

         std::vector<Int_t> emcData = m_hltInf->getEmcData();
         if( !emcData.empty() ) {
            cout << " EmcData= ";
            copy(emcData.begin(), emcData.end(),
                  ostream_iterator<int>(cout, " "));
            cout << endl;
         }

         std::vector<Int_t> mucData = m_hltInf->getMucData();
         if( !mucData.empty() ) {
            cout << " MucData= ";
            copy(mucData.begin(), mucData.end(),
                  ostream_iterator<int>(cout, " "));
            cout << endl;
         }

         std::vector<Int_t> conData = m_hltInf->getMdcData();
         if( !conData.empty() ) {
            cout << " ConData= ";
            copy(conData.begin(), conData.end(),
                  ostream_iterator<int>(cout, " "));
            cout << endl;
         }
      }
   }

   if( m_dstHltInf && m_dstHltInf->IsA() == TDstHltInf::Class() ) {
      cout << " L________> m_dstHltInf: " << endl;
      cout << " EventType= " << m_dstHltInf->getEventType()
         << " AlgProcess= " << m_dstHltInf->getAlgProcess()
         << " CriteriaTable= " << m_dstHltInf->getCriteriaTable()
         << " Version= " << m_dstHltInf->getVersion() << endl;
      cout << " TotalEnergy= " << m_dstHltInf->getTotalEnergy()
         << " Number= " << m_dstHltInf->getNumber()
         << " NCON= " << m_dstHltInf->getNCON() << endl;
   }
   cout << " ====================================================== "
      << endl;
}

#if BOSS_VER >= 661
//--------------------------------------------------------------------
void DstFormat::PrintEvtNavigator() const
//--------------------------------------------------------------------
{
   if( !m_TEvtNavigator ) return;
   m_TEvtNavigator->Print();
}
#endif


//--------------------------------------------------------------------
void DstFormat::ClearClasses()
//--------------------------------------------------------------------
{
   if( Verbose() ) {
      cout << " DstFormat::ClearClasses() " << endl;
   }

   // Clear() functions are empty!
   if( m_TEvtHeader ) {
      m_TEvtHeader->Clear();
   }

   if( m_TDstEvent ) {
      m_TDstEvent->Clear();

      auto ptr = m_TDstEvent;
      const_cast<TObjArray*>(ptr->getMdcTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getMdcDedxCol())->Delete();
      const_cast<TObjArray*>(ptr->getMdcKalTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getTofTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getExtTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getEmcTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getMucTrackCol())->Delete();
   }

   if( m_TEvtRecObject ) {
      m_TEvtRecObject->Clear();

      auto ptr = m_TEvtRecObject;
      const_cast<TObjArray*>(ptr->getEvtRecTrackCol())->Delete();
      const_cast<TObjArray*>(ptr->getEvtRecVeeVertexCol())->Delete();
      const_cast<TObjArray*>(ptr->getEvtRecDTagCol())->Delete();
      const_cast<TObjArray*>(ptr->getEvtRecPi0Col())->Delete();
      const_cast<TObjArray*>(ptr->getEvtRecEtaToGGCol())->Delete();
   }

   if( m_TMcEvent ) {
      m_TMcEvent->Clear();

      auto ptr = m_TMcEvent;
      const_cast<TObjArray*>(ptr->getMdcMcHitCol())->Delete();
      const_cast<TObjArray*>(ptr->getEmcMcHitCol())->Delete();
      const_cast<TObjArray*>(ptr->getTofMcHitCol())->Delete();
      const_cast<TObjArray*>(ptr->getMucMcHitCol())->Delete();
      const_cast<TObjArray*>(ptr->getMcParticleCol())->Delete();
   }

   if( m_TTrigEvent ) {
      m_TTrigEvent->Clear();
      m_TTrigEvent->clearTrigData();
   }

   if( m_TDigiEvent )  {
      m_TDigiEvent->Clear();

      const_cast<TObjArray*>(m_TDigiEvent->getMdcDigiCol())->Delete();
      const_cast<TObjArray*>(m_TDigiEvent->getEmcDigiCol())->Delete();
      const_cast<TObjArray*>(m_TDigiEvent->getTofDigiCol())->Delete();
      const_cast<TObjArray*>(m_TDigiEvent->getMucDigiCol())->Delete();
   }

   if( m_THltEvent ) {
      m_THltEvent->Clear();
      m_THltEvent->clearHltInf();
      m_THltEvent->clearDstHltInf();

      const_cast<TObjArray*>(m_THltEvent->getHltRawCol())->Delete();
   }
}
