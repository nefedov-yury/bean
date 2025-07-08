#ifndef _DstFormat_h
#define _DstFormat_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// DstFormat                                                        //
// * This is an abstract class that is the base class for ReadDst   //
// * It describes the data format in the DST file                   //
// * Functions for detailed printing of all fields of DST-events    //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <TSelector.h>

class TChain;
class TFile;
class TBranch;
class TObjArray;

class TEvtHeader;
class TDstEvent;
class TEvtRecObject;
class TMcEvent;
class TTrigEvent;
class TDigiEvent;
class THltEvent;
#if BOSS_VER >= 661
class TEvtNavigator;
#endif

class DstFormat : public TSelector
{
   public:
      DstFormat();
      virtual ~DstFormat();

      // inherited from TSelector ------------------------------------

      // We MUST define API version: version one uses Process(),
      // introduces the SlaveBegin()/SlaveTerminate() functions and
      // clarifies the role of Init() and Notify()
      // Init() call in TTreePlayer requires Version()  == 2
      virtual Int_t  Version() const { return 2;}

      virtual void   Begin(TTree* ) = 0;
      virtual void   SlaveBegin(TTree* ) = 0;
      virtual Bool_t Process(Long64_t entry) = 0;
      virtual void   SlaveTerminate() = 0;
      virtual void   Terminate() = 0;

      virtual void   Init(TTree* tree);
      virtual Bool_t Notify();
      virtual Int_t  GetEntry(Long64_t entry, Int_t getall = 0);
      // -------------------------------------------------------------

      virtual bool Verbose() const = 0;

      void   SetBranches(TTree *tree);
      void   Show(Long64_t entry = -1);
      TFile* GetCurrentFile() const;

      // print-functions
      void  PrintHeader() const;

      void PrintDst() const;
      void PrintMdcTracks(const TObjArray* m_mdcTrackCol) const;
      void PrintMdcDedx(const TObjArray* m_mdcDedxCol)    const;
      void PrintMdcKalTracks(const TObjArray* m_mdcKalTrackCol) const;
      void PrintTofTrack(const TObjArray* m_tofTrackCol)  const;
      void PrintExtTrack(const TObjArray* m_extTrackCol)  const;
      void PrintEmcTrack(const TObjArray* m_emcTrackCol)  const;
      void PrintMucTrack(const TObjArray* m_mucTrackCol)  const;

      void PrintEvtRec() const;
      void PrintRecTrack(const TObjArray* m_evtRecTrackCol) const;
      void PrintRecVeeVertex(
            const TObjArray* m_evtRecVeeVertexCol) const;
      void PrintRecPi0(const TObjArray* m_evtRecPi0Col)     const;
      void PrintRecEtaToGG(const TObjArray* m_evtRecEtaToGGCol) const;
      void PrintRecDTag(const TObjArray* m_evtRecDTagCol)   const;

      void PrintMcEvent() const;
      void PrintMdcMC(const TObjArray* m_mdcMcHitCol)       const;
      void PrintEmcMC(const TObjArray* m_emcMcHitCol)       const;
      void PrintTofMC(const TObjArray* m_tofMcHitCol)       const;
      void PrintMucMC(const TObjArray* m_mucMcHitCol)       const;
      void PrintMCParticle(const TObjArray* m_mcParticleCol)const;
      void PrintMcDecayTree(int root = -99, int shift = 0)  const;

      void PrintTrigEvent() const;
      void PrintDigiEvent() const;
      void PrintHltEvent() const;
#if BOSS_VER >= 661
      void PrintEvtNavigator() const;
#endif

      // access functions:
      const TEvtHeader* GetEvtHeader() const {return m_TEvtHeader;}
      const TDstEvent* GetDstEvent() const {return m_TDstEvent;}
      const TEvtRecObject* GetEvtRecObject() const {
         return m_TEvtRecObject;
      }
      const TMcEvent* GetMcEvent() const {return m_TMcEvent;}
      const TTrigEvent* GetTrigEvent() const {return m_TTrigEvent;}
      const TDigiEvent* GetDigiEvent() const {return m_TDigiEvent;}
      const THltEvent* GetHltEvent() const {return m_THltEvent;}
#if BOSS_VER >= 661
      const TEvtNavigator* GetEventNavigator() const {
         return m_TEvtNavigator;
      }
#endif

   protected:
      TTree*               fChain; //!pointer to the analyzed TChain

      // Declaration of leaf types
      TEvtHeader*          m_TEvtHeader;
      TDstEvent*           m_TDstEvent;
      TEvtRecObject*       m_TEvtRecObject;
      TMcEvent*            m_TMcEvent;
      TTrigEvent*          m_TTrigEvent;
      TDigiEvent*          m_TDigiEvent;
      THltEvent*           m_THltEvent;
#if BOSS_VER >= 661
      TEvtNavigator*       m_TEvtNavigator;
#endif

      // List of branches
      TBranch*             b_TEvtHeader;
      TBranch*             b_TDstEvent;
      TBranch*             b_TEvtRecObject;
      TBranch*             b_TMcEvent;
      TBranch*             b_TTrigEvent;
      TBranch*             b_TDigiEvent;
      TBranch*             b_THltEvent;
#if BOSS_VER >= 661
      TBranch*             b_TEvtNavigator;
#endif

      virtual void ClearClasses();

   // ClassVersionID=0 because we don't need object I/O
   ClassDef(DstFormat,0); // DST format description
   // if class definition use `override` keyword
   //    ClassDefOverride(DstFormat,0); // DST format description
};
#endif
