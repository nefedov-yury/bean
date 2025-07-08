#ifndef _ReadDst_h
#define _ReadDst_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// ReadDst                                                          //
//  * primary class to read DST files                               //
//  * overridden methods of the TSelect class are used              //
//  * there is basic EntryList functionality                        //
//  * user functions stored in the BEAN class are called            //
//  * histograms and selected DST events are saved                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "Bean.h"
#include "DstFormat.h"

// forward declaration of classes
class ReadDst;
class TNamed;
class TEntryList;
class TMergeableMap;

#if USE_PROOF != 0
class TProofOutputFile;
#endif

// pointers on user functions:
// for start/end job
typedef void (*ssfn) (ReadDst*);
// per each event
typedef bool (*pufn) (ReadDst*,
      TEvtHeader*,TDstEvent*,TEvtRecObject*,
      TMcEvent*,TTrigEvent*,TDigiEvent*,THltEvent*);
// function names
typedef std::vector<TNamed* > VecObj;

class ReadDst : public DstFormat
{
   public:
      ReadDst();
      virtual     ~ReadDst();

      bool        LoadConfig(Bean* bn=nullptr);
      bool        Verbose() const;
      void        SetVerbose() {bean->SetVerbose();}
      void        SetSilent()  {bean->SetSilent();}
      std::string GetBaseDir() const;
      std::string AbsPath(const std::string& rel_path) const;

      const TObjArray* GetEvtRecTrkCol() const{return m_evtRecTrkCol;}

      // functions to work with EntryList
      void     SetEntryList(TEntryList* el);
      Long64_t GetEntryNumber();
      void     SaveEntryInList(TEntryList* el);

      // TSelector functions
      void   Begin(TTree* );
      void   SlaveBegin(TTree* );
      void   Init(TTree* tree);
      Bool_t Notify();
      Bool_t Process(Long64_t entry);
      void   SlaveTerminate();
      void   Terminate();

      // call user functions
      void UserStartJob();
      bool UserEvent(TTree* T);
      void UserEndJob();
      void CreateEvtRecTrkCol();

      // histogramming
      void RegInDir(const VecObj& vhst, const char* dir=nullptr);
      void RegInDir(const VecObj* vhst, const char* dir=nullptr);
      void RegInDir(const TList* lhst, const char* dir=nullptr);
      void RegInDir(const TList& lhst, const char* dir=nullptr);
      void Save_histo() const;

   private:
      Bean*       bean; //-> all configuration parameters are here

      Long64_t    current_entry;
      TEntryList* selected_entries;

      // track linking
      TObjArray*  m_evtRecTrkCol;

      // for selected events
      TTree*      T_select;
      TFile*      f_select;
      size_t      n_select_events;

      // to save histograms to the directory specified in RegInDir
      // this is map(TObjString*,TObjString*) hst_name:dir_name
      TMergeableMap*       m_hst_dir;

      // time measuring
      time_t      start_time;
      size_t      n_events;

      // internal functions
      void CheckDupName(TObject* obj);
      void WriteJobInfo();

#if USE_PROOF != 0
      TProofOutputFile* fp_select;
#endif

   // ClassVersionID=0 because we don't need object I/O
   ClassDef(ReadDst,0); // Primary class to read DST
   // if class definition use `override` keyword
   // ClassDefOverride(ReadDst,0); // Primary class to read DST
};
#endif
