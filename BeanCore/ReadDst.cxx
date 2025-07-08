#include <iostream>
#include <cstdlib>
#include <csignal>
#include <ctime>

// I define 'bean_termination' here, as part of BeanCore library
// and use the external declaration in the main.cxx
volatile sig_atomic_t bean_termination = 0;  // signal handler

#include <RVersion.h> // Root version macros
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TEntryList.h>

#if USE_PROOF != 0
#include <TProofOutputFile.h>
#include <TProof.h>
#include <TMacro.h>
#include <TFileCollection.h>
// #include <TObjectTable.h>
#endif

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"

#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"
#include "RootEventData/TJobInfo.h"

#include "DstEvtRecTracks.h"
#include "TMergeableMap.h"
#include "ReadDst.h"

using namespace std;

//--------------------------------------------------------------------
ReadDst::ReadDst(): DstFormat()
//--------------------------------------------------------------------
{
   bean = nullptr;

   current_entry = -1;
   selected_entries = nullptr;

   m_evtRecTrkCol = new TObjArray();

   T_select = nullptr; // to save selected events
   f_select = nullptr;
   n_select_events = 0;

   n_events = 0;
   start_time = time(0);

   m_hst_dir = new TMergeableMap();
   m_hst_dir->SetName("DirectoryMap");

#if USE_PROOF != 0
   fp_select = 0;
#endif
}

//--------------------------------------------------------------------
ReadDst::~ReadDst()
//--------------------------------------------------------------------
{
   delete m_evtRecTrkCol;
}

//--------------------------------------------------------------------
bool ReadDst::LoadConfig(Bean* bn)
//--------------------------------------------------------------------
{
   if( bean ) {
      cout << "ReadDst::LoadConfig Warning: "
         "config is already loaded" << endl;
   } else if( bn ) {
      bean = bn;
   } else {
      // fInput is in Tselector or TProofPlayer
      bean = (Bean *) fInput->FindObject("Bean");
      if( !bean ) {
         cout << "ReadDst::LoadConfig ERROR: "
            "No Bean in fInput" << endl;
      }
   }
   return (bean != nullptr);
}

//--------------------------------------------------------------------
bool ReadDst::Verbose() const
//--------------------------------------------------------------------
{
   if( !bean ) {
      cout << "ReadDst::Verbose Warning: "
         "config is not yet loaded" << endl;
      return true;
   }
   return bean->Verbose();
}

//--------------------------------------------------------------------
std::string ReadDst::GetBaseDir() const
//--------------------------------------------------------------------
{
   if( !bean ) {
      cout << "ReadDst::GetBaseDir ERROR: "
         "config is not yet loaded" << endl;
      exit(EXIT_FAILURE);
   }
   return bean->GetBaseDir();
}

//--------------------------------------------------------------------
std::string ReadDst::AbsPath(const std::string& rel_path) const
//--------------------------------------------------------------------
{
   return (this->GetBaseDir() + "/" + rel_path);
}

//--------------------------------------------------------------------
void ReadDst::SetEntryList(TEntryList* el)
//--------------------------------------------------------------------
{
   if( selected_entries ) {
      cout << " ReadDst::SetEntryList WARNING: "
         "you will rewrite list of selected entries" << endl;
   }
   selected_entries = el;

   if( fChain ) {
      if( fChain->GetEntryList() ) {
         cout << " ReadDst::SetEntryList WARNING: "
            "list of selected entries of the fChain "
            "will be overridden" << endl;
      }
      fChain->SetEntryList(selected_entries);
   }
}

//--------------------------------------------------------------------
Long64_t ReadDst::GetEntryNumber()
//--------------------------------------------------------------------
{
   return fChain->GetChainEntryNumber(current_entry);
}

//--------------------------------------------------------------------
void ReadDst::SaveEntryInList(TEntryList* el)
//--------------------------------------------------------------------
{
   el->Enter(this->GetEntryNumber(),fChain);
}

//--------------------------------------------------------------------
void ReadDst::Begin(TTree* )
//--------------------------------------------------------------------
{
   // This method is called before looping on the events in the Tree.
   // The user can create his histograms in this function.
   // When using PROOF Begin() is called on the client only. Histogram
   // creation should preferable be done in SlaveBegin() in that case.

   if( !bean ) {
      LoadConfig();
   }
   if( Verbose() ) {
      cout << "ReadDst::Begin()" << endl;
   }
}

//--------------------------------------------------------------------
void ReadDst::SlaveBegin(TTree* )
//--------------------------------------------------------------------
{
   // The SlaveBegin() function is called after the Begin() function.
   // This method is called on each PROOF worker node.
   // The user can create his histograms in this method.
   // In local mode this method is called on the client too.

   // The Init() function is called after the SlaveBegin() function
   // therefore the fChain is not accessible here on PROOF-workers

   if( !bean && !LoadConfig() ) {
      exit(EXIT_FAILURE);
   }
   if( Verbose() ) {
      cout << "ReadDst::SlaveBegin()" << endl;
      // if(fChain ) fChain->Print();
      // cout << " ... " << endl;
   }

   bean->LoadUserFcns();

   string new_dst_file = bean->DstFile();
   if( new_dst_file.size() > 0 ) {
      if( !bean->IsProof() ) {
         f_select = TFile::Open(new_dst_file.c_str(),"RECREATE");
         if( !f_select ) {
            cout  << " ReadDst::SlaveBegin:: ERROR "
               "can not open file " << new_dst_file << endl;
            exit(EXIT_FAILURE);
         }
      }
#if USE_PROOF != 0
      if( bean->IsProof() ) {
         if( bean->DstFileIsDataset() ) {
            fp_select = new TProofOutputFile(
                  bean->DstFileName().c_str(),
                  TProofOutputFile::kDataset,
                  TProofOutputFile::kRegister |
                  TProofOutputFile::kVerify |
                  TProofOutputFile::kOverwrite );

            fp_select->GetFileCollection()
               ->SetDefaultTreeName("/Event");
         } else {
            // here we are using old constructor
            fp_select = new TProofOutputFile(
                  bean->DstFileName().c_str(), "M" );

            if( !bean->DstFileInSandbox() ) {
               fp_select->SetOutputFileName(bean->DstFile().c_str());
            }
         }

         f_select = fp_select->OpenFile("RECREATE");

         if( !f_select ) {
            cout << " ReadDst::SlaveBegin ERROR: "
               "can not open " << fp_select->GetDir() << "/"
               << fp_select->GetFileName() << endl;
            exit(EXIT_FAILURE);
         }
      }
#endif
   }

   UserStartJob();
}

//--------------------------------------------------------------------
void ReadDst::Init(TTree* tree)
//--------------------------------------------------------------------
{
   // The Init() function is called when the selector needs to
   // initialize a new tree or chain. Typically here the branch
   // addresses and branch pointers of the tree will be set.

   DstFormat::Init(tree);

   if( selected_entries ) {
      // file names (in fChain and in "selected_entries")
      // are expanded  to full path names
      fChain->SetEntryList(selected_entries);
   }
}

//--------------------------------------------------------------------
Bool_t ReadDst::Notify()
//--------------------------------------------------------------------
{
   // The Notify() function is called when a new file is opened.
   // This can be either for a new TTree in a TChain or when a new
   // TTree is started when using PROOF.
   // The return value is currently not used.

   // if( fChain ) {
   // fChain->SetParallelUnzip(kFALSE);
   // fChain->SetCacheSize(-1);
   // }

   DstFormat::Notify();

   if( Verbose() ) {
      cout << "ReadDst::Notify()" << endl;
   }

   if( f_select && T_select == nullptr ) {
      // Note that only active branches are copied !
      // MUST be called DstFormat::Init() before this line
      T_select = fChain->CloneTree(0);
      T_select->SetDirectory(0);
      T_select->SetDirectory(f_select);
      // T_select->SetAutoSave(); // use default !
      // T_select->SetBasketSize("*",1024);
      T_select->SetAutoFlush(100);
      T_select->AutoSave();
   }

#if USE_PROOF != 0
   // The cloned tree stays connected with this tree until this tree
   // is deleted. The ProofPlayer::Proccess(dataset,selector) creates
   // new tree for each element of the dataset. Therefore we MUST
   // reassign all branch addresses of T_select every new file
   if( bean->IsProof() && T_select ) {
      SetBranches(T_select);
   }
#endif

   return kTRUE;
}

//--------------------------------------------------------------------
Bool_t ReadDst::Process(Long64_t entry)
//--------------------------------------------------------------------
{
   // This method is called to process an event. It is the user's
   // responsibility to read the corresponding entry in memory
   // (may be just a partial read).
   // Once the entry is in memory one can apply a selection and if
   // the event is selected histograms can be filled.
   //
   // The processing can be stopped by calling Abort()
   // in this case (root-6.26) termination functions are not called

   current_entry = entry;

   ClearClasses(); // must be before "GetEntry"

   Int_t nbytes=0;

   // WARNING: entry is always the local entry number
   //          in the current tree!
   nbytes += fChain->GetTree()->GetEntry(entry);

   if( Verbose() ) {
      cout << "ReadDst::Process() Entry# " << entry
         << " nbytes= " << nbytes << endl;
      if( entry == 0 ) {
         Show(entry);
      }
   }

   if( bean->IsDump() ) {
      PrintHeader();
      PrintDst();
      PrintEvtRec();
      PrintMcEvent();
      PrintTrigEvent();
      PrintDigiEvent();
      PrintHltEvent();
#if BOSS_VER >= 661
      PrintEvtNavigator();
#endif
   }

   CreateEvtRecTrkCol();

   if( UserEvent(T_select) ) {
      T_select->Fill();
      n_select_events++;
      if( Verbose() ) {
         cout << " + n_select_events= " << n_select_events << endl;
      }
   }

   n_events++;
   if( bean_termination != 0 ) {
      printf("User signal '%s' had been received on event# %zu\n",
#if defined (__unix__) || defined (__APPLE__)
            strsignal(bean_termination),
#elif defined _WIN32
            (string("Interrupt# ") +
             to_string(bean_termination)).c_str(),
#endif
            n_events);
      this->Abort("Stop the event loop...");
   }

   return kTRUE;
}

//--------------------------------------------------------------------
void ReadDst::SlaveTerminate()
//--------------------------------------------------------------------
{
   // This method is called at the end of the loop on all PROOF worker
   // nodes. In local mode this method is called on the client too.

   if( Verbose() )  {
      cout << "ReadDst::SlaveTerminate()" << endl;
      // gObjectTable->Print();
   }

   UserEndJob();

#if USE_PROOF != 0
   if( fp_select ) { // Proof only
      Bool_t cleanup = kFALSE;
      TDirectory *savedir = gDirectory;
      if( Verbose() ) {
         cout << "ReadDst::SlaveTerminate n_select_events= "
            << n_select_events << endl;
      }
      if( n_select_events > 0 ) {
         f_select->cd();
         WriteJobInfo();
         T_select->Write();
         if( Verbose() ) {
            cout << "events saved in file: " << f_select->GetName()
               << endl;
         }
         fOutput->Add(fp_select);
      } else {
         cleanup = kTRUE;
      }
      // #warning  T_select->SetDirectory(0); leads to crash.

      gDirectory = savedir;
      f_select->Close();

      // Cleanup, if needed
      if( cleanup ) {
         TUrl uf( *(f_select->GetEndpointUrl()) );
         SafeDelete(f_select);
         gSystem->Unlink(uf.GetFile());
         SafeDelete(fp_select);
         if( Verbose() ) {
            cout << "cleanup file: " << uf.GetFile() << endl;
         }
      }
   }
#endif

   fOutput->Add(m_hst_dir);
}

//--------------------------------------------------------------------
void ReadDst::Terminate()
//--------------------------------------------------------------------
{
   // This method is called at the end of the loop on all events.
   // When using PROOF Terminate() is call on the client only.
   // Typically one performs the fits on the produced histograms
   // or write the histograms to file in this method.

   if( Verbose() ) {
      cout << "ReadDst::Terminate()" << endl;
      if( fChain ) {
         fChain->PrintCacheStats("a");
      }
   }

   if( !bean->IsProof() ) {
      if( T_select ) {
         f_select->cd();
         T_select->Write();
         cout << "ReadDst::Terminate " << n_select_events
            << " events saved in file: " << f_select->GetName()
            << endl;

         WriteJobInfo();
         f_select->Close();
         SafeDelete(f_select);
      }

      size_t time_taken = (long)(time(0) - start_time);
      double vel = double(n_events) / double(time_taken);
      printf("Performance summary: "
            "%zu events processed in %zus: ~%.1f events/sec\n",
            n_events, time_taken, vel);
   }

   Save_histo();

#if USE_PROOF != 0
   // retrieving output DST file from master
   if( bean->IsProof()
         && !bean->Proof()->IsLite()
         && bean->DstFileName().size() > 0
         && bean->DstFileInSandbox() ) {
      string workdir = bean->GetProofMasterWorkdir();

      // GetFile() return 0 on success, -1 on error
      int getfile_rc = -1;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,2)
      if( !bean->ProofXrdOutput() ) {
         // download file using GetFile (default)
         string relative_path = workdir + "/" + bean->DstFileName();
         getfile_rc = bean->Proof()->GetManager()
            ->GetFile(relative_path.c_str(),bean->DstFile().c_str());
      }
#endif
      if( (bean->ProofXrdOutput()) || (getfile_rc == -1) ) {
         // download file with Xrootd: construct URL using
         // master host, master workdir and Dst file name
         string serverDstLocation = "root://" + bean->ProofClr()
            + "/" + workdir + "/" + bean->DstFileName();
         if( Verbose() ) {
            cout << "serverDstLocation:" << serverDstLocation << endl;
         }
         TFile::Cp(serverDstLocation.c_str(),bean->DstFile().c_str());
      }
   }
#endif
}

//--------------------------------------------------------------------
void ReadDst::UserStartJob()
//--------------------------------------------------------------------
{
   // run all 'Start-Job' functions:
   const VecUF& Ufn_start = bean->GetStartJobFns();
   for( const auto& uf : Ufn_start ) {
      ssfn user_func = (ssfn) (uf);
      user_func(this);
   }
}

//--------------------------------------------------------------------
bool ReadDst::UserEvent(TTree* T)
//--------------------------------------------------------------------
{
   // run all 'Evens' functions:
   const VecUF& Ufn_event = bean->GetUserEventFns();
   bool save_this = false;
   for( const auto& uf : Ufn_event ) {
      pufn user_func = (pufn) (uf);
      bool ret = user_func(this,
            m_TEvtHeader,m_TDstEvent,m_TEvtRecObject,
            m_TMcEvent,m_TTrigEvent,m_TDigiEvent,m_THltEvent);
      if( T && ret ) {
        save_this = true;
      }
   }

   if( T && Ufn_event.empty() ) { // copy dst
      save_this = true;
   }
   return save_this;
}

//--------------------------------------------------------------------
void ReadDst::UserEndJob()
//--------------------------------------------------------------------
{
   // run all 'End-Job' functions:
   const VecUF& Ufn_end = bean->GetEndJobFns();
   for( const auto& uf : Ufn_end ) {
      ssfn user_func = (ssfn) (uf);
      user_func(this);
   }
}

//--------------------------------------------------------------------
void ReadDst::CreateEvtRecTrkCol()
//--------------------------------------------------------------------
{
   if( Verbose() ) {
      cout << " CreateEvtRecTrkCol "
         << m_evtRecTrkCol->GetEntriesFast() << endl;
   }
   m_evtRecTrkCol->Delete();
   if( m_TEvtRecObject ) {
      const TObjArray* m_evtRecTrackCol =
         m_TEvtRecObject->getEvtRecTrackCol();
      TIter evtRecTrackIter(m_evtRecTrackCol);
      TEvtRecTrack* evtRecTrack = nullptr;
      while( (evtRecTrack=(TEvtRecTrack*)evtRecTrackIter.Next()) ) {
         m_evtRecTrkCol->AddLast(
               new DstEvtRecTracks(evtRecTrack,m_TDstEvent)
               );
      }
   }
}

//--------------------------------------------------------------------
void ReadDst::RegInDir(const VecObj& vhst, const char* dir)
//--------------------------------------------------------------------
{
   if( Verbose() ) {
      cout << " ReadDst::RegInDir(VecObj) " << endl;
   }
   // add the prefix "dir_" to the name so that it is unique in ROOT
   for( const auto hst : vhst ) {
      if( hst ) {
         if( dir ) {
            string nn = string(dir) + "_" + hst->GetName();
            hst->SetName(nn.c_str());
            m_hst_dir->Add(
                  new TObjString(hst->GetName()),
                  new TObjString(dir)
                  );
         }
         CheckDupName(hst);
         fOutput->Add(hst);
         if( Verbose() ) {
            cout << " + " << hst->GetName() << endl;
         }
      }
   }
}

//--------------------------------------------------------------------
void ReadDst::RegInDir(const VecObj* vhst, const char* dir)
//--------------------------------------------------------------------
{
   RegInDir(*vhst,dir);
}

//--------------------------------------------------------------------
void ReadDst::RegInDir(const TList* lhst, const char* dir )
//--------------------------------------------------------------------
{
   if( Verbose() ) {
      cout << " ReadDst::RegInDir(TList) " << endl;
   }
   VecObj vhst;
   TIter next(lhst);
   while( TNamed* obj = (TNamed*)next() ) {
      vhst.push_back(obj);
   }
   RegInDir(vhst, dir);
}

//--------------------------------------------------------------------
void ReadDst::RegInDir(const TList& lhst, const char* dir )
//--------------------------------------------------------------------
{
   RegInDir(&lhst,dir);
}

//--------------------------------------------------------------------
void ReadDst::Save_histo() const
//--------------------------------------------------------------------
{
   if( fOutput->GetSize() == 0 ) {
      return;
   }

   string hst_file = bean->HstFile();
   TFile* c_out = TFile::Open(hst_file.c_str(),"RECREATE");
   if( !c_out ) {
      cout << " ReadDst::Save_histo ERROR: "
         "can not open " << hst_file << endl;
      return;
   }
   c_out->cd();

   int nhisto = 0;
   int ntrees = 0;

   TMergeableMap* directoryMap =
      (TMergeableMap*) fOutput->FindObject("DirectoryMap");
   if( Verbose() ) {
      cout << " Save_histo directoryMap:" << directoryMap << endl;
   }

   TIter next(fOutput);
   while( TObject* obj = next() ) {
      if( obj == directoryMap ) {
         ; // this object is for internal use
#if USE_PROOF != 0
      } else if( obj->IsA() == TProofOutputFile::Class() ) {
         ; // this object is for internal use
#endif
      } else if( obj->InheritsFrom(TObject::Class()) &&
            obj->IsA()->GetMethodWithPrototype("GetName","") ) {
         TNamed* hst  = (TNamed*) obj;
         const char* hname = hst->GetName();
         if( obj->InheritsFrom("TH1") ) { // this is histograms
            if( ((TH1*) obj) -> GetEntries() < 0.01 ) {
               if( Verbose() ) {
                  cout << " skip empty histogram: " << hname << endl;
               }
               continue;
            }
            nhisto++;
         } else if( obj->InheritsFrom("TTree") ) { // this is tree
            if( ((TTree*) obj) -> GetEntries() < 0.01 ) {
               if( Verbose() ) {
                  cout << " skip empty tree: " << hname << endl;
               }
               continue;
            }
            ntrees++;
         }

         // save objects in the directory if it was specified
         // in this case remove the prefix "dir_" from the names
         // see RegInDir(VecObj*)
         if( directoryMap ) {
            TObjString* dir =
               (TObjString*)directoryMap->GetValue(hname);
            if( dir ) {
               const char* dir_name = dir->GetName();
               string prefix(dir_name);
               string nn = string(hname).substr(prefix.size()+1);
               hst->SetName(nn.c_str());

               if( !c_out->GetDirectory( dir_name ) ) {
                  c_out->mkdir( dir_name );
               }
               c_out->cd( dir_name );
            } else {
               c_out->cd();
            }
         } else {
            c_out->cd();
         }

         // to store TMap/TCollection/etc as single object
         obj->Write(NULL, TObject::kSingleKey);
      } else {
         cout << "ReadDst::Save_histo Warning: "
            << "an object in fOutput doesn't have a name, typeid= "
            << typeid(*obj).name() << endl;
         TDirectory* dir = gDirectory;
         c_out->cd();
         obj->Write();
         dir->cd();
      }
   }

   c_out->Close();
   if( Verbose() ) {
      cout << "Save " << nhisto << " histograms " << endl
         << "  and " << ntrees << " trees in file " << hst_file
         << endl;
   }
}

//--------------------------------------------------------------------
void ReadDst::CheckDupName(TObject *obj)
//--------------------------------------------------------------------
{
   if( obj ) {
      TObject *org = fOutput->FindObject(obj->GetName());
      if( org && (org != obj) ) {
         cout << " ReadDst::CheckDupName:: FATAL ERROR " << endl
            << " object with name: " << obj->GetName()
            << " already in the fOutput list " << endl;
         exit(EXIT_FAILURE);
      }
   }
}

//--------------------------------------------------------------------
void ReadDst::WriteJobInfo()
//--------------------------------------------------------------------
{
   // copy&pasted from BOSS
   TJobInfo* jobInfo = new TJobInfo;
   TTree* m_jobInfoTree = new TTree("JobInfoTree","Job info");
   m_jobInfoTree->Branch("JobInfo",&jobInfo);

   string m_bossVer = "BEAN";
   jobInfo->setBossVer(m_bossVer);

   // fill all members of TJobInfo with something
   // jobInfo->setDecayOptions(string(10,'+'));
   // jobInfo->addJobOptions(string(10,'+'));
   // vector<int> m_totEvtNo(2,1);
   // jobInfo->setTotEvtNo(m_totEvtNo);

   m_jobInfoTree->Fill();
   f_select->cd();
   m_jobInfoTree->Write();
}

