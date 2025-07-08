#include <iostream>
#include <string>
#include <cstdlib>

#if USE_PROOF != 0
#include <TProof.h>
#include <TMacro.h>
#include <TObjString.h>
#include <TUrl.h>
#endif

#include "RootEventData/TEvtHeader.h"
#include "RootEventData/TDstEvent.h"
#include "RootEventData/TEvtRecObject.h"
#include "RootEventData/TMcEvent.h"
#include "RootEventData/TTrigEvent.h"
#include "RootEventData/TDigiEvent.h"
#include "RootEventData/THltEvent.h"

#include "Bean.h"

// These includes (espesially windows.h) MUST BE after
// "RootEventData/TMcEvent.h" include
// Predefined system-specific compiler macros:
// __unix__  for Unix
// __APPLE__ for MacOS
// _Win32    for Windows (32 or 64 bit)
#if defined (__unix__) || defined (__APPLE__)
#include <dlfcn.h>    // interface to dynamic linking loader
static void* usrlib_handle = 0; //! this is not-persistent member
#elif defined _WIN32
#include <windows.h>
static HINSTANCE usrlib_handle = 0; //! this is not-persistent member
#endif

ClassImp(Bean);

using namespace std;

//--------------------------------------------------------------------
Bean::Bean() : TObject()
//--------------------------------------------------------------------
{
   verbose = false;
   event_dump = false;
   hst_file = "bean_histo.root";
   max_number_events = 0; // no restriction on the number of events

#if USE_PROOF != 0
   proof = nullptr;
   proof_lite = false;
   dst_file_in_sandbox = true;
   dst_file_is_dataset = false;
   proof_xrd_output = false;
#endif
}

//--------------------------------------------------------------------
void Bean::SetDstFile(const char* opt)
//--------------------------------------------------------------------
{
   dst_file = opt;
   // save only dst_file name without path
   auto pos = dst_file.rfind('/'); // search from right to left
   if( pos == std::string::npos ) { // not found
      dst_file_name = dst_file;
   } else {
      dst_file_name = dst_file.substr(pos+1); // [pos+1,end)
   }

#if USE_PROOF != 0
   // set flags: _is_sandbox and _is_dataset

   // In case of PROOF Lite the output file is not in sandbox
   if( proof_lite ) {
      dst_file_in_sandbox = false;
   }

   std::string dataset_name = ParseDatasetName(opt);
   if( !dataset_name.empty() ) {
      dst_file_in_sandbox = false;
      dst_file_is_dataset = true;
   }
#endif
}

//--------------------------------------------------------------------
void Bean::PrintOptions() const
//--------------------------------------------------------------------
{
   if( Verbose() ) {
      cout << string(20,'=') << " Bean::Options "
         << string(20,'=') << endl;
      cout << "> Bean base dir is: " << GetBaseDir() << endl;
      if( MaxNumberEvents() ) {
         cout << "> Process not more than " << MaxNumberEvents()
            << " events in total " << endl;
      }
      if( IsDump() ) {
         cout << "> A detailed printout of DST events will be made"
            << endl;
      }
      cout << "> Histograms will be saved in '" << HstFile()
         << "' file" << endl;
      if( !DstFile().empty() ) {
         cout << "> Selected dst events will be saved in '"
            << DstFile() << "' file" << endl;
         if( dst_file_name != dst_file ) {
            cout << "      dst_file_name= " << DstFileName() << endl;
         }
      }
      cout << "> Bean will call the following user functions:"
         << endl;
      for( const auto& f : Ufn_names ) {
         cout << "     " << f << endl;
      }
      if( NUserFns() == 0 ) {
         cout << "     -- NO FUNCTIONS --" << endl;
      }
#if USE_PROOF != 0
      if( IsProof() ) {
         if( IsProofLite() ) {
            cout << "> Use PROOF-Lite: " << ProofClr() << endl;
         } else {
            cout << "> Use PROOF cluster: " << ProofClr() << endl;
         }
         cout << "> PROOF parameter are: " << ProofParam() << endl;
         if( DstFileInSandbox() ) {
            cout << "> DST file will saved in sandbox" << endl;
         }
         if( DstFileIsDataset() ) {
            cout << "> DST file will saved as dataset" << endl;
         }
         if( ProofXrdOutput() ) {
            cout << "> Xrootd will be used to get DST file "
               "from master node" << endl;
         }
      }
#endif
      cout << string(55,'=') << endl;
   }
}

//--------------------------------------------------------------------
void Bean::AddUserFcn(const char* name)
//--------------------------------------------------------------------
{
   Ufn_names.push_back(string(name));
}

//--------------------------------------------------------------------
void Bean::LoadUserLib()
//--------------------------------------------------------------------
{
   if( verbose ) {
      cout << "Bean::LoadUserLib()" << endl;
   }

#if defined (__unix__) || defined (__APPLE__)
   // string usrlib_name = "libUser.so";
   // if( IsProof() ) usrlib_name = "BeanUser/" + usrlib_name;
   // usrlib_handle =
   //           dlopen(usrlib_name.c_str(), RTLD_NOW | RTLD_LOCAL);
   // RTLD_NOW - all undefined symbols in the library are resolved
   //            before dlopen()  returns;
   // RTLD_LOCAL - symbols defined in this library are not made
   //              available to resolve references in subsequently
   //              loaded libraries.

   usrlib_handle = dlopen(0, RTLD_NOW | RTLD_GLOBAL);
#elif defined _WIN32
   // macro BEANLIBDIR defined in BeanCore/CMakeLists.txt
   if( !SetDllDirectory(BEANLIBDIR) ) {
      cout << "Bean::LoadUserLib() ==> cannot "
         "SetDllDirectory(" BEANLIBDIR "): " << GetLastError()<<endl;
      exit(EXIT_FAILURE);
   }
   usrlib_handle = LoadLibrary("BeanUser");
#endif

   if( !usrlib_handle ) {
      cout << "Bean::LoadUserLib() ==> cannot ";
#if defined (__unix__) || defined (__APPLE__)
      // cout << "dlopen(" << usrlib_name << ")"
      cout << "dlopen(0)" << endl;
      cout << dlerror() << endl;
#elif defined _WIN32
      cout << R"(LoadLibrary("BeanUser.dll"): )"
         << GetLastError() << endl;
#endif  // _WIN32
      exit(EXIT_FAILURE);
   }

   if( verbose ) {
      cout << "User library is loaded" << endl;
   }
}

//--------------------------------------------------------------------
void Bean::LoadUserFcn(const std::string& name)
//--------------------------------------------------------------------
{
   if( verbose ) {
      cout << "Bean::LoadUserFcn(" << name << ")" << endl;
   }

   if( !usrlib_handle ) {
      LoadUserLib();
   }

   string nameJobStart = name + "StartJob";
   string nameEvent = name + "Event";
   string nameJobEnd = name + "EndJob";
   string vout = "User functions: ";

#if defined (__unix__) || defined (__APPLE__)
   void* sym = dlsym(usrlib_handle, nameEvent.c_str());
#elif defined _WIN32
   // FARPROC is a generic function pointer
   FARPROC sym = GetProcAddress(usrlib_handle, nameEvent.c_str());
#endif
   if( sym ) {
      Ufn_event.push_back(sym);
      vout += nameEvent + " ";
   } else {
      cout << "ERROR: User function " << nameEvent
         << " does not exist!" << endl;
#if defined _WIN32
      if( verbose ) {
         cout << "check name with: 'dumpbin -exports BeanUser.dll'"
            << endl;
      }
#endif  // _WIN32
      cout << " This is mandatory function. Stop!" << endl;
      exit(EXIT_FAILURE);
   }

#if defined (__unix__) || defined (__APPLE__)
   sym = dlsym(usrlib_handle, nameJobStart.c_str());
#elif defined _WIN32
   sym = GetProcAddress(usrlib_handle, nameJobStart.c_str());
#endif
   if( sym ) {
      Ufn_start.push_back(sym);
      vout += nameJobStart + " ";
   }

#if defined (__unix__) || defined (__APPLE__)
   sym = dlsym(usrlib_handle, nameJobEnd.c_str());
#elif defined _WIN32
   sym = GetProcAddress(usrlib_handle, nameJobEnd.c_str());
#endif
   if( sym ) {
      Ufn_end.push_back(sym);
      vout += nameJobEnd + " ";
   }

   if( verbose ) {
      cout << vout << "are loaded" << endl;
   }
}

//--------------------------------------------------------------------
void Bean::LoadUserFcns()
//--------------------------------------------------------------------
{
   // load all user functions
   if( verbose ) {
      cout << "Bean::LoadUserFcns()" << endl;
   }

   Ufn_event.clear();
   Ufn_start.clear();
   Ufn_end.clear();

   for( const auto& ufn : Ufn_names ) {
      LoadUserFcn( ufn );
   }
}

//--------------------------------------------------------------------
std::string Bean::ParseDatasetName(std::string filename)
//--------------------------------------------------------------------
{
   if ( filename.substr(0,5) == "ds://" ) { // pos,len
      return filename.substr(5);            // pos,end()
   }
   return std::string();
}

#if ( USE_PROOF != 0 )
//--------------------------------------------------------------------
void Bean::SetProof(TProof* val)
//--------------------------------------------------------------------
{
   proof = val;
   if( proof ) {
      proof_lite = proof->IsLite();

      // In case of PROOF Lite the output file is not in sandbox
      if( proof_lite ) {
         dst_file_in_sandbox = false;
      }
   }
}

//--------------------------------------------------------------------
std::string Bean::GetProofMasterWorkdir()
//--------------------------------------------------------------------
{
   if( proof ) {
      RedirectHandle_t rh;
      TString tempFileName;
      gSystem->TempFileName(tempFileName);
      gSystem->RedirectOutput(tempFileName, "a", &rh);
      proof->Print();
      gSystem->RedirectOutput(0, 0, &rh);

      TMacro macro;
      macro.ReadFile(tempFileName);
      gSystem->Unlink(tempFileName);

      TObject* line = macro.GetLineWith("Working directory:");
      if( line ) {
         TString s = line->GetName();
         TString ss(s(s.First(":")+1, s.Length()));
         TString sss(ss.Strip(TString::kBoth));
         return sss.Data();
      }
   }
   return std::string();
}
#endif
