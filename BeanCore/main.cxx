#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <regex>

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <csignal>

// Predefined system-specific compiler macros:
// __unix__  for Unix
// __APPLE__ for MacOS
// _Win32    for Windows (32 or 64 bit)
#if defined (__unix__) || defined (__APPLE__)
#include <getopt.h>
#elif defined _WIN32
// do not define min() and max() macros in windows.h
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include "win_getopt/getopt.h"
#endif

#include <TSystem.h>
#include <TEnv.h>
#include <TChain.h>

// #include <TObjectTable.h> // to debug memory leak
// ## Put following in .rootrc file
// ## memstat, memcheck have been removed in root-6.26
// ## Root.MemStat: 1
// ## Root.MemStat.size: -1
// ## Root.MemStat.cnt: -1
// ## Activate TObject statistics
// Root.ObjectStat: 1

#if USE_PROOF != 0
#include <map>

#include <TROOT.h>
#include <TProof.h>
#include <TProofLog.h>

#include <TObjString.h>
#include <TObjArray.h>
#include <TMap.h>
#include <TFileCollection.h>

#include "DatabaseSvc/DatabaseSvc.h"
#include "unix_cp.h"
#endif

#include "Bean.h"
#include "ReadDst.h"

// Function prototypes
void Usage(int argc, char **argv, bool verbose=false);
static bool str2long(const char* str, long& val);
static std::pair<std::string,std::string> str2strs(const char* opt);

static void set_user_termination(bool verbose);
static void termination_handler(int isig);
#if defined _WIN32
BOOL CtrlHandler( DWORD fdwCtrlType );
void win_exit(){::ExitProcess(0);}
#endif

#if USE_PROOF != 0
void graceful_proof_exit();
static void set_proof_termination();
static void proof_segfault_handler(int isig);

static void prooflite_ld_library_path();
static void check_upload_enable(TProof* proof, const char* pname);
static std::string DatasetStr(
      const std::vector<std::string>& dataset_names);
#endif

// Global variables
static Bean* bean = nullptr;
#define TREE_CACHE_SIZE 10000000  //10 MB

// defined in BeanCore library (see ReadDst.cxx)
extern volatile sig_atomic_t bean_termination;

#if USE_PROOF != 0
static bool save_proof_logs;
#endif

using namespace std;

//--------------------------------------------------------------------
int main(int argc, char **argv)
//--------------------------------------------------------------------
{
   // Switch on synchronization with the standard C streams.
   ios_base::sync_with_stdio(true);

#if defined _WIN32
   // Enable two-digit exponent format: obsolete for VS>=15
   // _set_output_format(_TWO_DIGIT_EXPONENT);

   // suppress the abort message
   _set_abort_behavior( 0, _WRITE_ABORT_MSG );

   // Disable the message box for errors, unrecoverable problems, etc
   _CrtSetReportMode( _CRT_ERROR, 0 );

   // Disable assertion failures
   _CrtSetReportMode( _CRT_ASSERT, 0 );

   // function win_exit() will run at exit()
   // auto win_exit = []() -> void{::ExitProcess(0);};
   atexit( win_exit );
#endif

   bean = new Bean;
   // set the base directory; BEANBASE is defined by cmake
#ifdef BEANBASE
   bean->SetBaseDir(BEANBASE);
#else
#error "macro BEANBASE not defined"
#endif

   // ===================== PARSE COMMAND OPTIONS ====================
   bool is_error = false;
   bool verbose = bean->Verbose(); // use default verbose of the bean

   string list_fnames;
   string prefix;

#if USE_PROOF != 0
   string sqlite_db[2];
   multimap<string, string> proof_params;
   TString workers_disable_str;
   save_proof_logs = false;
   int proof_log_level = 0;
#endif

   bool print_opt = false; // change to true to check optarg strings
   auto vprt = [print_opt](const char* fmt, const char* arg) -> void {
      if( print_opt ) {
         printf(fmt,arg);
      }
   };

   string optstring = ":L:r:u:vDN:h:o:e:";
#if USE_PROOF != 0
   optstring += "S:a:A:d:E:glp:V:x";
#endif
   int oc = 0;
   while( (oc = getopt(argc,argv,optstring.data())) != -1 ) {
      switch( oc ) {

         case 'L': vprt("L filelist: %s\n",optarg);
                   list_fnames = optarg;
                   break;

         case 'r': vprt("r Prefix: %s\n",optarg);
                   prefix = optarg;
                   break;

         case 'u': vprt("u User function: %s\n",optarg);
                   bean->AddUserFcn(optarg);
                   break;

         case 'v': vprt("%s","v Verbose output\n");
                   verbose = true;
                   bean->SetVerbose();
                   break;

         case 'D': vprt("%s","D events dump\n");
                   bean->SetEventDump();
                   break;

         case 'N': vprt("N 'num' events %s\n",optarg);
                {
                   long val = 0;
                   if ( !str2long(optarg,val) || val <= 0 ) {
                      is_error = true;
                      printf("-N incorrect argument: %s "
                            "(must be a positive integer)\n",optarg);
                   } else {
                      bean->SetMaxNumberEvents(val);
                   }
                   break;
                }

         case 'h': vprt("h hst_file: %s\n",optarg);
                   bean->SetHstFile(optarg);
                   break;

         case 'o': vprt("o out_file: %s\n",optarg);
                   bean->SetDstFile(optarg);
                   break;

         case 'e': vprt("e Add to gEnv: %s\n",optarg);
                {
                   string name, val;
                   tie(name,val) = str2strs(optarg);
                   if ( name.empty() || val.empty() ) {
                      is_error = true;
                      printf("-e incorrect argument: %s "
                            "(must have form: Name=Val)\n",optarg);
                   } else {
                      gEnv->SetValue( name.data(), val.data() );
                   }
                   break;
                }

#if USE_PROOF != 0
         case 'S': vprt("S initialize databases: %s\n",optarg);
                   sqlite_db[0] = bean->GetBaseDir() +
                      "/Analysis/DatabaseSvc/dat";
                   if ( optarg ) {
                      sqlite_db[1] = string(optarg);
                   }
                   break;

         case 'a': vprt("a Set PROOF parameter: %s\n",optarg);
                   bean->SetProofParam(optarg);
                   break;

         case 'A': vprt("A set PROOF INPUT param: %s\n",optarg);
                {
                   auto par = str2strs(optarg);
                   if ( par.first.empty() || par.second.empty() ) {
                      is_error = true;
                      printf("-A incorrect argument: %s "
                            "(must have form: Name=Val)\n",optarg);
                   } else {
                      proof_params.insert(par);
                   }
                   break;
                }

         case 'd': vprt("d Disable some workers: %s\n",optarg);
                   workers_disable_str = optarg;
                   break;

         case 'E': vprt("E PROOF environment: %s\n",optarg);
                {
                   string name, val;
                   tie(name,val) = str2strs(optarg);
                   if ( name.empty() || val.empty() ) {
                      is_error = true;
                      printf("-E incorrect argument: %s "
                            "(must have form: Name=Val)\n",optarg);
                   } else {
                      TProof::AddEnvVar( name.data(), val.data() );
                   }
                   break;
                }

         case 'g': vprt("%s","g save PROOF error logs\n");
                   save_proof_logs = true;
                   break;

         case 'l': vprt("%s","l use PROOF-Lite\n");
                   bean->SetProofClr("lite://");
                   break;

         case 'p': vprt("p PROOF cluster: %s\n",optarg);
                   bean->SetProofClr(optarg);
                   break;

         case 'V': vprt("V 'log_level' for PROOF: %s\n",optarg);
                {
                   long val = 0;
                   if ( !str2long(optarg,val) || val < 0 ) {
                      is_error = true;
                      printf("-V incorrect argument: %s (must be"
                            " not negative integer)\n",optarg);
                   } else {
                      proof_log_level = val;
                   }
                   break;
                }

         case 'x': vprt("%s","x Use Xrootd\n");
                   bean->SetProofXrdOutput(true);
                   break;
#endif

         case ':': // no argument, if leadind ':' in optstring
                   is_error = true;
                   printf("option `-%c' requires an argument\n",
                         optopt);
                   break;

         case '?': // errors
         default:
                   is_error = true;
                   printf(" invalid option `-%c'\n",optopt);
                   break;
      }
   }

   if( is_error ) {
      Usage(argc,argv,verbose);
   }

#if USE_PROOF != 0
   // initialize databases
   if( !sqlite_db[0].empty() ) {
      vprt("sqlite_db[0] = %s\n",sqlite_db[0].data());
      // initialize db
      DatabaseSvc* dbs = DatabaseSvc::instance();
      if( !sqlite_db[1].empty() ) { // copy db to new path
         vprt("sqlite_db[1] = %s\n",sqlite_db[1].data());
         string db_new = copy_dir_temp(sqlite_db[0], sqlite_db[1]);
         dbs->SetDBFilePath(db_new);
      } else {                      //  without argument
         dbs->SetDBFilePath(sqlite_db[0]);
      }
   }
#endif

   // ============= GET LIST OF FILES TO PROCESS  ====================
   vector<string> file_names;

   // the remaining arguments are filenames
   file_names.reserve(32);
   for( int iarg = optind; iarg < argc; ++iarg ) {
      file_names.push_back( prefix+string(argv[iarg]) );
   }

   // read names of files from list_fnames
   if( !list_fnames.empty() ) {
      istream* list_s = nullptr;
      if( list_fnames == "-" ) {
         list_s = &std::cin;
      } else {
         list_s = new ifstream(list_fnames);
         if( !list_s ) {
            printf("Unable to open file '%s'\n",list_fnames.data());
            exit(EXIT_FAILURE);
         }
      }

      // https://en.cppreference.com/w/cpp/string/basic_string/getline
      for( string line; getline(*list_s,line); ) {
         // strip left and right white spaces
         line = regex_replace(line, regex{R"(^\s+|\s+$)"}, "");
         file_names.push_back( prefix+line );
      }

      if( list_s != &cin ) {
         delete list_s;
      }
   }

   if( file_names.empty() ) {
      // skip print if no options or only one option "-v"
      if( argc > 2 || (argc==2 && verbose != true ) ) {
         printf("ERROR: list of input ROOT files is required\n");
      }
      Usage(argc,argv,verbose);
   } else {
      cout << "-- Total " << file_names.size() <<
        " files in the list for processing --" << endl;
      if( verbose ) {
         for ( const auto& f : file_names ) {
            cout << f << endl;
         }
         cout << "-- End of list --" << endl;
      }
   }

   // CHECKING THE CORRECTNESS OF THE PARAMETERS AND PRINTING THEM ==
   // parse filenames and fill in dataset names
   vector<string> dataset_names;
   dataset_names.reserve(file_names.size());
   for( const auto& f : file_names ) {
      string dsn = bean->ParseDatasetName(f);
      if( !dsn.empty() ) {
         dataset_names.push_back(dsn);
      }
   }

   if( !bean->IsProof() ) {
      if( !dataset_names.empty() ) {
         printf("ERROR: Dataset input is not supported without PROOF"
               "\n       filenames in the format 'ds://DatasetName'"
               " are not allowed\n");
         is_error = true;
      }

#if USE_PROOF != 0
      if( bean->DstFileIsDataset() ) {
         printf("ERROR: dataset output "
               "is not allowed in no-PROOF mode\n");
         is_error = true;
      }

      if( !bean->DstFileInSandbox() ) {
         printf("ERROR: remote file output "
               "is not allowed in no-PROOF mode\n");
         is_error = true;
      }
#endif

   } else { // Proof mode
      if( dataset_names.size() != file_names.size() ) {
         printf("ERROR: Mixing datasets "
               "and other files is not supported\n");
         is_error = true;
      }
   }

   if( is_error ) {
      Usage(argc,argv,verbose);
   }
   bean->PrintOptions();

   // ================= HANDLE THE TERMINATION SIGNALS ===============
#if USE_PROOF != 0
   if( bean->IsProof() ) {
      set_proof_termination();
   } else {
      set_user_termination(verbose);
   }
#else
   set_user_termination(verbose);
#endif

   // =========== RUN A LOOP THROUGH EVENTS NON PROOF MODE ===========
   if( !bean->IsProof() ) {
      if( verbose ) {
         cout << "Run a loop through events NON-PROOF mode" << endl;
      }
      TChain chain("Event");
      Long64_t nentries = bean->MaxNumberEvents();

      for( const auto& f : file_names ) {
         chain.Add( f.data() );
      }

      ReadDst* selector = new ReadDst;

      TList* inputList  = new TList;
      inputList->Add(bean);
      selector->SetInputList(inputList);

      if( nentries == 0 ) { // no restriction on the number of events
         chain.Process(selector);
      } else {
         chain.Process(selector,"",nentries);
      }

      // if process was interrupted with Abort() call
      // Terminate functions by hand:
      if( selector->GetAbort() == TSelector::kAbortProcess ) {
         selector->SlaveTerminate();
         selector->Terminate();
      }
      delete selector;
   }

   // display the contents of the memory table
   // gObjectTable->Print();

   exit(EXIT_SUCCESS);

#if USE_PROOF != 0
   // =================== RUN WITH PROOF =============================
   TProof* proof = nullptr;
   gROOT->SetBatch(true); // Prevents ROOT from enabling graphics

   const string& clr = bean->ProofClr();
   const string& param = bean->ProofParam();

   // to force creation of a new session
   string clr_open = clr + "/?N";
   // TProof::Reset(clr.c_str());

   proof = TProof::Open(clr_open.c_str(), param.c_str());
   gEnv->SetValue("Proof.StatsHist",1);
   gEnv->SetValue("Proof.StatsTrace",1);
   gEnv->SetValue("Root.Stacktrace","no"); // do not print stack trace

   proof->SetLogLevel(proof_log_level);

   // set default proof parameters
   proof->SetParameter("PROOF_SavePartialResults", "1");

   // set user-specified proof parameters
   for( auto const& par : proof_params ) {
      proof->SetParameter(par.first.c_str(), par.second.c_str());
   }

   // TODO: do we need it? after or befor TProof::Open ?
   // prooflite_ld_library_path();

   if ( bean->DstFileIsDataset() ) {
      // TODO: make some option to force overwrite
      string dst = bean->DstFile();
      if( proof->GetDataSets(dst.c_str())->GetSize() != 0 ) {
         printf("WARNING: dataset '%s' already exists "
               "on cluster and will be overwriten\n",dst.c_str());
      }
   }

   // Disable some workers
   TIter next(workers_disable_str.Tokenize(","));
   TObjString* worker;
   while( (worker = (TObjString*) next()) ) {
      proof->DeactivateWorker(worker->GetName());
   }

   check_upload_enable(proof, "par/Analysis.par");
   check_upload_enable(proof, "par/BeanUser.par");
   string BeanLib = "par/BeanLib_" + to_string(BOSS_VER) + ".par";
   check_upload_enable(proof, BeanLib.c_str());

   bean->SetProof(proof);
   proof->AddInput(bean);

   Long64_t nentries = bean->MaxNumberEvents();
   if( dataset_names.empty() ) {
      // use chain of file_names
      TChain chain("Event");
      for( const auto& f : file_names ) {
         chain.Add( f.data() );
      }
      chain.SetProof();
      if( !nentries ) {
         chain.Process("ReadDst");
      } else {
         chain.Process("ReadDst","",nentries);
      }
   } else {
      // use datasets string
      auto dataset_string = DatasetStr(dataset_names);
      if( !nentries ) {
         proof->Process(dataset_string.c_str(),"ReadDst");
      } else {
         proof->Process(dataset_string.c_str(),"ReadDst","",nentries);
      }
   }

   // clean up after termination
   if( !sqlite_db[1].empty() ) { // remove copy of database
      DatabaseSvc* dbs = DatabaseSvc::instance();
      string new_dir = dbs->GetDBFilePath();
      rm_whole_dir(new_dir.c_str());
      printf("rm tmp sqlite bd directory: %s\n",new_dir.data());
   }

   graceful_proof_exit();
   exit(EXIT_SUCCESS);
#endif
}

//--------------------------------------------------------------------
void Usage(int argc, char **argv, bool verbose)
//--------------------------------------------------------------------
{
   cout << endl;
   cout << "Usage: " << argv[0] << " [ -option(s)] dst_file(s)\n";
   if( !verbose ) {
      cout << "  Note: This is a short note,"
         " use the '-v' option for more details\n";
   } else {
      cout << "  Note: The program is sutable for dst-files after"
         " the BOSS-" << BOSS_VER << endl;
#if USE_PROOF == 0
      cout << "        Running in PROOF environment"
         " is not supported" << endl;
#endif
   }

   cout << R"(
Arguments:
  dst_file(s)   input ROOT files are specified as:
     /path/to/local/file - path to the local file
     root://user@host/path/to/file - remote file via xrootd protocol
)";
#if USE_PROOF != 0
   if( verbose ) {
      cout << "     ds://DatasetName"
         " - use dataset registered in the PROOF cluster\n";
   }
#endif
   cout << "     Specifying files is mandatory, see also '-L' flag\n";

   cout << R"(
Options for dst_files:
  -L filelist   read dst filenames from given list: one file per line
                to read from standard input stream use '-L-'
  -r prefix     add prefix to each dst_file name
)";

   cout << R"(
Other options:
  -u Uname      use "BeanUser/Uname.cxx" function to process events
                option '-u' could be specified more than once
  -v            set more verbose output
  -D            detailed printout of dst-events
  -N num        process not more than "num" events in total
  -h hst_file   histograms filename, default: bean_histo.root
  -o out_file   file to save selected dst-events, default: no save
)";
#if USE_PROOF != 0
   if( verbose ) { // PROF-experts
      cout << string(R"(
                In PROOF mode the filename could be:
                1) out_file.root - the file will be merged on master
                   node and then fetched back to client
                2) root://user@host/path/to/file.root - to save
                   output at remote xrootd server
                3) ds://DatasetName - to register output as dataset
                   on PROOF cluster (experts prefer)
)").data()+1; // skip first empty line
   }
#endif

   if( verbose ) { // Advanced options
      cout << R"(
Advanced options:
  -e var=value  add variable to gEnv, for example: MyParam=1
)";

#if USE_PROOF != 0
      cout << R"(
  -S path       Copies the 'run' and 'offlinedb' databases to the
                'path/temporary_dir/' and initialize them for use
                in the BEAN; can be used to avoid NFS locking issues

  *** PROOF management ***
  -a params     set PROOF parameters: "valgrind", "workers=42", etc
  -A key=value  set PROOF INPUT parameter: see TProof::SetParameter()
  -d workers    disable specified workers, comma separated
  -E var=value  add variable to PROOF environment, for example:
                PROOF_WRAPPERCMD=valgrind_opts:--leak-check=full
  -g            save proof error logs to bean_proof.log file
  -l            use PROOF-Lite: one PC with multi-core processor
  -p proof_clr  use the specified PROOF cluster
  -V level      set PROOF log level
  -x            use xroot to get output DST-file specified with '-o'
                option from master node)";
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,2)
      cout << ", default: PROOF sandbox\n";
#else
      cout << R"(: only this method available,
                in order to use PROOF sandbox, use ROOT >= 5.25/2
)";
#endif
#endif
   } // end of Advanced options

#if USE_PROOF != 0
   graceful_proof_exit();
#endif
   exit(EXIT_SUCCESS);
}

//--------------------------------------------------------------------
static bool str2long(const char* str, long& val)
//--------------------------------------------------------------------
{
   // https://en.cppreference.com/w/cpp/string/byte/strtol
   char* endptr = nullptr;
   errno = 0;    // To distinguish success/failure after call
   val = strtol(str, &endptr, 10);
   // Check for errors
   if ( errno != 0  || endptr == str ) {
      return false;
   }
   return true;
}

//--------------------------------------------------------------------
static std::pair<std::string,std::string> str2strs(const char* opt)
//--------------------------------------------------------------------
{
   std::string s(opt);
   auto p = s.find('=');
   if ( p == std::string::npos ) {
      return std::make_pair("","");
   }
   return std::make_pair(s.substr(0,p),s.substr(p+1));
}

//--------------------------------------------------------------------
static void set_user_termination(bool verbose)
//--------------------------------------------------------------------
{
   if( verbose ) {
      printf("INFO: start %s()\n",__func__);
   }
   gSystem->ResetSignals(); // Reset 'ROOT-signals' to default

#if defined (__unix__) || defined (__APPLE__)
   struct sigaction new_action, old_action;

   // Set up the structure to specify the new action
   new_action.sa_handler = &termination_handler;
   int ret  = sigemptyset(&new_action.sa_mask);
   if ( ret != 0 ) {
      printf("ERROR: sigemptyset return %i\n", ret);
   }
   new_action.sa_flags = 0;

   // avoid handling signals previously set to be ignored
   sigaction(SIGINT, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
      ret = sigaction(SIGINT, &new_action, NULL);
      if( verbose ) {
         printf("INFO: set SIGINT: ret=  %i\n", ret);
      }
   } else {
      if( verbose ) {
         printf("INFO: SIGINT ignored\n");
      }
   }
   sigaction(SIGTERM, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
      ret = sigaction(SIGTERM, &new_action, NULL);
      if( verbose ) {
         printf("INFO: set SIGTERM: ret=  %i\n", ret);
      }
   } else {
      if( verbose ) {
         printf("INFO: SIGTERM ignored\n");
      }
   }
   sigaction(SIGHUP, NULL, &old_action);
   if( old_action.sa_handler != SIG_IGN ) {
      ret = sigaction(SIGHUP, &new_action, NULL);
      if( verbose ) {
         printf("INFO: set SIGHUP: ret=  %i\n", ret);
      }
   } else {
      if( verbose ) {
         printf("INFO: SIGHUP ignored\n");
      }
   }

   // ignoring SIGSEGV results in undefined behavior
   ret = sigaction(SIGSEGV, &new_action, NULL);
   if( verbose ) {
      printf("INFO: set SIGSEGV: ret=  %i\n", ret);
   }

#elif defined _WIN32
   // https://learn.microsoft.com/en-us/cpp/
   //       c-runtime-library/reference/signal
   signal(SIGINT,  termination_handler);
   signal(SIGTERM, termination_handler);
   signal(SIGSEGV, termination_handler);

   // https://learn.microsoft.com/en-us/
   // windows/console/registering-a-control-handler-function
   if( SetConsoleCtrlHandler(CtrlHandler,TRUE )) {
      if( verbose ) {
         printf("INFO: The Control Handler is installed\n");
      }
   } else {
      printf("ERROR: Could not set CtrlHandler\n");
      exit(EXIT_FAILURE);
   }
#endif
}

//--------------------------------------------------------------------
static void termination_handler(int isig)
//--------------------------------------------------------------------
{
   // POSIX 2008 edition says:
   // the behavior is undefined if the signal handler refers to any
   // object other than 'volatile sig_atomic_t',
   // or if the signal handler calls any function except one of the
   // functions listed in the table...
   // There are __no printf()__ functions in the list,
   // only _Exit() and abort().

   switch( isig ) {
      case SIGINT:  // "program interrupt" (the user types CTRL-C )
         if ( bean_termination != 0 ) { // second CTRL-C
            abort();
         }
      case SIGTERM: // politely ask a program to terminate.
#if defined (__unix__) || defined (__APPLE__)
      case SIGHUP:  // "hang up" - the user's terminal is disconnected
#endif
         bean_termination = isig;
         break;

         // SIGSEGV is the signal sent to a process when it makes an
         //         invalid memory reference, or segmentation fault.
      case SIGSEGV:
         abort();

      default:
         break;
   }
}

#if defined _WIN32
//--------------------------------------------------------------------
BOOL CtrlHandler( DWORD fdwCtrlType )
//--------------------------------------------------------------------
{
   printf(" Ctrl signal '%lu' had been received\n",fdwCtrlType);

   static DWORD CtrlSignal = -1;
   switch( fdwCtrlType ) {

      case CTRL_C_EVENT: // Handle the CTRL-C signal.
         printf("Ctrl-C event\n");
         if( CtrlSignal == CTRL_C_EVENT ) {
            printf("second CRTL-C had been detected. Abort!\n\n");
            abort();
         }
         break;

      case CTRL_CLOSE_EVENT: // confirm that the user wants to exit
         printf("Ctrl-Close event\n");
         break;

      default:
         return FALSE;
   }
   CtrlSignal = fdwCtrlType;

   // event loop interrupt:
   bean_termination = SIGINT;

   return( TRUE );
}
#endif


#if USE_PROOF != 0
//--------------------------------------------------------------------
void graceful_proof_exit()
//--------------------------------------------------------------------
{
   if( gProof ) {
      if( save_proof_logs ) {
         gProof->GetManager()->GetSessionLogs()
            ->Save("*", "bean_proof.log");
      }

      gProof->Close();
      delete gProof->GetManager();
   }
}

//--------------------------------------------------------------------
static void set_proof_termination()
//--------------------------------------------------------------------
{
#if defined (__unix__) || defined (__APPLE__)
   struct sigaction new_action;

   // Set up the structure to specify the new action
   new_action.sa_handler = &proof_segfault_handler;
   sigemptyset(&new_action.sa_mask);
   new_action.sa_flags = 0;

   // ignoring SIGSEGV results in undefined behavior
   sigaction(SIGSEGV, &new_action, NULL);
#elif defined _WIN32
   signal(SIGSEGV, proof_segfault_handler);
#endif
}

//--------------------------------------------------------------------
static void proof_segfault_handler(int isig)
//--------------------------------------------------------------------
{
   // Normal behavior is not guaranteed, but worth a try:
   bean->Proof()->GetManager()->GetSessionLogs()
      ->Save("*", "bean_proof.log");
   abort();
}

//--------------------------------------------------------------------
#if __cplusplus >= 201703L
[[maybe_unused]]
#endif
static void prooflite_ld_library_path()
//--------------------------------------------------------------------
{
   const string& clr = bean->ProofClr();

   // We use manager to check whether this will be Proof-Lite
   TProofMgr* manager = TProofMgr::Create( clr.c_str() );

   if( manager->IsLite() ) {
#if ROOT_VERSION_CODE < ROOT_VERSION(5,29,1)
      TProof::AddEnvVar("ROOTPROOFLITE", "1");
#endif

      // If BEAN is built with xlinked ROOT there is no ROOT
      // libraries in the ld.so search PATH. But proofserv.exe
      // need this libraries to work. So in case of ProofLite we
      // should set LD_LIBRARY_PATH to appropriate one
#ifdef ROOTLIBDIR
      string new_library_path;
      new_library_path += ROOTLIBDIR;
      new_library_path +=":$LD_LIBRARY_PATH";
      TProof::AddEnvVar("LD_LIBRARY_PATH",new_library_path.c_str());
#endif
   }
}

//--------------------------------------------------------------------
static void check_upload_enable(TProof* proof, const char* pname)
//--------------------------------------------------------------------
{
   if( proof->UploadPackage(pname) != 0 ) {
      printf("ERROR: Cannot upload package '%s' on PROOF\n",pname);
      graceful_proof_exit();
      exit(EXIT_FAILURE);
   }
   if( bean->Verbose() ) {
      printf("INFO: Package '%s' successfully uploaded on PROOF\n",
            pname);
   }

   Bool_t notOnClient = kTRUE; // to enable packages also on the client
   if( proof->EnablePackage(pname, notOnClient) != 0 ) {
      printf("ERROR: Cannot enable package '%s' on PROOF\n",pname);
      graceful_proof_exit();
      exit(EXIT_FAILURE);
   }
   if( bean->Verbose() ) {
      printf("INFO: Package '%s' successfully enabled on PROOF\n",
            pname);
   }
}

//--------------------------------------------------------------------
static std::string DatasetStr(
      const std::vector<std::string>& dataset_names)
//--------------------------------------------------------------------
{
   std::string dataset_string;
   for( size_t i = 0; i < dataset_names.size(); ++i ) {
      if( i != 0 ) {
         dataset_string += "|";
      }
      dataset_string += dataset_names[i];
      dataset_string += "#Event";
   }
   return dataset_string;
}
#endif
