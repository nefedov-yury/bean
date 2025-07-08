#ifndef _Bean_h
#define _Bean_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// Bean                                                             //
//  * BEAN parameters                                               //
//  * Dynamic load and storage of user functions                    //
//  * Basic parameters for working in PROOF mode                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <TObject.h>

#if USE_PROOF != 0
class TProof;
#endif

typedef std::vector<void* > VecUF;
typedef std::vector<std::string> VecName;

class Bean : public TObject
{
   public:
      Bean();
      virtual      ~Bean()=default;

      // -- set/get options functions:
      void SetBaseDir(const char* dir) {base_dir = dir;}
      void SetVerbose()                {verbose = true;}
      void SetSilent()                 {verbose = false;}
      void SetEventDump()              {event_dump = true;}
      void SetHstFile(const char* opt) {hst_file = opt;}
      void SetDstFile(const char* opt);
      void SetMaxNumberEvents(long val) {max_number_events = val;}

      std::string GetBaseDir() const      {return base_dir;}
      bool        Verbose() const         {return verbose;}
      bool        IsDump() const          {return event_dump;}
      std::string HstFile() const         {return hst_file;}
      std::string DstFile() const         {return dst_file;}
      std::string DstFileName() const     {return dst_file_name;}
      long        MaxNumberEvents() const {return max_number_events;}

      void PrintOptions() const;

      // -- user functions
      void         AddUserFcn(const char* name);
      size_t       NUserFns() const        {return Ufn_names.size();}
      const VecUF& GetStartJobFns() const  {return Ufn_start;}
      const VecUF& GetUserEventFns() const {return Ufn_event;}
      const VecUF& GetEndJobFns() const    {return Ufn_end;}

      void LoadUserLib();
      void LoadUserFcn(const std::string& name);
      void LoadUserFcns();

      // IsProof() should be available even if PROOF is not available
      bool IsProof() const {
#if USE_PROOF != 0
         return !proof_clr.empty();
#else
         return false;
#endif
      }
      // to check that the filename is the name of the dataset
      std::string ParseDatasetName(std::string filename);

   private:
      std::string base_dir; // directory with Analysis, BeanUser, etc.
      bool        verbose;    // more verbose output
      bool        event_dump; // detailed printout of dst-events
      std::string hst_file;   // histograms filename
      std::string dst_file;   // file to save selected dst-events
      std::string dst_file_name; // only dst_file name without path
      long        max_number_events;// max number of events to process

      VecUF       Ufn_start; //! vector of user functions StartJob
      VecUF       Ufn_event; //! vector of user Event functions
      VecUF       Ufn_end;   //! vector of user functions EndJob
      VecName     Ufn_names; // names of loaded functions

#if USE_PROOF != 0
   public:
      void SetProof(TProof* val);
      void SetProofClr(const char* opt)   {proof_clr=opt;}
      void SetProofParam(const char* opt) {proof_param=opt;}
      void SetProofXrdOutput(bool val)    {proof_xrd_output=val;}

      TProof*     Proof() const      {return proof;}
      std::string ProofClr() const   {return proof_clr;}
      std::string ProofParam() const {return proof_param;}

      bool IsProofLite() const      {return proof_lite;}
      bool DstFileInSandbox() const {return dst_file_in_sandbox;}
      bool DstFileIsDataset() const {return dst_file_is_dataset;}
      bool ProofXrdOutput() const   {return proof_xrd_output;}

      std::string GetProofMasterWorkdir();

   private:
      TProof*     proof;       //! this is not-persistent member
      std::string proof_clr;
      std::string proof_param;

      bool        proof_lite;
      bool        dst_file_in_sandbox; // store in 'local' sandbox
      bool        dst_file_is_dataset; // store in dataset
      bool        proof_xrd_output;    // use Xroot to get dst file
#endif

   // ClassVersionID=1 because we have need object I/O
   ClassDef(Bean,1); // Primary class to read DST
};
#endif
