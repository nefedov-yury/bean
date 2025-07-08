//////////////////////////////////////////////////////////////////////
//
// TestDb
//
// This is test example of user functions for SqliteInterface
//
// It does not need event loop. Use "bean.exe -N 1 -u TestDb ..."
// to test it.
//
// NOTE: You MUST specify the path to the database files.
//       Use the function AbsPath(file) which returns absolute path
//       in the form: "path_to_base_Bean_dir/" + file
//
//////////////////////////////////////////////////////////////////////

#include "DLLDefines.h"         // mandatory!

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "TBufferFile.h"
#include "TTree.h"

#include "RootEventData/TEvtHeader.h"
#include "ReadDst.h"
#include "DatabaseSvc/DatabaseSvc.h"
#include "VertexFit/VertexDbSvc.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

static DatabaseSvc* m_dbsvc = 0;

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestDbStartJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestDbStartJob() " << endl;

   // initialize access to Db:
   m_dbsvc = DatabaseSvc::instance();

   // set path to directory with databases:
   m_dbsvc->SetDBFilePath(
         selector->AbsPath("Analysis/DatabaseSvc/dat") );

   // ----------------------------------------------------------------
   // example of using VertexDbSvc::ReadOneTime()
   // ----------------------------------------------------------------
   VertexDbSvc* vtxsvc = VertexDbSvc::instance();
   // vtxsvc -> SetBossVer("6.6.3");
   // vtxsvc -> ReadOneTime(9947,10878);
   vtxsvc -> SetBossVer("7.0.5");
   vtxsvc -> ReadOneTime(64314,64320);
   vector<int> RunList = vtxsvc -> ListReadRuns();
   cout << " read " << RunList.size() << " runs" << endl;
   // for ( int r : RunList ) {
   // cout << " " << r;
   // }
   // cout << endl;

   // ----------------------------------------------------------------
   // this is example from DatabaseSvc/TestDbAlg
   // ----------------------------------------------------------------

   cout << "******************************************************\n";
   cout << "Test 1: Numbers" << std::endl;
   cout << "******************************************************\n";

   char stmt[255];
   snprintf( stmt, sizeof(stmt),
         "select Vx, Vy, Vz, SigmaVx, SigmaVy, SigmaVz, RunNo "
         "from BeamPar "
         "where RunNo > 9947 and RunNo < 10100 and SftVer='6.6.3'" );

   DatabaseRecordVector res;
   int row_no = 0;
   row_no = m_dbsvc->query("offlinedb",stmt,res);
   if( row_no < 0 ){
      std::cout << "Query \"" << stmt << "\" failed" << std::endl;
      // return;
   }

   double vx = 0;
   double svx = 0;
   for(int row = 0; row < row_no; row++) {
      DatabaseRecord& records = *res[row];
      sscanf(records["Vx"], "%lf", &vx);
      sscanf(records["SigmaVx"], "%lf", &svx);
      cout << "Read from DB: RunNo " << records["RunNo"]
         << " Vx= " << vx << " SigmaVx= " << svx << endl;
      cout << "                    " << records.GetLong("RunNo")
         << " Vx= " << records.GetDouble("Vx")
         << " SigmaVx= " << records.GetDouble("SigmaVx") << endl;
   }

   // Check BLOBs (you need full 'offlinedb' database for this check)
   cout << "******************************************************\n";
   cout << "Test 2: BLOBs" << std::endl;
   cout << "******************************************************\n";

   res.clear();
   snprintf( stmt, sizeof(stmt), "select EndTofPar,BarTofPar "
         "from TofCalConst where RunFrom <= 11000 and RunTo >= 11000 "
         "and SftVer='6.6.3'" );
   row_no = m_dbsvc->query("offlinedb", stmt, res);
   if( row_no < 0 ){
      cout << " Query \"" << stmt << "\" failed" << endl;
      cout << " Try to use full 'offlinedb' database for this check"
         << endl;
   }

   for(int row = 0; row < row_no; row++) {
      DatabaseRecord& records = *res[row];
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
      TBufferFile *buf1 = new TBufferFile(TBuffer::kRead);
#else
      TBuffer *buf1 = new TBuffer(TBuffer::kRead);
#endif
      buf1->SetBuffer(records["EndTofPar"],32768,kFALSE);
      TTree* curvetree = new TTree();
      curvetree->Streamer(*buf1);
      double cnvAtten[8];
      curvetree -> SetBranchAddress("Atten0", &cnvAtten[0]);
      curvetree -> SetBranchAddress("Atten1", &cnvAtten[1]);
      curvetree -> SetBranchAddress("Atten2", &cnvAtten[2]);
      curvetree -> SetBranchAddress("Atten3", &cnvAtten[3]);
      curvetree -> SetBranchAddress("Atten4", &cnvAtten[4]);
      Long64_t entries=curvetree->GetEntries();
      if(entries>10) entries = 10;
      for(int iiii = 0; iiii < entries; iiii++) {
         curvetree->GetEntry(iiii);
         for(int jjj = 0; jjj < 5; jjj++) {
            std::cout << "cnvAtten[" << jjj
               << "]=" << cnvAtten[jjj] << " ";
         }
         std::cout << std::endl;
      }
   }

   // Test strings (you need full 'offlinedb' database for this check)
   cout << "******************************************************\n";
   cout << "Test 3: Strings" << std::endl;
   cout << "******************************************************\n";

   res.clear();
   snprintf( stmt, sizeof(stmt),
         "select XtTree,QtTree,T0Tree,SdTree,RunFrom,RunTo,CalParVer,"
         "FileName from MdcCalConst where SftVer = '6.6.3'" );
   row_no = m_dbsvc->query("offlinedb",stmt,res);
   if( row_no < 0 ) {
      cout << " Query \"" << stmt << "\" failed" << endl;
      cout << " Try to use full 'offlinedb' database for this check"
         << endl;
      // return;
   }

   for(int row = 0; row < row_no; row++) {
      DatabaseRecord& records = *res[row];
      cout << "Read from DB: Runs " << records["RunFrom"] << " "
         << records["RunTo"]
         << " FileName = " << records["FileName"] << endl;
   }

}

//-----------------------------------------------------------------
BeanUserShared_EXPORT
bool TestDbEvent(ReadDst* selector,
                 TEvtHeader* m_TEvtHeader,
                 TDstEvent* m_TDstEvent,
                 TEvtRecObject* m_TEvtRecObject,
                 TMcEvent* m_TMcEvent,
                 TTrigEvent* m_TTrigEvent,
                 TDigiEvent* m_TDigiEvent,
                 THltEvent* m_THltEvent)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestDbEvent() " << endl;
   int run = abs( m_TEvtHeader->getRunId() );

   // ----------------------------------------------------------------
   // this is test for MagneticField/ConnectionDB
   // ----------------------------------------------------------------

   string stmt("select Magnet_Current,SCQL,SCQR "
         "from SC_magnet where run_number = ");
   stmt += to_string(run);
   DatabaseRecordVector res;
   int row_no = m_dbsvc->query("run",stmt,res);

   if( row_no != 1 ) {
      cout << " ERROR: can not found Magnet_Current information"
         << endl << "QUERY: " << stmt << endl;
   } else {
      // res.PrintRecords();
      DatabaseRecord& rec = *res[0];
      double ssm_curr = rec.GetDouble("Magnet_Current");
      double scql_curr = rec.GetDouble("SCQL");
      double scqr_curr = rec.GetDouble("SCQR");

      cout << " currents= " << ssm_curr << " "
         << scql_curr << " " << scqr_curr << endl;
   }

   // ----------------------------------------------------------------
   // this is an example to use the VertexFit/VertexDbSvc.h service
   // ----------------------------------------------------------------

   VertexDbSvc* vtxsvc = VertexDbSvc::instance();
   vtxsvc->handle(run); // MUST BE !
   cout << " status: " << vtxsvc->isVertexValid() << endl;
   double* a1 = vtxsvc->PrimaryVertex();
   double* a2 = vtxsvc->SigmaPrimaryVertex();
   cout << " vx: " << a1[0] << " vy: " << a1[1] << " vz: " << a1[2]
      << endl;
   cout << " vx sigma: " << a2[0] << " vy sigma: " << a2[1]
      << " vz sigma: " << a2[2] << endl;

   // ----------------------------------------------------------------
   // this is example from DatabaseSvc/TestDbAlg
   // Check BLOBs => this works only for full db !!!
   stmt =  "select EndTofPar,BarTofPar from TofCalConst where ";
   stmt += "RunFrom <= 11000 and RunTo >= 11000 and SftVer='6.6.3'";
   res.clear();
   row_no = m_dbsvc->query("offlinedb", stmt.c_str(), res);
   // cout << " row_no= " << row_no << endl;

   for(int row = 0; row < row_no; row++) {
      DatabaseRecord& records = *res[row];
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
      TBufferFile *buf1 = new TBufferFile(TBuffer::kRead);
#else
      TBuffer *buf1 = new TBuffer(TBuffer::kRead);
#endif
      buf1->SetBuffer(records["EndTofPar"],32768,kFALSE);
      TTree* curvetree = new TTree();
      curvetree->Streamer(*buf1);
      double cnvAtten[8];
      curvetree -> SetBranchAddress("Atten0", &cnvAtten[0]);
      curvetree -> SetBranchAddress("Atten1", &cnvAtten[1]);
      curvetree -> SetBranchAddress("Atten2", &cnvAtten[2]);
      curvetree -> SetBranchAddress("Atten3", &cnvAtten[3]);
      curvetree -> SetBranchAddress("Atten4", &cnvAtten[4]);
      int entries=curvetree->GetEntries();
      for(int iiii = 0; iiii < entries; iiii++) {
         curvetree->GetEntry(iiii);
         for(int jjj = 0; jjj < 5; jjj++) {
            cout << " cnvAtten[" << jjj
               << "]=" << cnvAtten[jjj] << " ";
         }
         cout << endl;
      }
   }
   // ----------------------------------------------------------------

   return false;
}

//--------------------------------------------------------------------
BeanUserShared_EXPORT
void TestDbEndJob(ReadDst* selector)
//--------------------------------------------------------------------
{
   if( selector->Verbose() ) cout << " TestDbEndJob() " << endl;
}

#ifdef __cplusplus
}
#endif
