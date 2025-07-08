#ifndef BEAN
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Bootstrap.h"

#include "facilities/Util.h"
#else
#include "DatabaseSvc/DatabaseSvc.h"
#endif

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "MagneticField/ConnectionDB.h"

using namespace std;

namespace FieldDBUtil {

#ifndef BEAN
  ConnectionDB::ConnectionDB() {
    StatusCode sc = Gaudi::svcLocator()->service("DatabaseSvc", m_dbsvc, true);
    if (sc .isFailure() ) {
      std::cout << "ERROR: In ConnectionDB()--> Unable to find DatabaseSvc " << std::endl;
    }
  }
#else
static DatabaseSvc* m_dbsvc = 0;

void ConnectionDB::SetDBFilePath(std::string new_DBFilePath){
    if( !m_dbsvc ) {
      m_dbsvc = DatabaseSvc::instance();
    }
    m_dbsvc->SetDBFilePath(new_DBFilePath);
}
#endif

  ConnectionDB::eRet ConnectionDB::getReadSC_MagnetInfo( std::vector<double>& current, int runNo) {

#ifdef BEAN
    if( !m_dbsvc ) {
      m_dbsvc = DatabaseSvc::instance();
    }
#endif

    //Read magnetic field map
    char stmt1[200];
    int run_No = std::abs(runNo);

    // sprintf(stmt1,"select Magnet_Current,SCQL,SCQR from SC_magnet where run_number = %d ",run_No);
    snprintf( stmt1, sizeof(stmt1), "select Magnet_Current,SCQL,SCQR "
          "from SC_magnet where run_number = %d ", run_No );
    //cout<<stmt1<<endl;
    DatabaseRecordVector results;
    results.clear();
    int status = m_dbsvc->query("run",stmt1,results);

    if(status<0){
      std::cout << "ERROR Read the SSM and SCQ current from the Database" << endl;
      exit(1);
    }

    int RowNumber = results.size();
    if(RowNumber!=1){
      std::cout<<"ERROR:error searching SC_Magnet Data in the database, check your selection criterions"<<std::endl;
      return RETMySQLError;
    }

    DatabaseRecord& rec = *results[0];
    double ssm_curr = rec.GetDouble("Magnet_Current");
    double scql_curr = rec.GetDouble("SCQL");
    double scqr_curr = rec.GetDouble("SCQR");

    // save results in vector
    current.resize(3);
    current[0] = ssm_curr;
    current[1] = scql_curr;
    current[2] = scqr_curr;
    //cout<<"run:"<<run_No<<"Magnet_Current:"<<current[0]<<"SCQL:"<<current[1]<<"SCQR:"<<current[2]<<endl;
    return RETOk;
  }

bool ConnectionDB::getBeamEnergy( std::map<int, std::vector<double> >& m_mapBeamEnergy, int runFrom, int runTo)
{
  char stmt1[200];
   int run_From = std::abs(runFrom);
   int run_To = std::abs(runTo);
   DatabaseRecordVector results;
   results.clear();

#ifndef BEAN
   IDatabaseSvc* m_dbsvc;
  StatusCode  sc = Gaudi::svcLocator()->service("DatabaseSvc",m_dbsvc,true);
  if (sc .isFailure() ) {
    std::cout<< "MSG::ERROR " << "Unable to find DatabaseSvc " << std::endl;
    exit(1) ;
  }
#else
    if( !m_dbsvc ) {
      m_dbsvc = DatabaseSvc::instance();
    }
#endif

  // sprintf(stmt1,"select run_number,BPR_PRB,BER_PRB from RunParams where run_number >= %d and run_number <= %d ",run_From,run_To);
  snprintf( stmt1, sizeof(stmt1), "select run_number,BPR_PRB,BER_PRB "
        "from RunParams where run_number >= %d and run_number <= %d ",
        run_From, run_To );
  //cout<<stmt1<<endl;
  int row_no = m_dbsvc->query("run",stmt1,results);
  if(row_no<=0){
      std::cout <<"ERROR:"<< "Run:"<<run_From<<" Can not read the beam energy from the Database" << endl;
      exit(1);                                                                                                             }

   /*int RowNumber = results.size();
   if(RowNumber == 0) {
       beamE.push_back(1.843); // for positron
       beamE.push_back(1.843); // for electron
       std::cout << "No beam energy, so set BPR_PRB=1.843, BER_PRB=1.843"<<std::endl;
      exit(1);
      }*/
   if( row_no > 0 )
     {
      for(int i=0;i<row_no;i++)
        {
         DatabaseRecord& dbrec = *results[i];
         int run_No = dbrec.GetInt("run_number");
         std::vector<double> beamEnergy;
         beamEnergy.push_back(dbrec.GetDouble("BPR_PRB"));
         beamEnergy.push_back(dbrec.GetDouble("BER_PRB"));
         m_mapBeamEnergy[run_No] = beamEnergy;
#ifndef BEAN
         float beam1,beam2;
         beam1=beamEnergy[0];
         beam2=beamEnergy[1];
         //cout<<"map of run:"<<run_No<<"BPR_PRB:"<<beam1<<"BER_PRB:"<<beam2<<endl;
#endif

        }
      return true;
    }
   return false;
}

  bool ConnectionDB::getReadSC_MagnetInfo(std::map<int, std::vector<double> >& m_mapMagnetInfo, int runFrom, int runTo)
  {
   char stmt1[200];
   int run_From = std::abs(runFrom);
   int run_To = std::abs(runTo);
   DatabaseRecordVector results;
   results.clear();

#ifndef BEAN
   IDatabaseSvc* m_dbsvc;
  StatusCode  sc = Gaudi::svcLocator()->service("DatabaseSvc",m_dbsvc,true);
  if (sc .isFailure() ) {
    std::cout<< "MSG::ERROR " << "Unable to find DatabaseSvc " << std::endl;
    exit(1) ;
  }
#else
    if( !m_dbsvc ) {
      m_dbsvc = DatabaseSvc::instance();
    }
#endif
   // sprintf(stmt1,"select run_number,Magnet_Current,SCQL,SCQR from SC_magnet where run_number >= %d and run_number<=%d ",run_From,run_To);
   snprintf( stmt1, sizeof(stmt1),
         "select run_number,Magnet_Current,SCQL,SCQR from SC_magnet "
         "where run_number >= %d and run_number<=%d ",
         run_From, run_To );
   //cout<<stmt1<<endl;
   int row_no = m_dbsvc->query("run",stmt1,results);
   if(row_no<=0){
      std::cout << "ERROR Read the SSM and SCQ current from the Database"<< " Run:"<<run_From << endl;                                          exit(1);                                                                                                             }
   if( row_no > 0 )
    {
     for(int i=0;i<row_no;i++)
      {
       DatabaseRecord& dbrec = *results[i];
       int run_No = dbrec.GetInt("run_number");
       std::vector<double> SCMagnet;
       SCMagnet.push_back(dbrec.GetDouble("Magnet_Current"));
       SCMagnet.push_back(dbrec.GetDouble("SCQL"));
       SCMagnet.push_back(dbrec.GetDouble("SCQR"));
       m_mapMagnetInfo[run_No] = SCMagnet;
       //cout<<"map of run:"<<run_No<<"Magnet_Current:"<<(m_mapMagnetInfo[run_No])[0]<<"SCQL:"<<(m_mapMagnetInfo[run_No])[1]<<"SCQR:"<<(m_mapMagnetInfo[run_No])[2]<<endl;
      }

     }
   return true;

   }

  ConnectionDB::eRet ConnectionDB::getBeamEnergy( std::vector<double>& beamE, int runNo) {

#ifdef BEAN
    if( !m_dbsvc ) {
      m_dbsvc = DatabaseSvc::instance();
    }
#endif

    //Read magnetic field map
     char stmt1[200];
     int run_No = std::abs(runNo);

     // sprintf(stmt1,"select BPR_PRB,BER_PRB from RunParams where run_number = %d ",run_No);
     snprintf( stmt1, sizeof(stmt1), "select BPR_PRB,BER_PRB "
           "from RunParams where run_number = %d ", run_No );
     //cout<<stmt1<<endl;
     DatabaseRecordVector results;
     results.clear();
     int status = m_dbsvc->query("run",stmt1,results);
     if(status<0){
       std::cout << "ERROR Read the beam energy from the Database" << std::endl;
       exit(1);
     }

     int RowNumber = results.size();

     if(RowNumber == 0) {
       beamE.push_back(1.843); // for positron
       beamE.push_back(1.843); // for electron

       return RETOk;
     }

     if(RowNumber!=1){
         std::cout<<"ERROR:error searching beam energy in the database, check your selection criterions"<<std::endl;
        return RETMySQLError;
     }

     beamE.push_back(atof((*results[0])["BPR_PRB"])); // for positron
     beamE.push_back(atof((*results[0])["BER_PRB"])); // for electron
#ifndef BEAN
     float beam1,beam2;
     beam1=atof((*results[0])["BPR_PRB"]);
     beam2=atof((*results[0])["BER_PRB"]);
     //cout<<"map of run:"<<run_No<<"BPR_PRB:"<<beam1<<"BER_PRB:"<<beam2<<endl;
#endif
     return RETOk;
  }

}
