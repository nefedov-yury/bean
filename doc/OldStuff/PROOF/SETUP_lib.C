Int_t SETUP()
{
   // This is automatically generated file.
   // Edit "cmake/SETUP_lib.C" file.

   // Set LD_LIBRARY_PATH
   char ld_lib_path[1024];
   sprintf(ld_lib_path,"%s:%s",
           gSystem->Getenv("LD_LIBRARY_PATH"),gSystem->WorkingDirectory());
   gSystem->Setenv("LD_LIBRARY_PATH",ld_lib_path);
   cout << " LD_LIBRARY_PATH= " << gSystem->Getenv("LD_LIBRARY_PATH") << endl;

   // this function should return -1 on failure and 0 on success
   // cout << " Load " << endl;
   // if ( gSystem->Load("") < 0 ) return -1;

   cout << " Load libTree" << endl;
   if ( gSystem->Load("libTree") < 0 ) return -1;

   cout << " Load libHist" << endl;
   if ( gSystem->Load("libHist") < 0 ) return -1;

   cout << " Load libPhysics" << endl;
   if ( gSystem->Load("libPhysics") < 0 ) return -1;

   cout << " Load libEG" << endl;
   if ( gSystem->Load("libEG") < 0 ) return -1;

   cout << " Load libMLP" << endl;
   if ( gSystem->Load("libMLP") < 0 ) return -1;

   cout << " Load libRootEventData" << endl;
   if ( gSystem->Load("libRootEventData") < 0 ) return -1;

   cout << " Load libBeanCore" << endl;
   if ( gSystem->Load("libBeanCore") < 0 ) return -1;

   cout << " Load libDatabaseSvc" << endl;
   if ( gSystem->Load("libDatabaseSvc") < 0 ) return -1;

   cout << " Load libMagneticField" << endl;
   if ( gSystem->Load("libMagneticField") < 0 ) return -1;

   cout << " Load libVertexFit" << endl;
   if ( gSystem->Load("libVertexFit") < 0 ) return -1;

   cout << " Load libParticleID" << endl;
   if ( gSystem->Load("libParticleID") < 0 ) return -1;

   cout << " Load libEventTag" << endl;
   if ( gSystem->Load("libEventTag") < 0 ) return -1;

   cout << " Load libAbsCor" << endl;
   if ( gSystem->Load("libAbsCor") < 0 ) return -1;

   cout << " Load libBeanUser" << endl;
   if ( gSystem->Load("libBeanUser") < 0 ) return -1;

   cout << " lib/PROOF-INF/SETUP() done" << endl;
   return 0;
}
