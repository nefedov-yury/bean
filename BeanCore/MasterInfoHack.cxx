//////////////////////////////////////////////////////////////////////
//                                                                  //
// MasterInfoHack                                                   //
// * The class to get workdir location on master when using PROOF   //
//   This is done using the 'fake' Merge method which executes on   //
//   master                                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include "MasterInfoHack.h"

ClassImp(MasterInfoHack);

MasterInfoHack::MasterInfoHack()
{
   // fDir = std:string();
   this->SetName("MasterInfoHack");
}

Long64_t MasterInfoHack::Merge(TCollection* list)
{
   fDir = gSystem->WorkingDirectory();
   return 0;
}
