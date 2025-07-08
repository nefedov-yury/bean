#ifndef _MasterInfoHack_h
#define _MasterInfoHack_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// MasterInfoHack                                                   //
// * The class to get workdir location on master when using PROOF   //
//   This is done using the 'fake' Merge method which executes on   //
//   master                                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <string>
#include <TNamed.h>
#include <TCollection.h>

class MasterInfoHack : public TNamed
{
   public:
      MasterInfoHack();
      ~MasterInfoHack() = default;

      Long64_t Merge(TCollection* list);
      std::string MasterWorkdir() const {return fDir;}

   private:
      std::string fDir;

   ClassDef(MasterInfoHack,1);
};
#endif
