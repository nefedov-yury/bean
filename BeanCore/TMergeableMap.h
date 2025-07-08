#ifndef _TMergeableMap_h
#define _TMergeableMap_h
//////////////////////////////////////////////////////////////////////
//                                                                  //
// TMergeableMap                                                    //
//  * Tmap with merging support: to add an object to fOutput,       //
//    the Merge() method MUST BE implemented                        //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <TMap.h>

class TMergeableMap : public TMap
{
   public:
      Long64_t Merge(TCollection* list);

      ClassDef(TMergeableMap,1); // Tmap with merging support
};
#endif
