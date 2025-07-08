//////////////////////////////////////////////////////////////////////
//                                                                  //
// TMergeableMap                                                    //
//  * Tmap with merging support: to add an object to fOutput,       //
//    the Merge() method MUST BE implemented                        //
//                                                                  //
//////////////////////////////////////////////////////////////////////
// #include <iostream>
#include <TObjString.h>

#include "TMergeableMap.h"

ClassImp(TMergeableMap);

Long64_t TMergeableMap::Merge(TCollection* li)
{
   if (!li) {
      return 0;
   }
   TIter next(li);
   while( TMap* map = (TMap*)next() ) {
      if( map == this ) {
         continue;
      }
      TMapIter next_key(map);
      while( TObjString* key = (TObjString*)next_key() ) {
         //check whether this key already exists in map
         TObjString* present_value =
            (TObjString*) GetValue( key->GetName() );
         TObjString* new_value = (TObjString*) map->GetValue( key );
         if( present_value ) {
            if(new_value->GetString() != present_value->GetString()) {
               Error(
                     "TMergeableMaps::Merge",
                     "Trying to merge TMergeableMaps with different"
                     " values(%s,%s) of the same key %s!",
                     new_value->GetString().Data(),
                     present_value->GetString().Data(),
                     key->GetString().Data()
                    );
            }
         } else {
            Add(key, new_value);
         }
      }
   }
   return 0;
}
