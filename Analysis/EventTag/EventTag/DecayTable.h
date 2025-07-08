#ifndef DecayTable_H
#define DecayTable_H
//====================================================================
//
// Table of decays:
// *) variables:
//    decays: table string with final particles -> number of decays
//            this is "global" table one per job
//    vecdes: vector with decay products (pdg-codes) of studied
//            particle; this variable may change for each event
//    ntot:   total number of decays
// *) functions:
//    Size():   size of table
//    Reset():  clear vecdec
//    StrDec(): convert vecdec to string with pdg names
//    Add():    fill in the table according to the decay stored
//              in the vecdec
//    Print(min_percent):
//              print table sorted by percentage of each decay
//              up to min_percent
//
//====================================================================

#include <map>
#include <vector>
#include <string>
#include <cstdio> // for size_t

struct DecayTable {
   std::map<std::string,size_t> decays; // table for the studied particle
   std::vector<long int>        vecdec; // pdg-codes; one per event
   size_t                       ntot;   // total number of decays

   DecayTable() { // constructor
      vecdec.reserve(8);
      ntot = 0;
   }

   size_t Size() { // size of table
      return decays.size();
   }

   void Reset() { // clear vecdec
      vecdec.clear();
   }

   std::string StrDec(); // convert vecdec to string with pdg names
   void Add();           // fill in the table
   void Print(double min_percent = 0); // print sorted by values;
};
#endif
