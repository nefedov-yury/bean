#include "EventTag/DecayTable.h"

#include <algorithm>
#include <TDatabasePDG.h>   // from ROOT for names of particle

using namespace std;

//====================================================================
string DecayTable::StrDec() {
//====================================================================
   // convert vecdec to string with pdg names
   string strdec;
   if ( !vecdec.empty() ) {
      // sort by absolute value but positive first
      sort( vecdec.begin(), vecdec.end(),
            [](long int x, long int y) {
                  long int ax = abs(x), ay = abs(y);
                  return (ax == ay) ? x>y : ax>ay;
            }
          );
      for ( const auto& pdg : vecdec ) {
         if ( !strdec.empty() ) {
            strdec += string(" ");
         }
         TParticlePDG* pdgParticle =
            TDatabasePDG::Instance()->GetParticle(pdg);
         strdec += ( pdgParticle ) ?
            string(pdgParticle->GetTitle()) : to_string(pdg);
      }
   } else {
      strdec = string("unknown");
   }
   return strdec;
}

//====================================================================
void DecayTable::Add() {
//====================================================================
   // fill the table
   ++ntot;
   decays[ StrDec() ] += 1;
}

//====================================================================
void DecayTable::Print(double min_percent) {
//====================================================================
   // print sorted by values (in decreasing order of decay frequency)
   // min_percent is the minimal % for printing of decays
   vector<pair<string,size_t>>
      vprt( decays.begin(),decays.end() );
   sort( vprt.begin(), vprt.end(),
         []( pair<string,size_t> elem1, pair<string,size_t> elem2 ) {
               return elem1.second > elem2.second;
         }
       );
   double R = 100./double( ntot );
   size_t sum_ev = 0;
   size_t vsize = vprt.size();
   for ( size_t i = 0; i < vsize; ++i ) {
      const auto& p = vprt[i];
      string name = p.first;
      size_t ev = p.second;
      double pc = R*ev;
      // there is no point to print 'others' if it is only one case
      bool print_other = (pc < min_percent) && (i != vsize-1);
      if( !print_other ) {
         sum_ev += ev;
         printf(" %-35.35s %9zu (%10.4f %%)\n", name.c_str(),ev,pc);
      } else {
         ev = ntot - sum_ev;
         pc = R*ev;
         printf(" %-35.35s %9zu (%10.4f %%)\n", "other decays",ev,pc);
         break;
      }
   }
}
