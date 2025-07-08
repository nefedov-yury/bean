#include "RootEventData/TMcHitTof.h"
#include <iostream>

ClassImp(TMcHitTof)

//Int_t TMcHitTof::s_count = 0;  
//************************************************
TMcHitTof::TMcHitTof() {
   //s_count++;
   //std::cout << "TMcHitTof count: " << s_count << std::endl;
   Clear();
}
//************************************************
TMcHitTof::~TMcHitTof (){
   //s_count--;
   //std::cout << "TMcHitTof count: " << s_count << std::endl;
   Clear();
}

