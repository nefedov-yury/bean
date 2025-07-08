#include "RootEventData/TTrigEvent.h"
#include <iostream>
#include "TCollection.h"  // Declares TIter

ClassImp(TTrigEvent)

// Allocate the TObjArray just once
TObject *TTrigEvent::s_staticTrigData = 0;

//***************************************************************
TTrigEvent::TTrigEvent() 
{
  if (! s_staticTrigData ) {
    s_staticTrigData = new TObject();
  }

  m_trigData = s_staticTrigData;

  Clear();
}

//*****************************************************************
TTrigEvent::~TTrigEvent() {
  if(m_trigData == s_staticTrigData ) s_staticTrigData = 0;
    //m_trigData->Delete();
    delete m_trigData;
    m_trigData = 0;
}

//*****************************************************************
void TTrigEvent::initialize(Bool_t fromMc){ 
    m_fromMc = fromMc;
}
  
//*****************************************************************
void TTrigEvent::Clear(Option_t *option) {

}

//*****************************************************************************
void TTrigEvent::Print(Option_t *option) const {
    TObject::Print(option);
}

///TrigData
void  TTrigEvent::addTrigData(TTrigData * trigData){
    m_trigData = trigData;
}


const TTrigData*  TTrigEvent::getTrigData() const {
        return (TTrigData*)m_trigData ;
}
