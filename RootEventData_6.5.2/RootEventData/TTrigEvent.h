#ifndef ROOT_TTrigEvent_H
#define ROOT_TTrigEvent_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TObjArray.h"

#include "TTrigData.h"

class TTrigEvent: public TObject {
public:

    TTrigEvent();
    virtual ~TTrigEvent();

    void initialize( Bool_t fromMc=true);

    void Clear(Option_t *option="");
 
    void Print(Option_t *option="") const;


    inline Bool_t getFromMc() { return m_fromMc; };

    //TrigData
   void   addTrigData(TTrigData * trigData);
   const  TTrigData*  getTrigData() const;
   void  clearTrigData() { m_trigData->Clear();}
     
private:

    
    /// Denote whether or not this data was simulated
    Bool_t m_fromMc;

    /// data members to store trigger data
    static TObject* s_staticTrigData;
    TObject* m_trigData;

    ClassDef(TTrigEvent,1) // Storage for trigger event and subsystem data
}; 
 
#endif





