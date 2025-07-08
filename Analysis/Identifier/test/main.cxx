#include <iostream>

#include "Identifier/Identifier.h"
#include "Identifier/EmcID.h"
#include "Identifier/HltID.h"

using namespace std;

int main()
{
   Identifier id;
   id=0;

   
   int x=1,y=2,z=3;

   cout<<id;

   cout.width(8);
   cout.fill('*');
   cout<<"Hello"<<endl;

   cout<<"==============EmcID Test================="<<endl;

   id=EmcID::crystal_id(x,y,z);

   cout<<"Identifier id="<<id<<endl;
   
   cout<<"barrel_ec    ="<<EmcID::barrel_ec(id)<<endl;
   cout<<"theta_module ="<<EmcID::theta_module(id)<<endl;
   cout<<"phi_module   ="<<EmcID::phi_module(id)<<endl;
   cout<<"is barrel?    "<<EmcID::is_barrel(id)<<endl;

   for ( int theta=-5; theta<=10; theta++ )  {
     cout<<"PHI_MAX(theta="<<theta<<"):"<<EmcID::getPHI_ENDCAP_MAX(theta)<<endl;
   }

   cout<<"ENDCAP_EAST  ="<<EmcID::getENDCAP_EAST()<<endl;
   cout<<"ENDCAP_WEST  ="<<EmcID::getENDCAP_WEST()<<endl;

   cout<<"===============HltID Test================="<<endl;

   x=0;y=2;
   id=HltID::data_type_id(x,y);
   
   cout<<"Identifier id="<<id<<endl;

   cout<<"detector     ="<<HltID::detector(id)<<endl;
   cout<<"id_in_sub    ="<<HltID::id_sub(id)<<endl;
   cout<<"is ef_result? "<<HltID::is_ef_result(id)<<endl;
   cout<<"is eventtype? "<<HltID::is_eventtype(id)<<endl;
   cout<<"is energy?    "<<HltID::is_energy(id)<<endl;
   cout<<"is algorithm? "<<HltID::is_algorithm(id)<<endl;
   cout<<"is mdc inf?   "<<HltID::is_mdc_inf(id)<<endl;

   cout<<"EMC          ="<<HltID::EMC<<endl;
   cout<<"DETECTOR_MAX ="<<HltID::getDETECTOR_MAX()<<endl;
   cout<<"this id's max="<<HltID::id_sub_max(id)<<endl;
   cout<<"ID_MDC_MAX   ="<<HltID::getID_MDC_MAX()<<endl;

}
