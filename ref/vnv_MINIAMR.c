///7070386492871253011
/// This file was automatically generated using the VnV-Matcher executable. 
/// The matcher allows for automatic registration of all VnV plugins and injection 
/// points. Building the matcher requires Clang. If Clang is not available on this machine,
/// Registration code should be written manually. 
/// 

//PACKAGENAME: MINIAMR

#include "VnV.h" 
const char* getFullRegistrationJson_MINIAMR(){
	 return "{\"Communicator\":{\"docs\":\"\",\"name\":\"mpi\",\"package\":\"VnV\"},\"Introduction\":\"\"}";}

INJECTION_REGISTRATION(MINIAMR){
	VnV_Declare_Communicator("MINIAMR","VnV","mpi");
	REGISTER_FULL_JSON(MINIAMR, getFullRegistrationJson_MINIAMR);
};



