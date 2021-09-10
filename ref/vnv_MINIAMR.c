///MINIAMR:13607191176472651629
/// This file was automatically generated using the VnV-Matcher executable. 
/// The matcher allows for automatic registration of all VnV plugins and injection 
/// points. Building the matcher requires Clang. If Clang is not available on this machine,
/// Registration code should be written manually. 
/// 

//PACKAGENAME: MINIAMR

#include "VnV.h" 
#include "version/version.h"


VNVVERSIONINFOCALLBACK(MINIAMR)

const char* getFullRegistrationJson_MINIAMR(){
	 return "{}";}

INJECTION_REGISTRATION(MINIAMR){
	REGISTER_FULL_JSON(MINIAMR, getFullRegistrationJson_MINIAMR);
};



