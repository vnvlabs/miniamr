///MINIAMR:1395143704222866919
/// This file was automatically generated using the VnV-Matcher executable. 
/// The matcher allows for automatic registration of all VnV plugins and injection 
/// points. Building the matcher requires Clang. If Clang is not available on this machine,
/// Registration code should be written manually. 
/// 

//PACKAGENAME: MINIAMR

#include "VnV.h" 
const char* getFullRegistrationJson_MINIAMR(){
	 return "{\"InjectionPoints\":{\"MY_FIRST_INJECTION_POINT\":{\"docs\":{\"description\":\"\",\"instructions\":\"\",\"params\":{},\"template\":\"\",\"title\":\"\"},\"name\":\"MY_FIRST_INJECTION_POINT\",\"packageName\":\"MINIAMR\",\"parameters\":{\"int main(int, char **)\":{\"argc\":\"int*\"}},\"stages\":{\"Begin\":{\"docs\":{\"description\":\"\",\"instructions\":\"\",\"params\":{},\"template\":\"\",\"title\":\"\"},\"info\":{\"Calling Function\":\"main\",\"Calling Function Column\":1,\"Calling Function Line\":43,\"filename\":\"main.c\",\"lineColumn\":5,\"lineNumber\":16}}}}}}";}

INJECTION_REGISTRATION(MINIAMR){
	Register_Injection_Point("MINIAMR","MY_FIRST_INJECTION_POINT","{\"int main(int, char **)\":{\"argc\":\"int*\"}}");
	REGISTER_FULL_JSON(MINIAMR, getFullRegistrationJson_MINIAMR);
}



