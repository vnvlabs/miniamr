//10294498724326657805
/// This file was automatically generated using the VnV-Matcher executable. 
/// The matcher allows for automatic registration of all VnV plugins and injection 
/// points. Building the matcher requires Clang. If Clang is not available on this machine,
/// Registration code should be written manually. 
/// 

//PACKAGENAME: MINIAMR

#include "VnV.h" 
const char* getFullRegistrationJson_MINIAMR(){
	 return "{\"Conclusion\":{\"configuration\":{},\"description\":\"\",\"instructions\":\"\",\"params\":{},\"shortTitle\":\"\",\"template\":\"\",\"title\":\"\"},\"Executables\":{\"configuration\":{},\"default\":{},\"description\":\"\",\"instructions\":\"\",\"lib\":\"executables\",\"params\":{},\"shortTitle\":\"\",\"template\":\"\",\"title\":\"\"},\"InjectionPoints\":{\"driver\":{\"docs\":{\"configuration\":{},\"description\":\"\",\"instructions\":\"\",\"params\":{},\"shortTitle\":\"\",\"template\":\"\",\"title\":\"\"},\"loop\":true,\"name\":\"driver\",\"packageName\":\"MINIAMR\",\"parameters\":{\"int main(int, char **)\":{\"init_x\":\"int*\",\"init_y\":\"int*\",\"init_z\":\"int*\"}},\"stages\":{\"Begin\":{\"docs\":{\"configuration\":{},\"description\":\"\",\"instructions\":\"\",\"params\":{},\"shortTitle\":\"\",\"template\":\"\",\"title\":\"\"},\"info\":{\"Calling Function\":\"main\",\"Calling Function Column\":1,\"Calling Function Line\":43,\"filename\":\"main.c\",\"lineColumn\":5,\"lineNumber\":18}},\"End\":{\"info\":{\"Calling Function\":\"main\",\"Calling Function Column\":1,\"Calling Function Line\":43,\"filename\":\"main.c\",\"lineColumn\":5,\"lineNumber\":22}}}}},\"Introduction\":{\"configuration\":{},\"description\":\"\",\"instructions\":\"\",\"params\":{},\"shortTitle\":\"\",\"template\":\"\",\"title\":\"\"}}";}

INJECTION_REGISTRATION(MINIAMR){
	Register_Injection_Point("MINIAMR","driver","{\"int main(int, char **)\":{\"init_x\":\"int*\",\"init_y\":\"int*\",\"init_z\":\"int*\"}}");
	REGISTER_FULL_JSON(MINIAMR, getFullRegistrationJson_MINIAMR);
}



FORTRAN_INJECTION_REGISTRATION(MINIAMR)

