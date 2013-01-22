#include "param.h"

namespace Rfold {
namespace Parameter {

double loghairpin[31];
double logtetra[30];
double logmismatchH[7][5][5];
double logmismatchI[7][5][5];
double logstack[7][7];
double logbulge[31];
double logTermAU;
double logint11[8][8][5][5];
double logint21[8][8][5][5][5];
double logint22[8][8][5][5][5][5];
double loginternal[31];
double logdangle5[8][5];
double logdangle3[8][5];
double logninio[MAXLOOP+1];
double logMLintern;    
double logMLclosing;    
bool initialized;




}
}