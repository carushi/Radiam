#include <string>
#include <iostream>
#include <cmath>
#include "Running_interface.h"
#define DEBUG
#define RADIAM

using namespace std;

std::string set(int count);

bool bind(char a, char b) 
{
	if (a == 'U') {
		if (b == 'A' || b == 'G') return true;
	} else if (a =='A' && b == 'U') return true;
	else if (a == 'G') {
        if (b == 'U' || b == 'C') return true;
	} else if (a == 'C' && b == 'G') return true;
	return false;
}

int main(int argc, char** argv)
{
    Rfold::Parameter::Init_ener_param();
    /*
    if (argc >= 2) {
        string str = argv[1];
        Rfold::Rfold_Lang model2;               
        model2.calculation(-1, str);
        model2.Write_bpp();
        return 0;
    }
    */
    for (int constraint = 10; constraint < 100; constraint += 30) {
        std::cerr << constraint << endl;        
         for (int count = 0; count < 1; count++) {
            std::string seq = set(count);
            //std::string seq = "CAUGGCCGAUUAAUAGCAAAGAAGGCCGAACCC";
            Rfold::Running_interface intf((constraint <= 0) ? (int)seq.length() : constraint, seq);
            //intf.Check_Mutation(seq);
            intf.Run_Radiam(seq, count == 0);
          	//intf.Run_BPP_Rfold_Model(seq);
            //intf.Raw_Compare_BPP_Rfold_Model(seq, count == 0);
        }
        exit(0);

    }
    return 0;
}
