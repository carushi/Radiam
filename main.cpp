#include <string>
#include <iostream>
#include <cmath>
#include <getopt.h>
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

struct option* option()
{
    struct option *options = new struct option[5];
    options[0].name = "constraint";
    options[0].has_arg = required_argument;
    options[0].flag = NULL; options[0].val = 1;
    options[1].name = "window";
    options[1].has_arg = required_argument;
    options[1].flag = NULL; options[1].val = 1;
    options[2].name = "threshold";
    options[2].has_arg = required_argument;
    options[2].flag = NULL; options[2].val = 1;
    options[3].name = 0; options[3].has_arg = 0;
    options[3].flag = 0; options[3].val = 0;
    return options;
}

void Run_interface(string seq, int mtype, int constraint, int window, int threshold)
{
    Rfold::Running_interface intf((constraint <= 0) ? (int)seq.length() : constraint, seq);
    intf.Run_Radiam(seq, mtype, threshold, window);
}

int main(int argc, char** argv)
{
    int window = 0;
    int mtype = 2;
    int threshold = DEF_PRE;
    int constraint = 50;    
    bool const_flag = false;    
    int option_index;    
    string seq;
    struct option *options = option();
    Rfold::Parameter::Init_ener_param();    
    while (1) {
        int opt = getopt_long(argc, argv, "dims:", options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 0: case 1:
                if (option_index == 0 && optarg) {
                    constraint = atoi(optarg); const_flag = true;
                } else if (option_index == 1 && optarg) window = atoi(optarg);
                else if (option_index == 2 && optarg) threshold = atoi(optarg);
                break;
            case 'd': mtype = 1; break;
            case 'i': mtype = 0; break;
            case 'm': mtype = 2; break;
            case 's': if (optarg != NULL) seq = string(optarg); break;
            default: break;
        }
    }
    if (const_flag && (int)seq.length() > 0) Run_interface(seq, mtype, constraint, window, threshold);
    else if (const_flag) {
        for (int count = 0; count < 5; count++) {
            seq = set(count);
            Run_interface(seq, mtype, constraint, window, threshold);
        }
    } else {
        for ( ; constraint <= 400; constraint += 50) {
            std::cerr << constraint << endl;
            if ((int)seq.length() > 0) Run_interface(seq, mtype, constraint, window, threshold);
            else {
                for (int count = 0; count < 5; count++) {
                    seq = set(count);
                    Run_interface(seq, mtype, constraint, window, threshold);
                    //intf.Check_Mutation(seq);
                    //intf.Run_BPP_Rfold_Model(seq);
                   //intf.Raw_Compare_BPP_Rfold_Model(seq, count == 0);
                }
            }
            break;
        }
    }
    delete[] options;
    return 0;
}
