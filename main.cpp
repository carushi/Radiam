#include <string>
#include <iostream>
#include <cmath>
#include <getopt.h>
#include "Running_interface.h"
#define DEBUG
#define OPTION (6)
#define ARRAY (15)

using namespace std;


string set(int count);
void Read_Rfam_file(Rfold::Arg& arg, string& filename, int array);

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
    struct option *options = new struct option[OPTION+1];
    options[0].name = "constraint";
    options[0].has_arg = required_argument;
    options[0].flag = NULL; options[0].val = 1;
    options[1].name = "window";
    options[1].has_arg = required_argument;
    options[1].flag = NULL; options[1].val = 1;
    options[2].name = "threshold";
    options[2].has_arg = required_argument;
    options[2].flag = NULL; options[2].val = 1;
    options[3].name = "output";
    options[3].has_arg = required_argument;
    options[3].flag = NULL; options[3].val = 1;
    options[4].name = "import";
    options[4].has_arg = required_argument;
    options[4].flag = NULL; options[4].val = 1;
    options[5].name = "importseq";
    options[5].has_arg = required_argument;
    options[5].flag = NULL; options[5].val = 1;
    options[OPTION].name = 0; options[OPTION].has_arg = 0;
    options[OPTION].flag = 0; options[OPTION].val = 0;
    return options;
}

void Read_raw_file(Rfold::Arg& arg, const string& filename)
{
    string str;
    ifstream ifs(filename.c_str());
    int constraint = (arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint;
    Rfold::Running_interface intf(constraint, arg.seq);
    if (arg.filename == "") arg.filename = "output.txt";
    for (int i = 0; getline(ifs, str); i++) {
        arg.seq = str;
        intf.Run_Radiam(arg, i == 0);
    }
}

void Run_interface(Rfold::Arg& arg, bool init = false)
{
    if (arg.filename == "") arg.filename = "output.txt";    
    Rfold::Running_interface intf((arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint, arg.seq);
    arg.get_range();
    if (arg.example) {
        intf.Run_BPP_Rfold_Model(arg.seq); //cut sequence test;
    } else if (arg.range != "") 
        intf.Run_Radiam(arg, init, arg.leftw, arg.rightw);
    else
        intf.Run_Radiam(arg, init);
    //intf.Check_Mutation(seq);
    //intf.Raw_Compare_BPP_Rfold_Model(seq, count == 0);
}

int main(int argc, char** argv)
{
    bool import_flag = false;
    bool rawseq_import_flag = false;
    int option_index;
    int array = 0;   
    struct option *options = option();
    string importfile;
    Rfold::Arg arg;
    Rfold::Parameter::Init_ener_param();    
    while (1) {
        int opt = getopt_long(argc, argv, "cledimas:r:p:", options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 0: case 1:
                if (option_index == 0 && optarg) {
                    arg.constraint = atoi(optarg); arg.const_flag = true;
                } else if (option_index == 1 && optarg) arg.window = atoi(optarg);
                else if (option_index == 2 && optarg) arg.threshold = atoi(optarg);
                else if (option_index == 3 && optarg) arg.filename = string(optarg);
                else if (option_index == 4 && optarg) {
                    importfile = string(optarg); import_flag = true;
                } else if (option_index == 5 && optarg) {
                    importfile = string(optarg); rawseq_import_flag = true;
                } break;
            case 'd': arg.mtype = 1; break;
            case 'i': arg.mtype = 0; break;
            case 'm': arg.mtype = 2; break;
            case 'e': arg.example = true; break;            
            case 's': if (optarg != NULL) arg.seq = string(optarg); break;
            case 'a': arg.acc_flag = true; break;
            case 'r': if (optarg != NULL) arg.range = string(optarg); break;
            case 'p': if (optarg != NULL) array = atoi(optarg); break;
            case 'l': arg.longer = true;
            case 'c': arg.dist = true;
            default: break;
        }
    }
    if (import_flag) {
        Read_Rfam_file(arg, importfile, array);
    } else if (rawseq_import_flag) {
        Read_raw_file(arg, importfile);
    } else if (arg.const_flag && (int)(arg.seq.length()) > 0) 
        Run_interface(arg, true);
    else {
        for (int count = 0; count < 1000; count++) {
            arg.seq = set(count);
            Run_interface(arg, count != 0);
        }
    }
    delete[] options;
    return 0;
}
