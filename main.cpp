#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include "file_explorer.h"
#include "Running_interface.h"
#define DEBUG
#define OPTION (8)


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
    options[6].name = "help";
    options[6].has_arg = 0;
    options[6].flag = NULL; options[6].val = 1;
    options[7].name = "snv";
    options[7].has_arg = 0;
    options[7].flag = NULL; options[6].val = 1;
    options[OPTION].name = 0; options[OPTION].has_arg = 0;
    options[OPTION].flag = 0; options[OPTION].val = 0;
    return options;
}

void Stem_probability(Rfold::Arg& arg)
{
    Rfold::Rfold_Lang model;
    Rfold::Running_interface::RNA_transform(arg.seq);
    cout << "#seq " << arg.seq << endl;
    model.calculation(arg.constraint, arg.seq, false);
    /*
    model.Write_stem(arg.acc_flag);
    model.Print_Mat(true);
    model.Print_Mat(false);
    */
    model.Print_dout(true);
    cout << endl;
    model.Print_dout(false);
}

string Get_filename(string& original, int array)
{
    ostringstream oss;
    oss << original << "_" << array << ".txt";
    return oss.str();
}

void Read_raw_file(Rfold::Arg& arg, const string& filename)
{
    string str;
    ifstream ifs(filename.c_str());
    int constraint = (arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint;
    Rfold::Running_interface intf(constraint, arg.seq);
    if (arg.filename == "") arg.filename = "output.txt";
    if (arg.array > 0) arg.filename = Get_filename(arg.filename, arg.array);
    for (int i = 0; getline(ifs, str); i++) {
        if (arg.array > 0 && i%ARRAY != arg.array-1) continue;
        arg.seq = str;
        if (arg.normal) Stem_probability(arg);
        else intf.Run_Radiam(arg, i/ARRAY == 0, i);
    }
}

void Run_interface(Rfold::Arg& arg, bool init = false)
{
    if (arg.filename == "") arg.filename = "output.txt";    
    if (arg.normal) {
        Stem_probability(arg);
        return;
    }
    Rfold::Running_interface intf((arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint, arg.seq);
    arg.get_range();
    if (arg.example) {
        intf.Run_BPP_Rfold_Model(arg.seq); //cut sequence test;
    } else if (arg.range != "") 
        intf.Run_Radiam(arg, init, arg.leftw, arg.rightw);
    else {
        intf.Run_Rfold_Model(arg.seq, arg.acc_flag);        
        //intf.Run_Radiam(arg, init);
        //intf.Check_Mutation(seq);
        //intf.Raw_Compare_BPP_Rfold_Model(seq, count == 0);
    }
}

void Print_option_help() 
{
    cout << "-d\tsimulate deletion\n"
         << "-i\tsimulate insertion\n"
         << "-m\tsimulate substitution\n"
         << "-a\tuse accessibility (default: base pairing probability)"
         << "-k [num]\tthe number of mutation\n"
         << "-e\tfor debug...\n"
         << "-s [sequence]\tdirect input of sequence\n"
         << "-l\tsimulate only long sequences\n"
         << "-r [start-end]\tcalculate correlation coefficient within range\n"
         << "-p [num]\tparallel computing (max 15)\n"
         << "-c\tprint correlation coefficient of distance (eucledian, eSDC)\n" 
         << "-b\tprint stem probability if sequence or import file is defined\n" 
         << "-o\tsimulate only one sample sequence\n"
         << "--help\tprint these sentences\n"
         << "--snv\tsnv analysis"
         << "--constraint [num]\tset the maximal span num\n"
         << "--window [num]\tset the window size for correlation coefficient num\n"         
         << "--threshold [num]\tset the precision num\n"         
         << "--output [filename]\toutput to filename\n"
         << "--import [filename]\tdeal filename as Rfam input\n"
         << "--importseq [filename]\tdeal filename as sequence file input\n"
         << endl;
}

int main(int argc, char** argv)
{
    bool once_flag = false;
    bool import_flag = false;
    bool rawseq_import_flag = false;
    bool help_flag = false;
    bool snv_flag = false;
    int option_index;
    struct option *options = option();
    string importfile;
    Rfold::Arg arg;
    Rfold::Parameter::Init_ener_param();    
    while (1) {
        int opt = getopt_long(argc, argv, "bcledimaos:r:p:k:", options, &option_index);
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
                } else if (option_index == 6) {
                    help_flag = true;
                } else if (option_index == 7) {
                    snv_flag = true;
                }
                break;
            case 'd': arg.mtype = 1; break;
            case 'i': arg.mtype = 0; break;
            case 'm': arg.mtype = 2; break;
            case 'e': arg.example = true; break;            
            case 's': if (optarg != NULL) arg.seq = string(optarg); break;
            case 'a': arg.acc_flag = true; break;
            case 'r': if (optarg != NULL) arg.range = string(optarg); break;
            case 'p': if (optarg != NULL) arg.array = atoi(optarg); break;
            case 'l': arg.longer = true; break;
            case 'c': arg.dist = true; break;
            case 'k': if (optarg != NULL) arg.num = atoi(optarg); break;
            case 'b': arg.normal = true; break;
            case 'o': once_flag = true; break;
            default: break;
        }
    }
    if (help_flag) {
        Print_option_help();
    } else if (snv_flag) {
        Transcript::Read_snp_file(importfile, arg);
    } else if (import_flag) { 
        Rfam::Read_Rfam_file(arg, importfile, arg.array);
    } else if (rawseq_import_flag) {
        Read_raw_file(arg, importfile);
    } else if (arg.const_flag && (int)(arg.seq.length()) > 0) 
        Run_interface(arg, true);
    else {
        for (int count = 0; count < 1000; count++) {
            if (count == 0)
                arg.seq = "atgttgtacctggaaaacaatgcccagtcccagtatagcgagccgcagtacacgaacctggggctcctgaacagcatggaccagcaggttcagaatggctcttcctccaccagcccctacaacacggagcacgcgcagaacagcgtcacggccccctcgccttacgcccagcccagctccacttttgatgccctctcgccctccccagccatcccttccaacacagactacccgggacctcacagcttcgacgtatcatttcaacaatccagcacagcaaagtctgcaacgtggacgtattccactgaactgaagaagctgtactgccagattgccaagacatgccccattcagatcaaagtgatgaccccaccaccccagggagctgtcatccgggctatgccagtctacaaaaaagcagggcacgtcaccgaagtggtcaaacgctgcccgaaccacgagctgagccgggagttcaatgaggggcagattgcacctcctagccacctgatcagagtggaaggaaacagccatgcccagtatgtggaagaccccatcactgggagacagagcgtgctggtcccatatgagccaccccaggttggtaccgagttcacaacagtcctgtacaacttcatgtgtaacagcagctgtgtaggagggatgaaccgtcgcccaattctcatcattgttacactggaaaccagagatgggcaagtcttgggccgccgatgttttgaagctcgcatttgcgcttgcccaggcagagatcgcaaagcagatgaggacagcatccgcaagcagcaagtctctgacagcacaaagaatggtgatgcttttcggcaaggaactcatggcatacagatgacatctatcaagaaaagacgttctccagatgatgagctcttgtacttgccggtgaggggacgagaaacatatgaaatgctactgaagatcaaagagtccctggaacctatgcagtaccttccccagcacacaattgagacttaccggcagcagcagcaacagcagcaccagcacttgctccagaagcagacctccattcagtcacagtcatcctatggctccaactcaccgccgctcagcaagatgaacagcatgaacaagctgccctcggtcagccagctcataaacccccagcagcgcaacgcactgaccccaaccaccatccctgacggcatgggaacaaacattcccatgatgggcactcacatggccatgaccggcgacatgaatgtcctcagccccacgcaggcgctgcctcctcccctctccatgccttcaacgtcccactgcactcctcctcctccataccccacagactgcagcattgtcagcttcttagcgaggttgggctgctcatcctgtgtggattatttcacgacccaagggctgaccaccatctatcatattgagcattactccatggatgatctggtgagcctcaagatcccggagcagttccgccacgccatctggaagggcatcctggaccaccggcagctccatgacttctcctctcctccccacctcctgcgtacccccagcggtgcctccaccgtcagcgtgggctccagcgaaacccggggggagcgggtcatcgatgcagtccgcttcactctccgccagaccatttccttcccgccccgcgacgagtggaacgatttcaacttcgacatggatgcccgccgcaacaaacagcaacgcatcaaggaggaaggggagtga";
            else
                arg.seq = set(count);
            Run_interface(arg, count != 0);
            if (once_flag) break;
        }
    }
    delete[] options;
    return 0;
}
