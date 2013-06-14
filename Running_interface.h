#ifndef _RUNNING_INTERFACE_H
#define _RUNNINGINTERFACE_H

#include "Radiam.h"
#include <sstream>
#include <fstream>


namespace Rfold {

typedef vector<double> Vec;
typedef vector<Vec> Mat;


class Arg
{
public:
    string seq;
    string range;
    string filename;
    int mtype;
    int constraint;
    int window;
    int threshold;
    int rightw;
    int leftw;
    bool acc_flag;
    bool const_flag;
    bool init;
    bool example;
    bool longer;
    bool dist;
    Arg() {
        window = 0;
        mtype = 2;
        threshold = DEF_PRE;
        acc_flag = false;
        example = false;
        longer = false;
        dist = false;
        constraint = 50;
    }
    void get_range() {
        if (range != "") {
            string::size_type ind = range.find_first_of('-');
            if (ind == string::npos) return;
            leftw = atoi(range.substr(0, ind).c_str());
            rightw = atoi(range.substr(ind).c_str()+1);
            cout << leftw << " " << rightw << endl;
        }
    }
};

class Running_interface {
private:
    int constraint;
    string sequence;
    static const int MAXTYPE = 7;
    static string base;
    vector<string> files;
    vector<string> filenames;
    void Set_Files();
    void Init_Files();
    void Same_Header(const char* header, int, char, char, const string&);    
    void compare(int, const Vec&, const Vec&);
    void compare(int, const Mat&, const Mat&);
    bool compare_same(int, const Vec&, const Vec&);    
    bool compare_same(int, const Mat&, const Mat&);
    bool Check_Difference(const class Matrix&, const class Matrix&);     
    void Output_Difference(const class Matrix& ori, const class Matrix& mut); 
    void Set_consensus_index(Radiam&, string, string, vector<int>&);
public:
    Running_interface() {}
    Running_interface(int Constraint, string Sequence) : constraint(Constraint), sequence(Sequence) {
        Set_Files();
    }
    ~Running_interface() {}
    void RNA_transform(string&);
    void Check_Mutation(string);
    void Raw_compare_BPP_Rfold_Model(string, bool);
    void Run_BPP_Rfold_Model(string);
    void Run_Radiam(Arg&, vector<string>&, vector<string>&, bool);
    void Run_Radiam(Arg& arg, vector<string>&, vector<string>&, string&, bool);    
    void Run_Radiam(Arg&, bool, int, int);
    void Run_Radiam(Arg&, bool);

};

}

#endif