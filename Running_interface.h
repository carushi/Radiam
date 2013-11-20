#ifndef _RUNNING_INTERFACE_H
#define _RUNNING_INTERFACE_H

#include "Radiam.h"
#include <sstream>
#include <fstream>
#define ARRAY (100)
#define TOOMUCH (700)
#define TOOLONG (500)


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
    int num;
    int array;
    bool acc_flag;
    bool const_flag;
    bool init;
    bool example;
    bool longer;
    bool dist;
    bool normal;
    Arg() {
        window = 0;
        mtype = 2;
        num = 1;
        threshold = DEF_PRE;
        array = 0;
        acc_flag = false;
        example = false;
        longer = false;
        dist = false;
        normal = false;
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
    void Set_consensus_index(Radiam&, string);    
    void Set_consensus_index(Radiam&, string, string);
public:
    Running_interface() {}
    Running_interface(int Constraint, string Sequence) : constraint(Constraint), sequence(Sequence) {
        //Set_Files();
    }
    ~Running_interface() {}
    static void RNA_transform(string&);
    /*
    void Check_Mutation(string);
    void Raw_compare_BPP_Rfold_Model(string, bool);
    */
    void Run_BPP_Rfold_Model(string);
    void Run_Rfold_Model(string, bool);
    void Run_Radiam(Arg&, vector<string>&, vector<string>&, bool);
    void Run_Radiam(Arg& arg, vector<string>&, vector<string>&, string&, bool);    
    void Run_Radiam(Arg&, bool, int, int);
    void Run_Radiam(Arg&, bool, int = 0);
    void Run_Radiam(Arg&, bool, int, vector<bool>&);    
    void Run_Radiam_for_mRNA(Arg&, string&, vector<bool>&);    

};

}

#endif