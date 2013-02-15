#ifndef _RUNNING_INTERFACE_H
#define _RUNNINGINTERFACE_H

#include "Radiam.h"
#include <sstream>
#include <fstream>


namespace Rfold {

typedef vector<double> Vec;
typedef vector<Vec> Mat;


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
    void Run_Radiam(string, int, int = 0, int = 0);

};

}

#endif