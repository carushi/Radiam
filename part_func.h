#ifndef _PART_FUNC_H
#define _PART_FUNC_H

#include <cstdio>
#include "matrix.h"
#include "param.h"
namespace Rfold {

class Rfold_Lang {
private:
public:
    Rfold_Lang() {}
    virtual ~Rfold_Lang() {}
    int _length;    // the length of original sequence
    int _constraint;
    Sequence seq;
    Matrix alpha;
    Matrix beta;
    static const bool noout = true;
    static const bool print = (false && !noout);
    static const bool debug = (false && print);
    double Calc_in_stem(int, int);
    double Calc_in_multiBif(int, int);
    double Calc_in_multi2(int, int);
    double Calc_in_multi1(int, int); 
    double Calc_in_multi(int, int);     
    double Calc_in_stemend(int, int);
    void Calc_in_outer(int);
    void Calc_inside_mat(int, int);    
    void Calc_inside();

    double Calc_out_stem(int, int);
    double Calc_out_multiBif(int, int);
    double Calc_out_multi2(int, int);    
    double Calc_out_multi1(int, int); 
    double Calc_out_multi(int, int);     
    double Calc_out_stemend(int, int);
    void Calc_out_outer();
    void Calc_outside_mat(int, int);
    void Calc_outside();

    double hairpin_acc(int, int);
    double interior_acc(int, int);
    double multi_acc(int, int);    
    double multi2_acc(int, int);
    void Calc_acc();
    double bpp(int, int);
    //void Write_bpp(Vec&);
    void Write_bpp(Mat&);    
    void Write_bpp();

    void Print_Mat(bool);
    void Initialize();
    void Set_sequence(const string&);
    void calculation(int, string);
    bool Is_range(int i, int j) { return (i >= j-_constraint-1); }
    void Set_Constraint(int constraint, int length) 
    {
        _length = length;
        _constraint = (constraint > 0) ? min(constraint, _length-1) : _length-1;
    }

};


}
#endif
