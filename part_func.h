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
    DOUBLE Calc_in_stem(int, int);
    DOUBLE Calc_in_multiBif(int, int);
    DOUBLE Calc_in_multi2(int, int);
    DOUBLE Calc_in_multi1(int, int); 
    DOUBLE Calc_in_multi(int, int);     
    DOUBLE Calc_in_stemend(int, int);
    void Calc_in_outer(int);
    void Calc_inside_mat(int, int);    
    void Calc_inside();

    DOUBLE Calc_out_stem(int, int);
    DOUBLE Calc_out_multiBif(int, int);
    DOUBLE Calc_out_multi2(int, int);    
    DOUBLE Calc_out_multi1(int, int); 
    DOUBLE Calc_out_multi(int, int);     
    DOUBLE Calc_out_stemend(int, int);
    void Calc_one_out_outer(int);    
    void Calc_out_outer();
    void Calc_outside_mat(int, int);
    void Calc_outside();

    DOUBLE hairpin_acc(int, int);
    DOUBLE interior_acc(int, int);
    DOUBLE multi_acc(int, int);    
    DOUBLE multi2_acc(int, int);
    void Calc_acc();
    DOUBLE bpp(int, int, bool deb = false);
    DOUBLE acc(int, int);
    void Write_acc(Mat&);
    void Write_bpp(Mat&);    
    void Write_stem(Vec&, const Mat&);
    void Write_stem(Vec&);
    void Write_stem(bool);
    void Write_bpp(bool);    
    void Write_bpp_part(bool, int, int);

    void Print_Mat(bool);
    void Print_dout(bool);
    void Initialize();
    void Set_sequence(const string&);
    void calculation(int, string, bool shrink = true);
    bool Is_range(int i, int j) { return (i >= j-_constraint-1); }
    void Set_Constraint(int constraint, int length) 
    {
        _length = length;
        _constraint = (constraint > 0) ? min(constraint, _length-1) : _length-1;
    }
    void Set_raw_Constraint(int constraint, int length) 
    {
        _length = length;
        _constraint = constraint;
    }
};


}
#endif
