#ifndef _RADIAM_H
#define _RADIAM_H
#define DEF_PRE 4
#include "part_func.h"
#include "param.h"

namespace Rfold {

typedef vector<double> Vec;
typedef vector<Vec> Mat;

using std::transform;
class Wobble
{
    friend class Radiam;
public:
    double min;
    double max;
    Vec a;
    Vec b;
    Vec imax;
    Vec imin;    
    Wobble() {}
    Wobble(int size) {
        a = Vec(size, 0.0); b = Vec(size, 0.0);
        imax = Vec(size, -INF); imin = Vec(size, INF);
    }
    ~Wobble() {}
};

class Radiam : public Rfold_Lang {
private:
	vector<int> _mpoint;
    vector<int> _mlist;
    vector<int> _index;
    vector<int> _right_limit;  // not constant limit;
	vector<int> _left_limit;   // not constant limit;
    Vec _constant;
    Vec bpp;
    Mat bppm;
    const int _precision;      // # of digits to keep correct;
    static const bool _omit = false;    
    static const bool _matrix = false;    
    static const bool rdebug = false;
    static const char* base;
    Mat& Get_inner(bool, int, bool);
    void Add_constant(int, int, double, bool);
    void Add_constant_inner(int, int, double);
    bool Set_constant_in(int, int, int, int);
    bool Set_constant_out(int, int, int, int);
    bool Set_constant_inner(int, int, int, int);
    double Get_min_inner(int, int);
    void Copy_const_inner(int, int, Matrix&, double);
    void Calc_inside();
    void Calc_out_outer();
    void Calc_outside_inner(int, int);    
    void Add_outside_inner(int, int&);
    void Calc_outside();

    void Initialize_seq(string&);
    void Calc_matrix(int, string&);
    void Copy_matrix(int);
    void Set_limit();
    void Set_mpoint(int);
    void Set_index();
    void Change_sequence(int, int, int, string);
    void All_calculation(int, int, string&);

    double Calc_bpp_cor(const Vec&, int);

public:
    int Mtype;
	Matrix ori_alpha;
	Matrix ori_beta;
    Wobble _wob;
    enum Type { In, Del, Mut, Stem, Stemend, Multi, Multi1, Multi2, Multibif };
    Radiam(int precision = DEF_PRE) : Rfold_Lang(), _precision(precision) {}
	virtual ~Radiam(){}
	void Mutation_calculation(int, string&);
    void Get_ori_matrix(const string&);
    void Correlation_of_bpp(int, int, int, string);
    /////////
    void Debug_confirm(int, string&);
    void Debug_output(int, int, bool, Rfold_Lang&);
    void compare(int, const Vec&, const Vec&);
    void compare(int, const Mat&, const Mat&);
    bool compare_same(int, const Vec&, const Vec&);    
    bool compare_same(int, const Mat&, const Mat&);
    int Check_Difference(const class Matrix&, const class Matrix&);     
    void Output_Difference(int, const class Matrix& ori, const class Matrix& mut);
    ////////
    bool Under_Prec(double max, double min, double value) {
        return (log10(fabs(max-min))-log10(fabs(value)) < -(_precision)-3);
    }
    int In_range(int j, int start) {
        if (Mtype == Mut) return j-(start+_constraint+2);
        else if (Mtype == Del) return j-(start+_constraint+2);
        else return j-(start+_constraint+3); /////////
    }
    int Out_range(int j, int start) {
        if (Mtype == Mut) return j-(start-_constraint-2);
        else return j-(start-_constraint-3);
    }
    bool Is_out_range(int j, int const_end, int mp) {
        if (Mtype == Mut) return (j >= const_end && Out_range(j, _mpoint[mp]) <= 0);
        else return (j >= const_end && Out_range(j, _mpoint[mp]) <= 0);        
    }
    bool Out_in_range(int i, int j, int mp) {
        if (Mtype == Mut) return (j > _mpoint[mp]+1 && i < _mpoint[mp]-1);
        else return (j > _mpoint[mp]+_constraint+1 && i < _mpoint[mp]-_constraint-1);
    }
};

}


#endif