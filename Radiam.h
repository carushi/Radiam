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
        imax = Vec(size, 0.0); imin = Vec(size, 0.0);
    }
    ~Wobble() {}
};

class Radiam : public Rfold_Lang {
private:
	vector<int> _mpoint;
    vector<int> _mlist;
    vector<int> _right_limit;  // not constant limit;
	vector<int> _left_limit;   // not constant limit;
    Vec _constant;
    Vec bpp;
    Mat bppm;
    const int _precision;      // # of digits to keep correct;
    //static const bool _reduced = false;
    static const bool rdebug = true;
    static const bool _matrix = false;    
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
    void Calc_outside();

    void Initialize_seq(string&);
    void Calc_matrix(int, string&);
    void Copy_matrix(int);
    void Set_limit();
    void Set_mpoint(int);
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
    bool Under_Prec(double range, double value) {
        return (log10(abs(range))-log10(value) < -(_precision)-1);
    }
    bool Is_Range(int i, int j, int mp) {
        return (j > _mpoint[mp]+1 && i < _mpoint[mp]);
    }
    int index(int j, int mp, bool inside) {
        if (!inside) mp--;
        if (Mtype == Mut) return j;
        else if (Mtype == In) return j-mp;
        else return j+mp;
    }
};

}


#endif