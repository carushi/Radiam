#ifndef _RADIAM_H
#define _RADIAM_H
#define DEF_PRE (4)
#define BAND (1)
#include <limits>
#include <ctime>
#include <cassert>
#include <sys/time.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "part_func.h"
#include "param.h"

namespace Rfold {

typedef vector<double> Vec;
typedef vector<Vec> Mat;

using std::transform;
using std::pair;

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
    vector<int> _inner_limit;
    Vec _constant;
    vector<pair<int, double> > _window_max;
    vector<pair<char, double> > _position_max;
    Mat bppm;
    bool _init;
    const int _precision;      // # of digits to keep correct;        
    static const char* base;
    Mat& Get_inner(bool, int, bool);
    void Add_constant(int, int, double, bool);
    void Add_constant_inner(int, int, double);
    bool Set_constant_in(int, int, int, int);
    bool Set_constant_out(int, int, int, int);
    bool Set_constant_inner(int, int, int, int);
    double Get_min_inner(int, int);
    void Copy_const_inner(int, int, Matrix&, double);
    //
    void Calc_inside();
    void Calc_out_outer();
    void Calc_outside_inner(int, int);
    double Calc_partition_function();
    void Calc_min_inside();
    void Calc_min_out_outer();
    void Calc_min_outside(double);
    void Add_outside_inner(int, int&);
    void Calc_outside();
    double Calc_bpp_cor(const Vec&, const Vec&);
    void Get_cor_vec(int, int, Vec&, Vec&);    
    void Get_cor_vec(int, int, Vec&, Vec&, const Mat&);
    void Calc_one_bpp_cor();
    void Calc_bpp_cor();
    //
    void Initialize_seq(string&);
    void Copy_matrix(int);
    void Set_limit();
    void Set_mpoint(int);
    void Set_index();
    void Print_mlist(int, string&);
    void Change_sequence(int, int, int, string);    
    void Calc_matrix(int, string&);    
    void All_calculation(int, int, string&);
    void Part_calculation(int, int, string&);
    void Output_common_data();
    void Output_delta_ene_line(double, int, int);
    void Output_delta_ene_line(double, double, int);
    int Calc_bpp_ene(Vec&);
    int Store_stem(Vec&, Vec&);
    int Store_bpp(Vec&, Vec&);
    void Out_header(int, ofstream&);
    void Out_value(double, double, int, ofstream&);
    void Output_delta_ene();
    double Euclidean_Distance(const Vec&, const Vec&);
    void Output_ED_CC();
    void Output_acc_ene();
    void Output_correlation(const double);    
    void Output_correlation(const Vec&);
    void Output_storage(const string&);
    void Storage_max(const Vec&);
    void Set_Output_file(int);
    void Set_Correlation(int, int, string&, bool);
    //
    void Calc_time(int, string&);
    void Calc_debug(int, string&);    
    double Get_diff(double, double);
    double Get_diff_at_same(int, int, const Mat&);
    bool Within_range(int, const Mat&);
    bool Within_range(int, int, const Mat&);
    void Write_bppm_fluc(const Mat&);
    void Write_bppm_fluc_abs(const Mat&);    
    void Write_accm_fluc();
    void Write_accm_fluc_abs(int);
    void Write_accm_fluc_abs();
    void Write_bppm_dif(const Mat&, const Mat&);
    void Debug_bppm(int, string&);
    void Debug_bppm(int, string&, int, int);
    void Debug_confirm(int, string&);
    void Debug_output(int, int, bool, Rfold_Lang&);
    void compare(int, const Vec&, const Vec&);
    void compare(int, const Mat&, const Mat&);
    bool compare_same(int, const Vec&, const Vec&);    
    bool compare_same(int, const Mat&, const Mat&);
    int Check_Difference(const class Matrix&, const class Matrix&);     
    void Output_Difference(int, const class Matrix& ori, const class Matrix& mut);

public:
    int Mtype;
    int window, rightw, leftw;
    vector<string> output_char;
	Matrix ori_alpha;
	Matrix ori_beta;
    Wobble _wob;
    string outf;
    int id;
    string ori_seq;
    string bp_seq;
    vector<int> cons_index;
    enum Type { In, Del, Mut, Stem, Stemend, Multi, Multi1, Multi2, Multibif };
    enum Out { Abs, Rel, Acc, RelAcc };
    enum Distance { StemDif, CC };
    enum Debug { Reg, Dif, Time, Bpp, Outer, Analyze };
    enum Calc { No, Conv, Omit, OmitOut };
    enum Cor { One, All, Spe };
    int _cor;
    int _output;
    int _debug;
    int _dist;
    static const int _calc = Calc::No;
    static const bool rdebug = false;
    Radiam(bool acc = false, int precision = DEF_PRE, int window = 0, bool init = false, bool dist = false) 
    : Rfold_Lang(), _precision(precision), window(window), _init(init)
    {
        _output = (acc) ? Out::Acc : Out::Abs;
        _dist = (dist) ? Distance::CC : Distance::StemDif;
        _cor = Cor::One;
        _debug = Debug::Bpp;
    }
    Radiam(string& filename, bool acc, int precision, int window, bool init, bool dist = false)
    : Rfold_Lang(), _precision(precision), window(window), _init(init)
    {
        if (filename != "") outf = filename;
        _output = (acc) ? Out::Acc : Out::Abs;
        _dist = (dist) ? Distance::CC : Distance::StemDif;
        _cor = Cor::One;
        _debug = Debug::Bpp; 
    }
	virtual ~Radiam(){}
	void Mutation_calculation(int, string&);
    void Get_ori_matrix(const string&);
    void Set_example();
    void Set_annotate(string&, string&, vector<int>&);
    void Initialize_before_calc(int, int, int&, string&);
    void Correlation_of_bpp(int, vector<int>&, int, string);
    void Correlation_of_bpp(int, int, int, string);
    void Write_bpp_binary();
    void Set_bpp_binary(int, int, string&);

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
    void Set_region(int l, int r) {
        if (l > r) swap(l, r);
        leftw = l; rightw = r;
        _cor = Cor::Spe;
    }
    double eSDC(double val) {
        return sqrt(seq.length)*(1.0-val);
    }
};

}


#endif