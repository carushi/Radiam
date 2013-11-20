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

class Wobble // for threshold cut;
{
    friend class Radiam;
public:
    double min;
    double max;
    Vec a;      // for threshold cut of inside outer calculation;
    Vec b;      // for threshold cut of outside outer calculation;
    Vec imax;   // for threshold cut of outside inner calculation;
    Vec imin;   // for threshold cut of outside inner calculation;
    Wobble() {}
    Wobble(int size) {
        a = Vec(size, 0.0); b = Vec(size, 0.0);
        imax = Vec(size, -INF); imin = Vec(size, INF);
    }
    ~Wobble() {}
};

class Radiam : public Rfold_Lang {
private:
	vector<int> _mpoint;        // mutation point to be calculated on matrix;
    vector<int> _mlist;         // mutation point on original sequence;
    vector<int> _index;         // position index of mutated seuqence on original sequence;
    vector<int> _right_limit;   // inside outer limit, which are not constant;
	vector<int> _left_limit;    // outside outer limit, which are not constant;
    vector<int> _inner_limit;   // outside inner limit, which are not constant;
    Vec _constant;              // difference of outer value for threshold cut;
    //vector<pair<int, double> > _window_max;     // max correlation coefficient within window;
    //vector<pair<char, double> > _position_max;  // max correlation coefficient at position;
    Mat bppm;                   // base pairing probability matrix of original sequence;
    bool _init;                 // flag of file initiation;
    bool _byte;                 // output to binary file
    const int _precision;       // # of digits to keep correct;
    static const char* base;    // 4 bases;
    //functions for threshold cut;
    Mat& Get_inner(bool, int, bool);
    void Add_constant(int, int, double, bool);
    void Add_constant_inner(int, int, double);
    bool Set_constant_in(int, int, int, int);
    bool Set_constant_out(int, int, int, int);
    bool Set_constant_inner(int, int, int, int);
    double Get_min_inner(int, int);
    void Copy_const_inner(int, int, Matrix&, double);
    //functions for calculation of mutated secondary structure;
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
    //functions for setting;
    void Initialize_seq(string&);
    void Copy_matrix(int);
    void Set_limit();
    void Set_mpoint(int);
    void Set_index();
    void Print_mlist(int, string&);
    void Change_sequence(int, int, int, string);    
    void Print_diff();
    void Print_neighborhood();
    void Calc_matrix(int, string&);    
    void All_calculation(int, int, string&);
    void Part_calculation(int, int, string&);

    void SNV_Change_sequence(int, int, int, string, vector<bool>&);
    void SNV_calculation(int, int, string&, vector<bool>&);

    //functions for output;
    void Output_common_data();
    void Output_delta_ene_line(double, int, int);
    void Output_delta_ene_line(double, double, int);
    int Calc_bpp_ene(Vec&);
    int Get_ori_start(int, int);
    int Get_ori_end(int, int);    
    int Get_original_stem(Vec&, int, int);
    void Get_mut_stem(Vec&, int, int);   
    int Store_stem(Vec&, Vec&, int, int);

    int Store_bpp(Vec&, Vec&);
    int Get_position_binary();    
    void Open_output_file(ofstream& ofs);
    void Out_header(int, ofstream&);
    void Out_value(double, double, int, ofstream&);
    void Output_delta_ene();
    void Check_stem(const Vec&, const Vec&);
    void Output_delta_around();    
    double Euclidean_Distance(const Vec&, const Vec&);
    void Output_ED_CC();
    void Output_acc_ene();
    void Output_acc_around();    
    void Output_correlation(const double);    
    void Output_correlation(const Vec&);
    void Output_storage(const string&);
    void Storage_max(const Vec&);
    void Set_Output_file(int);
    void Set_Correlation(int, int, string&, bool);
    //functions for debugging;
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
    void Get_mut_ori_mat(Mat&, Mat&, int, string);
    void Debug_bppm(int, string&);
    void Debug_bppm(int, string&, int, int);
    void Debug_multiple_mutations_stem();    
    void Debug_confirm(int, string&);
    void Debug_output(int, int, bool, Rfold_Lang&);
    void compare(int, const Vec&, const Vec&);
    void compare(int, const Mat&, const Mat&);
    bool compare_same(int, const Vec&, const Vec&);    
    bool compare_same(int, const Mat&, const Mat&);
    int Check_Difference(const class Matrix&, const class Matrix&);     
    void Output_Difference(int, const class Matrix& ori, const class Matrix& mut);

public:
    int Mtype;                  // mutation type;
    int window, rightw, leftw;  // range data for calculation of correlation coefficient;
	Matrix ori_alpha;          // inside matrix of original sequence;
	Matrix ori_beta;           // outside matrix of original sequence;
    Wobble _wob;               // matrix to store difference for threshold calculation;
    string outf;               // output filename;
    string binf;               // output filename for binary data;
    int id;                    // sequnce alignment id;
    int all_length;               // for multiple mutations;
    string ori_seq;            // original sequence;
    string bp_seq;             // annotation of base pair;
    vector<int> cons_index;    // consensus index of original sequence;
    enum Type { In, Del, Mut, Stem, Stemend, Multi, Multi1, Multi2, Multibif };
    enum Out { Abs, Rel, Acc, RelAcc };
    enum Distance { StemDif, CC };
    enum Debug { Reg, Dif, Time, Bpp, Outer, Analyze };
    enum Calc { No, Conv, Omit, OmitOut };
    enum Cor { One, All, Spe };
    // flags for calculation type and debug output;
    int _cor;
    int _output;
    int _debug;
    int _dist;
    static const int _calc = Calc::No;
    static const bool rdebug = false;
    Radiam(bool acc = false, int precision = DEF_PRE, int window = 0, bool init = false, bool dist = false) 
    : Rfold_Lang(), _init(init), _precision(precision), window(window)
    {
        _byte = false;
        id = 0;
        _output = (acc) ? Out::Acc : Out::Abs;
        _dist = (dist) ? Distance::CC : Distance::StemDif;
        _cor = Cor::One;
        _debug = Debug::Reg;
    }
    Radiam(string& filename, bool acc, int precision, int window, bool init, bool dist = false)
    : Rfold_Lang(), _init(init), _precision(precision), window(window)
    {
        _byte = false;
        if (filename != "") outf = filename;
        id = 0;                
        _output = (acc) ? Out::Acc : Out::Abs;
        _dist = (dist) ? Distance::CC : Distance::StemDif;
        _cor = Cor::One;
        _debug = Debug::Reg; 
    }
	virtual ~Radiam(){}

	void Mutation_calculation(int, string&);
    void Get_ori_matrix(const string&);
    void Set_example();
    void Set_annotate(string, string, vector<int>&, int length = 0, bool diff = true);
    void Initialize_before_calc(int, int, int&, string&);
    void Correlation_of_bpp(int, vector<int>&, int, string);
    void Correlation_of_bpp(int, int, int, string);
    void Mutation_analysis(int, int, int, string&, vector<bool>&);    
    void SNV_analysis(int, int, int, string&, vector<bool>&);        
    // functions to print base pairing probability of original sequence;
    void Write_bpp_binary();
    void Set_bpp_binary(int, int, string&);
    void Check_probability(double value) {
        if (value > 1.0 || value < 0.0) 
            cerr << "error " << value << endl;
    }
    bool Including(int i, int j) {
        if (_index[i-1] < 0 || _index[j-1] < 0) return false;
        else return true;
    }
    bool Both_including(int i) {
        if (_index[i-1] < 0 || (int)bppm.size() <= _index[i-1]) return false;
        else return true;
    }
    bool Both_including(int i, int j) {
        if (_index[j-1] < 0 || (int)bppm[_index[i-1]].size() < _index[j-1]-_index[i-1]) return false;
        else return true;
    }
    double bppm_get(int i, int j) {
        return bppm[_index[i-1]][_index[j-1]-_index[i-1]-1];
    }
    bool Under_Prec(double max, double min, double value) {
        return (log10(fabs(max-min))-log10(fabs(value)) < -(_precision)-3);
    }
    int In_range(int j, int start) {
        if (Mtype == Mut) return j-(start+_constraint+2);
        else if (Mtype == Del) return j-(start+_constraint+2);
        else return j-(start+_constraint+3); 
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