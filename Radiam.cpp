#include "Radiam.h"

namespace Rfold {

const char* Radiam::base = "ACGU";
struct Add {
    Add(double Constant) : constant(Constant) {}
    double constant;
    double operator() (double value) const {
        return Logsum(value, constant);
    }
};

Mat& Radiam::Get_inner(bool inside, int type, bool mut)
{
    Matrix& structure = (mut ? ((inside) ? alpha : beta) : ((inside) ? ori_alpha : ori_beta));
    switch(type) {
        case Stem: return structure.stem;
        case Stemend: return structure.stemend;
        case Multi: return structure.multi;
        case Multi1: return structure.multi1;
        case Multi2: return structure.multi2;
        case Multibif: default: return structure.multibif;
    }
}

void Radiam::Add_constant(int start, int end, double cons, bool inside)
{
    Add add(cons);
    if (start > end) return;
    (inside) ? transform(alpha.outer.begin()+start, alpha.outer.begin()+end+1, alpha.outer.begin()+start, add)
             : transform(beta.outer.begin()+start, beta.outer.begin()+end+1, beta.outer.begin()+start, add);
    if (_conv_out) cout << "! add cons " << inside << " " << " " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;
}

void Radiam::Add_constant_inner(int start, int end, double cons)
{
    Add add(cons);
    if (start > end) return;
    for (int type = Stem; type <= Multibif; type++) {
        Mat& structure = Get_inner(false, type, true);        
        for (int i = start; i <= end; i++)
            transform(structure[i].begin(), structure[i].end(), structure[i].begin(), add);
    }
    if (_conv_out) cout << "! add inner " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;    
}

bool Radiam::Set_constant_in(int j, int left, int right, int mp) 
{
    _wob.a[j] = alpha.outer[j]-ori_alpha.outer[_index[j]];
    if (In_range(j, left) >= 0) {
        double pre = _wob.a[j-_constraint-2];
        if (In_range(j, left) == 0 || pre == _wob.min || pre == _wob.max) {
            _wob.max = *max_element(_wob.a.begin()+j-_constraint, _wob.a.begin()+j+1);
            _wob.min = *min_element(_wob.a.begin()+j-_constraint, _wob.a.begin()+j+1);
        }
        if (_outer_fluc && _wob.max-_wob.min > 0.0) cout << j << " " << _wob.max-_wob.min << " " << alpha.outer[j] << endl;
        if (!_omit) return false;
        if (Under_Prec(_wob.max, _wob.min, alpha.outer[j])) {
            Add_constant(_right_limit[mp] = j+1, right, _constant[mp] = _wob.a[j], true);
            return true;
        }
    }
    return false;
}

bool Radiam::Set_constant_out(int j, int left, int right, int mp) 
{
    _wob.b[j] = beta.outer[j]-ori_beta.outer[_index[j]];
    if (Out_range(j, right) <= 0) {
        double pre = _wob.b[j+_constraint+2];
        if (Out_range(j, right) == 0 || pre == _wob.min || pre == _wob.max) {
            _wob.max = *max_element(_wob.b.begin()+j, _wob.b.begin()+j+_constraint+1);
            _wob.min = *min_element(_wob.b.begin()+j, _wob.b.begin()+j+_constraint+1);
        }
        if (_outer_fluc && _wob.max-_wob.min > 0.0) cout << j << " " << _wob.max-_wob.min << " " << beta.outer[j] << endl;
        if (!_omit) return false;
        if (Under_Prec(_wob.max, _wob.min, beta.outer[j])) {
            Add_constant(left, _left_limit[mp] = j-1, _wob.b[j], false);
            return true;
        }
    }
    return false;
}

bool Radiam::Set_constant_inner(int j, int left, int right, int mp) 
{
    Get_min_inner(j, mp);
    if (Out_range(j, right) <= 0) {
        double pre_max = _wob.imax[j+_constraint+2];
        double pre_min = _wob.imin[j+_constraint+2];
        if (Out_range(j, right) == 0 || pre_min == _wob.min || pre_max == _wob.max) {
            _wob.max = *max_element(_wob.imax.begin()+j, _wob.imax.begin()+min(right, j+_constraint+1));
            _wob.min = *min_element(_wob.imin.begin()+j, _wob.imin.begin()+min(right, j+_constraint+1));
        }
        if (Under_Prec(_wob.max, _wob.min, _wob.max/2.0+_wob.min/2.0)) {
            if (rdebug) cout << _wob.max << " " << _wob.min << endl;
            Add_constant_inner(left, j-1, _wob.max/2.0+_wob.min/2.0);
            return true;
        }
    }
    return false;
}

double Radiam::Get_min_inner(int j, int mp)
{
    double min_val;
    for (int type = Stem; type <= Multibif; type++) {
        Mat& mut = Get_inner(false, type, true);
        Mat& ori = Get_inner(false, type, false);
        for (int i = 0; i < (int)mut[j].size(); i++) {
            double value = mut[j][i]-ori[_index[j]][i];
            if (Is_INF(mut[j][i]) && Is_INF(ori[_index[j]][i])) continue;
            else if (_wob.imin[j] > value) _wob.imin[j] = value;
            else if (_wob.imax[j] < value) _wob.imax[j] = value;
        }
        min_val = min(min_val, *min_element(mut[j].begin(), mut[j].end()));
    }
    if (rdebug) cout << "Get min " << _wob.min << " " << _wob.max << endl;
    return min_val;
}

void Radiam::Copy_const_inner(int i, int j, Matrix& mut, double constant)
{
    mut.stem[j][j-i] = Logsum(mut.stem[j][j-i], constant);
    mut.stemend[j][j-i] = Logsum(mut.stemend[j][j-i], constant);
    mut.multi[j][j-i] = Logsum(mut.multi[j][j-i], constant);
    mut.multi1[j][j-i] = Logsum(mut.multi1[j][j-i], constant);
    mut.multi2[j][j-i] = Logsum(mut.multi2[j][j-i], constant);
    mut.multibif[j][j-i] = Logsum(mut.multibif[j][j-i], constant);
}

void Radiam::Calc_inside()
{
    for (int mp = 0; mp < (int)_mpoint.size(); mp++) 
    {
        int start = _mpoint[mp];
        int end = (mp == (int)_mpoint.size()-1) ? seq.length : _mpoint[mp+1]-1;
        for (int j = max(start, TURN+1); j <= end; j++) {
            if (rdebug) cout << "-------------------\n-j " << j << endl;    
            if (In_range(j, start) <= 0) {
                for (int i = min(j-TURN, start+1); i >= max(0, j-_constraint-1); i--) {
                    if (rdebug) cout << "--i " << i << endl;                
                    Calc_inside_mat(i, j);
                }
            }
            Calc_in_outer(j);
            if (!_general && Set_constant_in(j, start, end, mp)) break;  //(start+1)+_constraint+1 
        }
    }
}

void Radiam::Calc_out_outer()
{
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) 
    {
        int start = _mpoint[mp]+1;
        int end = (mp > 0) ? _mpoint[mp-1]+2 : 0; 
        for (int j = min(start, seq.length-1); j >= end; j--) {
            beta.outer[j] = beta.outer[j+1];
            for (int k = j+TURN; k <= min(j+_constraint+1, seq.length); k++) {
                double temp = Logsum(alpha.stem[k][k-j], beta.outer[k], 
                                     Parameter::Sum_Dangle(bp(j+1, k, seq.sequence), j, k+1, seq));
                beta.outer[j] = Logsumexp(beta.outer[j], temp);
            }
            if (!_general && Set_constant_out(j, end, start, mp)) break;  // (start-1)-(_constraint+1) 
        }
    }
}

void Radiam::Calc_outside_inner(int j, int mp)
{
    for (int i = max(0, j-_constraint-1); i < j-TURN-1; i++) 
    {
        if (rdebug) cout << "--i " << i << " " << _mpoint[mp] << endl;        
        if (i > 0 && j < seq.length) {
            if (!_general && !_omit && Out_in_range(i, j, mp)) {
                if (mp > 0) {
                    double value = _constant[(int)_mlist.size()-1]-(_constant[mp]-_constant[mp-1]);
                    Copy_const_inner(i, j, beta, value);
                }
            } else {
                Calc_outside_mat(i, j);
                if (rdebug) cout << "calc" << endl;
            }
        }
        beta.stem[j][j-i] = Calc_out_stem(i, j);
    } 
}


void Radiam::Calc_min_inside()
{
    for (int mp = 0; mp < (int)_mpoint.size(); mp++) 
    {
        int start = _mpoint[mp];
        int end = (mp == (int)_mpoint.size()-1) ? seq.length : _mpoint[mp+1]-1;
        for (int j = max(start, TURN+1); j <= end; j++) {
            if (rdebug) cout << "-------------------\n-j " << j << endl;    
            if (In_range(j, start) <= 0) {
                for (int i = min(j-TURN, start+1); i >= max(0, j-_constraint-1); i--) {
                    if (rdebug) cout << "--i " << i << endl;                
                    Calc_inside_mat(i, j);
                }
            }
            Calc_in_outer(j);
            if (!_general && Set_constant_in(j, start, end, mp)) break;  //(start+1)+_constraint+1 
        }
    }
}

void Radiam::Calc_min_outside()
{
    ;
}

void Radiam::Calc_min_out_outer()
{
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) 
    {
        int start = _mpoint[mp]+1;
        int end = (mp > 0) ? _mpoint[mp-1]+2 : 0; 
        for (int j = min(start, seq.length-1); j >= end; j--) {
            beta.outer[j] = beta.outer[j+1];
            for (int k = j+TURN; k <= min(j+_constraint+1, seq.length); k++) {
                double temp = Logsum(alpha.stem[k][k-j], beta.outer[k], 
                                     Parameter::Sum_Dangle(bp(j+1, k, seq.sequence), j, k+1, seq));
                beta.outer[j] = Logsumexp(beta.outer[j], temp);
            }
            if (!_general && Set_constant_out(j, end, start, mp)) break;  // (start-1)-(_constraint+1) 
        }
    }
}

void Radiam::Add_outside_inner(int mp, int& j) 
{
    if (In_range(j, _right_limit[mp]) >= 0) {
        Add_constant_inner(_right_limit[mp]+_constraint+1, j, _constant[(int)_mlist.size()-1]);
        j = _right_limit[mp]+_constraint;
    } else if (j == seq.length+1) j = seq.length;
}

void Radiam::Calc_outside()
{
    Calc_out_outer();
    int j = seq.length;
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) {
        if (!_general) Add_outside_inner(mp, j);
        int end = (mp > 0) ? _mpoint[mp-1]+1 : TURN+1;
        int const_end = (mp > 0) ? max(_left_limit[mp], _right_limit[mp-1]) : TURN+1;
        for (; j >= end; j--) {
            if (rdebug) cout << "----------------\n-j " << j << endl;
            Calc_outside_inner(j, mp);
        }
    }
}

double Radiam::Calc_bpp_cor(const Vec& bpp_mut, const Vec& bpp_ori) 
{
    double sum_sq_x = 0.0, sum_sq_y = 0.0, sum_coproduct = 0.0, mean_x = 0.0, mean_y = 0.0, n = 1.0;
    for (Vec::const_iterator it = bpp_ori.begin(), it2 = bpp_mut.begin(); it != bpp_ori.end() && it2 != bpp_mut.end(); it++, it2++) 
    {
        double x = *it, y = *it2;
        double delta_x = (x-mean_x), delta_y = (y-mean_y), sweep = (n-1.0)/n;
        sum_sq_x += delta_x * delta_x * sweep;
        sum_sq_y += delta_y * delta_y * sweep;
        sum_coproduct += delta_x * delta_y * sweep;
        mean_x += delta_x / n;
        mean_y += delta_y / n++;
    }
    if (sum_sq_x == 0.0 || sum_sq_y == 0.0) return numeric_limits<double>::quiet_NaN();
    else return sum_coproduct/sqrt(sum_sq_x*sum_sq_y);
}

void Radiam::Get_cor_vec(int start, int end, Vec& bpp_mut, Vec& bpp_ori)
{
    int limit = min(end, seq.length);
    for (int i = max(start, 1); i < min(limit, seq.length); i++) {
        if (_index[i-1] < 0) continue;
        for (int j = i+1; j <= min(limit, i+_constraint); j++) {
            if (_index[j-1] < 0 || _index[j-1] > _index[i-1]+_constraint) continue;
            bpp_mut.push_back((_acc) ? acc(i, j) : bpp(i, j));
            bpp_ori.push_back(bppm[_index[i-1]][_index[j-1]-_index[i-1]-1]);
        }
    }
}

void Radiam::Get_cor_vec(int start, int end, Vec& bpp_mut, Vec& bpp_ori, const Mat& bppm_mut)
{
    int limit = min(end, seq.length);
    for (int i = max(start, 1); i < min(limit, seq.length); i++) {
        if (_index[i-1] < 0) continue;
        for (int j = i+1; j <= min(limit, i+_constraint); j++) {
            if (_index[j-1] < 0 || _index[j-1] > _index[i-1]+_constraint) continue;
            bpp_mut.push_back(bppm_mut[i-1][j-i-1]);
            bpp_ori.push_back(bppm[_index[i-1]][_index[j-1]-_index[i-1]-1]);
        }
    }
}

void Radiam::Calc_one_bpp_cor()
{
    Vec bpp_mut, bpp_ori;
    Get_cor_vec(_mlist[0]+1-window, _mlist[0]+1+window, bpp_mut, bpp_ori);
    Output_correlation(Calc_bpp_cor(bpp_mut, bpp_ori));
}

void Radiam::Calc_bpp_cor() 
{
    Vec output;
    Mat bppm_mut;
    Write_bpp(bppm_mut);
    for (int start = window; start <= seq.length; start += window) {
        Vec bpp_mut, bpp_ori;
        int limit = min(seq.length, start+window);
        Get_cor_vec(start-window, limit, bpp_mut, bpp_ori, bppm);
        output.push_back(Calc_bpp_cor(bpp_mut, bpp_ori));
    }
    (_general) ? Storage_max(output) : Output_correlation(output);
}

}