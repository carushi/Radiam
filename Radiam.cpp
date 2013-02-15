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
    if (analyze) cout << "! add cons " << inside << " " << " " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;
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
    if (analyze) cout << "! add inner " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;    
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
        if (rdebug || _all_omit) cout << j << " " << _wob.max-_wob.min << " " << alpha.outer[j] << endl;
        if (_all_omit) return false;
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
        if (rdebug || (_all_omit && _wob.max-_wob.min > 0.0)) cout << j << " " << _wob.max-_wob.min << " " << beta.outer[j] << endl;
        if (_all_omit) return false;
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
        if (rdebug) cout << "innner" << _wob.max << " " << _wob.min << endl;
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
            if (Set_constant_in(j, start, end, mp)) break;  //(start+1)+_constraint+1 
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
            if (Set_constant_out(j, end, start, mp)) break;  // (start-1)-(_constraint+1) 
        }
    }
}

void Radiam::Calc_outside_inner(int j, int mp)
{
    for (int i = max(0, j-_constraint-1); i < j-TURN-1; i++) 
    {
        if (rdebug) cout << "--i " << i << " " << _mpoint[mp] << endl;        
        if (i > 0 && j < seq.length) {
            if (Out_in_range(i, j, mp)) {
                if (mp > 0) Copy_const_inner(i, j, beta, _constant[mp-1]);
            } else {
                Calc_outside_mat(i, j);
                if (rdebug) cout << "calc" << endl;
            }
        }
        beta.stem[j][j-i] = Calc_out_stem(i, j);
    } 
}

void Radiam::Add_outside_inner(int mp, int& j) 
{
    if (In_range(j, _right_limit[mp]) >= 0) {
        Add_constant_inner(_right_limit[mp]+_constraint+1, j, _constant[(int)_mlist.size()-1]);
        j = _right_limit[mp]+_constraint;
    } else if (j == seq.length+1) j = seq.length;
    /*
    if (In_range(j, _right_limit[mp]) >= 0) {
        Add_constant_inner(_right_limit[mp]+_constraint+1, j, _constant[mp]);
        j = _right_limit[mp]+_constraint;
    } else if (j == seq.length+1) j = seq.length;
    */
}

void Radiam::Calc_outside()
{
    Calc_out_outer();
    int j = seq.length;
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) {
        Add_outside_inner(mp, j);
        int end = (mp > 0) ? _mpoint[mp-1]+1 : TURN+1;
        int const_end = (mp > 0) ? max(_left_limit[mp], _right_limit[mp-1]) : TURN+1;
        for (; j >= end; j--) {
            if (rdebug) cout << "----------------\n-j " << j << endl;
            Calc_outside_inner(j, mp);
        }
    }
    /*    
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) {
        Add_outside_inner(mp, j);
        int end = (mp > 0) ? _mpoint[mp-1]+1 : TURN+1;
        int const_end = (mp > 0) ? _right_limit[mp-1] : TURN+1;
        for (; j >= end; j--) {
            if (rdebug) cout << "----------------\n-j " << j << endl;
            Calc_outside_inner(j, mp);
            if (!_omit&& !_all_omit && j >= const_end && j < _mpoint[mp]
                && Set_constant_inner(j, end, _mpoint[mp], mp)) break;   
        }
        j = end;
    }*/
}

//////////////////////////////

void Radiam::Initialize_seq(string& sequence)
{
    Set_sequence(sequence);
    _wob = Wobble(seq.length+1);
}


void Radiam::Copy_matrix(int type)
{
    Mtype = type;
    Set_mpoint(type);
    Set_limit();
    Set_index();
    alpha = ori_alpha; beta = ori_beta;
    if (type != Mut) {
        for (int i = (int)_mlist.size()-1; i >= 0; i--) {
            if (type == In) {
                alpha.insert(_mlist[i]); beta.insert(_mlist[i]);
            } else {
                alpha.erase(_mlist[i]);  beta.erase(_mlist[i]);
            }
        }
    }
    if (rdebug) cout << "size--" << alpha.stem.size() << endl;
}

void Radiam::Set_index()
{
    _index = vector<int>(seq.length+1);
    for (int ori = 0, mp = 0, mut = 0; mut < seq.length+1; ) {
        if (Mtype != Mut && _mpoint[mp] == ori) {
            if (Mtype == In) {
                _index[mut] = -1;
                mut++;
            } else ori++;
            mp++;
        } else {
            _index[mut] = ori;
            mut++; ori++;
        }
    }
}

void Radiam::Set_limit()
{
    int k = (int)_mpoint.size();
    _right_limit = vector<int>(k);
    for (int i = 0; i < k; i++) {
        if (i < k-1) _right_limit[i] = _mpoint[i+1]-1;
        else _right_limit[i] = seq.length+1;
    }
    _left_limit = vector<int>(k, -1);
    _inner_limit = vector<int>(k+1, -1);
    _constant = Vec(k, 0.0);
}

void Radiam::Set_mpoint(int type)
{
    _mpoint = vector<int>(_mlist.size()); 
    for (int i = (int)_mlist.size()-1; i >= 0; i--) {
        _mpoint[i] = _mlist[i];
        if (type == In && i > 0) 
            _mpoint[i] += i;
        else if (type == Del && i > 0) 
            _mpoint[i] -= i;
    }
    vector<int>::iterator end = unique(_mpoint.begin(), _mpoint.end());
    _mpoint.erase(end, _mpoint.end());
    if (rdebug) {
        cout << "mutation point : ";
        ostream_iterator<int> out_it(cout, ", ");
        copy(_mpoint.begin(), _mpoint.end(), out_it);        
        cout << endl;
    } 
}

void Radiam::Print_mlist(int type, string& sequence)
{
    cout << "* mlist: ";
    for (int i = (int)_mlist.size()-1; i >= 0; i--) {
        cout << _mlist[i] << " ";
        if (_mlist.size() == 1) {
            if (type != Del) cout << sequence[_mlist[0]];
            else cout << "X";
        } 
    }
    cout << endl;
}

void Radiam::Change_sequence(int type, int k, int i, string sequence)
{
    if (type == Del) {
        sequence.erase(sequence.begin()+i);
        All_calculation(type, k-1, sequence);
    } else {
        for (int j = 0; j < 4; j++) {
            if (type == Mut && base[j] == sequence[i]) continue;
            string str = sequence.substr(0, i)+base[j];
            if (i < (int)sequence.length()) {
                if (type == In) str += sequence[i];                
                if (i+1 < (int)sequence.length()) str += sequence.substr(i+1);
            }
            All_calculation(type, k-1, str);
        }
    }
}
void Radiam::Calc_matrix(int type, string& sequence)
{
    if (_time) { Calc_time(type, sequence); return; }
    if (_all_omit) cout << std::setprecision(10);
    Initialize_seq(sequence);
    Copy_matrix(type);
    Calc_inside();
    Calc_outside();
    if (_all_omit) {
        cout << "#inside" << endl; alpha.Print_Vec(_wob.a);
        cout << "#outside" << endl; alpha.Print_Vec(_wob.b);
        sleep(2);        
    } else if (analyze) {
        printf("* partition func %le\n", exp(alpha.outer[seq.length]));
        Debug_confirm(type, sequence);        
        Debug_bppm(type, sequence);
    } else if (window > 0) Calc_bpp_cor();
}

void Radiam::All_calculation(int type, int k, string& sequence)
{
    if (k > 0) {
        int max = (k < (int)_mlist.size()) ? _mlist[k] : (int)sequence.length();
        int min = k-1;
        if (type == In) { max++; min = 0; } // Insertion -> same position OK;
        for (int i = min; i < max; i++) {
            _mlist[k-1] = i;
            Change_sequence(type, k, i, sequence);
        }
    } else {
        if (analyze) Print_mlist(type, sequence);
        Calc_matrix(type, sequence);
    }
}
void Radiam::Part_calculation(int type, int k, string& sequence)
{
    if (k > 0) {
        Change_sequence(type, k, _mlist[k-1], sequence);
    } else {
        if (analyze) Print_mlist(type, sequence);
        Calc_matrix(type, sequence);
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

void Radiam::Calc_bpp_cor() 
{
    Vec output;
    Mat bppm_mut;
    Write_bpp(bppm_mut);
    for (int start = window; start <= seq.length; start += window) {
        Vec bpp_mut, bpp_ori;
        int limit = min(seq.length, start+window);
        for (int i = max(start-window, 1); i <= limit; i++) {
            if (_index[i-1] < 0) continue;
            for (int j = i+1; j <= min(limit, i+_constraint); j++) {
                if (_index[j-1] < 0 || _index[j-1] > _index[i-1]+_constraint) continue;
                bpp_mut.push_back(bppm_mut[i-1][j-i-1]);
                bpp_ori.push_back(bppm[_index[i-1]][_index[j-1]-_index[i-1]-1]);
            }
        }
        output.push_back(Calc_bpp_cor(bpp_mut, bpp_ori));
    }
    //if (!_general)
        Output_correlation(output);
    //else
        Storage_max(output);
}

void Radiam::Output_storage(const string& sequence)
{
    if (_mlist.size() == 1) {
        ofstream ofs("../window.txt");
        ofs << "s-e\tpos\tcor" << endl;
        for (int i = 0; i < _window_max.size(); i++)
            ofs << i*window << "-" << (i+2)*window << "\t" << _window_max[i].first << "\t" << _window_max[i].second << endl;
        ofs.close(); ofs.open("../position.txt");
        ofs << "pos\tbefore\tafter\tmin_cor" << endl;
        for (int i = 0; i < _position_max.size(); i++) {
            ofs << i << "\t";
            if (Mtype == Mut) ofs << sequence[i];
            else if (Mtype == Del) ofs << sequence[i];
            else {
                if (i > 0) ofs << sequence[i-1] << "[";
                if (i < sequence.length()) ofs << "]" << sequence[i];
            }
            ofs << "\t" << _position_max[i].first << "\t" << _position_max[i].second << endl;
        }
    }
}
void Radiam::Storage_max(const Vec& output)
{
    if (_mlist.size() == 1) {
        double min_elem = *min_element(output.begin(), output.end());
        if (min_elem <= _position_max[_mlist[0]].second) {
            _position_max[_mlist[0]].first = (Mtype != Del) ? seq.str[_mlist[0]+1] : 'X';
            _position_max[_mlist[0]].second = min_elem;
        }
        for (int i = 0; i < _window_max.size(); i++) {
            if (output[i] <= _window_max[i].second) 
                _window_max[i] = pair<int, double>(_mlist[0], output[i]);
        }
    }
}

void Radiam::Output_correlation(const Vec& output)
{
    if ((int)_mlist.size() > 1) return;
    char m = (Mtype == Del) ? 'X' : seq.str[_mpoint[0]+1];
    cout << Mtype << "\t" << m << "\t" << _constraint << "\t" << window << "\t";
    Print_Vec(_mpoint);
    cout << "\t";
    Print_Vec(_right_limit);
    cout << "\t";    
    Print_Vec(_left_limit);
    cout << "\t";        
    alpha.Print_Vec(output);
}
    
void Radiam::Mutation_calculation(int i, string& sequence) // special function
{
    _mlist = vector<int>(1, i);
    Calc_matrix(Mut, sequence);    
}

void Radiam::Get_ori_matrix(const string& sequence)
{
    Rfold_Lang model;
    model.calculation(_constraint, sequence);
    model.Write_bpp(bppm);
    ori_alpha = model.alpha;
    ori_beta = model.beta;
}
void Radiam::Initialize_before_calc(int type, int k, int constraint, string& sequence)
{
    cout << "#-----------------------\n# " << sequence << " dim: " << type << " thres: " 
         << _precision << " win:" << window << " cons:" << constraint << endl;
    int length = (int)sequence.length()+((type == In) ? k : 0);
    Set_Constraint(constraint, length);
    if (window > 0) {
        _window_max = vector<pair<int, double> >(length/window, pair<int, double>(0, 1.0));
        _position_max = vector<pair<char, double> >(length, pair<char, double>(0, 1.0));
    }
    Get_ori_matrix(sequence);
    if (!time && !_general) cout << "#type\tabase\tconst\twindow\tpos\trlim\tllim\tcor" << endl;
}

void Radiam::Correlation_of_bpp(int type, vector<int>& mlist, int constraint, string sequence)
{
    _mlist = mlist;
    Initialize_before_calc(type, (int)mlist.size(), constraint, sequence);
    Part_calculation(type, (int)_mlist.size(), sequence);
}

void Radiam::Correlation_of_bpp(int type, int k, int constraint, string sequence)
{
    _mlist = vector<int>(k);
    Initialize_before_calc(type, k, constraint, sequence);
    All_calculation(type, k, sequence);
    if (_general) Output_storage(sequence);
}

}
