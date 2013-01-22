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
    if (rdebug) cout << "! add cons " << inside << " " << cons << " " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;    
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
    if (rdebug) cout << "! add inner " << cons << " " << start << " " << end << " " << _constraint << " " << seq.length << " " << endl;    
}

bool Radiam::Set_constant_in(int j, int left, int right, int mp) 
{
    _wob.a[j] = alpha.outer[j]-ori_alpha.outer[index(j, mp, true)];
    double pre = _wob.a[j-_constraint-2];
    if (j == left+_constraint+2 || pre == _wob.min || pre == _wob.max) {
        _wob.max = *max_element(_wob.a.begin()+j-_constraint, _wob.a.begin()+j+1);
        _wob.min = *min_element(_wob.a.begin()+j-_constraint, _wob.a.begin()+j+1);
    }
    if (Under_Prec(_wob.max-_wob.min, alpha.outer[j])) {
        Add_constant(_right_limit[mp] = j+1, right, _constant[mp] = _wob.a[j], true);
        return true;
    }
    return false;
}

bool Radiam::Set_constant_out(int j, int left, int right, int mp) 
{
    _wob.b[j] = beta.outer[j]-ori_beta.outer[index(j, mp, false)];
    double pre = _wob.b[j+_constraint+2];
    if (j == right-_constraint-2 || pre == _wob.min || pre == _wob.max) {
        _wob.max = *max_element(_wob.b.begin()+j, _wob.b.begin()+j+_constraint+1);
        _wob.min = *min_element(_wob.b.begin()+j, _wob.b.begin()+j+_constraint+1);
    }
    if (Under_Prec(_wob.max-_wob.min, beta.outer[j])) {
        _constant[mp] = _wob.b[j];
        Add_constant(left, _left_limit[mp] = j-1, _constant[mp], false);
        return true;
    }
    return false;
}

bool Radiam::Set_constant_inner(int j, int left, int right, int mp) 
{
    double min_val = Get_min_inner(j, mp);
    double pre_max = _wob.imax[j+_constraint+2];
    double pre_min = _wob.imin[j+_constraint+2];
    if (j == right-_constraint-2 || pre_min == _wob.min || pre_max == _wob.max) {
        _wob.max = *max_element(_wob.imax.begin()+j, _wob.imax.begin()+j+_constraint+1);
        _wob.min = *min_element(_wob.imin.begin()+j, _wob.imin.begin()+j+_constraint+1);
    }
    if (Under_Prec(_wob.max-_wob.min, pre_max/2.0+pre_min/2.0)) {
        Add_constant_inner(left, j-1, pre_max/2.0+pre_min/2.0);
        return true;
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
            double value = mut[j][i]-ori[index(j, mp, false)][i];
            if (_wob.imin[j] > value) _wob.imin[j] = value;
            else if (_wob.imax[j] < value) _wob.imax[j] = value;
        }
        min_val = min(min_val, *min_element(mut[j].begin(), mut[j].end()));
    }
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
        for (int j = start; j <= end; j++) {
            if (rdebug) cout << "-------------------\n-j " << j << endl;    
            if (j <= start+_constraint+2) {
               for (int i = min(j-TURN, start+1); i >= max(0, j-_constraint-1); i--) {
                   if (rdebug) cout << "--i " << i << endl;                
                   Calc_inside_mat(i, j);
               }
            }
            Calc_in_outer(j);
            if (j >= start+_constraint+2 && Set_constant_in(j, start, end, mp)) break;  //(start+1)+_constraint+1 
        }
    }
    for (int j = TURN+1; j <= seq.length; j++) {
        if (debug) 
            cout << "-------------------\n-j " << j << endl;
        for (int i = j-TURN; i >= max(0, j-_constraint-1); i--) 
            Calc_inside_mat(i, j);
        Calc_in_outer(j);
    }

}

void Radiam::Calc_out_outer()
{
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) 
    {
        int start = _mpoint[mp]+1;
        int end = (mp > 0) ? _mpoint[mp-1]+2 : 0; 
        cout << "mp" << mp << " " << start << " " << end << endl;
        for (int j = min(start, seq.length-1); j >= end; j--) {
            beta.outer[j] = beta.outer[j+1];
            for (int k = j+TURN; k <= min(j+_constraint+1, seq.length); k++) {
                double temp = Logsum(alpha.stem[k][k-j], beta.outer[k], 
                                     Parameter::Sum_Dangle(bp(j+1, k, seq.sequence), j, k+1, seq));
                beta.outer[j] = Logsumexp(beta.outer[j], temp);
            }
            if (j <= start-_constraint-2 && Set_constant_out(j, end, start, mp)) break;  // (start-1)-(_constraint+1) 
        }

    }
}

void Radiam::Calc_outside_inner(int j, int mp)
{
    for (int i = max(0, j-_constraint-1); i < j-TURN-1; i++) {
        if (rdebug) cout << "--i " << i << " " << _mpoint[mp] << endl;        
        if (i != 0 && j < seq.length) {
            if (_mpoint.size() == 1) {
                if (Is_Range(i, j, mp)) Copy_const_inner(i, j, beta, _constant[mp]);
                else {
                    Calc_outside_mat(i, j);
                    if (rdebug) cout << "calc" << endl;
                }
            } else {
                Calc_outside_mat(i, j);
                if (rdebug) cout << "calc" << endl;
            }
        }
        beta.stem[j][j-i] = Calc_out_stem(i, j);
    } 
}

void Radiam::Calc_outside()
{
    Calc_out_outer();
    int preview = seq.length+1;
    for (int mp = (int)_mpoint.size()-1; mp >= 0; mp--) {
        if (mp > 0 && preview > _right_limit[mp]) {
            Add_constant_inner(_right_limit[mp], preview, _constant[mp]);
            preview = _right_limit[mp];
        }
        int j = preview-1;
        int end = (mp > 0) ? _mpoint[mp-1]+1 : TURN+1;
        int const_end = (mp > 0) ? _right_limit[mp-1] : TURN+1;
        if (rdebug) cout  << "start" << j << " " << end << " " << const_end << " " << _mpoint[mp] << endl;
        for ( ; j >= end; j--) {
            if (rdebug) cout << "----------------\n-j " << j << endl;
            Calc_outside_inner(j, mp);
            if (j >= const_end && j <= _mpoint[mp]-_constraint-2) {
                if (Set_constant_inner(j, end, _mpoint[mp], mp)) {
                    j = end-1; break;
                }
            }
        }
        if (rdebug) cout << "change " << end << " " << const_end << endl;
        preview = j+1;
    }
}

void Radiam::Initialize_seq(string& sequence)
{
    Set_sequence(sequence);
    _wob = Wobble(seq.length+1);
}

double Radiam::Calc_bpp_cor(const Vec& bpp_mut, int type) 
{
    double sum_sq_x = 0.0, sum_sq_y = 0.0, sum_coproduct = 0.0, mean_x = 0.0, mean_y = 0.0, n = 1.0;
    for (Vec::const_iterator it = bpp.begin(), it2 = bpp_mut.begin(); it != bpp.end() && it2 != bpp_mut.end(); it++, it2++) 
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


void Radiam::Calc_matrix(int type, string& sequence)
{
    Initialize_seq(sequence);    
    Copy_matrix(type);
    Calc_inside();
    Calc_outside();
    if (rdebug || true) Debug_confirm(type, sequence);
    return;
    if (_matrix) {
        Mat bppm_mut;
        Write_bpp(bppm_mut);
        //cout << "cor " << Calc_bpp_cor(bppm_mut) << endl;
    } else {
        Vec bpp_mut;
        Write_bpp(bpp_mut);
        cout << "cor " << Calc_bpp_cor(bpp_mut, type) << endl;
    }
}

void Radiam::Copy_matrix(int type)
{
    Mtype = type;
    Set_mpoint(type);
    Set_limit();
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

void Radiam::Set_limit()
{
    int k = (int)_mpoint.size();
    _right_limit = vector<int>(k);
    for (int i = 0; i < k; i++) {
        if (i < k-1) _right_limit[i] = _mpoint[i+1]-1;
        else _right_limit[i] = seq.length+1;
    }
    _left_limit = vector<int>(k, -1);
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

void Radiam::All_calculation(int type, int k, string& sequence)
{
    if (k > 0) {
        int max = (k < (int)_mlist.size()) ? _mlist[k] : (int)sequence.length();
        int min = k-1;        
        if (type == In) { max++; min = 0; } /* Insertion -> same position OK */
        for (int i = min; i < max; i++) {
            _mlist[k-1] = i;
            Change_sequence(type, k, i, sequence);
        }
    } else {
        if (rdebug) {
            cout << sequence.length() << " " << sequence << "\n! mlist: ";
            for (int i = (int)_mlist.size()-1; i >= 0; i--) 
                cout << _mlist[i] << " ";
            cout << " mpoint: ";
        }
        Calc_matrix(type, sequence);
    }
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
    (!_matrix) ? model.Write_bpp(bpp) : model.Write_bpp(bppm);
    ori_alpha = model.alpha;
    ori_beta = model.beta;

}

void Radiam::Correlation_of_bpp(int type, int k, int constraint, string sequence)
{
    Set_Constraint(constraint, (int)sequence.length()+((type == In) ? k : 0));
    _mlist = vector<int>(k); 
    Get_ori_matrix(sequence);
    cout << sequence << endl;
    All_calculation(type, k, sequence);
}

/*      For debug       */
void Radiam::Debug_output(int temp, int type, bool inside, Rfold_Lang& model)
{
    if (inside) {
        if (type == Mut) {
            Output_Difference(temp, alpha, model.alpha); Output_Difference(temp, model.alpha, ori_alpha);
        } else {
            Print_Mat(false); model.Print_Mat(false);
        }
    } else {
        if (type == Mut) {
            Output_Difference(temp, beta, model.beta); Output_Difference(temp, model.beta, ori_beta);
            Print_Mat(false); model.Print_Mat(false);            
        } else {
            Print_Mat(false); model.Print_Mat(false);
        }
    }
}

void Radiam::Debug_confirm(int type, string& sequence)
{
    printf("* partition func %le\n", exp(alpha.outer[seq.length]));        
    Rfold_Lang model2;                
    model2.calculation(_constraint, sequence);
    int temp = 0;
    if ((temp = Check_Difference(alpha, model2.alpha)) > 0) {
        cout << "Inside Differ!" << temp << endl;
        Debug_output(temp, type, true, model2);
        //exit(0);
    } else if ((temp = Check_Difference(beta, model2.beta)) > 0) {
        cout << "Outside Differ!" << temp << endl;
        if (temp == Multi) cout << "Oh no!" << endl;
        Debug_output(temp, type, false, model2);
        //exit(0);
    } else cout << "Same" << endl;
}

bool Radiam::compare_same(int type, const Vec& ori, const Vec& mut) 
{
    for (Vec::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) {
        double value = 0.0;
        if (Is_INF(*it) && Is_INF(*it2)) value = 0.0;
        else if (Is_INF(*it) || Is_INF(*it2)) value = INF;
        else value = abs(*it2-*it);
        if (log10(value)-log10(*it2) >= -DEF_PRE) {
            cout << "dif " << value << " mut " << *it << " ori " << *it2 << " " << log10(value) << " " << log10(*it) << endl;
            return true;
        }
    }
    return false;
}

bool Radiam::compare_same(int type, const Mat& ori, const Mat& mut) 
{
    for (Mat::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) 
        if (compare_same(type, *it, *it2)) return true;
    return false;
}

int Radiam::Check_Difference(const class Matrix& ori, const class Matrix& mut) 
{
    if (compare_same(0, ori.outer, mut.outer)) return 1;
    else if (compare_same(2, ori.stemend, mut.stemend)) return Stemend;    
    else if (compare_same(1, ori.stem, mut.stem)) return Stem;
    else if (compare_same(3, ori.multi, mut.multi)) return Multi;
    else if (compare_same(4, ori.multi2, mut.multi2)) return Multi2;
    else if (compare_same(5, ori.multibif, mut.multibif)) return Multibif;
    else if (compare_same(6, ori.multi1, mut.multi1)) return Multi1;
    else return 0;
}

void Radiam::compare(int type, const Vec& ori, const Vec& mut) 
{
    for (Vec::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) {
        double value = 0.0;
        if (Is_INF(*it) && Is_INF(*it2)) value = 0.0;
        else if (Is_INF(*it) || Is_INF(*it2)) value = INF;
        else value = abs(*it-*it2);
        cout << value << "\t";
    }
    cout << endl;
}

void Radiam::compare(int type, const Mat& ori, const Mat& mut) 
{
    for (Mat::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) 
        compare(type, *it, *it2);
    cout << endl;
}

void Radiam::Output_Difference(int temp, const Matrix& ori, const Matrix& mut) 
{
    switch(temp) {
        case 1: compare(0, ori.outer, mut.outer); break;
        case Stem: compare(1, ori.stem, mut.stem); break;
        case Stemend: compare(2, ori.stemend, mut.stemend); break;
        case Multi: compare(3, ori.multi, mut.multi); break;
        case Multi2: compare(4, ori.multi2, mut.multi2); break;
        case Multibif: compare(5, ori.multibif, mut.multibif); break;
        case Multi1: default: compare(6, ori.multi1, mut.multi1); break;
    }
}
/*                   */
}
