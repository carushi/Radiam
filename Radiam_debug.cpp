#include "Radiam.h"

namespace Rfold {

void Radiam::Calc_time(int type, string& sequence)
{
    clock_t t1, t2;
    t1 = clock();
    Initialize_seq(sequence); Copy_matrix(type); Calc_inside(); Calc_outside();
    t2 = clock();
    cout << "* time\t" << (double)(t2-t1)/CLOCKS_PER_SEC;
    Rfold_Lang model2;
    t1 = clock();
    model2.calculation(_constraint, sequence);
    t2 = clock();
    cout << "\t" << (double)(t2-t1)/CLOCKS_PER_SEC << endl;
    Mat bpp_mut, bpp_ori;
    Write_bpp(bpp_mut);
    model2.Write_bpp(bpp_ori);
    Write_bppm_dif(bpp_mut, bpp_ori);
}

void Radiam::Calc_debug(int type, string& sequence)
{
    switch(_debug) {
        case Debug::Outer:
            cout << "#inside" << endl; Print_Vec(_wob.a);
            cout << "#outside" << endl; Print_Vec(_wob.b); 
            break;
        case Debug::Bpp:
            if (_output == Out::RelAcc) Write_accm_fluc();
            else if (_output == Out::Acc) {
                for (int i = 1; i < 6; i++) {
                    Write_accm_fluc_abs(2*i-1);
                }
            } else {
                Mat bpp_mut;
                Write_bpp(bpp_mut);
                (_output == Out::Rel) ? Write_bppm_fluc(bpp_mut) : Write_bppm_fluc_abs(bpp_mut);
            }
            break;
        case Debug::Analyze:
            printf("* partition func %le\n", exp(alpha.outer[seq.length]));
            Debug_confirm(type, sequence); Debug_bppm(type, sequence);
            break;
        default: break;
    }
}

void Radiam::Debug_output(int temp, int type, bool inside, Rfold_Lang& model)
{
    if (inside) {
        if (type == Mut) {
            Output_Difference(temp, alpha, model.alpha); Output_Difference(temp, model.alpha, ori_alpha);
            //Print_Mat(true); model.Print_Mat(true);            
        } else {
            Output_Difference(temp, alpha, model.alpha);
            //Print_Mat(true); model.Print_Mat(true);
        }
    } else {
        if (type == Mut) {
            Output_Difference(temp, beta, model.beta); Output_Difference(temp, model.beta, ori_beta);
            //Print_Mat(false); model.Print_Mat(false);            
        } else {
            Output_Difference(temp, beta, model.beta); 
            //Print_Mat(false); model.Print_Mat(false);
        }
    }
}

double Radiam::Get_diff(double a, double b)
{
    double diff = fabs(a-b);
    if (diff != 0.0) {
        double relative_diff = log10(diff)-log10(a);
        return relative_diff;
    }
    return -INF;
}

double Radiam::Get_diff_at_same(int i, int j, const Mat& bpp_mut) 
{
    return bpp_mut[i][j-i-1]-bppm[_index[i]][_index[j]-_index[i]-1];
}

bool Radiam::Within_range(int i, const Mat& bpp_mut) 
{
    if (i < (int)bpp_mut.size() && _index[i] < (int)bppm.size()) return true;
    else return false;
}

bool Radiam::Within_range(int i, int j, const Mat& bpp_mut) 
{
    if (_index[i] < 0) return false;
    int max_ind = i+(int)bpp_mut[i].size();
    int max_oind = _index[i]+(int)bppm[_index[i]].size();
    if (j >= (int)bpp_mut.size() || _index[j] >= (int)bppm.size()) return false;
    else if (j <= max_ind && _index[j] <= max_oind) return true;
    else return false;
}

void Radiam::Write_bppm_fluc(const Mat& bpp_mut)
{
    vector<pair<double, double> > max_fluc((int)bpp_mut.size(), pair<double, double>(-INF, 0.0));
    for (int i = 0; Within_range(i, bpp_mut); i++) {
        for (int j = i+1; Within_range(i, j, bpp_mut); j++) {
            if (_index[j] < 0) continue;
            double diff = fabs(Get_diff_at_same(i, j, bpp_mut));
            if (diff == 0.0) continue;
            double relative_diff = log10(diff)-log10(bpp_mut[i][j-i-1]);
            if (max_fluc[i].first < relative_diff) { 
                max_fluc[i] = pair<double, double>(relative_diff, bpp_mut[i][j-i-1]);
            }
            if (max_fluc[j].first < relative_diff) { 
                max_fluc[j] = pair<double, double>(relative_diff, bpp_mut[i][j-i-1]); 
            }
        }
    }
    for (int i = 0; i < (int)max_fluc.size(); i++) {
        if (_index[i] < 0) continue;
        if (max_fluc[i].first > -INF)
            cout << "* bpp_fluc " << _index[i]+1 << " " << max_fluc[i].first << " " << max_fluc[i].second << endl;
    }
}



void Radiam::Write_bppm_fluc_abs(const Mat& bpp_mut)
{
    vector<pair<double, double> > max_fluc((int)bpp_mut.size(), pair<double, double>(0, 0.0));
    vector<pair<double, double> > min_fluc((int)bpp_mut.size(), pair<double, double>(0, 0.0));
    for (int i = 0; Within_range(i, bpp_mut); i++) {
        for (int j = i+1; Within_range(i, j, bpp_mut); j++) {
            if (_index[j] < 0) continue;
            double diff = Get_diff_at_same(i, j, bpp_mut);
            if (diff >= 0.0) {
                if (max_fluc[i].first < diff) max_fluc[i] = pair<double, double>(diff, bpp_mut[i][j-i-1]);
                if (max_fluc[j].first < diff) max_fluc[j] = pair<double, double>(diff, bpp_mut[i][j-i-1]); 
            } else {
                if (min_fluc[i].first > diff) min_fluc[i] = pair<double, double>(diff, bpp_mut[i][j-i-1]);
                if (min_fluc[j].first > diff) min_fluc[j] = pair<double, double>(diff, bpp_mut[i][j-i-1]); 
            }
        }
    }
    for (int i = 0; i < (int)max_fluc.size(); i++) {
        if (_index[i] < 0) continue;
        cout << "* bpp_abs "  << _index[i]+1 << " " << max_fluc[i].first << " " << max_fluc[i].second << " "
        << min_fluc[i].first << " " << min_fluc[i].second << endl;
    }
}

void Radiam::Write_accm_fluc()
{
    vector<pair<double, double> > max_fluc(seq.length, pair<double, double>(-INF, 0.0));
    for (int i = 1; i <= seq.length; i++) {
        double a = acc(i, i), relative_diff = Get_diff(a, bppm[i-1][0]);
        if (max_fluc[i-1].first < relative_diff) { 
            max_fluc[i-1] = pair<double, double>(relative_diff, a);
        }
        if (i == seq.length) break;
        a = acc(i, i+1); relative_diff = Get_diff(a, bppm[i-1][1]);
        if (max_fluc[i-1].first < relative_diff) { 
            max_fluc[i-1] = pair<double, double>(relative_diff, a);
        }
        if (max_fluc[i].first < relative_diff) { 
            max_fluc[i] = pair<double, double>(relative_diff, a);
        }
    }
    for (int i = 0; i < (int)max_fluc.size(); i++) {
        if (max_fluc[i].first > -INF)
            cout << "* acc_fluc " << i+1 << " " << max_fluc[i].first << " " << max_fluc[i].second << endl;
    }
}

void Radiam::Write_accm_fluc_abs(int band)
{
    vector<pair<double, double> > max_fluc(seq.length, pair<double, double>(-INF, 0.0));
    for (int i = 1; i+band <= seq.length; i++) {
        if (_index[i-1] < 0 || _index[i-1+band] < 0) continue;
        double acc_mut = acc(i, i+band);
        double acc_dif = acc_mut-bppm[i-1][band];
        cout << "* acc_abs\t" << band+1 << "\t" << i << "\t" << acc_dif << "\t" << acc_mut 
        << " \t" << -Parameter::kT*(log(acc_mut)-log(bppm[i-1][band]))/1000. << "\t" 
        << -Parameter::kT*log(acc_mut)/1000. << endl;
    }
}

void Radiam::Write_accm_fluc_abs()
{
    int start[] = {4, -50, -100, -200, -500};
    for (int band = 1; band <= 11; band += 2) {
        start[0] = band/2-1;
        for (int i = 0; i < 5; i++) {
            int left = _mlist[0]+1+start[i]-(band-1); int right = _mlist[0]+1+start[i];
            if (left < 1 || right > seq.length) continue;
            double acc_mut = acc(left, right), a = acc_mut-bppm[left-1][right-left];
            cout << "* acc_abs\t" << band << "\t" << left << "\t" << right << "\t" << a << "\t" << acc_mut 
            << " \t" << -Parameter::kT*(log(acc_mut)-log(bppm[left-1][right-left]))/1000. << "\t" 
            << -Parameter::kT*log(acc_mut)/1000. << endl;
        }
        for (int i = 1; i < 5; i++) {
            int left = _mlist[0]+1-start[i]; int right = _mlist[0]+1-start[i]+(band-1);
            if (left < 1 || right > seq.length) continue;
            double acc_mut = acc(left, right), a = acc_mut-bppm[left-1][right-left];
            cout << "* acc_abs\t" << band << "\t" << left << "\t" << right << "\t" << a << "\t" << acc_mut
            << "\t" << -Parameter::kT*(log(acc_mut)-log(bppm[left-1][right-left]))/1000. << "\t"
            << -Parameter::kT*log(acc_mut)/1000. << endl;
        }
    }
}

void Radiam::Write_bppm_dif(const Mat& bpp_mut, const Mat& bpp_ori)
{
    double tmax = -INF, t1 = 0.0;
    for (int i = 0; i < (int)bpp_mut.size(); i++) {
        for (int j = 0; j < (int)bpp_mut[i].size(); j++) {
            double diff = (_output == Out::Abs) ? fabs(bpp_ori[i][j]-bpp_mut[i][j])
                                   : Get_diff(bpp_ori[i][j], bpp_mut[i][j]);
            if (rdebug) cout << i << "," << i+1+j << ": " << diff << " " << bpp_mut[i][j] << endl;
            if (tmax < diff) {
                tmax = diff; t1 = bpp_ori[i][j];
            }
        }
    }
    cout << "* same_bpp_fluc " << tmax << " " << t1 << endl;
}

void Radiam::Write_bpp_binary()
{
    ofstream ofs;
    if (_init) ofs.open(outf.c_str(), ios::trunc | ios::binary);
    else ofs.open(outf.c_str(), ios::app | ios::binary);
    vector<double> stem = vector<double>(seq.length, 0.0);
    int count = Calc_bpp_ene(stem);    
    Out_header(count, ofs);
    for (int i = 1; i <= seq.length; i++) {
        if (_index[i-1] < 0) continue;
        Out_value(stem[i-1], 0.0, i, ofs);
    }
}    

void Radiam::Set_bpp_binary(int type, int constraint, string& sequence)
{
    int length = (int)sequence.length()+((type == In) ? 1 : 0);
    Set_Constraint(constraint, length);
    Get_ori_matrix(sequence);
    Mtype = type;
    _mlist.push_back(0);    
    Initialize_seq(sequence);
    Copy_matrix(type);
}

void Radiam::Get_mut_ori_mat(Mat& bpp_mut, Mat& bpp_ori, int constraint, string sequence)
{
    Rfold_Lang model2;
    model2.calculation(_constraint, sequence);
    if (_output >= Out::Acc) {    
        Write_acc(bpp_mut); model2.Write_acc(bpp_ori);
    } else {
        Write_bpp(bpp_mut); model2.Write_bpp(bpp_ori);
    }
}

void Radiam::Debug_bppm(int type, string& sequence)
{
    Mat bpp_mut, bpp_ori;
    Get_mut_ori_mat(bpp_mut, bpp_ori, _constraint, sequence);
    Write_bppm_dif(bpp_mut, bpp_ori);
}

void Radiam::Debug_bppm(int type, string& sequence, int l, int r)
{
    Mat bpp_mut, bpp_ori;
    Get_mut_ori_mat(bpp_mut, bpp_ori, _constraint, sequence);
    if (type != Mut) {
        cerr << "cannot compare" << endl;
    }
    cerr << _mpoint[0] << " " << l << " " << r << endl;
    double tmax = -INF, t1 = 0.0, t2 = 0.0;
    for (int i = 0; i < (int)bpp_mut.size(); i++) {
        for (int j = 0; j < min(_constraint, (int)bpp_mut[i].size()); j++) {
            if (l < i+1 || r > i+j+1) continue;
            //if (l > i+j+1 || i+1 > r) continue;
            double diff = (_output == Out::Abs) ? fabs(bpp_ori[i][j]-bpp_mut[i][j])
                                   : Get_diff(bpp_ori[i][j], bpp_mut[i][j]);
            Check_probability(bpp_ori[i][j]);
            Check_probability(bpp_mut[i][j]);            
            if (diff > 1e-6) {
                cerr << i+1 << " " << i+j+1 << " " << diff << " " << bpp_ori[i][j] << " " << bpp_mut[i][j] << endl;
            }
            if (tmax < diff) {
                tmax = diff; t1 = bpp_mut[i][j]; t2 = bpp_ori[i][j];
                if (diff > 0.00001)
                    cerr << "* diff " << tmax << " " << t1 << " " << t2 << endl;
            }
        }
    }
    cerr << "* same_bpp_fluc " << tmax << " " << t1 << " " << t2 << endl;    
}

void Radiam::Debug_multiple_mutations_stem()
{
    Mat bpp_mut, bpp_ori;
    Get_mut_ori_mat(bpp_mut, bpp_ori, _constraint, seq.str.substr(1));
    vector<double> stem_o = vector<double>(seq.length, 0.0);
    vector<double> stem_m = vector<double>(seq.length, 0.0);    
    for (int i = 0; i < (int)bpp_mut.size(); i++) {
        for (int j = 0; j < (int)bpp_mut[i].size(); j++) {
            stem_m[i] += bpp_mut[i][j];
            stem_o[i] += bpp_ori[i][j];
            stem_m[i+1+j] += bpp_mut[i][j];
            stem_o[i+1+j] += bpp_ori[i][j];
            assert(fabs(bpp_mut[i][j]-bpp_ori[i][j]) < 0.00000001);
        }
    }
}


void Radiam::Debug_confirm(int type, string& sequence)
{
    Rfold_Lang model2;                
    model2.calculation(_constraint, sequence);
    int temp = 0;
    if ((temp = Check_Difference(alpha, model2.alpha)) > 0) {
        cerr << "Dif" << endl;
        cout << "Inside Differ!" << temp << endl;
        Debug_output(temp, type, true, model2);
    } else if ((temp = Check_Difference(beta, model2.beta)) > 0) {
        cerr << "Dif" << endl;        
        cout << "Outside Differ!" << temp << endl;
        Debug_output(temp, type, false, model2);
    }
}

bool Radiam::compare_same(int type, const Vec& ori, const Vec& mut) 
{
    for (Vec::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) {
        double value = 0.0;
        if (Is_INF(*it) && Is_INF(*it2)) value = 0.0;
        else if (Is_INF(*it) || Is_INF(*it2)) value = INF;
        else value = *it2-*it;
        if (value != 0.0 && (*it2 != 0 || log10(value)-log10(*it2) >= -DEF_PRE)) {
            cout << "dif " << value << " mut " << *it << " ori " << *it2 << " " << log10(value) << " " << log10(*it2) << endl;
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
        else value = *it-*it2;
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

}
