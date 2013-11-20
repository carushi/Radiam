#include "Radiam.h"

namespace Rfold {

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
    for (int ori = 0, mp = 0, mut = 0; mut <= seq.length; ) {
        if (Mtype != Mut && (mp < (int)_mpoint.size() && _mpoint[mp] == ori)) {
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

void Radiam::Print_diff() 
{
    if (_output == Out::Acc) Output_acc_ene();
    else (_dist) ? Output_ED_CC() : Output_delta_ene();
}

void Radiam::Print_neighborhood()
{
    if (_output == Out::Acc) Output_acc_around();    
    else Output_delta_around();
}

void Radiam::Calc_matrix(int type, string& sequence)
{
    if (_debug == Debug::Outer) cout << std::setprecision(10);
    Initialize_seq(sequence);
    Copy_matrix(type);
    if (_mlist.size() == 1 && _debug == Debug::Reg) {
        if (_cor == Cor::One && window > 0) {
            leftw = max(1, _mpoint[0]+1-window);
            rightw = min(seq.length, _mpoint[0]+1+window);
        } else {
            leftw = max(1, _mpoint[0]-_constraint);
            rightw = min(seq.length, _mpoint[0]+_constraint);
        }
        /*
        Calc_inside(); 
        Calc_outside();
        Copy_matrix(type);        
        */
        Calc_min_inside();
        double diff = Calc_partition_function();
        Calc_min_outside(diff);
        Print_neighborhood();        
        //Debug_bppm(type, sequence, leftw, rightw);        
    } else {
        Calc_inside(); 
        Calc_outside();
        if (window > 0 && _debug == Debug::Reg) Calc_bpp_cor();
        else if (_debug != Debug::Reg) {
            if (_mlist.size() < 3 && _debug == Debug::Dif) Print_diff();
            else Calc_debug(type, sequence);
        }
    }
    
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
        if (_debug == Debug::Analyze) Print_mlist(type, sequence);
        if (_debug == Debug::Time) Calc_time(type, sequence);
        else Calc_matrix(type, sequence);
    }
}
void Radiam::Part_calculation(int type, int k, string& sequence)
{
    if (k > 0) {
        Change_sequence(type, k, _mlist[k-1], sequence);
    } else {
        if (_debug == Debug::Analyze) Print_mlist(type, sequence);
        if (_debug == Debug::Time) Calc_time(type, sequence);
        else Calc_matrix(type, sequence);
    }
}

void Radiam::Output_delta_ene_line(double diff, int i, int j)
{
    cout << cons_index[i-1] << "\t" << cons_index[j-1] << "\t" << diff << endl;
}
void Radiam::Output_delta_ene_line(double diff, double value, int i)
{
    cout << cons_index[i-1] << "\t" << diff << "\t" << value << endl;
}

int Radiam::Calc_bpp_ene(Vec& stem)
{
    int count = 0;
    for (int i = 1; i <= seq.length; i++) {
        if (!Both_including(i)) continue;
        for (int j = i+TURN; j <= min(seq.length, i+_constraint); j++) {
            if (!Both_including(i, j)) continue;
            double bpp_ori = bppm_get(i, j);
            stem[i-1] += bpp_ori; stem[j-1] += bpp_ori;
        }
        count++;
    }
    return count;
}


int Radiam::Get_ori_start(int start, int end)
{
    for (int i = start; i <= end; i++) {
        if (Both_including(i))  return i;
    }
    return end;
}

int Radiam::Get_ori_end(int start, int end)
{
    for (int i = end; i >= start; i--) {
        if (Both_including(i)) return i;
    }
    return start;
}

int Radiam::Get_original_stem(Vec& stem_o, int start, int end)
{
    int ostart = Get_ori_start(start, end);
    int oend = Get_ori_end(start, end);
    assert(ostart < oend);
    int count = 0;
    for (int i = ostart; i <= oend; i++) {
        if (!Both_including(i)) continue;
        int pos = _index[i-1]+1;
        for (int j = max(pos-_constraint, 1); j <= min(seq.length, pos+_constraint); j++) {
            if (abs(pos-j) <= TURN) continue;
            double bpp_ori = bppm[min(j, pos)-1][abs(j-pos)-1];
            stem_o[i-start] += bpp_ori;
        }
        count++;
    }
    return count;
}

void Radiam::Get_mut_stem(Vec& stem_m, int start, int end)
{
    for (int i = start; i <= end; i++) {
        if (!Both_including(i)) continue;
        int pos = i;
        for (int j = max(pos-_constraint, 1); j <= min(seq.length, pos+_constraint); j++) {
            if (abs(pos-j) <= TURN) continue;
            double bpp_ori = bpp(i, j, true);
            if (rdebug) cout << "bpp " << i << " " << j << " " << bpp_ori << endl;
            stem_m[i-start] += bpp_ori;
        }
    }
}


int Radiam::Store_stem(Vec& stem_m, Vec& stem_o, int start = -1, int end = -1)
{
    if (start < 0) start = 1;
    if (end < 0) end = seq.length;
    int count = Get_original_stem(stem_o, start, end);
    Get_mut_stem(stem_m, start, end);
    return (end-start+1);
}

int Radiam::Store_bpp(Vec& bpp_m, Vec& bpp_o)
{
    int count = 0;
    for (int i = 1; i <= seq.length; i++) {
        if (!Both_including(i)) continue;
        for (int j = i+TURN; j <= min(seq.length, i+_constraint); j++) {
            if (!Both_including(i, j)) continue;
            double bpp_mut = bpp(i, j);
            double bpp_ori = bppm_get(i, j);
            bpp_m.push_back(bpp_mut);
            bpp_o.push_back(bpp_ori);
            count++;            
        }
    }
    return count;
}

int Radiam::Get_position_binary() 
{
    int position = -1;
    for (int i = 0; i < (int)_mlist.size(); i++) {
        int tpos = (_mlist[i] >= (int)cons_index.size()) 
            ? cons_index[(int)cons_index.size()-1]+1 : cons_index[_mlist[i]];
        position = (position+1)*(all_length+1)+tpos;
    }
    return position;
}

void Radiam::Open_output_file(ofstream& ofs)
{
    if (_init) {
        if (_byte) {
            ofs.open(outf.c_str(), ios::trunc | ios::binary);
        } else {
            ofs.open(outf.c_str(), ios::trunc);
        }
    } else {
        if (_byte) {
            ofs.open(outf.c_str(), ios::app | ios::binary);
        } else {
            ofs.open(outf.c_str(), ios::app | ios::binary);        
        }
    }
}    

void Radiam::Out_header(int count, ofstream& ofs)
{
    int position = Get_position_binary();
    char after = 'X';
    if (Mtype != Del) 
        after = seq.str[_mlist[0]+1];
    char before = 'X';
    if (_mlist[0] < (int)_index.size() && _index[_mlist[0]] >= 0) 
        before = ori_seq[_mlist[0]];
    if (_byte) {
        ofs.write((char*)&id, sizeof(int));
        ofs.write((char*)&(position), sizeof(int));
        ofs.write((char*)&(after), sizeof(char));
        ofs.write((char*)&(before), sizeof(char));
        ofs.write((char*)&count, sizeof(int));
    } else {
        ofs << id << "\t" << position << "\t" << after << "\t"
            << before << "\t" << count;
    }
}

void Radiam::Out_value(double diff, double value, int i, ofstream& ofs)
{
    if (_byte) {
        ofs.write((char*)&(cons_index[_index[i-1]]), sizeof(int));
        ofs.write((char*)&(diff), sizeof(double));    
        ofs.write((char*)&(value), sizeof(double));
    } else {
        ofs << cons_index[_index[i-1]] << "\t" << diff << "\t" << value << endl;
    }
}

void Radiam::Output_delta_ene()
{
    ofstream ofs;
    Open_output_file(ofs);
    if (!ofs) return;
    _init = false;        
    Vec stem_m = Vec(seq.length, 0.0);
    Vec stem_o = Vec(seq.length, 0.0);
    int count = Store_stem(stem_m, stem_o);
    Out_header(count, ofs);
    for (int i = 1; i <= seq.length; i++) {
        if (!Both_including(i)) continue;
        double diff = stem_o[i-1]-stem_m[i-1];
        Out_value(diff, stem_m[i-1], i, ofs);
        //double diff_ene = -Parameter::kT*(log(stem_o[i-1])-log(stem_m[i-1]))/1000.;
        //Output_delta_ene_line(diff, stem_m[i-1], i, ofs);
    }
    ofs.close();   
}

void Radiam::Check_stem(const Vec& stem_m, const Vec& stem_o)
{
    cerr << "Check_ori" << endl;
    for (int i = leftw; i <= rightw; i++) {
        Check_probability(stem_o[i-leftw]);
    }
    cerr << "Check_mut" << endl;    
    for (int i = leftw; i <= rightw; i++) {
        Check_probability(stem_m[i-leftw]);    
    }
}

void Radiam::Output_delta_around()
{
    ofstream ofs;
    Open_output_file(ofs);
    if (!ofs) return;
    _init = false;
    Vec stem_m = Vec(rightw-leftw+1, 0.0);
    Vec stem_o = Vec(rightw-leftw+1, 0.0);
    int count = Store_stem(stem_m, stem_o, leftw, rightw);
    Check_stem(stem_m, stem_o);
    Out_header(rightw-leftw+1, ofs);
    double diff_max = 0;
    for (int i = leftw; i <= rightw; i++) {
        if (!Both_including(i)) continue;
        diff_max = max(fabs(stem_o[i-leftw]-stem_m[i-leftw]), diff_max);
    }
    double CC_stem = Calc_bpp_cor(stem_m, stem_o);    
    Out_value(diff_max, CC_stem, 0, ofs);    
    ofs.close();   
}

double Radiam::Euclidean_Distance(const Vec& mut, const Vec& ori) 
{
    double distance = 0.0;
    for (int i = 0; i < (int)mut.size(); i++) {
        distance += (mut[i]-ori[i])*(mut[i]-ori[i]);
    }
    return sqrt(distance);
}

void Radiam::Output_ED_CC()
{
    ofstream ofs;
    if (_init) ofs.open(outf.c_str(), ios::trunc | ios::binary);
    else ofs.open(outf.c_str(), ios::app | ios::binary);
    _init = false;
    Vec stem_m = Vec(seq.length, 0.0), stem_o = stem_m;
    Vec bpp_m = Vec(seq.length, 0.0), bpp_o = bpp_m;
    Store_stem(stem_m, stem_o);
    Store_bpp(bpp_m, bpp_o);
    Out_header(3, ofs);
    double ED_stem = Euclidean_Distance(stem_m, stem_o);
    double ED_bpp = Euclidean_Distance(bpp_m, bpp_o);
    double CC_stem = Calc_bpp_cor(stem_m, stem_o);
    double CC_bpp = Calc_bpp_cor(bpp_m, bpp_o);
    Out_value(ED_stem, ED_bpp, 0, ofs);
    Out_value(CC_stem, CC_bpp, 0, ofs);
    Out_value(eSDC(CC_stem), eSDC(CC_bpp), 0, ofs);
    /*cout << ED_stem << "\t" << ED_bpp << "\t" << CC_stem << "\t" << CC_bpp 
    << "\t" << eSDC(CC_stem) << "\t" << eSDC(CC_bpp) << endl;*/
    ofs.close();  
}

void Radiam::Output_acc_around()
{
    cout << "null" << endl;
}

void Radiam::Output_acc_ene()
{
    for (int i = 1; i <= seq.length; i++) {
        for (int j = i; j <= min(seq.length, i+_constraint); j++) {
            if (!Including(i, j)) continue;
            double acc_mut = acc(i, j);
            double acc_ori = bppm_get(i, j);
            //double diff = acc_ori-acc_mut;
            double diff_ene = -Parameter::kT*(log(acc_ori)-log(acc_mut))/1000.;
            Output_delta_ene_line(diff_ene, i, j);
        }
    }
}
void Radiam::Output_common_data()
{
    char m = (Mtype == Del) ? 'X' : seq.str[_mpoint[0]+1];
    cout << Mtype << "\t" << m << "\t" << _constraint << "\t" << window << "\t";
    Print_Vec(_mpoint, false);
    cout << "\t";
    Print_Vec(_right_limit, false);
    cout << "\t";    
    Print_Vec(_left_limit, false);
    cout << "\t";        
}

void Radiam::Output_correlation(const double output)
{
    if ((int)_mlist.size() > 1) return;
    Output_common_data();
    cout << output << endl;
}

void Radiam::Output_correlation(const Vec& output)
{
    if ((int)_mlist.size() > 1) return;
    Output_common_data();
    Print_Vec(output);
}
    
void Radiam::Mutation_calculation(int i, string& sequence) // special function
{
    _mlist = vector<int>(1, i);
    Calc_matrix(Mut, sequence);    
}

void Radiam::Get_ori_matrix(const string& sequence)
{
    Rfold_Lang model;
    model.calculation(_constraint, sequence, false);
    (_output >= Out::Acc) ? model.Write_acc(bppm) : model.Write_bpp(bppm);
    ori_alpha = model.alpha;    
    ori_beta = model.beta;
}

void Radiam::Set_example()
{
    _debug = Debug::Bpp;
}

void Radiam::Set_annotate(string bp, string ori, vector<int>& cons, int length, bool diff)
{
    bp_seq = bp;
    ori_seq = ori;
    all_length = max(length, (int)ori.length());
    cons_index = cons;
    if (diff)
        _debug = Debug::Dif;
}

void Radiam::Set_Output_file(int type)
{
    if (outf != "") return;    
    stringstream ss;
    if (_debug == Debug::Reg) {
        if (_cor == Cor::All)
            ss << "../correlation/all_" << window;
        else if (_cor == Cor::Spe) {
            if (_output == Out::Acc) ss << "../correlation/acc_" << leftw << "_" << rightw;
            else ss << "../correlation/bpp_" << leftw << "_" << rightw;
        } else ss << "../correlation/one_" << window;
        ss << "_cons_" << _constraint << "_type_" << type << ".txt";                        
        outf = ss.str();
    } 
}

void Radiam::Initialize_before_calc(int type, int k, int& constraint, string& sequence)
{
    int length = (int)sequence.length();
    if (type == In) length += k;
    if (constraint > length) constraint = length;
    Set_Constraint(constraint, length);
    Set_Output_file(type);
    if (cons_index.size() == 0) {
        vector<int> cons;
        for (int i = 0; i < length; i++) cons.push_back(i);
        Set_annotate(string(""), sequence, cons, length, false);
    }
    Get_ori_matrix(sequence);
}

void Radiam::Set_Correlation(int type, int constraint, string& sequence, bool all) 
{
    cout << "#-----------------------\n# " << sequence << " dim: " << type << " thres: " 
         << _precision << " win:" << window << " cons:" << constraint << " id:" << id << endl;
    cout << "# ori_len: " << all_length << " ori: " << ori_seq << " bp: " << bp_seq << endl;
    cout << "# " << outf << endl;
    if (_debug == Debug::Reg) cout << "#type\tabase\tconst\twindow\tpos\trlim\tllim\tcor" << endl;         
    (all) ? All_calculation(type, (int)_mlist.size(), sequence) : Part_calculation(type, (int)_mlist.size(), sequence);
}

void Radiam::Correlation_of_bpp(int type, vector<int>& mlist, int constraint, string sequence)
{
    _mlist = mlist;    
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);    
    if (_debug == Debug::Reg) {
        streambuf* last = cout.rdbuf(); 
        ofstream ofs(outf.c_str(), (_init) ? ios::trunc : ios::app);
        cout.rdbuf(ofs.rdbuf());
        Set_Correlation(type, constraint, sequence, false);
        Part_calculation(type, (int)_mlist.size(), sequence);
        ofs.close(); cout.rdbuf(last);    
    } else Set_Correlation(type, constraint, sequence, false);
}

void Radiam::Correlation_of_bpp(int type, int k, int constraint, string sequence)
{
    _mlist = vector<int>(k);
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);
    if (_debug == Debug::Reg) {
        streambuf* last = cout.rdbuf();
        ofstream ofs(outf.c_str(), (_init) ? ios::trunc : ios::app);
        cout.rdbuf(ofs.rdbuf());
        Set_Correlation(type, constraint, sequence, true);
        All_calculation(type, k, sequence);
        ofs.close(); cout.rdbuf(last);            
    } else Set_Correlation(type, constraint, sequence, true);
}

void Radiam::SNV_Change_sequence(int type, int k, int i, string sequence, vector<bool>& flag)
{
    if (type == Del) {
        sequence.erase(sequence.begin()+i);
        SNV_calculation(type, k-1, sequence, flag);
    } else {
        for (int j = 0; j < 4; j++) {
            if (type == Mut && base[j] == sequence[i]) continue;
            string str = sequence.substr(0, i)+base[j];
            if (i < (int)sequence.length()) {
                if (type == In) str += sequence[i];                
                if (i+1 < (int)sequence.length()) str += sequence.substr(i+1);
            }
            SNV_calculation(type, k-1, str, flag);
        }
    }
}

void Radiam::SNV_calculation(int type, int k, string& sequence, vector<bool>& flag)
{
    if (k > 0) {
        int max = (k < (int)_mlist.size()) ? _mlist[k] : (int)sequence.length();
        int min = k-1;
        if (type == In) { max++; min = 0; } // Insertion -> same position OK;
        for (int i = min; i < max; i++) {
            if (flag[i]) {
                _mlist[k-1] = i;
                SNV_Change_sequence(type, k, i, sequence, flag);
            }
        }
    } else {
        if (_debug == Debug::Analyze) Print_mlist(type, sequence);
        if (_debug == Debug::Time) Calc_time(type, sequence);
        else {
            Calc_matrix(type, sequence);
        }
    }
}

void Radiam::Mutation_analysis(int type, int k, int constraint, string& sequence, vector<bool>& flag) 
{
    _mlist = vector<int>(k);
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);
    if (_debug == Debug::Reg) {
        streambuf* last = cout.rdbuf();
        ofstream ofs(outf.c_str(), (_init) ? ios::trunc : ios::app);
        cout.rdbuf(ofs.rdbuf());
        SNV_calculation(type, k, sequence, flag);
        ofs.close(); cout.rdbuf(last);            
    } else Set_Correlation(type, constraint, sequence, true);
}

void Radiam::SNV_analysis(int type, int k, int constraint, string& sequence, vector<bool>& flag) 
{
    _mlist = vector<int>(k);
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);
    SNV_calculation(type, k, sequence, flag);
}

}
