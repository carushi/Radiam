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
    if (_outer_fluc) cout << std::setprecision(10);
    Initialize_seq(sequence);
    Copy_matrix(type);
    Calc_inside();
    Calc_outside();
    if (_outer_fluc) {
        cout << "#inside" << endl; Print_Vec(_wob.a);
        cout << "#outside" << endl; Print_Vec(_wob.b);
    } else if (_analyze) {
        printf("* partition func %le\n", exp(alpha.outer[seq.length]));
        Debug_confirm(type, sequence);
        Debug_bppm(type, sequence);
    } else if (_bpp_fluc && Mtype == Mut) {
        if (_acc) Write_accm_fluc();
        else {
            Mat bpp_mut; Write_bpp(bpp_mut);
            Write_bppm_fluc(bpp_mut);            
        }
    } else if (window > 0) {
        (_single) ? Calc_one_bpp_cor() : Calc_bpp_cor();
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
        if (_analyze) Print_mlist(type, sequence);
        Calc_matrix(type, sequence);
    }
}
void Radiam::Part_calculation(int type, int k, string& sequence)
{
    if (k > 0) {
        Change_sequence(type, k, _mlist[k-1], sequence);
    } else {
        if (_analyze) Print_mlist(type, sequence);
        Calc_matrix(type, sequence);
    }
}


void Radiam::Output_storage(const string& sequence)
{
    if (_mlist.size() == 1) {
        ofstream ofs(winf.c_str(), (_init) ? ios::app : ios::trunc);
        ofs << "*s-e\tpos\tcor" << endl;
        for (int i = 0; i < _window_max.size(); i++)
            ofs << i*window << "-" << (i+2)*window << "\t" << _window_max[i].first << "\t" << _window_max[i].second << endl;
        ofs << "************" << endl;
        ofs.close(); 
        ofs.open(posif.c_str(), (_init) ? ios::app : ios::trunc);
        ofs << "*pos\tbefore\tafter\tmin_cor" << endl;
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
        ofs << "************" << endl; ofs.close();    
        _init = true;
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
    model.calculation(_constraint, sequence);
    (_acc) ? model.Write_acc(bppm) : model.Write_bpp(bppm);
    ori_alpha = model.alpha;
    ori_beta = model.beta;
}

void Radiam::Set_Output_file()
{
    stringstream ss;
    if (_single) {
        ss << "../correlation/one_" << window << "_cons_" << _constraint << ".txt";
        onef = ss.str();
    } else {
        ss << "../correlation/window_" << window << "_cons_" << _constraint << ".txt";
        winf = ss.str(); ss.str("");
        ss << "../correlation/position_" << window << "_cons_" << _constraint << ".txt";
        posif = ss.str();
    }
}

void Radiam::Initialize_before_calc(int type, int k, int constraint, string& sequence)
{
    int length = (int)sequence.length()+((type == In) ? k : 0);
    Set_Constraint(constraint, length);
    Set_Output_file();
    cout << onef << endl;
    if (window > 0) {
        _window_max = vector<pair<int, double> >(length/window, pair<int, double>(0, 1.0));
        _position_max = vector<pair<char, double> >(length, pair<char, double>(0, 1.0));
    }
    Get_ori_matrix(sequence);
}

void Radiam::Set_Correlation(int type, int constraint, string& sequence, bool all) 
{
    cout << "#-----------------------\n# " << sequence << " dim: " << type << " thres: " 
         << _precision << " win:" << window << " cons:" << constraint << endl;
    if (!_time && (_single || !_single && !_general)) cout << "#type\tabase\tconst\twindow\tpos\trlim\tllim\tcor" << endl;         
    (all) ? All_calculation(type, (int)_mlist.size(), sequence) : Part_calculation(type, (int)_mlist.size(), sequence);
}

void Radiam::Correlation_of_bpp(int type, vector<int>& mlist, int constraint, string sequence)
{
    _mlist = mlist;    
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);    
    if (General_Output() && _single && window > 0) {
        streambuf* last = cout.rdbuf(); 
        ofstream ofs(onef.c_str(), (_init) ? ios::app : ios::trunc);
        cout.rdbuf(ofs.rdbuf());
        Set_Correlation(type, constraint, sequence, false);
        ofs.close(); cout.rdbuf(last);    
    } else Set_Correlation(type, constraint, sequence, false);
}

void Radiam::Correlation_of_bpp(int type, int k, int constraint, string sequence)
{
    _mlist = vector<int>(k);
    Initialize_before_calc(type, (int)_mlist.size(), constraint, sequence);
    if (General_Output() && _single && window > 0) {
        streambuf* last = cout.rdbuf();     
        ofstream ofs(onef.c_str(), (_init) ? ios::app : ios::trunc);
        cout.rdbuf(ofs.rdbuf());
        Set_Correlation(type, constraint, sequence, true);
        All_calculation(type, k, sequence);
        ofs.close(); cout.rdbuf(last);            
    } else Set_Correlation(type, constraint, sequence, true);
}

}
