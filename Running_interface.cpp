#include "Running_interface.h"

namespace Rfold{

void Running_interface::RNA_transform(string& str) 
{
    string::iterator it = str.begin();
    while (it != str.end()) {
        switch(*it) {
            case 'A': case 'C': case 'G': case 'U': break;
            case 'a': case 'c': case 'g': case 'u': *it = toupper(*it); break;
            case 'T': case 't': *it = 'U'; break;
            default: *it = 'A';
        }
        ++it;
    }
}

void Running_interface::Set_Files()
{
    ostringstream ss;
	for (int i = 0; i < MAXTYPE; i++) {
        switch(i) {
            case 0: ss << "outer"; break;
            case 1: ss << "stem"; break;
            case 2: ss << "stemend"; break;
            case 3: ss << "multi"; break;
            case 4: ss << "multi"; break;
            case 5: ss << "multi2"; break;
            case 6: ss << "multibif"; break;
            default: ss << "multi1";
        }
        ss << '_' << (int)sequence.length() << '_' << constraint << ".txt";
        filenames.push_back(ss.str());
        ss.str("");
	}
}

void Running_interface::Init_Files()
{
    ofstream ofs;
    for (int i = 0; i < MAXTYPE; i++) {
        ofs.open(filenames[i].c_str(), ios::trunc);
        ofs.close();
    }
}

void Running_interface::Same_Header(const char* header, int point, char ori, char mut, const string& seq)
{
    ostringstream oss;
    oss << header << "\t" << point << "\t" << ori << "\t" << mut << "\t" << seq;
    ofstream ofs;
    for (int i = 0; i < MAXTYPE; i++) {
        ofs.open(filenames[i].c_str(), ios::app);
        ofs << oss.str() << endl;
        ofs.close();
    }
}

bool Running_interface::compare_same(int type, const Vec& ori, const Vec& mut) 
{
    for (Vec::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) {
        double value = 0.0;
        if (Is_INF(*it) && Is_INF(*it2)) value = 0.0;
        else if (Is_INF(*it) || Is_INF(*it2)) value = INF;
        else value = exp(*it)-exp(*it2);
        if (log10(value)-log10(exp(*it)) > -DEF_PRE) {
            cout << value << " " << *it << " " << log10(value) << " " << log10(*it) << " " << DEF_PRE<< endl;
            return true;
        }
    }
    return false;
}

void Running_interface::compare(int type, const Vec& ori, const Vec& mut) 
{
    ofstream ofs(filenames[type].c_str(), ios::app);
    ofs << "* ";
    for (Vec::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) {
        double value = 0.0;
        if (Is_INF(*it) && Is_INF(*it2)) value = 0.0;
        else if (Is_INF(*it) || Is_INF(*it2)) value = INF;
        else value = *it-*it2;
        ofs << value << "\t";
    }
    ofs << endl;
    ofs.close();
}

bool Running_interface::compare_same(int type, const Mat& ori, const Mat& mut) 
{
    for (Mat::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) 
        if (compare_same(type, *it, *it2)) return true;
    return false;
}

void Running_interface::compare(int type, const Mat& ori, const Mat& mut) 
{
    for (Mat::const_iterator it = ori.begin(), it2 = mut.begin(); it != ori.end() && it2 != mut.end(); it++, it2++) 
        compare(type, *it, *it2);
    ofstream ofs(filenames[type].c_str(), ios::app);    
    ofs << endl;    
}

bool Running_interface::Check_Difference(const class Matrix& ori, const class Matrix& mut) 
{
    if (compare_same(0, ori.outer, mut.outer) || compare_same(1, ori.stem, mut.stem) ||
        compare_same(2, ori.stemend, mut.stemend) || compare_same(3, ori.multi, mut.multi) || 
        compare_same(4, ori.multi2, mut.multi2) || compare_same(5, ori.multibif, mut.multibif) || 
        compare_same(6, ori.multi1, mut.multi1)) return true;
}

void Running_interface::Output_Difference(const Matrix& ori, const Matrix& mut) 
{
    compare(0, ori.outer, mut.outer);
    compare(1, ori.stem, mut.stem);
    compare(2, ori.stemend, mut.stemend);
    compare(3, ori.multi, mut.multi);
    compare(4, ori.multi2, mut.multi2);
    compare(5, ori.multibif, mut.multibif);
    compare(6, ori.multi1, mut.multi1);
}

void Running_interface::Check_Mutation(string sequence) 
{
    Vec bpp;
    const char* base = "ACGU";        
    int length = sequence.length();
    Radiam model1;
    model1.Set_Constraint(constraint, length);
    model1.Get_ori_matrix(sequence);
    Init_Files();    
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < 4; j++) {
            if (sequence[i] != base[j]) {
                string str = sequence;
                str[i] = base[j];
                model1.Mutation_calculation(i, str);
                Rfold_Lang model2;                
                model2.calculation(constraint, str);
                if (Check_Difference(model1.alpha, model2.alpha)) {
                    cout << "Inside Differ!" << endl;
                    model1.Print_Mat(true);
                    model2.Print_Mat(true);
                    Output_Difference(model1.alpha, model2.alpha);
                } else if (Check_Difference(model1.beta, model2.beta)) {
                    cout << "Outside Differ!" << endl;
                    model1.Print_Mat(false);
                    model2.Print_Mat(false);
                    Output_Difference(model1.beta, model2.beta);
                }
            }
        }
    }
}

void Running_interface::Raw_compare_BPP_Rfold_Model(string str, bool first) 
{
    const char* base = "ACGU";
    int length = str.length();
    Rfold_Lang model1, model2;
    model1.calculation(constraint, str);
    if (first) Init_Files();
    int i = length/2;
    for (int j = 0; j < 4; j++) {
        if (str[i] != base[j]) {
            string temp = str;
            temp[i] = base[j];
            model2.calculation(constraint, temp);
            Same_Header("> alpha", i, sequence[i], base[j], str);
            Output_Difference(model1.alpha, model2.alpha);
            Same_Header("> beta", i, sequence[i], base[j], str);    
            Output_Difference(model1.beta, model2.beta);
        }
    }
}

void Running_interface::Run_BPP_Rfold_Model(string str)
{
    int band = constraint;
    if (str.length() < band*2) return;
    Rfold_Lang model;
    model.calculation(constraint, str);
    cout << "# all" << endl;
    model.Write_bpp_part(false, str.length()/2-band+1, str.length()/2+band);
    for (int i = band; i < str.length()/2 && i*2 <= constraint*10; i += band) {
        string seq = str.substr(str.length()/2-i, i*2);
        model.calculation(constraint, seq);
        cout << "# " << i << endl;
        model.Write_bpp_part(false, seq.length()/2-band+1, seq.length()/2+band);
    }
}

void Running_interface::Set_consensus_index(Radiam& radiam, string sequence, string bp_seq, vector<int>& index)
{
    if (bp_seq.length() != sequence.length()) return;
    for (int i = (int)sequence.length()-1; i >= 0; i--) {
        if (sequence[i] == '.') {
            sequence.erase(i, 1);
            bp_seq.erase(i, 1);
        } else index.push_back(i);
    }
    reverse(index.begin(), index.end());
    //Print_Vec(index, true);
    RNA_transform(sequence);
    radiam.Set_annotate(bp_seq, sequence, index);
}


void Running_interface::Run_Radiam(Arg& arg, vector<string>& seqlist, vector<string>& namelist, bool bpp)
{
    vector<string>::iterator it = find(namelist.begin(), namelist.end(), "SS_cons");
    if (it == namelist.end()) return;
    int cons = (int)distance(namelist.begin(), it);
    for (string::iterator it1 = seqlist[cons].begin(); it1 != seqlist[cons].end(); it1++)
        if (isalpha(*it1)) return;
    cout << seqlist[cons] << endl;
    for (int i = 0; i < cons; i++) {
        vector<int> index;
        Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, 0, i == 0, arg.dist);
        radiam.id = i;
        radiam.Set_region(arg.rightw, arg.  leftw);
        Set_consensus_index(radiam, seqlist[i], seqlist[cons], index);
        if (radiam.ori_seq.length() == 0) {
            cout << "error" << endl;
            continue;
        }
        if (bpp) {
            radiam.Correlation_of_bpp(arg.mtype, 1, constraint, radiam.ori_seq);
        } else {
            radiam.Set_bpp_binary(arg.mtype, constraint, radiam.ori_seq);
            radiam.Write_bpp_binary();
        }
    }
}

void Running_interface::Run_Radiam(Arg& arg, vector<string>& seqlist, vector<string>& namelist, string& bp_seq, bool bpp)
{
    cout << bp_seq << endl;
    for (int i = 0; i < (int)seqlist.size(); i++) {
        vector<int> index;
        Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, 0, i == 0, arg.dist);
        radiam.id = i;
        radiam.Set_region(arg.rightw, arg.  leftw);
        Set_consensus_index(radiam, seqlist[i], bp_seq, index);
        if (radiam.ori_seq.length() == 0) {
            cout << "error" << endl;
            continue;
        }
        if (bpp) {
            radiam.Correlation_of_bpp(arg.mtype, 1, constraint, radiam.ori_seq);
        } else {
            radiam.Set_bpp_binary(arg.mtype, constraint, radiam.ori_seq);
            radiam.Write_bpp_binary();
        }
    }
}

void Running_interface::Run_Radiam(Arg& arg, bool init, int rightw, int leftw)
{
    if (arg.seq.length() == 0) return;
    Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, 0, init, arg.dist);
    RNA_transform(arg.seq);
    radiam.Set_region(rightw, leftw);
    radiam.Correlation_of_bpp(arg.mtype, 1, constraint, arg.seq);
}

void Running_interface::Run_Radiam(Arg& arg, bool init)
{
    if (arg.seq.length() == 0) return;    
    Radiam radiam(arg.acc_flag, arg.threshold, arg.window, init, arg.dist);
    RNA_transform(arg.seq);
    if (arg.example) radiam.Set_example();
    if (radiam._debug == Radiam::Debug::Outer || radiam._debug == Radiam::Debug::Bpp) {
        vector<int> mlist = vector<int>(1, (int)(arg.seq.length())/2);
        radiam.Correlation_of_bpp(arg.mtype, mlist, constraint, arg.seq);
    } else {
        radiam.Correlation_of_bpp(arg.mtype, 1, constraint, arg.seq);
    }
}

}