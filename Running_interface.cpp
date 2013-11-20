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

void Running_interface::Run_BPP_Rfold_Model(string str)
{
    int band = constraint;
    if ((int)str.length() < band*2) return;
    Rfold_Lang model;
    model.calculation(constraint, str);
    cout << "# all" << endl;
    model.Write_bpp_part(false, str.length()/2-band+1, str.length()/2+band);
    for (int i = band; i < (int)str.length()/2 && i*2 <= constraint*10; i += band) {
        string seq = str.substr(str.length()/2-i, i*2);
        model.calculation(constraint, seq);
        cout << "# " << i << endl;
        model.Write_bpp_part(false, seq.length()/2-band+1, seq.length()/2+band);
    }
}

void Running_interface::Run_Rfold_Model(string str, bool acc)
{
    if ((int)str.length() == 0) return;
    cout << "#" << str << endl;    
    Rfold_Lang model;
    model.calculation(constraint, str);
    model.Write_bpp(acc);
}

void Running_interface::Set_consensus_index(Radiam& radiam, string sequence)
{
    vector<int> index;
    int all_length = (int)sequence.length();
    for (int i = (int)sequence.length()-1; i >= 0; i--) {
        if (sequence[i] == '.' || sequence[i] == '-') {
            sequence.erase(i, 1);
        } else index.push_back(i);
    }
    reverse(index.begin(), index.end());
    RNA_transform(sequence);
    radiam.Set_annotate(string(""), sequence, index, all_length);
}

void Running_interface::Set_consensus_index(Radiam& radiam, string sequence, string bp_seq)
{
    vector<int> index;
    int all_length = (int)sequence.length();
    if (bp_seq.length() != sequence.length()) return;
    for (int i = (int)sequence.length()-1; i >= 0; i--) {
        if (sequence[i] == '.' || sequence[i] == '-') {
            sequence.erase(i, 1);
            bp_seq.erase(i, 1);
        } else index.push_back(i);
    }
    reverse(index.begin(), index.end());
    RNA_transform(sequence);
    radiam.Set_annotate(bp_seq, sequence, index, all_length);
}


void Running_interface::Run_Radiam(Arg& arg, vector<string>& seqlist, vector<string>& namelist, bool not_bpp)
{
    vector<string>::iterator it = find(namelist.begin(), namelist.end(), "SS_cons");
    if (it == namelist.end()) return;
    int cons = (int)distance(namelist.begin(), it);
    for (string::iterator it1 = seqlist[cons].begin(); it1 != seqlist[cons].end(); it1++)
        if (isalpha(*it1)) return;
    cout << seqlist[cons] << endl;
    for (int i = 0; i < cons; i++) {
        Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, 0, i == 0, arg.dist);
        radiam.id = i;
        radiam.Set_region(arg.rightw, arg.  leftw);
        Set_consensus_index(radiam, seqlist[i], seqlist[cons]);
        if (radiam.ori_seq.length() == 0) {
            cout << "error" << endl;
            continue;
        }
        if (not_bpp) {
            radiam.Correlation_of_bpp(arg.mtype, arg.num, constraint, radiam.ori_seq);
        } else {
            radiam.Set_bpp_binary(arg.mtype, constraint, radiam.ori_seq);
            radiam.Write_bpp_binary();
        }
    }
}

void Running_interface::Run_Radiam(Arg& arg, vector<string>& seqlist, vector<string>& namelist, string& bp_seq, bool not_bpp)
{
    cout << bp_seq << endl;
    for (int i = 0; i < (int)seqlist.size(); i++) {
        Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, 0, i == 0, arg.dist);
        radiam.id = i;
        radiam.Set_region(arg.rightw, arg.  leftw);
        Set_consensus_index(radiam, seqlist[i], bp_seq);
        if (radiam.ori_seq.length() == 0) {
            cout << "error" << endl;
            continue;
        }
        if (not_bpp) {
            radiam.Correlation_of_bpp(arg.mtype, arg.num, constraint, radiam.ori_seq);
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
    radiam.Correlation_of_bpp(arg.mtype, arg.num, constraint, arg.seq);
}

void Running_interface::Run_Radiam(Arg& arg, bool init, int id)
{
    if (arg.seq.length() == 0) return;
    Radiam radiam(arg.filename, arg.acc_flag, arg.threshold, arg.window, init, arg.dist);
    if (arg.example) radiam.Set_example();
    if (radiam._debug == Radiam::Debug::Outer || radiam._debug == Radiam::Debug::Bpp) {
        if (arg.num == 1) {
            RNA_transform(arg.seq);                    
            vector<int> mlist = vector<int>(1, (int)(arg.seq.length())/2);
            radiam.Correlation_of_bpp(arg.mtype, mlist, constraint, arg.seq);
            return;
        }
    } 
    radiam.id = id;        
    Set_consensus_index(radiam, arg.seq);        
    radiam.Correlation_of_bpp(arg.mtype, arg.num, constraint, radiam.ori_seq);
}

void Running_interface::Run_Radiam_for_mRNA(Arg& arg, string& filename, vector<bool>& flag)
{
    if (arg.seq.length() == 0) return;
    Radiam radiam(filename, arg.acc_flag, arg.threshold, 0, true, arg.dist);
    RNA_transform(arg.seq);
    if (arg.example) radiam.Set_example();
    if (arg.num == 1) {
        radiam.id = 0;
        radiam.SNV_analysis(arg.mtype, arg.num, constraint, arg.seq, flag);
    }
}


}