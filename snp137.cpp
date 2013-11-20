#include "file_explorer.h"

using namespace std;

char upper(char c) 
{
    c = toupper(c);
    if (c == 'U') c = 'T';
    else if (!isbase(c)) c = 'N';
    return c;
}

char complementary(char c)
{
    c = upper(c);
    if (c == 'A' || c == 'N') return 'T';
    else if (c == 'C') return 'G';
    else if (c == 'G') return 'C';
    else if (c == 'T') return 'A';
    else return 'T';
}

void get_comp_seq(string& str) 
{
    string seq;
    seq = str;
    str = "";
    transform(seq.begin(), seq.end(), back_inserter(str), complementary);
    reverse(str.begin(), str.end());
}

string Transcript::get_plus_strand(const string& chr, int start, int end)
{
    int count = 1;
    string str, seq;
    string filename = chr_dir(chr);
    ifstream ifs(filename.c_str());
    getline(ifs, str);  str.clear();
    while (getline(ifs, str) && count+(int)str.length() < start) {
        count += (int)str.length();
    }
    while (count < end && str.length() > 0) {
        int tstart = max(count, start);
        int tend = min(count+(int)str.length(), end);
        seq += str.substr(tstart-count, tend-tstart+1);
        count += (int)str.length();
        if (!getline(ifs, str)) return seq;
    }
    return seq;
}

string Transcript::get_pre_seq(const string& chr, int start, int end, bool mstrand)
{
    string seq = get_plus_strand(chr, start+1, end);
    if (mstrand) get_comp_seq(seq);
    return seq;
}

void Transcript::get_mature_seq(const string& chr, Words& data, string& seq)
{
	;
}


int Transcript::indel_type(const string& type_data)
{
    if (type_data.find("insertion") != string::npos) {
        cout << "* mut_ins\t";
        return 0;
    } else if (type_data.find("deletion") != string::npos) {
        cout << "* mut_del\t";
        return 1;
    } else {
        cout << "* mut_sub\t";
        return 2;
    }
}

string Transcript::Get_filename(string& name, int num, string& file)
{
    stringstream ss;
    if (file != "") {
        ss << file;
        if (file[(int)file.length()-1] != '/') ss << "/";
    } else ss << "/home/cawatchm/";
    ss << name;
    if (num > 0) ss << "_" << num;
    ss << ((arg.acc_flag) ? "_acc_" : "_bpp_") << arg.constraint;
    if (arg.mtype == 1) ss << "_del";
    else if (arg.mtype == 0) ss << "_in";
    ss << ".txt";
    cout << ss.str() << endl;
    return ss.str();
}

void Transcript::fill_flag(vector<bool>& flag, int dist, int constraint) 
{
    for (int i = max(0, dist-constraint); i < min(dist+constraint, (int)flag.size()); i++) {
        flag[i] = true;
    }
}

int Transcript::dist(int pos, string& snp) 
{
    Words snp_data;
    Parse(snp, snp_data);
    return dist(pos, atoi(snp_data[0].c_str()));
}

void Transcript::Set_flags(Words& snp, Flags& flags, int constraint, int pos)
{
    for (Words::iterator it = snp.begin(); it != snp.end(); it++) {
        int type = indel_type(*it);
        cout << *it << "\t" << arg.seq[dist(pos, *it)] << "\t" << dist(pos, *it) << endl;
        fill_flag(flags[type], dist(pos, *it), constraint/DIV);            
    }
}

void Transcript::Set_all_flag(Flags& flags, int dist, int range)
{
    for (int i = 0; i < 3; i++) {    
        fill_flag(flags[i], dist, range);
    }
}

void Transcript::Set_splicing(Words& gene_data, Flags& flags, int constraint, int pos)
{
    Words spl_data;
    int tss = atoi(gene_data[6].c_str());
    int tes = atoi(gene_data[7].c_str());
    if (tss == tes) return;
    Set_all_flag(flags, dist(pos, tss), constraint/DIV);
    Set_all_flag(flags, dist(pos, tes-1), constraint/DIV);
    Parse_csv(gene_data[9], spl_data);
    for (Words::iterator it = spl_data.begin(); it != spl_data.end(); it++) {
        Set_all_flag(flags, dist(pos, atoi(it->c_str())), constraint/DIV);
    }
    Parse_csv(gene_data[10], spl_data);
    for (Words::iterator it = spl_data.begin(); it != spl_data.end(); it++) {
        Set_all_flag(flags, dist(pos, atoi(it->c_str())-1), constraint/DIV);
    }
}

void Transcript::Set_and_calc(Words& gene_data, Words& snp, int constraint, int pos)
{
    Flags flags = Flags(3, vector<bool>((int)arg.seq.length()));
    Set_flags(snp, flags, constraint, pos);
    Set_splicing(gene_data, flags, constraint, pos);
    Rfold::Running_interface intf(constraint, arg.seq);
    string name = gene_data[2]+"_"+gene_data[1];
    string filename = Get_filename(name, 0, arg.filename);
    intf.Run_Radiam_for_mRNA(arg, filename, flags[arg.mtype]);
}

void Transcript::Calculate_snp(string& gene, Words& snp)
{
    if (gene == "") return;
	Words gene_data;
	Parse(gene, gene_data);
    int start = atoi(gene_data[4].c_str());
    int end = atoi(gene_data[5].c_str());
	arg.seq = get_pre_seq(gene_data[2], start, end, gene_data[3] == "-");
    cout << gene << " " << arg.seq.length() << endl;
    if (arg.seq.length() > 11000) return;
    int constraint = (arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint;    
    if (gene_data[3] == "-")    
        Set_and_calc(gene_data, snp, constraint, end-1);
    else
        Set_and_calc(gene_data, snp, constraint, start);        
}

void Transcript::Read_snp_file(string& filename, Rfold::Arg& targ)
{
    Transcript trans(targ);
    string str;
    string gene;
    Words snp;
    ifstream ifs(filename.c_str());
    for (int i = 0; getline(ifs, str); i++) {
        if ((int)str.length() == 0) continue;
        else if (str[0] != '*') {
            snp.push_back(str);
            i--;
        } else {
            if (i > 0 && (i-1)%ARRAY == targ.array-1) {
                trans.Calculate_snp(gene, snp);
            }
            refresh(gene, snp);
            gene = str;
        }
    }
    trans.Calculate_snp(gene, snp);
}
