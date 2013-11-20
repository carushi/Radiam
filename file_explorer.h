#ifndef _FILE_EXPLORE_H
#define _FILE_EXPLORE_H
#define DIV (10)
#include "Running_interface.h"
typedef vector<string> Words;
typedef vector<vector<bool> > Flags;


static void Parse(string& str, Words& v)
{
    v.clear();
    istringstream iss(str);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(v));
}

static void Parse_csv(string& str, Words& v)
{
    v.clear();
    string temp;
    istringstream iss(str);
    while (getline(iss, temp, ',')) {
        v.push_back(temp);
    }
}

char upper(char);
static char complementary(char); 
static bool isbase(char a) {
   return (a == 'A' || a == 'C' || a == 'G' || a == 'T' || a == 'U');
}
static void get_comp_seq(string&);

class Rfam
{
private:
    string bp_seq;
    Words seqlist;
    Words namelist;
    Rfold::Arg arg;
    ifstream ifs;
    void Get_name_seq(string&);
    string Get_filename(string&, int, string&);
    void Get_data_string(int&);
    bool Get_seq_string(int);
    int Scan_all_volume(string&);
    bool Check_pararell(int, int);
    bool Read_one_family(string&, string&, int&);
    bool Already_exist(string);
    void Run_single_calc(Rfold::Running_interface&);
    void Run_array_calc(int, Rfold::Running_interface&);
    void Cut_seq_list(int, int);
    void Run_too_much_calc(int, Rfold::Running_interface&);
    void init_seq() {
        bp_seq = "";
        seqlist.clear();
        namelist.clear();
    }
public:
    Rfam(Rfold::Arg targ) : arg(targ) {}
    virtual ~Rfam() {}
    static bool is_Long(string&, Words&);    
    void run(Rfold::Running_interface&, string&, int);
    static void Read_Rfam_file(Rfold::Arg&, string&, int);
};


class Transcript
{
private:
    Rfold::Arg arg;
    string get_plus_strand(const string&, int, int);
    string get_pre_seq(const string&, int, int, bool);
    void get_mature_seq(const string&, Words&, string&);
    int indel_type(const string&);
    int dist(int, string&);
    string Get_filename(string&, int, string&);
    void Set_flags(Words&, Flags&, int, int);
    void Set_all_flag(Flags&, int, int);   
    void Set_splicing(Words&, Flags&, int, int);    
    void Set_and_calc(Words&, Words&, int, int);   
    void Calculate_snp(string&, Words&);
    void fill_flag(vector<bool>&, int, int);    
    string chr_dir(const string& chr) {
        return "/home/cawatchm/chr/"+chr+".fa";
    }
    static void refresh(string& gene, Words& snp) {
        gene = "";
        snp.clear();
    }
    int dist(int pos, int tpos) {
        return abs(pos-tpos);
    }


public:
    Transcript(Rfold::Arg targ) : arg(targ) {}
    virtual ~Transcript(){}
    static void Read_snp_file(string&, Rfold::Arg&);


};



#endif