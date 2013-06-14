#include "Running_interface.h"

#define ARRAY (15)
#define TOOMUCH (700)
#define TOOLONG (500)

using namespace std;

void Parse(string& str, vector<string>& v)
{
    v.clear();
    istringstream iss(str);
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(v));
}

void Get_name_seq(string& str, vector<string>& seqlist, vector<string>& namelist, string& bp_seq)
{
    istringstream iss(str);
    int num;
    string name; string seq;         
    iss >> num >> name >> seq;
    if (name == "SS_cons") {
    	bp_seq = seq;
    } else if (name != "RF") {
	    seqlist.push_back(seq); namelist.push_back(name);
	}
}

string Get_filename(string& name, int num, Rfold::Arg& arg, string& file)
{
    stringstream ss;
    if (file != "") {
	    ss << file;
	    if (file[(int)file.length()-1] != '/') ss << "/";
	} else ss << "/gluster/cawatchm/correlation";
	//else ss << "../correlation";	
    if (arg.example) ss << "original_";
    ss << name;
    if (num > 0) ss << "_" << num;
    ss << ((arg.acc_flag) ? "_acc_" : "_bpp_") << arg.constraint;
    if (!arg.example) {
        if (arg.mtype == 1) ss << "_del";
        else if (arg.mtype == 0) ss << "_in";
    }
    if (arg.dist)
        ss << "_cor";
    ss << ".txt";
    cout << ss.str() << endl;
    return ss.str();
}

void Get_data_string(ifstream& ifs, int& size, Rfold::Arg& arg)
{
    string str, temp;
    getline(ifs, str); getline(ifs, str);
    istringstream iss(str);
    iss >> temp >> size;
    getline(ifs, str);
}

bool Get_seq_string(ifstream& ifs,
					int size,
					vector<string>& seqlist,
					vector<string>& namelist,
					string& bp_seq) 
{
	string str;
    for (int j = 0; j < size && getline(ifs, str); j++) 
        Get_name_seq(str, seqlist, namelist, bp_seq);
    assert(seqlist.size() == namelist.size() && (int)seqlist.size() >= 1);    
    for (string::iterator it1 = bp_seq.begin(); it1 != bp_seq.end(); it1++)
	    if (isalpha(*it1)) return false;
    return true;  
}

int Scan_all_volume(string& filename)
{
    int tsize = 0;
    string str;
    vector<string> v;
    ifstream ifs(filename.c_str());
    for (int i = 0; getline(ifs, str); i++) {
        Parse(str, v);        
        if (v.size() < 2 || v[0] != "size") continue;
        tsize += atoi(v[1].c_str());
    }
    return tsize;
}

bool Check_pararell(int target, int array)
{
    int max_list[] = {-1, 10, 20, 40, 60, 80, 100, 120, 160, 220, 260, 280, 300, 320, 370};
    if (array <= 0 || array > 15) return true;
    else if (target <= max_list[array-1]) return false;
    else {
        if (array > 14) return true;
        else if (target > max_list[array]) return false;
        else return true;
    }
}

void Run_array_calc(ifstream& ifs, int array, Rfold::Arg& arg, Rfold::Running_interface& intf)
{
    string str;
	string file = arg.filename;    
    vector<string> v;
    bool end = false;
    for (int i = 0; getline(ifs, str); i++) {
        Parse(str, v);
        if ((int)v.size() < 2 || v[0] != "id:") {
        	i--; continue;
        }
	    int size = 0;
	    string bp_seq;
        vector<string> seqlist, namelist;
        Get_data_string(ifs, size, arg);
        if (!Get_seq_string(ifs, size, seqlist, namelist, bp_seq)) continue;
        if (Check_pararell(i, array)) {
		    cout << ">\t" << v[1] << endl;
		    if (array > 0 && !arg.example && (seqlist.size() > TOOMUCH || seqlist[0].length() > TOOLONG)) {
		    	cout << "Too much... or Too long..." << endl; continue;
		    }
		    arg.filename = Get_filename(v[1], 0, arg, file);
		    intf.Run_Radiam(arg, seqlist, namelist, bp_seq, !arg.example);
		    if (array > 0 && Check_pararell(i, array+1)) break;
		} 
    }
}

void Cut_seq_list(int array, int all, vector<string>& seqlist, vector<string>& namelist)
{
	int start = (int)seqlist.size()*(array-1)/all;
	int end = (int)seqlist.size()*(array)/all;
	int size = (int)seqlist.size();
	if (start > 0) {
		seqlist.erase(seqlist.begin(), seqlist.begin()+start);
		namelist.erase(namelist.begin(), namelist.begin()+start);
	}
	if (end < size) {
		seqlist.erase(seqlist.begin()+(end-start), seqlist.end());
		namelist.erase(namelist.begin()+(end-start), namelist.end());
	}
}

void Run_too_much_calc(ifstream& ifs, 
					   int array,
					   Rfold::Arg& arg,
					   Rfold::Running_interface& intf)
{
    string str;
	string file = arg.filename;    
    vector<string> v;
    bool end = false;
    for (int i = 0; getline(ifs, str); i++) {
        Parse(str, v);
        if ((int)v.size() < 2 || v[0] != "id:") {
        	i--; continue;
        }
	    int size = 0;
	    string bp_seq;
        vector<string> seqlist, namelist;
        Get_data_string(ifs, size, arg);
        if (!Get_seq_string(ifs, size, seqlist, namelist, bp_seq)) continue;
	    cout << ">\t" << v[1] << endl;   
	    if (!arg.example && (seqlist.size() > TOOMUCH || seqlist[0].length() > TOOLONG)) {
	    	if ((int)seqlist.size() < array) continue;
	    	Cut_seq_list(array, ARRAY, seqlist, namelist);
		    arg.filename = Get_filename(v[1], array, arg, file);
		    intf.Run_Radiam(arg, seqlist, namelist, bp_seq, !arg.example);
		}
    }
}

void Read_Rfam_file(Rfold::Arg& arg, string& filename, int array)
{
    int constraint = (arg.constraint <= 0) ? (int)(arg.seq.length()) : arg.constraint;        
    ifstream ifs(filename.c_str());
    Rfold::Running_interface intf(constraint, arg.seq);
    if (arg.longer) {
    	if (array > 0) Run_too_much_calc(ifs, array, arg, intf);
    } else
		Run_array_calc(ifs, array, arg, intf);
}

