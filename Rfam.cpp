#include "file_explorer.h"

using namespace std;



void Rfam::Get_name_seq(string& str)
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

string Rfam::Get_filename(string& name, int num, string& file)
{
    stringstream ss;
    if (file != "") {
	    ss << file;
	    if (file[(int)file.length()-1] != '/') ss << "/";
	} else ss << "/home/cawatchm/";
	//else ss << "../correlation/";	
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

void Rfam::Get_data_string(int& size)
{
    string str, temp;
    getline(ifs, str); getline(ifs, str);
    istringstream iss(str);
    iss >> temp >> size;
    getline(ifs, str);
}

bool Rfam::Get_seq_string(int size)
{
    init_seq();
	string str;
    for (int j = 0; j < size && getline(ifs, str); j++) 
        Get_name_seq(str);
    assert(seqlist.size() == namelist.size() && (int)seqlist.size() >= 1);    
    for (string::iterator it1 = bp_seq.begin(); it1 != bp_seq.end(); it1++)
	    if (isalpha(*it1)) return false;
    return true;  
}

int Rfam::Scan_all_volume(string& filename)
{
    int tsize = 0;
    string str;
    Words v;
    ifstream tifs(filename.c_str());
    for (int i = 0; getline(tifs, str); i++) {
        Parse(str, v);        
        if (v.size() < 2 || v[0] != "size") continue;
        tsize += atoi(v[1].c_str());
    }
    return tsize;
}

bool Rfam::Check_pararell(int target, int array)
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

bool Rfam::is_Long(string& AC, Words& seqlist)
{
    if (AC.substr(0,5) == "full_") return true;
    else if (seqlist.size() > TOOMUCH) return true;
    else if (seqlist[0].length() > TOOLONG) return true;
    else return false;
}

bool Rfam::Read_one_family(string& str, string& AC, int& i)
{
    int size = 0;    
    Words v;    
    Parse(str, v);
    if ((int)v.size() < 2 || v[0] != "id:") {
        i--;
        return false;
    }
    AC = v[1];
    Get_data_string(size);
    return Get_seq_string(size);
}

bool Rfam::Already_exist(string file)
{
    FILE  *fp;    
    size_t index = file.find_last_of('/');
    if (index == string::npos)  index = 0;
    else index++;
    cout << "index: " << file.substr(index) << " " << index << endl;
    if ((fp = fopen(("/home/cawatchm/Radiam/file/mutation_full/"+file.substr(index)).c_str(), "r")) == NULL)
        return false;
    fclose(fp);
    return true;
}

void Rfam::Run_single_calc(Rfold::Running_interface& intf)
{
    string str;
    string file = arg.filename;    
    if (arg.example) return;
    for (int i = 0; getline(ifs, str); i++) {
        string AC;
        if (!Read_one_family(str, AC, i))
            continue;
        if (is_Long(AC, seqlist)) {
            cout << "Too much... or Too long..." << endl;
        } else {
            for (int j = 0; j < 3; j++) {
                arg.mtype = j;
                cout << ">\t" << AC << endl;
                arg.filename = Get_filename(AC, 0, file);
                if (Already_exist(arg.filename))    continue;                    
                cout << arg.filename << endl;
                intf.Run_Radiam(arg, seqlist, namelist, bp_seq, !arg.example);
            }
        }
     }
}

void Rfam::Run_array_calc(int array, Rfold::Running_interface& intf)
{
    string str;
	string file = arg.filename;    
    if (arg.example) return;
    for (int i = 0; getline(ifs, str); i++) {
        string AC;
        if (!Read_one_family(str, AC, i))
            continue;
        if (i%ARRAY == array-1) {
        //else if (Check_pararell(i, array)) {
		    if (is_Long(AC, seqlist)) {
		    	cout << "Too much... or Too long..." << endl;
		    } else {
                for (int j = 0; j < 3; j++) {
                    arg.mtype = j;
                    cout << ">\t" << AC << endl;
                    arg.filename = Get_filename(AC, 0, file);
                    if (Already_exist(arg.filename))    continue;                    
                    cout << arg.filename << endl;
                    intf.Run_Radiam(arg, seqlist, namelist, bp_seq, !arg.example);
                }
    		    //if (Check_pararell(i, array+1)) break;
            }
		} 
    }
}

void Rfam::Cut_seq_list(int array, int all)
{
	int start = (all >= (int)seqlist.size()) ? array-1 : (int)seqlist.size()*(array-1)/all;
	int end = (all >= (int)seqlist.size()) ? array : (int)seqlist.size()*(array)/all;
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

void Rfam::Run_too_much_calc(int array, Rfold::Running_interface& intf)
{
    string str;
	string file = arg.filename;    
    for (int i = 0; getline(ifs, str); i++) {
        string AC;
        if (!Read_one_family(str, AC, i))
            continue;
        if (AC != "full_RF00001" && AC != "full_RF00002" && AC != "full_RF00005")
            continue;
        //if (is_Long(AC, seqlist)) {
            cout << AC << endl;               
	    	if ((int)seqlist.size() < array) continue;
	    	Cut_seq_list(array, ARRAY);
            cout << "#num " << seqlist.size() << endl;
            for (int j = 0; j < 3; j++) {
                arg.mtype = j;
                arg.filename = Get_filename(AC, array, file);
                if (Already_exist(arg.filename))    continue;                    
                cout << arg.filename << endl;
                intf.Run_Radiam(arg, seqlist, namelist, bp_seq, !arg.example);
            }
		//}
    }
}

void Rfam::run(Rfold::Running_interface& intf, string& filename, int array)
{
    ifs.open(filename.c_str());
    if (arg.longer) {
        if (arg.example) return;        
        else if (array > 0) Run_too_much_calc(array, intf);
    } else {
        if (array > 0)
            Run_array_calc(array, intf);
        else
            Run_single_calc(intf);
    }
}

void Rfam::Read_Rfam_file(Rfold::Arg& targ, string& filename, int array)
{
    Rfam rfam(targ);
    int constraint = (targ.constraint <= 0) ? (int)(targ.seq.length()) : targ.constraint;        
    Rfold::Running_interface intf(constraint, targ.seq);
    rfam.run(intf, filename, array);
}
