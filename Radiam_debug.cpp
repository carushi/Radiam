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
    sleep(1);
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

void Radiam::Write_bppm_dif(const Mat& bpp_mut, const Mat& bpp_ori)
{
    double tmax = 0.0, t1 = 0.0, t2 = 0.0;
    for (int i = 0; i < (int)bpp_mut.size(); i++) {
        for (int j = i+1; j < (int)bpp_mut[i].size(); j++) {
            if (rdebug) cout << i << "," << j << ": " << bpp_mut[i][j] << endl;
            double diff = fabs(bpp_ori[i][j]-bpp_mut[i][j]);
            if (tmax < diff) { tmax = diff; t1 = bpp_ori[i][j]; t2 = bpp_mut[i][j]; }
        }
    }
    cout << "* same_bpp_diff_max " << tmax << " " << t1 << " " << t2 << endl;
}

void Radiam::Debug_bppm(int type, string& sequence)
{
    Rfold_Lang model2;
    model2.calculation(_constraint, sequence);
    Mat bpp_mut, bpp_ori;
    Write_bpp(bpp_mut);
    model2.Write_bpp(bpp_ori);
    Write_bppm_dif(bpp_mut, bpp_ori);
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
        if (log10(value)-log10(*it2) >= -DEF_PRE) {
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
