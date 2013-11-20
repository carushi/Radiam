#include "matrix.h"

namespace Rfold {

void Matrix::Initialize() {
    stem = stemend = multi = multibif = multi1 
    = multi2 = Mat(_length+1, Vec(_constraint+2, -INF));    
}

void Matrix::Print_Mat(const Mat& mat, const string& str) 
{
    for (int i = 0; i < (int)str.length(); i++)
        cout << str[i] << " ";
    cout << endl;
    for (Mat::const_iterator it = mat.begin(); it != mat.end(); it++)  
        Print_Vec(*it, true);
}

void Matrix::Print(const string& str)
{
    cout << "---stem" << endl;
    Print_Mat(stem, str);
    cout << "---stemend" << endl; 
    Print_Mat(stemend, str);        
    cout << "---multi" << endl;
    Print_Mat(multi, str);        
    cout << "---multiBif" << endl;
    Print_Mat(multibif, str);        
    cout << "---multi1" << endl;
    Print_Mat(multi1, str);        
    cout << "---multi2" << endl;
    Print_Mat(multi2, str);        
    cout << "---outer" << endl;
    Print_Vec(outer, true);

}

void Matrix::erase(int j) 
{
    stem.erase(stem.begin()+j);
    stemend.erase(stemend.begin()+j);
    multi.erase(multi.begin()+j);
    multibif.erase(multibif.begin()+j);
    multi1.erase(multi1.begin()+j);
    multi2.erase(multi2.begin()+j);
    outer.erase(outer.begin()+j);
}

void Matrix::insert(int j) 
{
    stem.insert(stem.begin()+j, Vec(_constraint+2, -INF));
    stemend.insert(stemend.begin()+j, Vec(_constraint+2, -INF));
    multi.insert(multi.begin()+j, Vec(_constraint+2, -INF));
    multibif.insert(multibif.begin()+j, Vec(_constraint+2, -INF));
    multi1.insert(multi1.begin()+j, Vec(_constraint+2, -INF));
    multi2.insert(multi2.begin()+j, Vec(_constraint+2, -INF));
    outer.insert(outer.begin()+j, 0.0);
}

void Matrix_reduced::Initialize() 
{
    stem = stemend = multi = multibif = multi1 
    = multi2 = Mat(_constraint*2, Vec(_constraint+2, -INF));
}

void Matrix_reduced::Pop_front() {
    //_start++;
}
void Matrix_reduced::Push_front() {
    //_start--;
}

}