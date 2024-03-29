
#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "energy_const.h"
#include "pair_mat.h"
#include "energy_par.h"

namespace Rfold {

using std::vector;
using std::deque;
using std::min;
using std::max;
using std::transform;
using std::ifstream;
using std::istringstream;
using std::ios;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::back_inserter;
using std::distance;
using std::swap;

#define DOUBLE double

typedef vector<DOUBLE> Vec;
typedef vector<Vec> Mat;


struct Base {
    int operator()(char c) {
        switch(c) {
            case 'A': case 'a': return 1;
            case 'C': case 'c': return 2;
            case 'G': case 'g': return 3;
            case 'T': case 'U': case 't': case 'u': return 4;
            default: return 1;
        }
    }
};


inline static bool Is_INF(const DOUBLE value) {
        return (value <= -INF);    
}
inline static DOUBLE Logsumexp(DOUBLE x, DOUBLE y) {
    if (x == -INF) return y;
    else if (y == -INF) return x;
    else if (x > y) return ((x + log(exp(y-x) + 1.0)));
    else return ((y + log(exp(x-y) + 1.0)));
}
inline static DOUBLE Logsum(DOUBLE a, DOUBLE b) {
    if (Is_INF(a) || Is_INF(b)) return -INF;
    else return a+b;
}
inline static DOUBLE Logsum(DOUBLE a, DOUBLE b, DOUBLE c) {
    return Logsum(a, Logsum(b, c));
}
inline static DOUBLE Logsum(DOUBLE a, DOUBLE b, DOUBLE c, DOUBLE d) {
    return Logsum(Logsum(a, b), Logsum(c, d));
}
inline static bool Can_bind(int type) {
    return (type > 0);    
}
inline static int bp(int i, int j, const vector<int>& _sequence) {
    return BP_pair[_sequence[i]][_sequence[j]];    
}
inline static int rbp(int i, int j, const vector<int>& _sequence) {
    return rtype[bp(i, j, _sequence)];    
}

class Sequence {
    friend class Matrix;
public:
    int length;    
    string str;
    vector<int> sequence;

    Sequence() {}
    Sequence(const string& str, const vector<int>& sequence, int length)
             : length(length), str(str), sequence(sequence) {}
    ~Sequence() {}
};
template<class Data>
void Print_Vec(const vector<Data>& elem, bool end = true)
{
    ostream_iterator<Data> out_it(cout, ", ");
    if ((int)elem.size() > 1) 
        copy(elem.begin(), elem.begin()+(int)elem.size()-1, out_it);    
    if ((int)elem.size() > 0)
        cout << *(elem.begin()+(int)elem.size()-1);
    if (end) cout << endl;
}

class Matrix {
public:
    int _length; // width
    int _constraint; // height 
    bool _inside;
	Vec outer;
	Mat stem;
	Mat stemend;
    Mat multi;
	Mat multi1;
	Mat multi2;
    Mat multibif;    
    Matrix() {}
    Matrix(int length, int constraint, bool inside)
     : _length(length), _constraint(constraint), _inside(inside) {
        Initialize();
        outer = Vec(_length+1, 0.0);
    }
    virtual ~Matrix(){}
    void operator=(const Matrix& right) {
        outer = right.outer;
        stem = right.stem;
        stemend = right.stemend;
        multi = right.multi;
        multi1 = right.multi1;
        multi2 = right.multi2;
        multibif = right.multibif;
        _length = right._length;
        _constraint = right._constraint;
        _inside = right._inside;
    }
    void Initialize(); 
    static void Print_Mat(const Mat&, const string&);
    void Print(const string&);
    int index(int j) const { return j; }
    void erase(int);
    void insert(int);

};


class Matrix_reduced : Matrix {
public:
    int _endp;
    int _endi;
    static const int _max = 2;
    Matrix_reduced() {}
    Matrix_reduced(int length, int constraint, bool inside) : Matrix(length, constraint, inside) {}
    virtual ~Matrix_reduced(){}
    void Print(const string&);
    void Pop_front();
    void Push_front();
    void Initialize();
    int index(int j) const {
        int temp_index = (_inside) ? j-_endp+_endi : _endp-j+_endi;
        if (temp_index < _max*_constraint) return temp_index;
        else return temp_index-_max*_constraint;
    }
};
}
#endif
