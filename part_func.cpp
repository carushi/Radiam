#include "part_func.h"

namespace Rfold {

double Rfold_Lang::Calc_in_stem(int i, int j)
{
    int bp1 = bp(i+1, j, seq.sequence);
    double temp = -INF;
	if (Can_bind(bp1)) {
        if (debug) cout << seq.str[i+1] << "-" << seq.str[j] << endl;
        temp = Logsumexp(temp, Logsum(alpha.stem[j-1][j-1-(i+1)], Parameter::LogLoopEnergy(i+1, j, i+2, j-1, seq)));
        temp = Logsumexp(temp, alpha.stemend[j-1][j-1-(i+1)]);
  	}
    return temp;
}

double Rfold_Lang::Calc_in_multiBif(int i, int j)
{
    double temp = -INF;
	for (int k = i+1; k < j; k++) 
        temp = Logsumexp(temp, Logsum(alpha.multi1[k][k-i], alpha.multi2[j][j-k]));
	return temp;
}

double Rfold_Lang::Calc_in_multi2(int i, int j) 
{
    int type = bp(i+1, j, seq.sequence);
	double temp = Logsum(alpha.stem[j][j-i], Parameter::Sum_Dangle(type, i, j+1, seq), Parameter::logMLintern);
    temp = Logsumexp(temp, alpha.multi2[j-1][j-1-i]+Parameter::logML_BASE());
    return temp;
}

double Rfold_Lang::Calc_in_multi1(int i, int j) {
    return Logsumexp(alpha.multi2[j][j-i], alpha.multibif[j][j-i]);
}

double Rfold_Lang::Calc_in_multi(int i, int j) { 
    double temp = Logsumexp(alpha.multibif[j][j-i], 
                            Logsum(alpha.multi[j][j-(i+1)], Parameter::logML_BASE()));
    return temp;
}

double Rfold_Lang::Calc_in_stemend(int i, int j) 
{
    int p = i, q = j+1;
    double temp = -INF;
    if (q <= seq.length && Can_bind(bp(p, q, seq.sequence))) {
        temp = Parameter::LogHairpinEnergy(p, q, seq);
        for (int ip = i; ip < j-TURN-1; ip++) {
            for (int jp = max(j-MAXLOOP+(ip-i), ip+TURN+2); jp <= j && (j-jp)+(ip-i) > 0; jp++) {
                if (debug) cout << ip << " " << jp << endl;
                if (Can_bind(bp(ip+1, jp, seq.sequence))) 
                    temp = Logsumexp(temp, Logsum(alpha.stem[jp][jp-ip], Parameter::LogLoopEnergy(p, q, ip+1, jp, seq)));   
            }
        }
        temp = Logsumexp(temp, Logsum(alpha.multi[j][j-i], Parameter::Sum_Dangle(rbp(p, q, seq.sequence), j, i+1, seq), 
                                      Parameter::logMLclosing+Parameter::logMLintern));
    }
    return temp;
}

void Rfold_Lang::Calc_in_outer(int j)
{
    if (j == 0) alpha.outer[0] = 0.0;
    else alpha.outer[j] = alpha.outer[j-1];
    for (int k = max(0, j-_constraint-1); k < j; k++) {
        double temp = Logsum(alpha.outer[k], alpha.stem[j][j-k], Parameter::Sum_Dangle(bp(k+1, j, seq.sequence), k, j+1, seq));
        alpha.outer[j] = Logsumexp(alpha.outer[j], temp);
    }
}

void Rfold_Lang::Calc_inside_mat(int i, int j) 
{
    alpha.stem[j][j-i] = Calc_in_stem(i, j);
    alpha.multibif[j][j-i] = Calc_in_multiBif(i, j);
    alpha.multi2[j][j-i] = Calc_in_multi2(i, j);
    alpha.multi1[j][j-i] = Calc_in_multi1(i, j);
    alpha.multi[j][j-i] = Calc_in_multi(i, j);
    alpha.stemend[j][j-i] = Calc_in_stemend(i, j);
}

void Rfold_Lang::Calc_inside()
{
    for (int j = TURN+1; j <= seq.length; j++) {
        if (debug) 
            cout << "-------------------\n-j " << j << endl;
        for (int i = j-TURN; i >= max(0, j-_constraint-1); i--) 
            Calc_inside_mat(i, j);
        Calc_in_outer(j);
    }
}

/* ///////////////////////////////////////////// */

void Rfold_Lang::Calc_out_outer()
{ 
    beta.outer[seq.length] = 0.0;
    for (int j = seq.length-1; j >= 0; j--) {
        beta.outer[j] = beta.outer[j+1];
        for (int k = j+TURN; k <= min(j+_constraint+1, seq.length); k++) {
            double temp = Logsum(alpha.stem[k][k-j], beta.outer[k], Parameter::Sum_Dangle(bp(j+1, k, seq.sequence), j, k+1, seq));
            beta.outer[j] = Logsumexp(beta.outer[j], temp);
        }
    }
}

double Rfold_Lang::Calc_out_stemend(int i, int j) {
    return (Is_range(i-1, j+1) ? beta.stem[j+1][j+1-(i-1)] : -INF);
}

double Rfold_Lang::Calc_out_multi(int i, int j) 
{
    double temp = -INF;
    if (Is_range(i-1, j)) {
        temp = Logsum(beta.multi[j][j-(i-1)], Parameter::logML_BASE());
        temp = Logsumexp(temp, Logsum(beta.stemend[j][j-i], Parameter::Sum_Dangle(rbp(i, j+1, seq.sequence), j, i+1, seq),
                                      Parameter::logMLclosing+Parameter::logMLintern));
    }
    return temp;
}

double Rfold_Lang::Calc_out_multi1(int i, int j) 
{
    double temp = -INF;
    for (int k = j+1 ; k <= min(seq.length, i+_constraint+1); k++) 
        temp = Logsumexp(temp, Logsum(beta.multibif[k][k-i], alpha.multi2[k][k-j]));
	return temp;
}

double Rfold_Lang::Calc_out_multi2(int i, int j) 
{
    double temp = beta.multi1[j][j-i];
    if (Is_range(i, j+1)) 
        temp = Logsumexp(temp, Logsum(beta.multi2[j+1][j+1-i], Parameter::logML_BASE()));
    for (int k = max(0, j-_constraint-1); k < i; k++) 
        temp = Logsumexp(temp, Logsum(beta.multibif[j][j-k], alpha.multi1[i][i-k]));
    return temp;
}

double Rfold_Lang::Calc_out_multiBif(int i, int j) {
    return Logsumexp(beta.multi1[j][j-i], beta.multi[j][j-i]);
}

double Rfold_Lang::Calc_out_stem(int i, int j) 
{
    int p = i+1, q = j, type = bp(p, q, seq.sequence);
    double temp = -INF;
    if (Can_bind(type)) {
        temp = Logsum(alpha.outer[i], beta.outer[j], Parameter::Sum_Dangle(type, i, j+1, seq));
        for (int ip = i; ip > 0; ip--) {
            for (int jp = min(ip+_constraint-1, min(j+MAXLOOP-(i-ip), seq.length-1)); jp >= j && (i-ip)+(jp-j) > 0; jp--) {
                if (Can_bind(bp(ip, jp+1, seq.sequence))) 
                    temp = Logsumexp(temp, Logsum(beta.stemend[jp][jp-ip], Parameter::LogLoopEnergy(ip, jp+1, p, q, seq)));
            }
        }
        temp = Logsumexp(temp, Logsum(beta.multi2[j][j-i], Parameter::Sum_Dangle(type, i, j+1, seq), Parameter::logMLintern));
        if (i >= 1 && j < seq.length && Is_range(i-1, j+1)) 
            temp = Logsumexp(temp, Logsum(beta.stem[j+1][j+1-(i-1)], Parameter::LogLoopEnergy(i, j+1, p, q, seq)));
    }
    return temp;
}

void Rfold_Lang::Calc_outside_mat(int i, int j)
{
    beta.stemend[j][j-i] = Calc_out_stemend(i, j);
    beta.multi[j][j-i] = Calc_out_multi(i, j);
    beta.multi1[j][j-i] = Calc_out_multi1(i, j);
    beta.multi2[j][j-i] = Calc_out_multi2(i, j);
    beta.multibif[j][j-i] = Calc_out_multiBif(i, j);
}

void Rfold_Lang::Calc_outside()
{
    Calc_out_outer();
    for (int j = seq.length; j >= TURN+1; j--) {
        if (debug) cout << "----------------\n-j " << j << endl;        
        for (int i = max(0, j-_constraint-1); i < j-TURN-1; i++) {        
            if (debug) cout << "--i " << i << endl;            
            if (i != 0 && j < seq.length) Calc_outside_mat(i, j);
            beta.stem[j][j-i] = Calc_out_stem(i, j);
        }
    }
    if (!noout) printf("* partition func %le\n", exp(alpha.outer[seq.length]));
    if (print) {
        Print_Mat(true);
        Print_Mat(false);
    }
}

double Rfold_Lang::hairpin_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = 0; i < x1; i++) 
        for (int j = x2; j <= min(seq.length, i+_constraint+1); j++) {
            value = Logsumexp(value, Logsum(beta.stem[j][j-i], Parameter::LogHairpinEnergy(i, j, seq)));
        }
    return value;
}

double Rfold_Lang::interior_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = 0; i < x1; i++) {
        for (int j = x2; j <= min(seq.length, i+_constraint+1); j++) {
            for (int k = x2; k < j; k++) 
                for (int l = k; l < j; l++)
                    value = Logsumexp(value, Logsum(beta.stemend[j][j-i], Parameter::LogLoopEnergy(i, k, l, j, seq), alpha.stem[l][l-k]));
            for (int k = x2; k < j; k++) 
                for (int l = k; l < j; l++)
                    value = Logsumexp(value, Logsum(beta.stemend[j][j-i], Parameter::LogLoopEnergy(i, k, l, j, seq), alpha.stem[l][l-k]));
        }
    }
    return value;
}

double Rfold_Lang::multi_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int j = x2; j <= min(seq.length, x1+_constraint+1); j++) {
        value = Logsumexp(value, Logsum(beta.multi[j][j-(x1-1)], (x2-x1+1)*Parameter::logML_BASE(), alpha.multi[j][j-x2]));
    }
    return value;
}

double Rfold_Lang::multi2_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = min(0, x2-_constraint-1); i <= x1; i++) {
        value = Logsumexp(value, Logsum(beta.multi2[x2][x2-i], (x2-x1+1)*Parameter::logML_BASE(), alpha.multi2[x1-1][x1-1-i]));
    }
    return value;
}

void Rfold_Lang::Calc_acc()
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);    
    for (int i = 1; i < seq.length; i++) {
        for (int j = max(2, i-_constraint); j <= min(seq.length, i+_constraint); j++) {
            double value = hairpin_acc(i, j);
            value = Logsumexp(value, interior_acc(i, j));
            value = Logsumexp(value, multi_acc(i, j));  
            value = Logsumexp(value, multi2_acc(i, j));
            value = Logsumexp(value, Logsum(beta.outer[j], alpha.outer[i-1]));
            if (!noout) cout << "* " << i << " " << j << " " << value << endl;        
        }
    }
}

double Rfold_Lang::bpp(int i, int j) 
{
    if (i > j) swap(i, j);
     double stack = Logsum(alpha.stem[j-1][j-1-i], beta.stem[j][j-(i-1)],
                           Parameter::LogLoopEnergy(i, j, i+1, j-1, seq));
    double stemend = Logsum(alpha.stemend[j-1][j-1-i], beta.stem[j][j-(i-1)]);
    return exp(Logsum(Logsumexp(stack, stemend), -alpha.outer[seq.length]));
 }

void Rfold_Lang::Write_bpp(Mat& data)
{
    data.clear();
    data = Mat(seq.length, Vec(_constraint, 0.0)); 
    for (int i = 1; i <= seq.length; i++) {
        Vec temp;
        for (int j = i+1; j <= min(seq.length, i+_constraint); j++) {
            data[i-1][j-i-1] = bpp(i, j);
            if (data[i-1][j-i-1] > 1.0) {
                cerr << "error? " << i << " " << j << " " << data[i-1][j-i-1] << endl;
                data[i-1][j-i-1] = 1.0;
            }
        }
    }
}

void Rfold_Lang::Write_bpp()
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);    
    for (int i = 1; i <= seq.length; i++) {
        double P[3] = { 1.0, 0.0, 0.0 };
        for (int j = max(1, i-_constraint); j <= min(seq.length, i+_constraint); j++) {
            if (j == i) continue;
            double value = bpp(i, j);
            P[0] -= value;            
            (j < i) ? P[2] += value : P[1] += value;
        }
        cout << "* " << i << " " << P[1] << " " << P[2] << endl;        
    }
}

void Rfold_Lang::Print_Mat(bool is_alpha) 
{
    if (is_alpha) {
        cout << "--alpha" << endl;
        alpha.Print(seq.str);
    } else {
        cout << "--beta" << endl;       
        beta.Print(seq.str);
    }
}

void Rfold_Lang::Initialize()
{
    alpha = Matrix(seq.length, _constraint, true);
    beta = Matrix(seq.length, _constraint, false);
}

void Rfold_Lang::Set_sequence(const string& sequence) 
{
    if (!noout)    cout << "* " << sequence << endl;
    string str = "$" + sequence;
    vector<int> num_seq;
    transform(str.begin(), str.end(), back_inserter(num_seq), Base());
    seq = Sequence(str, num_seq, (int)sequence.length());
}

void Rfold_Lang::calculation(int constraint, string sequence) 
{
    Set_Constraint(constraint, (int)sequence.length());
    Set_sequence(sequence);    
    Initialize();
    Calc_inside();
    Calc_outside();
}
	
}