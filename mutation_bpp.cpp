#include "part_func.h"

namespace Rfold {

void Rstate_model::Calc_bpp()
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);    
    for (int i = 1; i <= _length; i++) {
        double P[3] = { 1.0, 0.0, 0.0 };
        for (int j = 1; j <= _length; j++) {
            if (j == i) continue;
            double value = bpp(i, j);
            if (j < i) {
                P[2] += value;
                P[0] -= value;
            } else {
                P[1] += value;
                P[0] -= value;
            }
        }
        if (!noout) cout << "* " << i << " " << P[1] << " " << P[2] << endl;        
    }
}

/* ///////////////////////////////////////////// */

double Rstate_model::hairpin_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = 0; i < x1; i++) 
        for (int j = x2; j <= _length; j++) 
            value = Logsumexp(value, Logsum(beta.stem[i][j-i], LogHairpinEnergy(i, j)));
    return value;
}

double Rstate_model::interior_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = 0; i < x1; i++) {
        for (int j = x2; j <= _length; j++) {
            for (int k = x2; k < j; k++) 
                for (int l = k; l < j; l++)
                    value = Logsumexp(value, Logsum(beta.stemend[i][j-i], LogLoopEnergy(i, k, l, j), alpha.stem[k][l-k]));
            for (int k = x2; k < j; k++) 
                for (int l = k; l < j; l++)
                    value = Logsumexp(value, Logsum(beta.stemend[i][j-i], LogLoopEnergy(i, k, l, j), alpha.stem[k][l-k]));
        }
    }
    return value;
}

double Rstate_model::multi_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int j = x2; j < _length; j++) 
        value = Logsumexp(value, Logsum(beta.multi[x1-1][j-x1+1], (x2-x1+1)*ML_BASE, alpha.multi[x2][j-x2]));
    return value;
}

double Rstate_model::multi2_acc(int x1, int x2) 
{
    double value = 0.0;
    for (int i = 0; i <= x1; i++) 
        value = Logsumexp(value, Logsum(beta.multi2[i][x2-i], (x2-x1+1)*ML_BASE, alpha.multi2[i][x1-1-i]));
    return value;
}

void Rstate_model::Calc_acc()
{
    cout.setf(std::ios_base::fixed, std::ios_base::floatfield);    
    for (int i = 1; i < _length; i++) {
        for (int j = i+1; j <= _length; j++) {
            double value = hairpin_acc(i, j);
            value = Logsumexp(value, interior_acc(i, j));
            value = Logsumexp(value, multi_acc(i, j));  
            value = Logsumexp(value, multi2_acc(i, j));
            value = Logsumexp(value, Logsum(beta.outer[j], alpha.outer[i-1]));
            if (!noout) cout << "* " << i << " " << j << " " << value << endl;        
        }
    }
}

/* ///////////////////////////////////////////// */

double Rstate_model::Sum_Dangle(int type, int five, int three) 
{
    double temp = 0.0;
    if (five > 0) temp = Logsum(temp, logdangle5[type][_sequence[five]]);
    if (three <= _length) temp = Logsum(temp, logdangle3[type][_sequence[three]]);
    else if (Is_AU(type)) temp = Logsum(temp, logTermAU);
    return temp;
}
 
double Rstate_model::LogHairpinEnergy(int i, int j) 
{
    int type = bp(i, j);
    double q = -INF;
    int d = j-i-1;
    q = (d <= 30) ? loghairpin[d] : loghairpin[30]-lxc37*log(d/30.0)*10.0/kT;
    if (tetra && d==4) {
  	    string sub_seq = _str.substr(i,d+2);			
   	    size_t tel = Tetraloops.find(sub_seq);
	    if (tel != string::npos)
		    q = Logsum(q, logtetra[tel/7]);
    }
    if (d == 3) {
        if (Is_AU(type))
            q = Logsum(q, logTermAU);                
    } else 
        q = Logsum(q, logmismatchH[type][_sequence[i+1]][_sequence[j-1]]);
    if (!Is_INF(q) && debug) {
        cout << "hairpin " << i+1 << _str.substr(i+1, j-1-i) << j-1 << " " << exp(q) << endl;
        cout << "        " << i << _str[i] << "-" << _str[j] << j << endl;
        if (d <= 30) cout << d << " " << loghairpin[d] << " ";
        else cout << d << " " << loghairpin[30]-lxc37*log(d/30.0)*10.0/kT << " ";
        if (d != 3) cout << "Mis" << type << _sequence[i+1] << _sequence[j-1] << 
            _str[i+1] << _str[j-1] << " " << logmismatchH[type][_sequence[i+1]][_sequence[j-1]] << " ";
        else if (Is_AU(type)) cout << "AU" << logTermAU;
        cout << endl;
    }
    return q;
}
double Rstate_model::LogBulge(int u, int type, int type2) 
{
    double z = (u <= 30) ? logbulge[u] : logbulge[30]-static_cast<double>(lxc37)*log(u/30.0)*10.0/kT;  
    if (u == 1) z = Logsum(z, logstack[type][type2]);
    else {
        if (Is_AU(type)) z = Logsum(z, logTermAU);
        if (Is_AU(type2)) z = Logsum(z, logTermAU);
    }
    if (debug) { 
        if (u <= 30) cout << "bulge " << logbulge[u] << endl;
        else cout << "bulge " << logbulge[30]-static_cast<double>(lxc37)*log(u/30.0)*10.0/kT << endl;
    }
    return z;
}

double Rstate_model::LogLoopEnergy(int i, int j, int p, int q) 
{
    int type = bp(i, j), type2 = rbp(p, q);
    double z = -INF;
    int u1 = p-i-1, u2 = j-q-1, u = max(u1, u2); 
    if (u1 == 0 && u2 == 0) return logstack[type][type2];
    else if (no_closingGU && (Is_closeGU(type) || Is_closeGU(type2))) return z;
    if (debug) cout << i << _str.substr(i, p-i+1) << p << " loop\n" << q << _str.substr(q, j-q+1) << j << endl;
	if ((u1 == 0) || (u2 == 0)) { /* bulge */
        return LogBulge(u, type, type2);
	} else {
        if (u <= 2) {             /* short internal */
	        if (u1+u2 == 2) 
                z = logint11[type][type2][_sequence[i+1]][_sequence[j-1]];
	        else if (u1 == 1 && u2 == 2) 
                z = logint21[type][type2][_sequence[i+1]][_sequence[q+1]][_sequence[j-1]];                
            else if (u1 == 2 && u2 == 1)
                z = logint21[type2][type][_sequence[q+1]][_sequence[i+1]][_sequence[p-1]];
            else                 
                z = logint22[type][type2][_sequence[i+1]][_sequence[p-1]][_sequence[q+1]][_sequence[j-1]];
		} else {                 /* long internal */
		    z = loginternal[u1+u2];
		    double temp1 = logmismatchI[type][_sequence[i+1]][_sequence[j-1]];
		    double temp2 = logmismatchI[type2][_sequence[q+1]][_sequence[p-1]];
		    double temp3 = logninio[abs(u1-u2)];
		    z = Logsum(z, temp1, temp2, temp3);
   		}
        if (debug) cout << "internal" << u1 << " " << u2 << " " << type << type2 << " " << z*kT/10.0 << endl;
	}
    return z;
}
/* ///////////////////////////////////////////// */

void Rstate_model::calculation(int constraint, string& sequence) 
{
    if (_str == "") Initiallize(constraint, sequence);
    else {
        Set_sequence(sequence);
        beta.multi = alpha.multi = Mat(_length+1, vector<double>(_constraint+2, 0.0));
        for (int i = 0; i <= _length;i++)
            for (int j = 0; j <= _constraint+1; j++)
                alpha.multi[i][j] = beta.multi[i][j] = -INF;
    }
    alpha.Set_INF();
    beta.Set_INF();
    Calc_inside();
    Calc_outside();
}

}
