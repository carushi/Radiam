#ifndef _PARAM_H
#define _PARAM_H
#include "matrix.h"

namespace Rfold {
namespace Parameter {


extern double loghairpin[31];
extern double logtetra[30];
extern double logmismatchH[7][5][5];
extern double logmismatchI[7][5][5];
extern double logstack[7][7];
extern double logbulge[31];
extern double logTermAU;
extern double logint11[8][8][5][5];
extern double logint21[8][8][5][5][5];
extern double logint22[8][8][5][5][5][5];
extern double loginternal[31];
extern double logdangle5[8][5];
extern double logdangle3[8][5];
extern double logninio[MAXLOOP+1];
extern bool initialized;
extern double logMLintern;    
extern double logMLclosing;    

static const bool no_closingGU = false;
static const bool tetra = true;    
static const int temperature = 37;
static const double ML_BASE = 0;
static const double kT = (temperature+K0)*GASCONST;
static const double lxc37 = 107.856; 

inline bool Is_AU(const int type) {
    return (type > 2);    
}
inline bool Is_closeGU(const int type) {
    return (type == 3 || type == 4);    
}

static double LogEnergy(const int value) {
    if (value == INF) return -INF;
    else return -(static_cast<double>(value))*10.0/Parameter::kT;
}

inline double logML_BASE() { return LogEnergy(ML_BASE); }

static void Init_mis_stack()
{
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 5; k++) {
                logmismatchI[i][j][k] = LogEnergy(mismatchI37[i][j][k]);
                logmismatchH[i][j][k] = LogEnergy(mismatchH37[i][j][k]);
            }
        }
        for (int j = 0; j < 7; j++) 
            logstack[i][j] = LogEnergy(stack37[i][j]);
    }
}

static void Init_internal() 
{
    for (int i = 0; i < 31; i++) 
        loginternal[i] = LogEnergy(internal_loop37[i]);
    for (int i = 0; i <= 7; i++) {
        for (int j = 0; j <= 7; j++) {
            for (int k = 0; k <= 4; k++) {
                for (int l = 0; l <= 4; l++) {
                    logint11[i][j][k][l] = LogEnergy(int11_37[i][j][k][l]);
                    for (int m = 0; m <= 4; m++) {
                        logint21[i][j][k][l][m] = LogEnergy(int21_37[i][j][k][l][m]);
                        for (int n = 0; n <= 4; n++) {
                            logint22[i][j][k][l][m][n] = LogEnergy(int22_37[i][j][k][l][m][n]);
                        }
                    }
                }
            }
        }
    }
}

static void Init_dangle()
{
    logTermAU = LogEnergy(TerminalAU);  
    for (int i = 0; i <= 7; i++)
        for (int j = 0; j <= 4; j++) {
            logdangle5[i][j] = LogEnergy(dangle5_37[i][j]);
            logdangle3[i][j] = (Is_AU(i)) ? LogEnergy(dangle3_37[i][j]+TerminalAU) 
                                          : LogEnergy(dangle3_37[i][j]);
        }
}

static void Init_loop()
{
    for (int i = 0; i < 31; i++) {
        loghairpin[i] = LogEnergy(hairpin37[i]);
        logbulge[i] = LogEnergy(bulge37[i]);
    }
    for (int i = 0; i < 30; i++) 
        logtetra[i] = LogEnergy(tetra_energy37[i]);    
    for (int i = 0; i <= MAXLOOP; i++)
        logninio[i]=LogEnergy(min(MAX_NINIO,i*F_ninio37));
}

static void Init_ener_param()
{
    Init_loop();
    Init_mis_stack();
    logMLclosing = LogEnergy(ML_closing37);
    logMLintern = LogEnergy(ML_intern37);
    Init_internal();
    Init_dangle();
    initialized = true;
}


static double Sum_Dangle(int type, int five, int three, const Sequence& seq)
{
    double temp = 0.0;
    if (five > 0) temp = Logsum(temp, logdangle5[type][seq.sequence[five]]);
    if (three <= (int)seq.length) temp = Logsum(temp, logdangle3[type][seq.sequence[three]]);
    else if (Is_AU(type)) temp = Logsum(temp, logTermAU);
    return temp;
}
 
static double LogHairpinEnergy(int i, int j, const Sequence& seq)
{
    int type = bp(i, j, seq.sequence);
    double q = -INF;
    int d = j-i-1;
    q = (d <= 30) ? loghairpin[d] : loghairpin[30]-lxc37*log(d/30.0)*10.0/kT;
    if (tetra && d==4) {
        string sub_seq = seq.str.substr(i,d+2);            
        size_t tel = Tetraloops.find(sub_seq);
        if (tel != string::npos)
            q = Logsum(q, logtetra[tel/7]);
    }
    if (d == 3) {
        if (Is_AU(type))
            q = Logsum(q, logTermAU);                
    } else 
        q = Logsum(q, logmismatchH[type][seq.sequence[i+1]][seq.sequence[j-1]]);
    return q;
}

static double LogBulge(int u, int type, int type2) 
{
    double z = (u <= 30) ? logbulge[u] : logbulge[30]-static_cast<double>(lxc37)*log(u/30.0)*10.0/kT;  
    if (u == 1) z = Logsum(z, logstack[type][type2]);
    else {
        if (Is_AU(type)) z = Logsum(z, logTermAU);
        if (Is_AU(type2)) z = Logsum(z, logTermAU);
    }
    return z;
}

static double LogLoopEnergy(int i, int j, int p, int q, const Sequence&  seq)
{
    int type = bp(i, j, seq.sequence), type2 = rbp(p, q, seq.sequence);
    double z = -INF;
    int u1 = p-i-1, u2 = j-q-1, u = max(u1, u2); 
    if (u1 == 0 && u2 == 0) return logstack[type][type2];
    else if (no_closingGU && (Is_closeGU(type) || Is_closeGU(type2))) return z;
    if ((u1 == 0) || (u2 == 0)) { /* bulge */
        return LogBulge(u, type, type2);
    } else {
        if (u <= 2) {             /* short internal */
            if (u1+u2 == 2) 
                z = logint11[type][type2][seq.sequence[i+1]][seq.sequence[j-1]];
            else if (u1 == 1 && u2 == 2) 
                z = logint21[type][type2][seq.sequence[i+1]][seq.sequence[q+1]][seq.sequence[j-1]];
            else if (u1 == 2 && u2 == 1)
                z = logint21[type2][type][seq.sequence[q+1]][seq.sequence[i+1]][seq.sequence[p-1]];
            else                 
                z = logint22[type][type2][seq.sequence[i+1]][seq.sequence[p-1]][seq.sequence[q+1]][seq.sequence[j-1]];
        } else {                 /* long internal */
            z = loginternal[u1+u2];
            double temp1 = logmismatchI[type][seq.sequence[i+1]][seq.sequence[j-1]];
            double temp2 = logmismatchI[type2][seq.sequence[q+1]][seq.sequence[p-1]];
            double temp3 = logninio[abs(u1-u2)];
            z = Logsum(z, temp1, temp2, temp3);
        }
    }
    return z;
}
}
}

#endif