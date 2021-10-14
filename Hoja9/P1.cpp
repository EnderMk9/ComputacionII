#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
#include "libs/rwlib.hpp"
#include "libs/linearalgebra.hpp"

int main(){
    int n = 1000;
    Vector b  = fillvec(4,n);
    Vector ac = fillvec(-1,n-1);
    Vector f  = fillvec(200,n); f[0]=100; f[n-1]=100;
    //Matrix TD = TrDiag(ac,b,ac);
    //coutmat(TD); coutvec(f);
    //Vector SLU = LUSolve(TD,f);
    Vector STr = TrDiagSolve(ac,b,ac,f,1);
    double STrw[n]{};
    copy(STr.begin(), STr.end(), STrw);
    d_w_file_ln("STrRes.dat", STrw,n);
}
