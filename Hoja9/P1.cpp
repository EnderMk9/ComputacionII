//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

#include "libs/header.hpp" // my setup

int main(){
    int n = 1000;
    Vector b  = VecFull(4,n);
    Vector ac = VecFull(-1,n-1);
    Vector f  = VecFull(200,n); f[0]=100; f[n-1]=100;
    //Matrix TD = TrDiag(ac,b,ac);
    //coutmat(TD); coutvec(f);
    //Vector SLU = LUSolve(TD,f);
    Vector STr = TrDiagSolve(ac,b,ac,f,1);
    write_col_double("STrOut.dat", STr);
    Matrix M = read_matrix_double("rdata.txt",3);
    coutmat(M);
}
