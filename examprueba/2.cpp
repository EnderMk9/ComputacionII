//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

#include "libs/header.hpp" // my setup

Matrix A = {{12,-2,1,4,-0.5},
            {1,15,4,6,7},
            {9,0.5,12,3,1},
            {-3,-6,9,24,1},
            {1,-1,3,2,8}};

Vector b = {-1,21,43,32,55};
int main(){
    cout << "A = ";
    coutmat(A);
    cout << "b = ";
    coutvec(b);
    Vector xLU = LUSolve(A,b);
    cout << "xLU = ";
    coutvec(xLU);
    Vector x0 = VecFull(0,5);
    double emin = 1e-3;double emax = 1e-12; double de = 0.1;
    double imax = log10(emin/emax)+1; double e{};
    IVector evec(imax,0);
    Matrix xJacobi(imax,Vector(5,0));
    IVector iterJ(imax,0);
    for (int i = 0; i < imax; i++){
        e = emin*(pow(de,i));
        evec[i] = log10(e);
        xJacobi[i] = JacobiSolve(A,b,x0,e,iterJ[i],0,1);
    }
    Matrix xGauss(imax,Vector(5,0));
    IVector iterG(imax,0);
    for (int i = 0; i < imax; i++){
        e = emin*(pow(de,i));
        xGauss[i] = GaussSeidelSolve(A,b,x0,e,iterG[i],0,1);
    }
    cwrite_row_int("comp.dat", evec);
    cwrite_row_int("comp.dat", iterJ);
    cwrite_row_int("comp.dat", iterG);
}