//----------------------------------------------------------
// Abel Rosado Peinado - 5265 - 26 de Noviembre de 2121
//----------------------------------------------------------

#include "ARP_libs/header.hpp"

Vector b = {0,1,2,3,4,5,6,7,8,9};

Matrix A = {{10,1,1,1,1,1,1,1,1,0},
            {1,10,1,1,1,1,1,1,0,1},
            {1,1,10,1,1,1,1,0,1,1},
            {1,1,1,10,1,1,0,1,1,1},
            {1,1,1,1,10,0,1,1,1,1},
            {1,1,1,1,0,10,1,1,1,1},
            {1,1,1,0,1,1,10,1,1,1},
            {1,1,0,1,1,1,1,10,1,1},
            {1,0,1,1,1,1,1,1,10,1},
            {0,1,1,1,1,1,1,1,1,10}};

int main(){
    //----------------------------------------------------------
    // CHECK DETERMINANT
    Matrix U; int k; Matrix D =  JacobiDiag(A, U, k, 1e-8);
    Vector lambda = MatrixDiag(A); cout << "Det(A) = " << VecProd(lambda) << endl;
    //----------------------------------------------------------
    // SOLVE LU AND CHECK
    vector xLU = LUSolve(A,b,0); cout << "xLU = ";coutvec(xLU,12);
        //Matrix S = Vec2Mat(xLU,1);
        //Matrix j = matprod(A,S);
        //coutmat(j,5);
    //----------------------------------------------------------
    // JACOBI AND GAUSS SEIDEL
    Vector x0 = VecFull(0,10);
    double emin = 1e-3;double emax = 1e-12; double de = 0.1;
    double imax = log10(emin/emax)+1; double e{};
    IVector evec(imax,0);
    // JACOBI
    Matrix xJacobi(imax,Vector(10,0));
    IVector iterJ(imax,0);
    for (int i = 0; i < imax; i++){
        e = emin*(pow(de,i));
        evec[i] = log10(e);
        xJacobi[i] = JacobiSolve(A,b,x0,e,iterJ[i],0,1);
    }
    // GAUSS SEIDEL
    Matrix xGauss(imax,Vector(10,0));
    IVector iterG(imax,0);
    for (int i = 0; i < imax; i++){
        e = emin*(pow(de,i));
        xGauss[i] = GaussSeidelSolve(A,b,x0,e,iterG[i],0,1);
    }
    // EXPORT SOLUTIONS
    write_mat_double("Jdata.dat", xJacobi);
    write_mat_double("Gdata.dat", xGauss);
    // SHOW SOLUTION
    cout << "xJacobi = ";coutvec(xJacobi[imax-1],12);
    cout << "xGauss = ";coutvec(xGauss[imax-1],12);
    // EXPORT ITERATION DATA
    IMatrix Sdata(3,IVector(imax,0));
    Sdata[0] = evec; Sdata[1] = iterJ; Sdata[2] = iterG;
    Sdata = Itranspose(Sdata);
    write_mat_int("sdata.dat", Sdata);
}
