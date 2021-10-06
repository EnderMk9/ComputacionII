#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#include "matrixtoolbox.h"

int main(){
    Matrix A = {{2,2,8},
                {2,2,4},
                {2,2,6}}; // We define matrix A
    Matrix B = {{4,6},
                   {3,9},
                   {1,8}};   // We define matrix B
    cout << "A = "; coutmat(A);
    cout << "B = "; coutmat(B);
    Matrix C = matprod(A,B);
    cout << "C = "; coutmat(C);
    Matrix Ap = gramschmidt(A);
    cout << "Ap = "; coutmat(Ap);
    Matrix Z = transpose(Ap);
    Matrix I = matprod(Ap,Z);
    cout << "I = "; coutmat(I);
}