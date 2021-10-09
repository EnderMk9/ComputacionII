#include <iostream>
using namespace std;
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include "libs/rwlib.hpp"
#include "libs/linearalgebra.hpp"
// i, x0_1, x0_2, x0_3, x1_1, x1_2, x1_3, Δx_1, Δx_2, Δx_3, ϵ
int main(){
    Matrix A = {{-6,0,3,1},
                {1,4,1,-1},
                {3,-2,8,1},
                {1,2,0,3}}; // We define matrix A
    Vector b = {3,-4,5,-2};
    Vector AD = MatrixDiag(A);
    Vector x0 = VecDiv(b,AD);
    Vector sJ = JacobiSolve(A, b, x0, 1e-5,1);
    coutvec(sJ);
    Vector sG = GaussSeidelSolve(A, b, x0, 1e-5,1);
    Vector sLU = LUSolve(A,b);
    coutvec(sG); coutvec(sLU);
}
