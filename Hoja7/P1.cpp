#include <iostream>
using namespace std;
#include <math.h>
#include <vector>
#include "libs/linearalgebra.h"

int main(){
    Matrix A = {{1,1, 1},
                {4,3,-1},
                {3,5,3}}; // We define matrix A
    Vector b = {1,6,4};
    Vector x = LUSolve(A,b);
    coutvec(x);
}
