#include "libs/header.hpp"

int main(){
    Matrix A = {{1,1,1},
                {1,2,3},
                {1,1,2}}; // We define matrix A
    Vector b = {1,2,3};
    Vector x = LUSolve(A,b,1);
    coutvec(x);
    cout << determinant(A) << endl;
}
