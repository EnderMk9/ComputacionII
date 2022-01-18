#include "libs/header.hpp"

Matrix A = {{1,2,3,4,5},
            {2,6,4,3,2},
            {3,4,3,1,0},
            {4,3,1,2,9},
            {5,2,0,9,1}};

int main(){
  Matrix U = MatFull(0,5,5); int k{};
  Matrix D = JacobiDiag(A, U, k, 1e-10);
  Vector lambda = MatrixDiag(D);
  coutvec(lambda,5); coutmat(U,5); cout << k << endl;
  
}
