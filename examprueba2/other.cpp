#include "libs/header.hpp"

Matrix A ={{1,0.5,0.3},{0.5,2,0.25},{0.3,0.25,3}};

int main(){
    Matrix U; int k;
    Matrix D  = JacobiDiag(A,U,k,1e-8);
    coutmat(D,4); coutmat(U,4);
}
