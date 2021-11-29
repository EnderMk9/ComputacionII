#include "libs/header.hpp"
Matrix A={{7,6,5,7,6},{0,2,3,4,5},{0,0,3,2,1},{0,0,0,1,2},{0,0,0,0,6}}; 
Vector d={23,14,6,3,6};
int main(){
    coutmat(A);
    Vector s = LUSolve(A,d,1);
    coutvec(s);
    Matrix S = Vec2Mat(s,1);
    Matrix j = matprod(A,S);
    coutmat(j);
}