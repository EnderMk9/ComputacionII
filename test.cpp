#include "libs/header.hpp"

double sen(double x){
    return sin(x);
}

int main(){
    double root = secant(sen, 1, 0.5, 1e-5);
    cout << root << endl;
}