//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

#include "libs/header.hpp" // my setup

Vector f(Vector& x){
    Vector T(3,0);
    T[0] = 6*x[0]- 2*cos(x[1]*x[2]) - 1;
    T[1] = 9*x[1]+sqrt(x[0]*x[0]+sin(x[2])+1.06)+0.9;
    T[2] = 60*x[2]+3*exp(-1*x[0]*x[1])+10*M_PI-3;
    return T;
}

int main(){
    Vector x0 = {1,1,1};
    Vector f0 = NewtonNumSys(3,x0,f,0.001,10e-8,1);
    coutvec(f0);
}

