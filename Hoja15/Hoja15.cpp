#include "libs/header.hpp"

double m1 = 2; double m2 = 3.5;
double k1 = 2.5; double k2 = 3.5;

Vector g(double t, Vector x){
    Vector F = VecFull(0,4);
    F[0] = x[1]; // xdot
    F[1] = -(k1/m1)*x[0]-(k2/m1)*(x[0]-x[2]); // xddot
    F[2] = x[3]; // ydot
    F[3] = (k2/m2)*(x[0]-x[2]); // yddot
    return F;
}

int main(){
    int n = 1001; double t0 = 0; double tf = 100; Vector x0 = {3,0,4,0};
    Matrix x = RungeKutta4nO(g,t0,tf,x0,n);
    write_mat_double("funcX.dat",x);
    Matrix D = {{-(k1/m1)-(k2/m1),(k2/m1)},{(k2/m2),-(k2/m2)}};
    coutmat(D);
}