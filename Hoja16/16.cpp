#include "libs/header.hpp"

double f(double r,double V, double U){
    return U;
}
double g(double r,double V, double U){
    return -(2./r)*U;
}

int main(){
    double p = ShootSO(f,g,0.05, 0.1,110, 0, 1001, 1e-10);
    cout.precision(16); cout << p << endl;
    Matrix data = RungeKutta4SO(f,g,0.05,0.15,110,p,1001); data = transpose(data);
    write_mat_double("rk.dat",data);
}