#include "libs/header.hpp"
#include "libs/gnuplot.hpp"

double g = 9.81;

double a(double t, double v, double x, double gamma){
    return g-gamma*v*v;
}
double xdot(double t, double v, double x){
    return -v;
}

int main(){
    int n = 10e4; double t0 = 0; double tf = 25; double v0 = 0;
    double vlim = 57;
    function<double (double,double,double) > aB = bind(a,_1,_2,_3,g/(vlim*vlim));
    //Matrix DataE = EulerPO(aB,t0,tf,v0,n);
    //Matrix DataRK = RungeKutta4PO(aB,t0,tf,v0,n);
    Matrix DataRK = RungeKutta4SO(aB,xdot,t0,tf,v0,100,n);
    //DataE = transpose(DataE);
    //write_mat_double("funcE.dat",DataE);
    DataRK = transpose(DataRK);
    write_mat_double("funcRK.dat",DataRK);
}