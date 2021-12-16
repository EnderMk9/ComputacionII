#include "libs/header.hpp"
#include "libs/gnuplot.hpp"

double g = 9.81;
double vlim = 57;
double gamm = g/(vlim*vlim);

double a(double t, double v, double gamma){
    return g-gamma*v*v;
}
double xdot(double t, double v, double x){
    return -v;
}
Vector F(double t, Vector x){
    Vector F = VecFull(0,2);
    F[0] = -x[1];
    F[1] = g - gamm*x[1]*x[1];
    return F;
}

int main(){
    int n = 10e4; double t0 = 0; double tf = 25; Vector x0 = {100,0};
    //function<double (double,double) > aB = bind(a,_1,_2,gamm);
    //Matrix DataE = EulerPO(aB,t0,tf,v0,n);
    //Matrix DataRK = RungeKutta4FO(aB,t0,tf,x0[1],n);
    //Matrix DataRK = RungeKutta4SO(aB,xdot,t0,tf,x0[1],x0[0],n);
    //DataE = transpose(DataE);
    //write_mat_double("funcE.dat",DataE);
    //DataRK = transpose(DataRK);
    //write_mat_double("funcRK.dat",DataRK);
    Matrix x = RungeKutta4nO(F,t0,tf,x0,n);
    write_mat_double("funcX.dat",x);
}