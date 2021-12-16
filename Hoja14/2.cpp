#include "libs/header.hpp"
#include "libs/gnuplot.hpp"

double omega = 0.1;

Vector F(double t, Vector x){
    Vector F = VecFull(0,4);
    F[0] = x[1]; //  xdot = vx
    F[1] = 2*omega*x[3] + omega*omega*x[0]; // xddot = 2omega vy + omega² x
    F[2] = x[3]; // ydot = vy 
    F[3] = -2*omega*x[1] + omega*omega*x[2]; // yddot = 2omega vx + omega² y
    return F;
}

int main(){
    int n = 10e4; double t0 = 0; double tf = 300; Vector x0 = {1,0,1,0};
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
    GnuplotPipe gp; gp.sendLine("plot 'funcX.dat'");
}