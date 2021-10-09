#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "libs/gnuplot.hpp"
#include "libs/rwlib.hpp"
#include "libs/funlib.hpp"
#include "libs/calculus.hpp"
#include "libs/nonlinsolv.hpp"

double h(double t){
    double h0 = 300; double g  = 32.17; double m  = 0.25; double k  = 0.1;
    return h0 - (m*g/k)*t + (m*m*g/(k*k)) * (1-exp(-k*t/m));
}

double hp(double t){
    double h0 = 300; double g  = 32.17; double m  = 0.25; double k  = 0.1;
    return (exp(-k*t/m)-1)*(m*g/k);
}

int main(){
    double t0 = 0; double tf = 10; int d = 100; string wname = "eval.dat";
    double t [d+1] {}; double y [d+1] {};
    eval(h, t0, tf, d, t, y);
    d_w_file_2cols(wname, t, y,d+1);
    //GnuplotPipe gp; gp.sendLine("plot '" + wname + "'");
    cout.precision(11);
    cout << bisection(h,6.0,6.1,1e-8) << endl;
    cout << secant(h,6.2,6.1,1e-8) << endl;
    cout << NewtonRhapsonNum(h, 6.5, 1e-8, 0.01, 1) << endl;
}
