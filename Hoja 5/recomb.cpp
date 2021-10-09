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

const double eta = 6e-10; const double Ebind = 13.6; const double me = 0.510998928e6;
const double z3 = 1.202056903159594285399;

double f(double T){
    double Xe = 0.1;
    return (eta*((2*z3)/(M_PI*M_PI))*pow((2*M_PI*T)/me,3./2)*exp(Ebind/T) + (Xe-1)/(Xe*Xe));
}

int main(){
    double T0 = 0.25; double Tf = 0.5; int d = 200; string wname = "recomb.dat";
    double T [d+1] {}; double v [d+1] {};
    eval(f, T0, Tf, d, T, v);
    d_w_file_2cols(wname, T, v,d+1);
    GnuplotPipe gp; gp.sendLine("plot '" + wname + "'");
    cout.precision(15);
    cout << bisection(f, T0, Tf, 1e-5) << " eV" << endl;
    cout << NewtonRhapsonNum(f, 0.25, 1e-5, 0.001,0) << " eV" << endl;
}
