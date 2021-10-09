#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "libs/gnuplot.hpp"    // For plotting inside c++
#include "libs/rwlib.hpp"      // For r-w operations
#include "libs/funlib.hpp"     // For evaluating functions
#include "libs/calculus.hpp"   // For numerical derivatives
#include "libs/nonlinsolv.hpp" // For solving zeroes of non-lineal functions

const double V0 = 1; const double sigma = 1; // sigma = 1 nm

double V(double r){ //the variable r is in nm
    return V0*(pow(sigma/r,6)-exp(-r/sigma)); // Buckingham's Potential
}
double Vp(double r){ //the variable r is in nm
    return V0*(-6*pow(1/r,7)*pow(sigma,6)+exp(-r/sigma)/sigma); // The (central) force associated with V
}
double logfun(double r){
    return log10(V(r));
}

int main(){
    double r0 = 0.1; double rf = 2; int d = 200; // Parameters for plotting V
    string wname = "buckingham.dat";              // Savefile
    double r [d+1] {}; double v [d+1] {};         // arrays of domain and image
    eval(logfun, r0, rf, d, r, v);                     // Evaluate V
    d_w_file_2cols(wname, r, v,d+1);              // Write V
    GnuplotPipe gp; gp.sendLine("plot '" + wname + "'"); // Plot V using GNUPlot
    for (int i = 3; i < 13; i++){ // Test from 1e-3 to 1e-12 tolerande
        cout.precision(i+1);
        cout << "1e-" << i << endl;
        cout << bisection(Vp, 0.1, 2, pow(10.,-i)) << " nm" << endl; // Bisection
        cout << secant(Vp, 1.5, 2, pow(10.,-i)) << " nm" << endl;    // Secant
        cout << NewtonRhapsonNum(Vp, 1.5, pow(10.,-i), 1e-4,1) << " nm" << endl; // Newton
        cout << endl;
    }
    // Podemos observar que el método de Newton es el más rápido, pero solo si elegimos un punto cercano 
    // a la raiz, por otro lado la secante es mucho más sensible a los puntos iniciales debido a que la
    // precisión utilizada puede no ser suficiente para sumar dos números con muchos órdenes de diferencia 
    // y en ese caso no converge.
}
