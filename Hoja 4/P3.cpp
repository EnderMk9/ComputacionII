#include <iostream>
#include <math.h>
#include <array>
using namespace std;

double F(double (*g)(double), double x){
    double gv = g(x);
    if (gv < 1)          return 2.5;
    else if (1 < gv < 2) return 2.0;
    else if (2 < gv < 2) return 1.5;
    else if (3 < gv)     return 1.0;
}

double g(double x){
    return sin(x);
}

double h(double x){
    return exp(5*x);
}

double i(double x){
    return log(0.01+x);
}

void eval(double (*f)(double), double x0, double xf, int d, double x[]){
    double step = (xf-x0)/d;
    for (int i = 0; i <= d; i++){
        x[i] = f(x0+i*step);
        cout << "f(" << x0+i*step << ") = " << x[i] << endl;
    }
}

int main(){
    int d = 100; double x0 = 0; double xf = 1;
    double x [d+1] {};
    eval(g,x0,xf,d,x);
    cout << F(g,1.5) <<endl;
}