#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <functional>
using namespace std::placeholders;
#include "libs/gnuplot.hpp"
using namespace std;
#include "libs/funlib.hpp"
#include "libs/rwlib.hpp"

double a0 = 1; double a1 = 1;

// y''+y'-xy=0
double a(int n, double an_1, double an_3){
  return (an_3-an_1*(n-1))/(n*(n-1));
}

double series(double x, double terms[],int N, int n0){
  double y{};
  for (int i = n0; i < N; i++) {
    y += terms[i]*pow(x,i-n0);
    if (x < -4) {
    ofstream wfile;
    wfile.open("debug.dat", ios_base::app);
    wfile << "x: " << x << " i: " << i << " y: " << y << endl; }
  }
  return y;
}

int main(){
  const int N = 250;
  double as[N+1]{};
  as[1] = a0; as[2] = a1;
  for (int i = 2; i < N; i++) {
    as[i+1] = a(i+1,as[i],as[i-2]);
  }
  d_w_file_ln("terms.csv", as, N+1,50);
  function<double (double) > serie = bind(series,_1,as,N,1);
  double x [1001]{}; double y[1001]{};
  eval(serie, -8, 4, 1000, x,y);
  d_w_file_2cols("func.dat", x, y,1001);
  GnuplotPipe gp; gp.sendLine("plot 'func.dat'");
}
