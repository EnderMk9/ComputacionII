#include "libs/header.hpp"
#include "libs/gnuplot.hpp"

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
  function<double (double) > serie = bind(series,_1,as,N,1);
  Vector x(1001,0); Vector y(1001,0);
  eval(serie, -8, 4, 1000, x,y);
  Matrix out(2,Vector(1001,0));
  out[0] = x; out[1] = y;
  out = transpose(out);
  write_mat_double("func.dat", out);
  GnuplotPipe gp; gp.sendLine("plot 'func.dat'");
}
