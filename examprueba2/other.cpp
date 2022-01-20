#include "libs/header.hpp"

double gauss(double x){
  return exp(-(x*x)/2.)/sqrt(2*M_PI);
}

double h = 0.01; double xmin = -4; double xmax = 4;

int main(){
  int n = ((xmax-xmin)/h)+1;
  Matrix xy = IntIndSmp13(gauss,xmin,xmax,100,n);
  Vector dF = derivativeCDArr(xy[0], xy[1]);
  Vector f = eval(gauss, xy[0]);
  Vector e = VecDiff(dF, f);
  e = VecDiv(e, f); e = VecAbs(e);
  Matrix Data = MatFull(0,5,81);
  Data[0] = xy[0]; Data[1] = xy[1]; Data[2] = dF; Data[3] = f; Data[4] = e;
  Data = transpose(Data);
  write_mat_double("gauss0.01.dat",  Data);
}
