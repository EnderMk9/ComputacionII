#include "libs/header.hpp"

double g = 9.8; double l = 9.8; double m = 1; double F0 = 5; double b = 1e-3; double wf = 0.8;
double theta0 = 0; double thetadot0 = 0; double dt = 0.1; double t0 = 0; double tf = 160;

double f1(double t, double theta, double thetadot){
  return thetadot;
}
double f2(double t, double theta, double thetadot){
  return -(b/m)*thetadot -(g/l)*theta +(F0/(m*l))*cos(wf*t);
}

int main(){
  int n = tf/dt+1;
  Vector t = linspace(t0, tf, n);
  Matrix Data = RK4_2O(f1,f2,t0,tf,theta0,thetadot0, n);
  Data = transpose(Data);
  write_mat_double("data.dat", Data,0,5);
}
