#include "libs/header.hpp"

double f (double t, double theta, double omega){
  return omega;
}

double g (double t, double theta, double omega){
  return -sin(theta);
}

int main () {
  Matrix Data = RK4_2O(f,g,0,10,3,0,1000);
  Data = transpose(Data);
  write_mat_double("data.dat",Data);
}
