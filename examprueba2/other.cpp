#include "libs/header.hpp"

double q (double t){
  return -1;
}
double p (double t){
  return 1;
}
double r (double t){
  return 0;
}

int main(){
  Vector y = FiniteDiff(p,q,r,100,0,4,0,-2);
  coutvec(y);
}
