#include "libs/header.hpp"

double x2 (double x){
  return cos(x);
}

int main(){
  Vector x = linspace(0, 10, 100);
  Vector y = eval(x2, x); cout.precision(9);
  cout << IntSmpLoop(x2,0,10,1e-9,10,1) << endl;
  cout << sin(10) << endl;
}
