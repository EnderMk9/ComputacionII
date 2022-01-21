#include "libs/header.hpp"

double sin2 (double x){
  return sin(2*x);
}

int main(){
  cout.precision(10);
  cout << IntTr(sin2,0,2,12) << endl; 
  cout << IntSmp13(sin2,0,2,12) << endl; 
  cout << IntSmp38(sin2,0,2,12) << endl; 
  cout << IntQuadGL(sin2,0,2,5) << endl; 
}
