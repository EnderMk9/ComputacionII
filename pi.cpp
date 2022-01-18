#include "libs/header.hpp"

long double pi2;
long double pi3;
long double pi4;

int main(){
  pi2 = pow(M_PIl,M_PIl);
  cout.precision(256);
  cout << pi2 << endl;
  pi3 = pow(pi2,pi2);
  cout << pi3 << endl;
  pi4 = pow(pi3,pi3);
  cout << pi4 << endl;
}
