#include "libs/header.hpp"

double e = 1e-10;

double oneoverx(double x){
  return 1./x;
}

int main(){
  cout.precision(-log10(e)+1);
  for (int i = 2; i <= 9; i++) {
    int n = ceil(pow(6*pow(i-1,5.)/(180*e), 1./4));
    if (n % 2 != 0) n++;
    cout << "ln(" << i << ")  ";
    cout << "n = " << n << "  ";
    double lg = DefIntSimp13(oneoverx,1,i,n);
    cout << "numsol = " << lg << "  ";
    cout << "c++sol = " << log(i) << "  ";
    cout << "error = " << abs(lg-log(i)) << endl;
  }
}
