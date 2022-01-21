//----------------------------------------------------------
// Abel Rosado Peinado - 5265 - 21 de enero de 2022
//----------------------------------------------------------
#include "ARPLibs/header.hpp"

double e = 1e-10;

double oneoverw(double w){
  return 1./w;
}

int main(){
  cout.precision(-log10(e)+1);
  cout << "CALCULATING OPTIMAL N" << endl;
  for (int i = 1; i <= 9; i++) {
    int n = ceil(pow(24*pow(i-1,5.)/(180*e), 1./4));
    if (n % 2 != 0) n++;
    cout << "ln(" << i << ")  ";
    cout << "n = " << n << "  ";
    double lg = IntSmp13(oneoverw,1,i,n);
    cout << "numsol = " << lg << "  ";
    cout << "c++sol = " << log(i) << "  ";
    cout << "error = " << abs(lg-log(i)) << endl;
  }
  cout << "ITERATING AND CHECKING DIFFERENCE BETWEEN EACH ITERATION" << endl;
  for (int i = 1; i <= 9; i++) {
    cout << "ln(" << i << ")  ";
    double lgtol = IntSmp13Tol(oneoverw, 1, i, e, 10);
    cout << "numsol = " << lgtol<< "  ";  
    cout << "c++sol = " << log(i) << "  ";
    cout << "error = " << abs(lgtol-log(i)) << endl;
  }
}

