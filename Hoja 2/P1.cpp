#include <iostream>
#include "exp.h"
using namespace std;

double x;

int main(){
    cout << "Introduce the value x for which you want to calculate the exponential : ";
    cin >> x;
    double ex = exp(x);
    cout.precision(6);
    cout << "exp(" << x << ") = "<< ex << endl;
}
