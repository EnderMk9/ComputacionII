#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#include "linearalgebra.hpp"

int main(){
    Vector A = {1,2,3};
    Vector B = {1,1,1};
    cout << "Producto escalar de A y B : " << dotprod(A,B) << endl;
    cout << "Norma de A : " << norm(A) << endl;
    cout << "Norma de B : " << norm(B) << endl;
    cout << "Ángulo entre A y B : " << angle(A,B)*180/M_PI << "º" << endl;
}
