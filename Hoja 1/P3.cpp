#include <iostream>
#include <fstream> 
#include <gmp.h>
using namespace std;

unsigned int fint       = 1;
unsigned long int flint = 1;
float ffloat            = 1;
double fdouble          = 1;
long double fldouble    = 1;

int main(){
    ofstream wfile ("fact.csv");  // set write file
    wfile.precision(18);
    wfile << "it," << "int," << " long int," << " float," << " double," << " long double" << endl;
    for (int i = 1; i <= 1800; i++){
        fint     = fint    *i;
        flint    = flint   *i;
        ffloat   = ffloat  *i;
        fdouble  = fdouble *i;
        fldouble = fldouble*i;
        wfile << i << "," <<fint << "," << flint << "," << ffloat << "," << fdouble << ",'" << fldouble << "'"<<endl;
    }
}