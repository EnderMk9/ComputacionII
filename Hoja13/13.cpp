#include "libs/header.hpp"

double IEIFK(double phi, double k){
    return 1./sqrt(1-pow(k,2)*pow(sin(phi),2));
}

int main(){
    Matrix Data = MatFull(0,9,31);
    double a = 0; double b = M_PI_2; double k{};
    for (int i = 0; i < 31; i++){
        k = sin((i*M_PI)/72.);
        function<double (double) > TQ = bind(IEIFK,_1,k);
        Data[0][i] = (i*M_PI)/36.;
        Data[1][i] = (2./M_PI)*DefIntTrap(TQ, a, b, 4);
        Data[2][i] = (2./M_PI)*DefIntTrap(TQ, a, b, 5);
        Data[3][i] = (2./M_PI)*DefIntTrap(TQ, a, b, 20);
        Data[4][i] = (2./M_PI)*DefIntSimp(TQ, a, b, 4);
        Data[5][i] = (2./M_PI)*DefIntSimp(TQ, a, b, 5);
        Data[6][i] = (2./M_PI)*DefIntSimp(TQ, a, b, 20);
        Data[7][i] = (2./M_PI)*DefIntGaussLegendre(TQ, a, b, 2);
        Data[8][i] = (2./M_PI)*DefIntGaussLegendre(TQ, a, b, 5);
    }
    Data = transpose(Data);
    write_mat_double("data.dat", Data);
}