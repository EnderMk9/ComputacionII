//----------------------------------------------------------
// Abel Rosado Peinado - 5265 - 26 de Noviembre de 2021
//----------------------------------------------------------

#include "libs/header.hpp"

const double g = 9.81;

double h(double D, double L, double u, double lambda){
    return (lambda*L*u*u)/(2*D*g);
}

// Function to find root
double lambdaroot(double lambda, double D,double K, double Re){
    return pow(lambda,-1./2)+2*log10((2.51)/(pow(lambda,1./2)*Re)+(K)/(3.71*D));
}

int main(){
    cout.precision(9); // cout precision
    //----------------------------------------------------------
    // Define the function fixing the parameters
    double K = 0.25e-3; double D = 0.3; double Re = 200000;
    function<double (double) > lambdarootbind = bind(lambdaroot,_1,D,K,Re);
    //----------------------------------------------------------
    // BISECTION
    int bisi; double bisroot = bisection(lambdarootbind, 0.0001, 1.0, 1e-9, bisi);
    cout << "Bisection, i = " << bisi << ", lambda = " << bisroot << endl;
    //----------------------------------------------------------
    // NEWTON-RHAPSON
    int newi; double newroot = NewtonNum(lambdarootbind, 0.0001, 1e-9, 1e-9, 0, newi);
    cout << "Newton-Rhapson, i = " << newi << ", lambda = " << newroot << endl;
    //----------------------------------------------------------
    // SECANT (it displays x^{i+1} each iteration)
    int seci; double secroot = secant(lambdarootbind, 0.0001, 1.0, 1e-9, seci);
    cout << "Secant, i = " << seci << ", lambda = " << secroot << endl;
    //----------------------------------------------------------
    // Evaluate f(\lambda) and export data
    Vector lambda = linspace(0.0001, 1, 1000);
    vector f = eval(lambdarootbind, lambda); int nf = f.size();
    Matrix fdata(2,Vector(nf,0)); fdata[0] = lambda; fdata[1] = f;
    fdata = transpose(fdata); write_mat_double("fdata.dat", fdata);
}