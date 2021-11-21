//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

#include "libs/header.hpp" // my setup

double f(double x){
    return exp(x) - x*(x+5) -2;
}

double fp(double x){
    return exp(x) -2*x -5;
}

int main(){
    double tol = 1e-5; int i{};
    double xmin = -6; double xmax = 5;
    double dx = 0.1;
    double jmax = (xmax - xmin)/dx +1;
    IVector iter(jmax,0);
    Vector xs    = VecFull(0,jmax);
    Vector roots = VecFull(0,jmax);
    IVector index(jmax,0);
    for (int j = 0; j < jmax; j++){
        i = 0; xs[j] = xmin+j*dx;
        roots[j] = NewtonAnal(f,fp, xs[j], tol, i);
        iter[j]  = i;
    }
    Icoutvec(iter);
    cout << endl;
    coutvec(xs,2);
    cout << endl;
    roots = countdiffval(roots, tol/10, index);
    Icoutvec(index);
    coutvec(roots,-log10(tol));
    int n0{}; int n1{}; int n2{};
    for (int j = 0; j < jmax; j++){
        if(index[j] == 0){
            n0++;
        }
        if(index[j] == 1){
            n1++;
        }
        if(index[j] == 2){
            n2++;
        }
    }
    cout << n0;
    Vector xs0 = VecFull(0,n0);
    Vector xs1 = VecFull(0,n1);
    Vector xs2 = VecFull(0,n2);
    Vector it0 = VecFull(0,n0);
    Vector it1 = VecFull(0,n1);
    Vector it2 = VecFull(0,n2);
    int i0{};  int i1{};  int i2{};
    for (int j = 0; j < jmax; j++){
        if(index[j] == 0){
            xs0[i0] = xs[j];
            it0[i0] = iter[j];
            i0++;
        }
        if(index[j] == 1){
            xs1[i1] = xs[j];
            it1[i1] = iter[j];
            i1++;
        }
        if(index[j] == 2){
            xs2[i2] = xs[j];
            it2[i2] = iter[j];
            i2++;
        }
    }
    cwrite_row_double("root0.dat", xs0, 2);
    cwrite_row_double("root0.dat", it0);
    cwrite_row_double("root1.dat", xs1, 2);
    cwrite_row_double("root1.dat", it1);
    cwrite_row_double("root2.dat", xs2, 2);
    cwrite_row_double("root2.dat", it2);
}

