#include "libs/header.hpp"
#include "libs/gnuplot.hpp"

double p(double r){
    return -2./r;
}
double qr(double x){
    return 0;
}

int main(){
    double x0 = 0.05; double xf = 0.1;
    double y0 = 110;  double yf = 0;
    Vector y1 = FiniteDifferences(p,qr,qr,10,x0,xf,y0,yf);
    Vector x1 = linspace(x0,xf,10+2,1);
    write_2col_double("dat1.dat",x1,y1);
    Vector y2 = FiniteDifferences(p,qr,qr,100,x0,xf,y0,yf);
    Vector x2 = linspace(x0,xf,100+2,1);
    write_2col_double("dat2.dat",x2,y2);
    Vector y3 = FiniteDifferences(p,qr,qr,1000,x0,xf,y0,yf);
    Vector x3 = linspace(x0,xf,1000+2,1);
    write_2col_double("dat3.dat",x3,y3);
    //GnuplotPipe gp; gp.sendLine("plot 'dat3.dat'");
}
