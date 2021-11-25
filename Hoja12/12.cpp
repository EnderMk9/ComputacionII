#include "libs/header.hpp"

double f(double x){
    return 70*exp(-x/10.0)+7*x-70;
}

double df(double x){
    return 7-7*exp(-x/10.0);
}

double ddf(double x){
    return (7./10)*exp(-x/10.0);
}

int main(){
    Matrix rdata = read_matrix_double("posicion.txt",2);
    rdata = transpose(rdata);
    Vector t = rdata[0]; Vector x = rdata[1];
    Vector dx = derivativeCDArr(t,x);
    Vector ddx = derivative2CDArr(t,x);
    Vector xa = eval(f,t);
    Vector dxa = eval(df,t);
    Vector ddxa = eval(ddf,t);
    Vector ddiff = VecDiff(dxa,dx); ddiff = absVec(ddiff);
    Vector adxa = absVec(dxa); ddiff = VecDiv(ddiff,adxa);
    Vector dddiff = VecDiff(ddxa,ddx); dddiff = absVec(dddiff);
    Vector addxa = absVec(ddxa); dddiff = VecDiv(dddiff,addxa);
    Matrix wdata(9,Vector(t.size(),0));
    wdata[0] = t; wdata[1] = x; wdata[2] = dx; wdata[3] = ddx; wdata[4] = xa;
    wdata[5] = dxa; wdata[6] = ddxa; wdata[7] = ddiff; wdata[8] = dddiff;
    wdata = transpose(wdata);
    write_mat_double("out.dat", wdata);
}