    //----------------------------------------------------------
// Abel Rosado Peinado - 5265 - 21 de enero de 2022
//----------------------------------------------------------
#include "ARPLibs/header.hpp"

double omega = 7.27e-5; double psi = M_PI/4.; double l = 20; double g = 9.8; //double k = sqrt(g/l);
double t0 = 0; double tf = 300; Vector x0 = {5,0,0,0}; double P = 9;
// The position-velocity vector consists of [x,dx/dt,y,dy/dt]

Vector F(double t, Vector x){
  Vector F = VecFull(0, 4);
  F[0] = x[1];
  F[1] = 2*omega*sin(psi)*x[3]-(g/l)*x[0];
  F[2] = x[3];
  F[3] = -2*omega*sin(psi)*x[1]-(g/l)*x[2];
  return F;
} // Differential Equation, derivative of the position-velocity vector

int main(){
  int nf; Matrix X = RK_n0tol(F, t0,tf, x0, nf, 1e-6,100);
  Vector t = linspace(t0,tf,nf,1);
  write_col_double("ftime.dat", t);
  write_mat_double("foucault.dat", X, 0, 7);
  Vector asn = VecFull(0,8832); Vector tP = asn;
  for (int i = 0; i < 8832; i++) {
    tP[i] = t[i]/P;
    asn[i] = asin(X[i][2]/X[i][0]);
  }
  double theta = IntTrVec(tP, asn);
  cout << theta/3 << endl;
}
