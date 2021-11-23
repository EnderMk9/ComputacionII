#include "libs/header.hpp"

double m = 4000;
double k = 5000;

/*
m_i ddot(x_i) = -k_i (x_i-x_{i-1})+k_{i+1} (x_{i+1}-x_i)

m ddot(x_1) +2 k x_1 -k x_2 = 0
m ddot(x_2) +2 k x_2 -k x_1 - k x_3 = 0
m ddot(x_3) +2 k x_3 -k x_2 - k x_4 = 0
m ddot(x_4) +k (x_4-x_3) = 0

 */
Matrix K = {{ 2, -1, 0, 0},
            { -1, 2, -1, 0},
            { 0, -1, 2, -1},
            { 0, 0, -1, 1}};

Matrix M = {{ 1, 0, 0, 0},
            { 0, 1, 0, 0},
            { 0, 0, 1, 0},
            { 0, 0, 0, 1}};

int main(){
    K = ScalMultMat(K,k);
    Matrix MI = ScalMultMat(M,1./m);
    Matrix Kr = matprod(MI,K);
    coutmat(Kr);
    Matrix U; int k;
    Matrix D = JacobiDiag(Kr,U,k,1e-8);
    Vector l = MatrixDiag(D);
    coutvec(l);
    coutmat(U);
    Matrix UT = transpose(U);
    Matrix I = matprod(U,UT);
    coutmat(I);
}

// Si U X' = X
// ddot{x'_i} + \lambda_i x'_i = 0 --> x'_i = a e^{sqrt{-\lambda_i}t} + b e^{-sqrt{-\lambda_i}t}
// como \lambda_i > 0 --> x'_i = a cos(\omega_i t) + b sin(\omega_i t) tq \omega_i^2 = \lambda_i
// entonces X' = [a_i cos(\omega_i t)]_i + [b_i sin(\omega_i t)]_i
// y X = U ([a_i cos(\omega_i t)]_i +[b_i sin(\omega_i t)]_i)