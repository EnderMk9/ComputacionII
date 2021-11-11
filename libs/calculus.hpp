// Library for derivatives and integrals
// requires "linearalgebra.hpp" for matrices in multivariable calculus

// Calculates derivative using a loop to have an error less than the tolerance
double derivativeloop(double (*f)(double), double x0, double h0, double tolerance){
    double h   = h0;            // initial variation of x
    double err = 2*tolerance;   // initual value of error for the while to work
    double dp{}; double d;      // each term in the succesion
    while (err > tolerance){    // while the difference between terms is greater that the tolerance
        d = (f(x0+h)-f(x0-h))/(2*h);    // approximate value of the derivative using second order formula
        err = abs(d-dp);                // calculate the difference between the previous iteration
        h = h/2; dp = d;                // make h smaller and update dp
    }
    return d;                           
}

// Calculates the derivative using second order formula for a fixed variation of x
double derivative(double (*f)(double), double x0, double h){
    return (f(x0+h)-f(x0-h))/(2*h); // approximate value of the derivative using second order formula
}

// partial of a R^N to R^N function with respect to its j-th coordinate, at a point x_0
Vector partialVV(int n, int j, std::function<Vector (Vector&)> f, Vector& x0,double h0, double tolerance){
    Vector hv(n,0);              // h variation vector
    double h = h0;               // magnitude of hv
    hv[j] = h;                   // initial variation of x_j
    double err = 2*tolerance;    // initual value of error for the while to work
    Vector d(n,0); Vector dp(n,0);
    Vector dsum(n,0); Vector ddiff(n,0); double hnorm; double error; Vector deltad(n,0);
    Vector f1(n,0); Vector f2(n,0); Vector df(n,0);
    // definition of variables used
    while (err > tolerance){    // while the difference between terms is greater that the tolerance
        dsum = VecSum(x0,hv); ddiff = VecDiff(x0,hv); // x0+h and x0-h                         
        f1 = f(dsum); f2 = f(ddiff); df = VecDiff(f1,f2);
        d = ScalMult(df,1/(2*h));    // approximate value of the derivative using second order formula
        deltad = VecDiff(d, dp);
        err = norm(deltad);     // calculate the difference between the previous iteration
        h = h/2; hv[j] = h;     // make h smaller
        dp = d;                 // update dp
    }
    return d;
}


// Calculates the Jacobian using numerical methods for an R^n to R^n function
Matrix JacobianNum(int n, Vector& x0, std::function<Vector (Vector&)> f, double h0, double tolerance){
    Matrix J( n,vector<double>(n,0)); // Declaration of the variable
    for (int j = 0; j < n; j++){
        J[j] = partialVV(n,j,f,x0,h0,tolerance); // Each row is the partial derivative of f with
                                                 // respect to the j-th coordinate
    }
    J = transposeSqr(J);       // Transpose the matrix because rows must vary the coordenate
    return J;
}