//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// Library for derivatives and integrals
// requires "linearalgebra.hpp" for matrices in multivariable calculus

//-------------------------------------------------------------------------------------
// Derivatives
//-------------------------------------------------------------------------------------

// DERIVADA PRIMERA
// Usamos una fórmula de diferencias centradas f(x_0+h)-f(x_0-h) / 2h que tiene un error
// del orden O(h²), que puede demostrarse expandiendo en taylor f(x_0+h) y f(x_0-h) para
// h = 0, restar ambas expresiones y despejar f'(x_0).
// El error de truncado es siempre menor que h²/6 max_{x_0-h,x_0+h}|f⁽³⁾(x_0)|.
// El error de redondeo es siempre menor que e/h  max_{x_0-h,x_0+h}|f   (x_0)|.
// e es el error maximo que se puede producir el redondeo, que para el tamaño double
// tiene un valor de 5*2^-53 ≈ 6e-16, donde 53 es el numero de cifras significativas del tamaño.
// El valor óptimo de hace es la derivada respecto a h de la suma de los errores anteriores,
// teniendo en cuenta que si h es muy pequeño entonces
// max_{x_0-h,x_0+h}|f⁽³⁾(x_0)| ≈ |f⁽³⁾(x_0)| y max_{x_0-h,x_0+h}|f(x_0)| ≈ |f(x_0)| 
// Tal que entonces h* = \3th-root(3e|f(x_0)| / |f⁽³⁾(x_0)|)

// Calculates derivative using a loop to have an error less than the tolerance
double derivativeLoopCD(double (*f)(double), double x0, double h0, double tol){
    double h   = h0;            // initial variation of x
    double err = 2*tol;         // initual value of error for the while to work
    double dp{}; double d;      // each term in the succesion
    while (err > tol){          // while the difference between terms is greater that the tolerance
        d = (f(x0+h)-f(x0-h))/(2*h);
        err = abs(d-dp);                // calculate the difference between the previous iteration
        h = h/2; dp = d;                // make h smaller and update dp
    }
    return d;                           
}

// best used when knowing the order of magnitude of f⁽³⁾(x_0), because then you can
// estimate the error and optimize h and this function requieres less iterations
double derivativeCD(double (*f)(double), double x0, double h){
    return (f(x0+h)-f(x0-h))/(2*h); // approximate value of the derivative using second order formula
}

// DERIVADA SEGUNDA
// Usamos una fórmula de diferencias centradas f(x_0+h)+f(x_0-h) +2f(x_0) / h² que tiene un error
// del orden O(h²), que puede demostrarse expandiendo en taylor f(x_0+h) y f(x_0-h) para
// h = 0, sumar ambas expresiones y despejar f''(x_0)

// Calculates derivative using a loop to have an error less than the tolerance
double derivative2LoopCD(double (*f)(double), double x0, double h0, double tol){
    double h   = h0;            // initial variation of x
    double err = 2*tol;         // initual value of error for the while to work
    double dp{}; double d;      // each term in the succesion
    while (err > tol){          // while the difference between terms is greater that the tolerance
        d = (f(x0+h)+f(x0-h)-2*f(x0))/(h*h);
        err = abs(d-dp);                // calculate the difference between the previous iteration
        h = h/2; dp = d;                // make h smaller and update dp
    }
    return d;                           
}

// best used when knowing the order of magnitude of f⁽⁴⁾(x_0), because then you can
// estimate the error and optimize h and this function requieres less iterations
double derivative2CD(double (*f)(double), double x0, double h){
    return (f(x0+h)+f(x0-h)-2*f(x0))/(h*h);
}

// HIGHER DIMENSIONS

// Partial derivative of a scalar function f with respect to the i-th coordenate
double partialLoopCD(int n, double (*f)(Vector&), int i, Vector& x0, double h0, double tol){
    Vector hv(n,0);              // h variation vector
    double h = h0;               // magnitude of hv
    hv[i] = h;                   // initial variation of x_i
    double err = 2*tol;          // initual value of error for the while to work
    double dp{}; double d; Vector dsum(n,0); Vector ddif(n,0); double ff; double fb; // declaration
    while (err > tol){           // while the difference between terms is greater that the tolerance
        dsum = VecSum(x0,hv);           //x0+h
        ddif = VecDiff(x0,hv);          //x0-h                         
        ff = f(dsum); fb = f(ddif);     // eval in f
        d = (ff-fb)/(2*h);              // central difference formula
        err = abs(d-dp);                // calculate the difference between the previous iteration
        h = h/2; hv[i] = h;             // make h smaller
        dp = d;                         // update dp
    }
    return d;                           
}

// best used when knowing the order of magnitude of f⁽³⁾_i(x_0), because then you can
// estimate the error and optimize h and this function requieres less iterations
double partialCD(int n, double (*f)(Vector&), int i, Vector& x0, double h){
    Vector hv(n,0); hv[i] = h;                 // Variation
    Vector dsum = VecSum(x0,hv);               // x0+h
    Vector ddif = VecDiff(x0,hv);              // x0-h                         
    double ff = f(dsum); double fb = f(ddif);  // eval in f
    return (ff-fb)/(2*h); // approximate value of the derivative using second order formula
}

// partial of a R^N to R^N function with respect to its j-th coordinate, at a point x_0, using the
// central difference formula
Vector partialLoopCDVec(int n,int i,std::function<Vector (Vector&)> f,Vector& x0,double h0,double tol){
    Vector hv(n,0);              // h variation vector
    double h = h0;               // magnitude of hv
    hv[i] = h;                   // initial variation of x_i
    double err = 2*tol;    // initual value of error for the while to work
    Vector d(n,0); Vector dp(n,0);
    Vector dsum(n,0); Vector ddif(n,0); double hnorm; double error; Vector deltad(n,0);
    Vector ff(n,0); Vector fb(n,0); Vector df(n,0);
    // definition of variables used
    while (err > tol){    // while the step is greater that the tolerance
        dsum = VecSum(x0,hv);       // x0+h
        ddif = VecDiff(x0,hv);      // x0-h                           
        ff = f(dsum); fb = f(ddif); // evaluate
        df = VecDiff(ff,fb);        // difference of the functions
        d = ScalMult(df,1/(2*h));   // central difference formula
        deltad = VecDiff(d, dp);    // step taken
        err = norm(deltad);         // calculate the magnitud of the step
        h = h/2; hv[i] = h;         // make h smaller
        dp = d;                     // update dp
    }
    return d;
}

Vector partialCDVec(int n,int i,std::function<Vector (Vector&)> f,Vector& x0,double h){
    Vector hv(n,0);  hv[i] = h;               // Variation
    Vector dsum = VecSum(x0,hv);              // x0+h
    Vector ddif = VecDiff(x0,hv);             // x0-h                           
    Vector ff = f(dsum); Vector fb = f(ddif); // evaluate
    Vector df = VecDiff(ff,fb);               // difference of the functions
    Vector d = ScalMult(df,1/(2*h));          // central difference formula
    return d;
}

// Calculates the Jacobian using numerical methods for an R^n to R^n function
Matrix JacobianNum(int n,Vector& x0,std::function<Vector(Vector&)> f,double h,double tol=0,bool loop=0){
    Matrix J( n,vector<double>(n,0)); // Declaration of the variable
    if (loop){
        for (int j = 0; j < n; j++){
            J[j] = partialLoopCDVec(n,j,f,x0,h,tol);
        }
    }else if (not loop){
        for (int j = 0; j < n; j++){
            J[j] = partialCDVec(n,j,f,x0,h);
        }
    }
    J = transposeSqr(J);       // Transpose the matrix because rows must vary the coordenate
    return J;
}