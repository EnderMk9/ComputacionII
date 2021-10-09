// Library for derivatives and integrals

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