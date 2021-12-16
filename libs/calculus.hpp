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
double derivativeLoopCD(std::function<double(double)> f, double x0, double h0, double tol){
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
double derivativeCD(std::function<double(double)> f, double x0, double h){
    return (f(x0+h)-f(x0-h))/(2*h); // approximate value of the derivative using second order formula
}

// MUST BE EQUISPACED
Vector derivativeCDArr(Vector& x, Vector& y){
    int xs = x.size(); int ys = y.size();
    if (xs != ys){
        cout << "ERROR SIZES NOT COMPATIBLE" << endl;
        return{};
    }
    Vector dy(xs,0);
    dy[0] = (4*y[1]-3*y[0]-y[2])/(2*(x[1]-x[0]));
    for (int i = 1; i < xs-1; i++){
        dy[i] = (y[i+1]-y[i-1])/(2*(x[i+1]-x[i]));
    }
    dy[xs-1] = (y[xs-1]-y[xs-2])/(x[xs-1]-x[xs-2]);
    return dy;
}

// DERIVADA SEGUNDA
// Usamos una fórmula de diferencias centradas f(x_0+h)+f(x_0-h) +2f(x_0) / h² que tiene un error
// del orden O(h²), que puede demostrarse expandiendo en taylor f(x_0+h) y f(x_0-h) para
// h = 0, sumar ambas expresiones y despejar f''(x_0)

// Calculates derivative using a loop to have an error less than the tolerance
double derivative2LoopCD(std::function<double(double)> f, double x0, double h0, double tol){
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
double derivative2CD(std::function<double(double)> f, double x0, double h){
    return (f(x0+h)+f(x0-h)-2*f(x0))/(h*h);
}

// MUST BE EQUISPACED
Vector derivative2CDArr(Vector& x, Vector& y){
    int xs = x.size(); int ys = y.size();
    if (xs != ys){
        cout << "ERROR SIZES NOT COMPATIBLE" << endl;
        return{};
    }
    Vector ddy(xs,0);
    ddy[0] = (y[2]-2*y[1]+y[0])/((x[1]-x[0])*(x[1]-x[0]));
    for (int i = 1; i < xs-1; i++){
        ddy[i] = (y[i+1]+y[i-1]-2*y[i])/((x[i+1]-x[i])*(x[i+1]-x[i]));
    }
    ddy[xs-1] = (y[xs-1]-2*y[xs-2]+y[xs-3])/((x[xs-1]-x[xs-2])*(x[xs-1]-x[xs-2]));
    return ddy;
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

//-------------------------------------------------------------------------------------
// Integrals
//-------------------------------------------------------------------------------------

double DefIntTrap(std::function<double(double)> f, double a, double b, int n){
    double I; double h = (b-a)/n;
    I = (h/2)*(f(a)+f(b));
    for (int i = 1; i < n; i++){
        I += h*f(a+i*h);
    }
    return I;
}

// n even
double DefIntSimp13(std::function<double(double)> f, double a, double b, int n){
    double I; double h = (b-a)/n;
    I = (h/3)*(f(a)+f(b));
    for (int i = 1; i < n; i = i+2){
        I += (4./3)*h*f(a+i*h);
    }
    for (int i = 2; i < n-1; i = i+2){
        I += (2./3)*h*f(a+i*h);
    }
    return I;
}

// n multiple of 3
double DefIntSimp38(std::function<double(double)> f, double a, double b, int n){
    double I; double h = (b-a)/n;
    I = (3*h/8)*(f(a)+f(b));
    for (int i = 1; i < n; i = i+3){
        I += (9./8)*h*f(a+i*h);
    }
    for (int i = 2; i < n; i = i+3){
        I += (9./8)*h*f(a+i*h);
    }
    for (int i = 3; i < n; i = i+3){
        I += (6./8)*h*f(a+i*h);
    }
    return I;
}

//
double DefIntSimp(std::function<double(double)> f, double a, double b, int n){
    if (n % 2 == 0){
        return DefIntSimp13(f,a,b,n);
    } else if (n % 2 != 0){
        double I; double h = (b-a)/n;
        I = DefIntSimp13(f,a,b-3*h,n-3);
        I+= (3*h/8.)*(f(b-3*h)+f(b)+3*f((3*b-6*h)/3.)+3*f((3*b-3*h)/3.));
        return I;
    }
    return {};
}

Matrix Lw = {{2,0,0,0,0},
            {1,1,0,0,0},
            {8./9,5./9,5./9,0,0},
            {(18+sqrt(30))/(36.),(18+sqrt(30))/(36.),(18-sqrt(30))/(36.),(18-sqrt(30))/(36.),0},
            {128./225,(322+13*sqrt(70))/(900.),(322+13*sqrt(70))/(900.),(322-13*sqrt(70))/(900.),(322-13*sqrt(70))/(900.)}};
Matrix Lxi = {{0,0,0,0,0},
             {1./sqrt(3),-1./sqrt(3),0,0,0},
             {0,sqrt(3./5),-sqrt(3./5),0,0},
             {sqrt(3./7 - (2./7)*sqrt(6./5)),-sqrt(3./7 - (2./7)*sqrt(6./5)),sqrt(3./7 + (2./7)*sqrt(6./5)),-sqrt(3./7 + (2./7)*sqrt(6./5)),0},
             {0,(1./3)*sqrt(5-2*sqrt(10./7)),-(1./3)*sqrt(5-2*sqrt(10./7)),(1./3)*sqrt(5+2*sqrt(10./7)),-(1./3)*sqrt(5+2*sqrt(10./7))}};

double DefIntGaussLegendre(std::function<double(double)> f, double a, double b, int n = 4){
    double I{}; double p = (b-a)/2; double q = (b+a)/2;
    for (int i = 0; i < n; i++){
        I += p*Lw[n-1][i]*f(p*Lxi[n-1][i]+q);
    }
    return I;
}

//-------------------------------------------------------------------------------------
// EDO's
//-------------------------------------------------------------------------------------

Matrix EulerFO(std::function<double(double,double)> yp, double x0, double xf,double y0, int n){
    Vector x = linspace(x0,xf,n);
    double h = (xf-x0)/(n-1);
    Vector y = VecFull(0,n); y[0] = y0;
    for (int i = 1; i < n; i++){
        y[i] = y[i-1]+h*yp(x[i-1],y[i-1]);
    }
    Matrix xy = MatFull(0,2,n);
    xy[0] = x; xy[1] = y;
    return xy;
}

Matrix RungeKutta4FO(std::function<double(double,double)> yp, double x0, double xf,double y0, int n){
    Vector x = linspace(x0,xf,n);
    double h = (xf-x0)/(n-1);
    Vector y = VecFull(0,n); y[0] = y0;
    double k1; double k2; double k3; double k4;
    for (int i = 1; i < n; i++){
        k1 = yp(x[i-1],y[i-1]);
        k2 = yp(x[i-1]+h/2.,y[i-1]+(h*k1)/2.);
        k3 = yp(x[i-1]+h/2.,y[i-1]+(h*k2)/2.);
        k4 = yp(x[i-1]+h,y[i-1]+h*k3);
        y[i] = y[i-1]+h*(k1+2*k2+2*k3+k4)/6.;
    }
    Matrix xy = MatFull(0,2,n);
    xy[0] = x; xy[1] = y;
    return xy;
}

Matrix RungeKutta4SO(std::function<double(double,double,double)> f,std::function<double(double,double,double)> g, double x0, double xf, double y0, double z0, int n){
    Vector x = linspace(x0,xf,n);
    double h = (xf-x0)/(n-1);
    Vector y = VecFull(0,n); y[0] = y0;
    Vector z = VecFull(0,n); z[0] = z0;
    double k1; double k2; double k3; double k4;
    double l1; double l2; double l3; double l4;
    for (int i = 1; i < n; i++){
        k1 = f(x[i-1],y[i-1],z[i-1]);
        l1 = g(x[i-1],y[i-1],z[i-1]);
        k2 = f(x[i-1]+h/2.,y[i-1]+(h*k1)/2.,z[i-1]+(h*l1)/2.);
        l2 = g(x[i-1]+h/2.,y[i-1]+(h*k1)/2.,z[i-1]+(h*l1)/2.);
        k3 = f(x[i-1]+h/2.,y[i-1]+(h*k2)/2.,z[i-1]+(h*l2)/2.);
        l3 = g(x[i-1]+h/2.,y[i-1]+(h*k2)/2.,z[i-1]+(h*l2)/2.);
        k4 = f(x[i-1]+h,y[i-1]+h*k3, z[i-1]+h*l3);
        l4 = g(x[i-1]+h,y[i-1]+h*k3, z[i-1]+h*l3);
        y[i] = y[i-1]+h*(k1+2*k2+2*k3+k4)/6.;
        z[i] = z[i-1]+h*(l1+2*l2+2*l3+l4)/6.;
    }
    Matrix xyz = MatFull(0,3,n);
    xyz[0] = x; xyz[1] = y; xyz[2] = z;
    return xyz;
}

Matrix RungeKutta4nO(std::function<Vector(double,Vector&)> F, double t0, double tf, Vector x0, int n){
    Vector t = linspace(t0,tf,n);
    double h = (tf-t0)/(n-1); int m = x0.size(); double prevt; Vector prevx;
    Matrix x = MatFull(0,n,m); x[0] = x0;
    Matrix k = MatFull(0,4,m); Vector sum = x0; Vector prod = x0;
    for (int i = 1; i < n; i++){
        prevt = t[i-1];
        prevx = x[i-1];
        k[0] = F(prevt,prevx);
        prod = ScalMult(k[0],h/2.);
        sum  = VecSum(prevx,prod);
        k[1] = F(prevt+h/2.,sum);
        prod = ScalMult(k[1],h/2.);
        sum  = VecSum(prevx,prod);
        k[2] = F(prevt+h/2.,sum);
        prod = ScalMult(k[1],h);
        sum  = VecSum(prevx,prod);
        k[3] = F(prevt+h,sum);
        prod = ScalMult(k[1],2);
        sum  = VecSum(k[0],prod);
        prod = ScalMult(k[2],2);
        sum  = VecSum(sum,prod);
        sum  = VecSum(sum,k[3]);
        prod = ScalMult(sum,h/6.);
        x[i] = VecSum(prevx,prod);
    }
    return x;
}

double ShootRK4SO(std::function<double(double,double,double)> f,std::function<double(double,double,double)> g, double x0, double xf, double y0, double z0, int n){
    Vector x = linspace(x0,xf,n);
    double h = (xf-x0)/(n-1);
    Vector y = VecFull(0,n); y[0] = y0;
    Vector z = VecFull(0,n); z[0] = z0;
    double k1; double k2; double k3; double k4;
    double l1; double l2; double l3; double l4;
    for (int i = 1; i < n; i++){
        k1 = f(x[i-1],y[i-1],z[i-1]);
        l1 = g(x[i-1],y[i-1],z[i-1]);
        k2 = f(x[i-1]+h/2.,y[i-1]+(h*k1)/2.,z[i-1]+(h*l1)/2.);
        l2 = g(x[i-1]+h/2.,y[i-1]+(h*k1)/2.,z[i-1]+(h*l1)/2.);
        k3 = f(x[i-1]+h/2.,y[i-1]+(h*k2)/2.,z[i-1]+(h*l2)/2.);
        l3 = g(x[i-1]+h/2.,y[i-1]+(h*k2)/2.,z[i-1]+(h*l2)/2.);
        k4 = f(x[i-1]+h,y[i-1]+h*k3, z[i-1]+h*l3);
        l4 = g(x[i-1]+h,y[i-1]+h*k3, z[i-1]+h*l3);
        y[i] = y[i-1]+h*(k1+2*k2+2*k3+k4)/6.;
        z[i] = z[i-1]+h*(l1+2*l2+2*l3+l4)/6.;
    }
    return y[n-1];
}

double ShootSO(std::function<double(double,double,double)> f,std::function<double(double,double,double)> g, double x0, double xf, double y0, double yf, int n, double tol){
    double p0 = 0.5*(yf-y0)/(xf-x0); double p1; double p2;
    double err0; double err1; double err2; double y;
    y = ShootRK4SO(f,g,x0,xf,y0,p0,n); err0 = y-yf;
    if (p0*err0 > 0){
        p1 = p0/2.;
    } else if(p0*err0 < 0){
        p1 = 2*p0;
    } else if(p0 == 0 && err0 > 0){
        p1 = -1;
    } else if(p0 == 0 && err0 < 0){
        p1 = 1;
    }
    y = ShootRK4SO(f,g,x0,xf,y0,p1,n); err1 = y-yf;
    while (err1 > tol){
        p2 = p1 - err1*(p1-p0)/(err1-err0);
        y = ShootRK4SO(f,g,x0,xf,y0,p2,n); err2 = y-yf;
        p0 = p1; p1 = p2; err0 = err1; err1 = err2;
    }
    return p1;
}