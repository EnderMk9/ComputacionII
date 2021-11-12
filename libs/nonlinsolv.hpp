//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// Library to solve f=0 for non-linear functions
// requires calculus.hpp for derivative
// requires linealgebra.hpp for multivariable functions

double bisection(double (*f)(double), double x1, double x2, double tol){
    if (f(x1)*f(x2)>0){
        cout << "ERROR POINTS NOT VALID" << endl; return NAN;}
    double x3; double err = 2*tol; //int i{};
    while (err > tol){
        x3 = (x1+x2)/2;
        if (f(x1)*f(x3) <0){
            x2 = x3;
        } else if (f(x2)*f(x3) < 0){
            x1 = x3;
        }
        err = abs(x2-x1);
        //i++;
    }
    //cout << "bisection iterations : " << i << endl;
    return x3;
}

double secant(double (*f)(double), double x1, double x2, double tol){
    double x3; double err = 2*tol; //int i{};
    while (err > tol){
        x3 = x2-f(x2)*((x2-x1)/(f(x2)-f(x1)));
        x1 = x2; x2 = x3; err = abs(x1-x2);
        //i++;
    }
    //cout << "secant iterations : " << i << endl;
    return x3;
}

double NewtonNum(double (*f)(double), double x0, double tol, double h, bool loop){
    double xp = x0; double x; double d; double y;
    double err = 2*tol; //int i{};
    if (loop){
        while (err > tol){
            y = f(xp); d = derivativeLoopCD(f, xp, h, tol);
            x = xp - y/d;
            err = abs(xp-x); xp = x;
            //i++;
        }
        //cout << "Newton devloop iterations : " << i << endl;
        return x;
    } else if (!loop){
        while (err > tol){
            y = f(xp); d = derivativeCD(f, xp, h);
            x = xp - y/d;
            err = abs(xp-x); xp = x;
            //i++;
        }
        //cout << "Newton single diff iterations : " << i << endl;
        return x;
    }
    return 0;
}

double NewtonAnal(double (*f)(double),double (*fp)(double), double x0, double tol){
    double xp = x0; double x; double d; double y;
    double err = 2*tol; //int i{};
    while (err > tol){
        y = f(xp); d = fp(xp);
        x = xp - y/d;
        err = abs(xp-x); xp = x;
        //i++;
    }
    //cout << "Newton analytic iterations : " << i << endl;
    return x;
}

Vector NewtonNumSys(int n,Vector& x0,std::function<Vector(Vector&)> f,double h, double tol,bool loop=0){
    Matrix J( n,vector<double>(n,0));
    double err = 2*tol;
    Vector x = x0; Vector xp = x0;
    Vector b(n,0); Vector d(n,0);
    Vector dx(n,0);
    if (loop){                                  // PartialLoop option
        while (err > tol){
            J = JacobianNum(n,x, f, h, tol,1);  // Calculate Jacobian
            b = f(x); b = ScalMult(b,-1);       // independent term
            d = LUSolve(J, b); x = VecSum(x,d); // Solution
            dx = VecDiff(x,xp); err = norm(dx); // error
            xp = x;                             // update xp
        }
    }else if (not loop){                        // Partial option
        while (err > tol){
            J = JacobianNum(n,x, f, h);         // Calculate Jacobian
            b = f(x); b = ScalMult(b,-1);       // independent term
            d = LUSolve(J, b); x = VecSum(x,d); // Solution
            dx = VecDiff(x,xp); err = norm(dx); // error
            xp = x;                             // update xp
        }
    }
    //coutmat(J);
    return x;
}

Vector NewtonAnalSys(int n,Vector& x0,std::function<Vector(Vector&)> f,std::function<Matrix(Vector&)> J, double tol){
    Matrix J0( n,vector<double>(n,0));
    double err = 2*tol;
    Vector x = x0; Vector xp = x0;
    Vector b(n,0); Vector d(n,0);
    Vector dx(n,0);
    while (err > tol){
        J0 = J(x);                           // Calculate Jacobian
        b = f(x); b = ScalMult(b,-1);        // independent term
        d = LUSolve(J0, b); x = VecSum(x,d); // Solution
        dx = VecDiff(x,xp); err = norm(dx);  // error
        xp = x;                              // update xp
    }
    return x;
}