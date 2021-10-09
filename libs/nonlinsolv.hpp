// Library to solve f=0 for non-linear functions
// requires calculus.h for derivative

double bisection(double (*f)(double), double x1, double x2, double tolerance){
    if (f(x1)*f(x2)>0){
        cout << "ERROR POINTS NOT VALID" << endl; return NAN;}
    double x3;
    double err = 2*tolerance;
    int i{};
    while (err > tolerance){
        x3 = (x1+x2)/2;
        if (f(x1)*f(x3) <0){
            x2 = x3;
        } else if (f(x2)*f(x3) < 0){
            x1 = x3;
        }
        err = abs(x2-x1);
        i++;
    }
    cout << "bisection iterations : " << i << endl;
    return x3;
}

double secant(double (*f)(double), double x1, double x2, double tolerance){
    double x3; double err = 2*tolerance;
    int i{};
    while (err > tolerance){
        x3 = x2-f(x2)*((x2-x1)/(f(x2)-f(x1)));
        x1 = x2; x2 = x3; err = abs(x1-x2);
        i++;
    }
    cout << "secant iterations : " << i << endl;
    return x3;
}

double NewtonRhapsonNum(double (*f)(double), double x0, double tolerance, double h, bool loop){
    double xp = x0; double x; double d; double y;
    double err = 2*tolerance; int i{};
    if (loop){
        while (err > tolerance){
            y = f(xp); d = derivativeloop(f, xp, h, tolerance);
            x = xp - y/d;
            err = abs(xp-x); xp = x;
            i++;
        }
        cout << "Newton devloop iterations : " << i << endl;
        return x;
    } else if (!loop){
        while (err > tolerance){
            y = f(xp); d = derivative(f, xp, h);
            x = xp - y/d;
            err = abs(xp-x); xp = x;
            i++;
        }
        cout << "Newton single diff iterations : " << i << endl;
        return x;
    }
    return 0;
}

double NewtonRhapsonAnly(double (*f)(double),double (*fp)(double), double x0, double tolerance){
    double xp = x0; double x; double d; double y;
    double err = 2*tolerance; int i{};
    while (err > tolerance){
        y = f(xp); d = fp(xp);
        x = xp - y/d;
        err = abs(xp-x); xp = x;
        i++;
    }
    cout << "Newton analytic iterations : " << i << endl;
    return x;
}