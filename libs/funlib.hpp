// library for dealing with functions, evaluations, definitions, etc

// function to evaluate a single variable function in an interval with d subintervals, d+1 divisions
// double t [d+1] {}; double y [d+1] {}; are requiered to be defined previously and passed as an imput
// there is no return because arrays are directly modified in memory inside functions
void eval(std::function<double(double)> f, double x0, double xf, int d,double x[], double y[]){
    double step = (xf-x0)/d;
    double xi; double yi;
    for (int i = 0; i <= d; i++){
        xi = x0+i*step; yi = f(x0+i*step);
        x[i] = xi;
        y[i] = yi;
    }
}
void longeval(std::function<long double(long double)> f, long double x0, long double xf, int d, long double x[], long double y[]){
    long double step = (xf-x0)/d;
    long double xi; long double yi;
    for (int i = 0; i <= d; i++){
        xi = x0+i*step; yi = f(x0+i*step);
        x[i] = xi;
        y[i] = yi;
    }
}
