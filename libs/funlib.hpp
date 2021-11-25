//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// library for dealing with functions, evaluations, definitions, etc

// include <functional>
// use using namespace std::placeholders;

// function to evaluate a single variable function in an interval with d subintervals, d+1 divisions
// double t [d+1] {}; double y [d+1] {}; are requiered to be defined previously and passed as an imput
// there is no return because arrays are directly modified in memory inside functions
Vector eval(std::function<double(double)> f, Vector x){
    int n = x.size();
    Vector y(n,0);
    for (int i = 0; i < n; i++){
        y[i] = f(x[i]);
    }
    return y;
}
void longeval(std::function<long double(long double)> f, long double x0, long double xf, int d, Vector& x, Vector& y){
    long double step = (xf-x0)/d;
    long double xi; long double yi;
    for (int i = 0; i <= d; i++){
        xi = x0+i*step; yi = f(x0+i*step);
        x[i] = xi;
        y[i] = yi;
    }
}
