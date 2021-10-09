// calculates the exponential function up to a certain tolerance
double exp (double x){
    int i = 0;          //  begin iteration
    double sum = 1.0;   //  the first term of the sum
    double ex  = 1.0;   //  assuming we begin with exp(x)=1
    while (sum > 0.00001){
        sum = sum*x/i;  //  iterative calculation of the sum, based on the previous one
        ex += sum;      //  sum it to ex
        i++;
    }
    return ex;
}