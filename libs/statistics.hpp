//------------------------------------------
// Abel Rosado - 2021
//------------------------------------------

// requires <math.h>

double cuadraticmean(double data[], int lines){
    double csum{};
    for (int i = 0; i < lines; i++){
        csum += data[i]*data[i];
    }
    double cm = sqrt(csum/lines);
    return cm;
}

double harmonicmean(double data[], int lines){
    double isum{};
    for (int i = 0; i < lines; i++){
        isum += 1/data[i];
    }
    double hm = lines/isum;
    return hm;
}
