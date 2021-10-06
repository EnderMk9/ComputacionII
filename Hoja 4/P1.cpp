#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;
#include "rwlib.h"
#include "statistics.h"

int main(){
    string rname = "rdata.txt";
    string wname = "wdata.txt";
    int ln = file_ln(rname);
    double data[ln];
    d_r_file_ln(rname,data);
    double cm = cuadraticmean(data,ln);
    double hm = harmonicmean(data,ln);
    cout << "Media cuadrática : " << cm << endl;
    cout << "Media armónica : " << hm << endl;
}