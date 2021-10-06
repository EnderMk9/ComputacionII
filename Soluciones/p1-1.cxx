// Cerdeno 21/08/2020
//
// Sumas parciales para estudiar propiedades
// de aritmética con coma flotante

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    // Define all variables locally
    float suma1as=0.0, suma1bs=0.0, suma2as=0.0, suma2bs=0.0;
    double suma1ad=0.0, suma1bd=0.0, suma2ad=0.0, suma2bd=0.0;
    float PI=4.0*atan(1.0);
    double DPI=4.0*atan(1.0);
    const int n=2;

    // Sumas hasta 10000 términos (suma 1)
    // Hacia adelante (a)
    for (int i=1; i<=10000; i++) {
        // Simple precision
      suma1as=suma1as+(1.0/pow(float(i),n));
        // Double precision
      suma1ad=suma1ad+(1.0/pow(double(i),n));
    }
    // Hacia atras (b)
    for (int i=10000; i>=1; i--) {
        // Simple precision
      suma1bs=suma1bs+(1.0/pow(float(i),n));
        // Double precision
      suma1bd=suma1bd+(1.0/pow(double(i),n));
    }
    
    
    // Sumas hasta 1000000 términos (suma 2)
    for (int i=1; i<=1000000; i++) {
      // Simple precision
      suma2as=suma2as+(1.0/pow(float(i),n));
      // Double precision
      suma2ad=suma2ad+(1.0/pow(double(i),n));
    }
    // Hacia atras (b)
    for (int i=1000000; i>=1; i--) {
        // Simple precision
      suma2bs=suma2bs+(1.0/pow(float(i),n));
        // Double precision
      suma2bd=suma2bd+(1.0/pow(double(i),n));
    }
    
    // Setting the precision to print on screen to 16 digits
    cout << setprecision(16);
    
    //
    cout << "SUMA 100000 TERMINOS" << endl;
    cout << "SOLUCION EXACTA= " << pow(PI,n)/6.0  << endl;
    cout << "- Hacia adelante simple= " << suma1as << " -- Error= " << pow(PI,n)/6.0-suma1as << endl;
    cout << "- Hacia adelante doble = " << suma1ad << " -- Error= " << pow(PI,n)/6.0-suma1ad << endl;
    cout << "- Hacia atras simple= " << suma1bs << " -- Error= " << pow(PI,n)/6.0-suma1bs << endl;

    cout << "- Hacia atras doble = " << suma1bd << " -- Error= " << pow(PI,n)/6.0-suma1bd << endl;

    cout << "SUMA 10000 TERMINOS" << endl;
    cout << "SOLUCION EXACTA= " << pow(DPI,n)/6.0 << endl;
    cout << "- Hacia adelante simple= " << suma2as << " -- Error= " << pow(PI,n)/6.0-suma2as << endl;
    cout << "- Hacia adelante doble = " << suma2ad << " -- Error= " << pow(PI,n)/6.0-suma2ad << endl;
    cout << "- Hacia atras simple= " << suma2bs << " -- Error= " << pow(PI,n)/6.0-suma2bs << endl;
    cout << "- Hacia atras doble = " << suma2bd << " -- Error= " << pow(PI,n)/6.0-suma2bd << endl;


    // This is one of the surprising aspects of floating-point arithmetic: it actually matters what order you do things like addition in. (Formally, we say that floating-point addition is not commutative.)
    
    // Another example, easier to visualize
    // single precision = 7 digits of precision
    // 1000000.0 + 0.1 + 0.1 + 0.1 .....
    // 0.1 + 0.1 + 0.1 ... + 1000000.0
    // Are not the same

    float f1 = 1000000.0, f2 = 0.0;
    int i;
    for(i = 0; i < 100; i++) {
        f1 += 0.1;
        f2 += 0.1;
    }
    f2 += 1000000.0;
    printf("%.1f %.1f\n", f1, f2);

    
    return 0;

}
