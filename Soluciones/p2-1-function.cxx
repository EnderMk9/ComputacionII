// Cerdeno 07/09/2020

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;


// Apparently we need to define the functions before they are used in main
// Otherwise we need to declare them first, then write main, then the function

double exponencial(double);



int main()
{
    double x;
    
    cout << "introduce valor de x: ";
    cin >> x;
    
    cout << exponencial(x);
    cout << "error=" << fabs(exponencial(x)-exp(x)) << endl;

}


double exponencial(double x)
{
    // Iterative solution to the exponential function
    double tolerancia, suma, termino;
    int contador;
    
    tolerancia = 0.00001;
    suma = 0.0;
    termino = 1;
    
    contador = 1;
    while (termino > tolerancia) {
        suma = suma + termino;
        termino = termino*(x/contador);
        contador = contador + 1;
    }

    return suma;
}
