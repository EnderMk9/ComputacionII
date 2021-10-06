#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

int main () {
    double x;
    int m=0;

    // Crear un flujo de entrada y llamarlo fentrada
    ifstream fentrada("datos.dat");
    if (!fentrada) {
        cout << "No se puede abrir el archivo de entrada, ";
        return 1;
    }
    
    // Crear un flujo de salida y llamarlo fsalida
    ofstream fsalida ("salida.dat");
    if (!fsalida)
    {
        cout << "No se puede abrir el archivo de salida  ";
        return 1;
    }

    // utilizacion de eof
    while (!fentrada.eof())
    {
        fentrada >> x;
        cout << x << " " << sqrt(x) << endl;
        if( x > 1){fsalida << sqrt(x) << endl;
            cout << " se graba" << endl;}
        m=m+1;
    }
    cout << "Se han leido " << m << " datos " << endl;
    fentrada.close(); // cerrar el archivo
    fsalida.close(); // cerrar el archivo
    return 0;
    }
