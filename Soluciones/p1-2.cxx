// Bubble algorithm to reorder a series
// Cerdeno 21/08/2020

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
const int n=7;

int main()
{
    int datos[n]={34,70,5,89,6,3220,1};
    int swap;
    
    for (int j=0; j<n-1; j++) {
        for (int i=0; i<n-j-1; i++) {
            if (datos[i]>datos[i+1]) {
                // intercambiamos los datos
                swap=datos[i];
                datos[i]=datos[i+1];
                datos[i+1]=swap;
                cout << "Intercambiados" << datos[i] << datos[i+1];
            }
        }
    }
    
    cout << "Despues de ordenar" << endl;
    for (int i=0;i<=n-1;i++) {
        cout << datos[i] << endl;
    }

}
