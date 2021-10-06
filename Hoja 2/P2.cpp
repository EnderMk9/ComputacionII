#include <iostream>
using namespace std;

int main(){
    for (int i = 1; i <= 10; i++){
        cout << i*10 << " ";
    } cout << endl;
    int j = 1;
    while (j <=10){
        cout << j*10 << " "; j++;
    } cout << endl;
    int k = 1;
    do{
        cout << k*10 << " "; k++;
    } while (k <= 10); cout << endl;
}