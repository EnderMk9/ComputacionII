#include <iostream>
using namespace std;

int arr [7] = {34,70,5,89,6,3220,1};

int main(){
    int n = *(&arr + 1) - arr;
    for(int i = 0; i < n ; i++){
        for(int j = 0; j < n - i; j++){
            if (arr[j]>arr[j+1]){
                swap(arr[j],arr[j+1]);
            };
        }
    }
    for (int i = 0; i < n; ++i){
            cout << arr[i] << " ";
        } cout << endl;
}
