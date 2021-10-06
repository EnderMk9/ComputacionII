#include <iostream>
#include <fstream> // Library for reading and writing files
#include <math.h>
using namespace std;

int main(){
    ifstream rdata ("rdata.txt");  // set read file
    ofstream wdata ("wdata.txt");  // set write file
    float line;
    while (rdata >> line){   // read every line
      if (line > 1){         // if the number is greater than 1
          wdata << sqrt(line) << endl; // write the sqrt of the number to wdata
      }  
    }
}