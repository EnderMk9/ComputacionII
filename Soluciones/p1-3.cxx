
#include <iostream>
#include <iomanip>
using namespace std;
int main()
{
  int fact1=1;
  float fact2=float(1.0);
  double fact3=double(1.0);
  cout << setprecision(17);
  cout << scientific;
  
  for (int i=1; i<=180; i++) {
    fact1=fact1*(i);
    fact2=fact2*float(i);
    fact3=fact3*double(i);
    cout << i << "\t" << fact1 << "\t" << fact2;
    cout << "\t" << fact3 << endl;
  }
  return 0;
}

