#include <iostream> 
#include <math.h>

using namespace std;

double pi = 3.14159265358979;  // Real value for comparison
float   sum1as = 0, sum1bs = 0, sum2as = 0, sum2bs = 0;  // simple precision variables
double  sum1ad = 0, sum1bd = 0, sum2ad = 0, sum2bd = 0;  // double precision variables
// a means forward sum and b backward, 1 mean 10⁴ terms and 2 means 10⁶ terms
int n1 = 10000, n2 = 1000000; // 10⁴ and 10⁶

int main()
{
    for (long int i = 1; i <= n1; i++){
        sum1as += 1.0f/(i*i);}   // Simple precision, forward,  10⁴
    for (long int i = n1; i > 0; i--){
        sum1bs += 1.0f/(i*i);}   // Simple precision, backward, 10⁴
    for (long int i = 1; i <= n2; i++){
        sum2as += 1.0f/(i*i);}   // Simple precision, forward,  10⁶
    for (long int i = n2; i > 0; i--){
        sum2bs += 1.0f/(i*i);}   // Simple precision, backward, 10⁶
    for (long int i = 1; i <= n1; i++){
        sum1ad += 1.0/(i*i);}    // Double precision, forward,  10⁴
    for (long int i = n1; i > 0; i--){
        sum1bd += 1.0/(i*i);}    // Double precision, backward, 10⁴
    for (long int i = 1; i <= n2; i++){
        sum2ad += 1.0/(i*i);}    // Double precision, forward,  10⁶
    for (long int i = n2; i > 0; i--){
        sum2bd += 1.0/(i*i);}    // Double precision, backward, 10⁶
    float pi1as = sqrt(6*sum1as);
    float pi1bs = sqrt(6*sum1bs);
    float pi2as = sqrt(6*sum2as);
    float pi2bs = sqrt(6*sum2bs);
    float pi1ad = sqrt(6*sum1ad);
    float pi1bd = sqrt(6*sum1bd);
    float pi2ad = sqrt(6*sum2ad);
    float pi2bd = sqrt(6*sum2bd);
    // Calculate pi
    cout.precision(7);  // Float precision
    cout << "Simple, forward,  n = 10⁴; π is equal to: " << pi1as << "; error = " << 100*(pi1as-pi)/pi << " %" << endl;
    cout << "Simple, backward, n = 10⁴; π is equal to: " << pi1bs << "; error = " << 100*(pi1bs-pi)/pi << " %" << endl;
    cout << "Simple, forward,  n = 10⁶; π is equal to: " << pi2as << "; error = " << 100*(pi2as-pi)/pi << " %" << endl;
    cout << "Simple, backward, n = 10⁶; π is equal to: " << pi2bs << "; error = " << 100*(pi2bs-pi)/pi << " %" << endl;
    cout.precision(16); // Double precision
    cout << "Double, forward,  n = 10⁴; π is equal to: " << pi1ad << "; error = " << 100*(pi1ad-pi)/pi << " %" << endl;
    cout << "Double, backward, n = 10⁴; π is equal to: " << pi1bd << "; error = " << 100*(pi1bd-pi)/pi << " %" << endl;
    cout << "Double, forward,  n = 10⁶; π is equal to: " << pi2ad << "; error = " << 100*(pi2ad-pi)/pi << " %" << endl;
    cout << "Double, backward, n = 10⁶; π is equal to: " << pi2bd << "; error = " << 100*(pi2bd-pi)/pi << " %" << endl;
}