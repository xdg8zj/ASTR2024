//here's a slight edit to use a function where you modify the input variable using pass-by-reference:
 
#include <iostream>
#include <cmath>
#include <list>
#include <vector>
using namespace std;
 
// \int_0^1 dx*2*x = 1
 
void calculate_int(int n, double &integral) {
  vector <double> x(n);
  vector <double> y(n);
  double dx = 1.0/n;
  for (int i = 0; i < n; i++) {
    x[i] = dx*i;
    y[i] = pow(x[i],2);
  }
  double running_total = dx*0.5*(y[0]+y[n-1]);
  for(int i = 1; i < n-1; i++ ) {
    running_total = running_total + dx*y[i];
  }
  integral = running_total;
}
 
int main() {
  double exact=1.0/3.0;
  double integral=9.0;
  cout << "starting value =" << integral << endl;
  for (int n=2; n<1000; n=2*n){
    calculate_int(n,integral);
    cout << "n,exact,calc,err = " << n << ", " << exact << ", " << integral << ", " << integral-exact << endl;
  }
    return 0;
}

