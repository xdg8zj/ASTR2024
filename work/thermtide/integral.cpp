// Online C++ compiler to run C++ program online
//uses trapezoid rule to approximate integral
#include <iostream>
#include <cmath>
#include <list>
#include <vector>
using namespace std;

double calculate_int() {
     // Write C++ code here
    
    //set x coordinate vector:
    vector <double> x  = {0,2,4,6,8,10};
    int x_size = x.size();
    cout << "Number of x values:" << x_size << endl;
    //make y vectors:
    vector <double> y = {};
    y.resize(x.size());
    // cout << x[2] << endl;
    // make loop that assigns function of x to y vector
    for (int i = 0; i < x.size(); i++) {
         y[i] = 2*x[i];
        cout << "Y output:" << y[i] << endl;
       
        
    }
    
    //first x value
    double a = x[0];
    // cout << a;
    //last x value
    double b = x[5];
    // cout << b;
    //evenly spaced sections
    double n = 5;
    //make deltax or h depending on where you look
    double deltax = (b-a)/n;
    cout << "Delta x:" << deltax<< endl;
    
    //actual integral:
    double running_total = 0;
    for(int i = 0; i < y.size() ; i++ ) { //iterating through y values
        if(i == 0 || i == y.size()-1) { //if first or last:
           running_total = running_total + y[i]; //add it to the total
        }
        else { //for all values in between
            running_total = running_total + 2*(y[i]); //multiply them by two and sum them
        }
    }
    double integral = 0.5*deltax*running_total; //using trapezoid method: 1/2(deltax)*(y(0)+2(y(1)+y(2)...y(n-1))+y(n))
    
    return(integral);
    
}

double& integral_ref() {
    double& ref_int = calculate_int();
    return ref_int;
}
    


int main() {
    double integral_calc = calculate_int();
    cout << "Area under the curve:" << integral_calc << endl;
   
    
    


    return 0;
}



