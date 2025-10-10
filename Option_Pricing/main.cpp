

double K= 1.05;
double ttm = 1.0/3.0;
double r = 0.025;
double d = 0.02;
double vol = 0.21;
double S_0 =1.00;

#include <iostream>
#include <cmath>

int main(){
    double d1 = (log(S_0/K)+r-d+0.5*vol*vol*ttm)/(vol*sqrt(ttm));
    double d2=d1-vol*sqrt(ttm);
    double Nd1=0.5*(1+erf(d1/sqrt(2)));
    double Nd2=0.5*(1+erf(d2/sqrt(2)));
    double call_price = S_0*exp(-d*ttm)*Nd1-K*exp(-r*ttm)*Nd2;
    std::cout<<"The call price is: "<<call_price<<std::endl;
}
