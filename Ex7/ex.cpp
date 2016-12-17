#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "fstream"
using namespace std;

double f(double x) {
  double calc = 7./(pow(x,2)+1);
  return calc;
}

double trapezoid(double a, double b) {
  double calc = (b-a)/2. * (f(a)+f(b));
  return calc;
}

double trapezoid_err(double a,double b) {
  double calc = -1./12 *pow(b-a,3);
  return calc;
}

double simpson(double a, double b) {
  double calc = (b-a)/6. * (f(a)+4*f((b+a)*0.5)+f(b));
  return calc;
}

double simpson_err(double a, double b){
  double calc = -1./90 * pow((b-a)/2,5);
  return calc;
}

double simpson2(double a, double b) {
  double calc = (b-a)/90. * (7*f(a)+32*f(a+(b-a)/4)+12*f((b+a)/2)+32*f(a+3*(b-a)/4)+7*f(b));
  return calc;
}

double simpson2_err(double a, double b) {
  double calc = -3./80 * pow((b-a)/3.,5);
  return calc;
}

double gauss_leg_quad(double a, double b) {
  double n1 = (b-a)/2;
  double n2 = (b+a)/2;
  double x[5];
  x[0] = 0;
  x[1] = 1./3*sqrt(5-2*sqrt(10./7));
  x[2] = - 1./3*sqrt(5-2*sqrt(10./7));
  x[3] = 1./3*sqrt(5+2*sqrt(10./7));
  x[4] = - 1./3*sqrt(5+2*sqrt(10./7));
  double w[5];
  w[0] = 128./225;
  w[1] = (322.+13.*sqrt(70.))/900.;
  w[2] = w[1];
  w[3] = (322.-13.*sqrt(70.))/900.;
  w[4] = w[3];
  double calc = 0;
  for (size_t i = 0; i < 5; i++) {
    calc += w[i]*f(n1*x[i] + n2);
  }
  return n1*calc;

}

int main(int argc, char **argv)
{
  // if(argc != 2)
  //   {
  //     cout<<"Usage: ./ex <NUMBER>"<<endl;
  //     return -1;
  //   }
  //create a window
  double trap = trapezoid(0,5.);
  double a = 0;
  double b = 5;
  double trap_seg = 0;
  double sim_seg = 0;
  double sim2_seg = 0;
  size_t n = 1000;
  double h = abs((b-a)/n);
  for (size_t i = 0; i < n; i++) {
    double x0 = i*h;
    double x1 = h + i*h;
    trap_seg += trapezoid(x0,x1);
    sim_seg+= simpson(x0,x1);
    sim2_seg+= simpson2(x0,x1);
  }

  double trap_err = trapezoid_err(0,5);
  double sim = simpson(0,.5);
  double sim_err = simpson_err(0,5);
  double sim2 = simpson2(0,5.);
  double sim2_err = simpson2_err(0,5);
  std::cout << "1. Newton-Cotes formulas" << std::endl;
  std::cout << "trapezoid = "<< trap << " +- "<< trap_err << std::endl;
  std::cout << "trapezoid_seg = "<< trap_seg <<" for "<<n<<" steps" << std::endl;
  std::cout << "simpson = "<< sim<< " +- "<< sim_err << std::endl;
  std::cout << "sim_seg = "<< sim_seg <<" for "<<n<<" steps" << std::endl;
  std::cout << "simpson 3/8 = "<< sim2 << " +- "<< sim2_err << std::endl;
  std::cout << "sim2_seg = "<< sim2_seg <<" for "<<n<<" steps" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "to use gauss quadrature one needs to change the interval to [-1,1]" << std::endl;
  double gaus_leg = gauss_leg_quad(0,5.);
  std::cout << "gauss_leg_quad = "<< gaus_leg << std::endl;



  return 1;
}
