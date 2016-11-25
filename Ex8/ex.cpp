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
  double calc = 1/(1+ sin(pow(x,2)));
  return calc;
}

// double f(double x) {
//   double calc = sqrt(1.-pow(x,2))/(1.1-x);
//   return calc;
// }

double trapezoid(double a, double b) {
  double calc = (b-a)/2. * (f(a)+f(b));
  return calc;
}

double trapezoid_err(double a,double b) {
  double calc = -1./12 *pow(b-a,3);
  return calc;
}

double simpson(double a, double b) {
  double calc = (b-a)/6. * (f(a)+4.*f((b-a)/2.)+f(b));
  return calc;
}

double simpson_err(double a, double b){
  double calc = -1./90 * pow((b-a)/2,5);
  return calc;
}

double simpson2(double a, double b) {
  double calc = (b-a)/90. * (7*f(a)+32*f((b-a)/4)+12*f((b-a)/2)+32*f(3*(b-a)/4)+7*f(b));
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
  //     cout<<"Usage: ./ex <Ex NUMBER>"<<endl;
  //     return -1;
  //   }
  // unsigned int exNumber = atoi(argv[1]);



  double under = 0;
  double upper = 1;
  size_t k = 10;
  double h = abs(upper - under)/k;

  double I1 = 0;
  double I3 = 0;
  for (size_t i = 0; i < k; i++) {
    double x1 = under + i*h;
    double x2 = under + h + i*h;
    I1 += trapezoid(x1,x2);
    I3 += simpson(x1,x2);
  }

  k = 20;
  h = abs(upper - under)/k;
  double I2 = 0;
  for (size_t i = 0; i < k; i++) {
    double x1 = under + i*h;
    double x2 = under + h + i*h;
    I2 += trapezoid(x1,x2);
  }

  double I = (4.*I2 - I1)/3.;
  std::cout << "I trap = "<< I1 << std::endl;
  std::cout << "I Richard = "<< I << std::endl;
  std::cout << "I simpson = "<< I3 << std::endl;


  // double under = 0;
  // double upper = 1;
  //
  // size_t k = 5;
  //
  // TMatrixD A(2*k,k);
  //
  //
  // for (size_t j = 0; j < 2*k; j++) {
  //   size_t subdiv = pow(2,j);
  //   double h = abs(upper - under)/subdiv;
  //   double I = 0;
  //   for (size_t i = 0; i < subdiv; i++) {
  //     double x1 = under + i*h;
  //     double x2 = under + h + i*h;
  //     I += trapezoid(x1,x2);
  //   }
  //   A[j][0] = I;
  // }
  //
  // for (size_t n = 1; n < k; n++) {
  //   for (size_t j = 0; j < k-n; j++) {
  //     A[j][n] = 1/(pow(4,n)-1) * (pow(4,n)*A[j+1][n-1] - A[j][n-1] );
  //   }
  // }
  //  A.Print();
  //std::cout << "trap = "<< A[0][k-1]<< std::endl;
  //std::cout << "I = "<< A[k-1][0] << std::endl;





  return 1;
}
