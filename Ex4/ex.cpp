#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TApplication.h>
using namespace std;

double func1(double x) {
  double calc = x * exp(-(x*x));
  return calc;
}
double func2(double x) {
  double calc = cos(3*x) + 1./3 * sin(x);
  return calc;
}
double func3(double x) {
  double calc = func2(x)/(x - 0.585338);
  return calc;
}

int main(int argc, char **argv)
{
  // if(argc != 2)
  //   {
  //     cout<<"Usage: ./ex <NUMBER>"<<endl;
  //     return -1;
  //   }
  //create a window
  TApplication theApp("demoApplication",&argc,argv);
  // create a canvas
  TCanvas c1("c1","c1",1,1,1024,768);

  std::cout << "Numercial derivatives" << std::endl;
  std::cout << "a) foward difference quotient" << std::endl;
  double x = 1;
  double dxa[3];
  double dxb[3];
  double h[3];
  h[0] = 1e-1;
  h[1] = 1e-2;
  h[2] = 1e-3;
  for (size_t i = 0; i < 3; i++) {
    double dx = (func1(x+h[i]) - func1(x))/(h[i]);
    std::cout << "h = "<<h[i]<<" ; f'(x=1) = "<< dx << std::endl;
    dxa[i] = dx;
  }
  std::cout << "b) central difference quotient" << std::endl;
  for (size_t i = 0; i < 3; i++) {
    double dx = (func1(x+h[i]) - func1(x-h[i]))/(2*h[i]);
    std::cout << "h = "<<h[i]<<" ; f'(x=1) = "<< dx << std::endl;
    dxb[i] = dx;
  }
  std::cout << "c) comparison" << std::endl;
  for (size_t i = 0; i < 3; i++) {
    double diff = abs(dxa[i] - dxb[i]);
    std::cout << "for h"<<i<<" ; a-b = "<< diff << std::endl;
  }
  std::cout << "one can see that the with smaller h the precision of the diff quotient rises. Also the precision of b is higher than in a" << std::endl;

  std::cout << "-------------------" << std::endl;
  std::cout << "Newton-Raphson method in D=1" << std::endl;
  std::cout << "a) x_start = 1 ; t = 0.001" << std::endl;
  double diff = 1e32;
  size_t n = 0;
  double x0 = 1;
  while (diff > 0.001) {
    double h = 0.001;
    double dx = (func2(x0+h) - func2(x0-h))/(2*h);
    double x1 = x0 - func2(x0)/dx;
    diff = abs(func2(x1));
    x0 = x1;
    n++;
    if (n > 100000) {
      break;
    }
  }
  std::cout << "root of f(x) = "<< x0 << std::endl;
  std::cout << "b) x_start = 0.5 ; t = 1e-8" << std::endl;
  diff = 1e32;
  n = 0;
  x0 = 0.5;
  while (diff > 1e-8) {
    double h = 0.001;
    double dx = (func2(x0+h) - func2(x0-h))/(2*h);
    double x1 = x0 - func2(x0)/dx;
    diff = abs(func2(x1));
    x0 = x1;
    n++;
    if (n > 100000) {
      break;
    }
  }
  std::cout << "root of f(x) = "<< x0 << std::endl;

  std::cout << "c) x_start = 0.8 ; t = 1e-8" << std::endl;
  diff = 1e32;
  n = 0;
  x0 = 0.8;
  while (diff > 1e-8) {
    double h = 0.001;
    double dx = (func3(x0+h) - func3(x0-h))/(2*h);
    double x1 = x0 - func3(x0)/dx;
    diff = abs(func3(x1));
    x0 = x1;
    n++;
    if (n > 100000) {
      break;
    }
  }
  std::cout << "root of f(x) = "<< x0 << std::endl;

  std::cout << "d) x_start = 0.8 ; t = 1e-8" << std::endl;
  diff = 1e32;
  n = 0;
  x0 = 0.8;
  while (diff > 1e-8) {
    double h = 0.001;
    double dx = (func2(x0+h) - func2(x0-h))/(2*h);
    double x1 = x0 - func2(x0)/dx;
    diff = abs(func2(x1));
    x0 = x1;
    n++;
    if (n > 100000) {
      break;
    }
  }
  std::cout << "root of f(x) = "<< x0 << std::endl;
  std::cout << "It finds the same root as in b, since it is the nearest root." << std::endl;
  std::cout << "d) one could bound the the x values, so that the method would give out, 'no roots found' if it hits the boundary or take the boundary point as new x or take a random new starting point" << std::endl;


  return 1;
}
