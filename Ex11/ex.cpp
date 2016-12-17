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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
using namespace std;

double f_bisec(double x){
  double calc = pow(x,3)+exp(-x) -x;
  return calc;
}


double f1_bisec(double x, double h) {
  double calc = (f_bisec(x+h) - f_bisec(x-h))/(2*h);
  return calc;
}

int sign(double x){
  bool s;
  if (x >= 0) {
    s = 1;
  } else {
    s = -1;
  }
  return s;
}

double f_gold(double x){
  double calc = pow(x,4) + 2*pow(x,2) + 4*x + 1;
  return calc;
}

double f_ber(double x){
  double calc = -2*x + tan(x);
  return calc;
}

double bisection(double under, double upper, double prec) {
  double a;
  double b;
  if (under < upper) {
    a = under;
    b = upper;
  }else{
    a = upper;
    b = under;
  }
  double c = 0.5*(a+b);
  size_t n = 0;
  while (abs(b-a) > prec) {
    if (f_bisec(c) < f_bisec(a)) {
      a = c;
    } else {
      b = c;
    }
    c = 0.5*(a+b);
    n++;
  }
  std::cout << "minimum at x = "<<c<< " after "<< n<<" steps with "<<prec<<" precision" << std::endl;
  return c;
}

double goldenRule(double under, double upper, double prec) {
  TRandom3 *rnd = new TRandom3(under+upper);
  bool start = true;
  double a = rnd->Uniform(under,upper);
  double b = rnd->Uniform(under,upper);
  double c = rnd->Uniform(under,upper);
  while (start) {
    a = rnd->Uniform(under,upper);
    c = rnd->Uniform(a,upper);
    b = rnd->Uniform(a,c);
    if ( f_gold(a) > f_gold(b) && f_gold(c) > f_gold(b)) {
      start = false;
    }
  }
  double d = rnd->Uniform(a,c);
  double w = (3-sqrt(5.))/2.;
  size_t n = 0;
  while (abs(a-c) > prec) {
    if (f_gold(d) < f_gold(b)) {
      if (d < b) {
        c = b;
        b = d;
      } else {
        a = b;
        b = d;
      }
    } else {
      if (d < b) {
        a = d;
      } else {
        c = d;
      }
    }
    if (abs(b-a) > abs(c-b)) {
      d = a + w*abs(b-a);
    } else {
      d = b + w*abs(c-b);
    }
    n++;
  }

  double calc = (a+c)*0.5;
  std::cout << "minimum at x = "<<calc<< " after "<< n<<" steps with "<<prec<<" precision" << std::endl;
  return calc;
}


double brenet(double under, double upper, double prec) {

  TRandom3 *rnd = new TRandom3(under+upper);
  bool start = true;
  double a = rnd->Uniform(under,upper);
  double b = rnd->Uniform(under,upper);
  double c = rnd->Uniform(under,upper);
  while (start) {
    a = rnd->Uniform(under,upper);
    c = rnd->Uniform(a,upper);
    b = rnd->Uniform(a,c);
    if ( f_ber(a) > f_ber(b) && f_ber(c) > f_ber(b)) {
      start = false;
    }
  }
  double d = rnd->Uniform(a,c);
  size_t n = 0;
  while (abs(a-c) > prec) {
    // std::cout << "a = "<< a << std::endl;
    // std::cout << "c = "<< c << std::endl;
    // std::cout << "b = "<< b << std::endl;
    // std::cout << "d = "<< d << std::endl;
    // std::cout << "f_ber(d) = "<< f_ber(d) << std::endl;
    if (f_ber(d) < f_ber(b)) {
      if (d < b) {
        c = b;
        b = d;
      } else {
        a = b;
        b = d;
      }
    } else {
      if (d < b) {
        a = d;
      } else {
        c = d;
      }
    }
    d = 0.5*((pow(a,2)*(f_ber(c)-f_ber(b)) + pow(b,2)*(f_ber(a)+f_ber(c)) + pow(c,2)*(f_ber(b)-f_ber(a)))/(a*(f_ber(c)-f_ber(b)) + b*(f_ber(a)+f_ber(c)) + c*(f_ber(b)-f_ber(a))));

    if (a < d && d < c) {
      b = d;
    } else {
      d = 0.5*(a+c);
    }
    n++;
  }
  double calc = (a+c)*0.5;
  std::cout << "minimum at x = "<<calc<< " after "<< n<<" steps with "<<prec<<" precision" << std::endl;
  return calc;
}

int main(int argc, char **argv)
{
  if (argc != 2)
    {
      cout<<"Usage: ./ex <method> "<<endl;
      return -1;
    }
  const char* method = argv[1];

 // TApplication theApp("demoApplication",&argc,argv);
 // TCanvas c1("c1","c1",1,1,1024,768);

  double epsilon = 1e-4;
  if (strcmp("b", method) == 0) {
    double calc = brenet(0,1,epsilon);
  }
  if (strcmp("g", method) == 0) {
    double calc = goldenRule(-1,0,epsilon);
  }
  if (strcmp("bi", method) == 0) {
    double calc = bisection(0,1,epsilon);
  }


  //turns off the program with mous clic
 // theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
 //  //starts the canvas
 // theApp.Run();

  return 1;
}
