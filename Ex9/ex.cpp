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

double f(double x, double y) {
  double calc = 8.*(1.-2*x)*y;
  return calc;
}

TMatrixD* get_tableaus(const char* in) {
  TMatrixD *A = new TMatrixD(5,5);

  ifstream inputstream (in,ios::in);
  size_t row = 0;
  if(inputstream.is_open()) {
    double val1;
    double val2;
    double val3;
    double val4;
    double val5;
    string buffer;
    getline(inputstream,buffer);
    while(inputstream >> val1 >> val2 >> val3 >> val4 >> val5  || row != 5 ) {
      (*A)[row][0] = val1;
      (*A)[row][1] = val2;
      (*A)[row][2] = val3;
      (*A)[row][3] = val4;
      (*A)[row][4] = val5;
      row++;
    }
      inputstream.close();
    }
  else {
    cout << "Unable to open file" << endl;
  }

  return A;
}

double get_next_y(TMatrixD *bu, TF1 *fa, double y0, double x0, double h) {
  unsigned int stage = bu->GetNrows() ;
  double k[stage];
  k[0] = fa->Eval(x0,y0);
  double xn = x0;
  double yn = y0;
  for (size_t i = 1; i < stage ; i++) {
    double ytemp = yn;
    for (size_t j = 1; j < i; j++) {
      ytemp += h* (*bu)[i-1][j] *k[j];
//      std::cout << "a"<<i-1<<j<<" = "<< (*bu)[i-1][j] << std::endl;
    }
    k[i] = fa->Eval(xn + h*(*bu)[i][0], ytemp);
//    std::cout << "c"<<i<<" = "<< (*bu)[0][i] << std::endl;
    yn += h*(*bu)[stage-1][i]*k[i];
//    std::cout << "b"<<i<<" = "<< (*bu)[stage-1][i] << std::endl;
  }
  return yn;
}

TGraphErrors* solve_diff(double under, double upper, double yn, size_t n = 0){
  TMatrixD *A = get_tableaus("butcherA.txt");
  TMatrixD *B = get_tableaus("butcherB.txt");
  double dy = 0;
  double x00 = under;
  double y00 = yn;

  string str;
  ifstream myfile ("function.txt",ios::in);
  if(myfile.is_open()) {
    getline(myfile,str);
    myfile.close();
    }
  else {
    cout << "Unable to open file" << endl;
  }
  const char *c = str.c_str();
  TF1 *fa = new TF1("fa",c,under,upper);


  if (n == 0) {
    std::vector<double> xAd;
    std::vector<double> yAd;
    std::vector<double> exAd;
    std::vector<double> eyAd;
    std::vector<double> hAd;
    double h = 0.0001;
    size_t i =0;
    double xn = under;
    double x0 = under;
    while (xn < upper) {
      yn = get_next_y(B,fa,yn,x0,h);
      dy = abs(yn - get_next_y(A,fa,yn,x0,h));
      xAd.push_back(x0 + h);
      x0 += h;
      hAd.push_back(h);
      yAd.push_back(yn);
      exAd.push_back(0);
      eyAd.push_back(dy);
      if (dy > pow(10,-4)) {
        h = 0.9*pow(pow(10,-4)/dy,1./3.)*h;
      }
      if (dy < pow(10,-5)) {
        h = 1.1*pow(pow(10,-5)/dy,1./3.)*h;
      }
      i++;
      xn = x0;
      std::cout << "i,xn = "<< i<<" , "<<xn << std::endl;
    }
    double x[i+1];
    x[0] = x00;
    double y[i+1];
    y[0] = y00;
    double ex[i+1];
    ex[0]=0;
    double ey[i+1];
    ey[0]=0;
    double h_plot[i+1];
    h_plot[0] = 0.0001;
    for (size_t j = 1; j < i+1; j++) {
      x[j] = xAd[j];
      y[j] = yAd[j];
      ex[j] = exAd[j];
      ey[j] = eyAd[j];
      h_plot[j] = hAd[j];
    }
    TGraphErrors *gr = new TGraphErrors(i,x,y,ex,ey);
    //TGraphErrors *gr = new TGraphErrors(i,x,h_plot,ex,ex);
    return gr;
  } else {
    TVectorD *y = new TVectorD(n+1);
    TVectorD *x = new TVectorD(n+1);
    TVectorD *ey = new TVectorD(n+1);
    TVectorD *ex = new TVectorD(n+1);
    (*x)[0] = under;
    (*y)[0] = yn;
    (*ex)[0] = 0;
    (*ey)[0] = 0;
    for (size_t i = 0; i < n; i++) {
      double h = abs(upper - under)/n;
      double x0 = under +  h*i;
      yn = get_next_y(B,fa,yn,x0,h);
      dy = abs(yn - get_next_y(A,fa,yn,x0,h));
      (*x)[i+1] = x0 + h;
      (*y)[i+1] = yn;
      std::cout << yn << std::endl;
      (*ex)[i+1] = 0;
      (*ey)[i+1] = dy;
    }
    TGraphErrors *gr = new TGraphErrors(*x,*y,*ex,*ey);

    for (size_t i = 0; i < 10; i++) {
      double y0 = exp(-2);
      double x0 = 0;
      double h = 0.01 + 0.01*i;
      yn = get_next_y(B,fa,yn,x0,h);
      dy = abs(yn - get_next_y(A,fa,yn,x0,h));
      std::cout << "h = "<< h << std::endl;
      std::cout << "dy = "<< dy << std::endl;
      std::cout << "------------------" << std::endl;
    }
    return gr;
  }
}

int main(int argc, char **argv)
{
  // if(argc != 2)
  //   {
  //     cout<<"Usage: ./ex <Ex NUMBER>"<<endl;
  //     return -1;
  //   }
  // unsigned int exNumber = atoi(argv[1]);

  TApplication theApp("demoApplication",&argc,argv);
  TCanvas c1("c1","c1",1,1,1024,768);
  TCanvas c2("c2","c2",1,1,1024,768);
  TGraphErrors *grfixed = solve_diff(0,2,exp(-2),100);
  TGraphErrors *gradapt = solve_diff(0,2,exp(-2),0);
  c1.cd();
  grfixed->Draw("ALP");
  c2.cd();
  gradapt->Draw("ALP");

  //turns off the program with mous clic
  theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
  //starts the canvas
  theApp.Run();

  return 1;
}
