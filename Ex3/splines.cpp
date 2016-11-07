
#include "TCanvas.h"
#include "TF1.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <TH3F.h>
#include <time.h>
#include <math.h>
#include "fstream"
#include "TVectorD.h"

using namespace std;

#include <stdio.h>

double fak(int n)
{
  int i;
  double fak = 1;

  for ( i=2; i<=n; i++ )
    fak *= i;

  return fak;
}

int evensign(int n ) {
  if (n % 2 == 0) {
    return 1;
  }else{
    return -1;
  }
}

int main(int argc,  char **argv) {

  if(argc != 2)
  {
    cout<<"Usage: ./splines.x <number of points>"<<endl;
    return -1;
  }

  unsigned int interpolationN = atoi(argv[1]);

  //create a window
  TApplication theApp("demoApplication",&argc,argv);

  TCanvas *c1 = new TCanvas("c1","Example1",1920,1080);
  // TCanvas *c2 = new TCanvas("c2","Example2",1920,1080);
  // TCanvas *c3 = new TCanvas("c3","Example3",1920,1080);
  // TCanvas *c4 = new TCanvas("c4","Example4",1920,1080);
  // TCanvas *c5 = new TCanvas("c5","Example5",1920,1080);
  // TCanvas *c6 = new TCanvas("c6","Example6",1920,1080);

  TMultiGraph *graph = new TMultiGraph("plot","Plot of interpolation");

  ifstream input("input.txt");
  double val1;
  double val2;
  int dim = 2;
  std::vector<RooArgList*> pointList;

  while (input >> val1 >> val2) {
    RooRealVar *roovar1 = new RooRealVar("x","xvalue",val1);
    RooRealVar *roovar2 = new RooRealVar("y","yvalue",val2);
    RooArgList *roolist = new RooArgList(*roovar1,*roovar2);
    pointList.push_back(roolist);
  }
  unsigned int number = pointList.size();
  input.close();
  TVectorD xpoints(number);
  TVectorD ypoints(number);

  double xmax = -1e32;
  double xmin = 1e32;
  for (size_t i = 0; i < number; i++) {
      RooArgList *point = pointList[i];
      RooRealVar *xrooval = (RooRealVar*) point->at(0);
      RooRealVar *yrooval = (RooRealVar*) point->at(1);
      xpoints[i] = xrooval->getVal();
      ypoints[i] = yrooval->getVal();
      double xval = xrooval->getVal();
      if (xmax < xval) {
        xmax = xval;
      }
      if (xmin > xval) {
        xmin = xval;
      }
  }

  //first calculate m
  double d[number-1];
  double h[number-1];
  double m[number];
  double b[number];
  unsigned int N = interpolationN;
  m[0] = 0;
  for (size_t k = 0; k < number-1; k++) {
    m[k+1] = 0;
    h[k] = xpoints[k+1] - xpoints[k];
    d[k] = (ypoints[k+1] - ypoints[k])/h[k];
  }
  //create Matrix A in A x = b
  double A[number][number];
  for (size_t i = 0; i < number; i++) {
    for (size_t j = 0; j < number; j++) {
      A[i][j] = 0;
    }
  }

  for (size_t i = 1; i < number-1; i++) {
    for (size_t j = 1; j < number-1; j++) {
      if (i == j) {
        b[i] = d[j] - d[j-1];
        A[i][j-1] = 1./6 * h[j-1];
        A[i][j] = 1./3 *(h[j]+h[j-1]);
        A[i][j+1] = 1./6 * h[j];
      }
    }
  }

  for (size_t i = 0; i < number-1; i++) {
    for (size_t j = 0; j < number-1; j++) {
      std::cout << A[i][j]<< " , " ;
    }
    std::cout << " "<< std::endl;
  }

  double D[number][number];
  double L[number][number];
  double U[number][number];
  for (size_t i = 0; i < number; i++) {
    for (size_t j = 0; j < number; j++) {
      D[i][j] = 0;
      L[i][j] = 0;
      U[i][j] = 0;
      if (i>j) {
        L[i][j] = A[i][j];
      }
      if (i==j) {
        if (A[i][j] == 0) {
          D[i][j] = 0;
        }else{
          D[i][j] = 1./A[i][j];
        }
      }
      if (i<j) {
        U[i][j] = A[i][j];
      }
    }
  }

  for (size_t x = 0; x < 100000; x++) {
    double mk[number];
    unsigned int delta = 0;
    for (size_t i = 0; i < number; i++) {
      double val = 0;
      for (size_t j = 0; j < number; j++) {
        val += -D[i][j]*(L[i][j] + U[i][j])*m[j] + D[i][j]*b[j];
      }
      mk[i] = val;
    }
    for (size_t i = 0; i < number; i++) {
      if (abs(m[i] - mk[i]) < 0.001) {
        delta++;
      }
      m[i] = mk[i];
    }
    if (delta == number) {
      break;
    }
  }

  // for (size_t i = 0; i < number; i++) {
  //   std::cout << "m"<<i<<" = "<<m[i] << std::endl;
  // }

  //interpolation
  TVectorD xinter(interpolationN*(number-1));
  TVectorD yinter(interpolationN*(number-1));
  size_t n = 0;
  for (size_t k = 0; k < number-1; k++) {
    double xk = xpoints[k];
    double delta = abs((xpoints[k] - xpoints[k+1])/interpolationN);
    double a[4];
    a[0] = ypoints[k];
    a[1] = d[k] - h[k]/6. * (2*m[k] + m[k+1]);
    a[2] = m[k]/2.;
    a[3] = (m[k+1] - m[k])/(6.*h[k]);
    for (size_t i = 0; i < interpolationN; i++) {
      double x = xk + delta*i;
      double s = a[0] + a[1]*(x-xk) + a[2]*pow(x-xk,2) + a[3]*pow(x-xk,3);
      xinter[n+i] = x;
      yinter[n+i] = s;
    }
    n += interpolationN;
  }


  TGraph *gr1 = new TGraph(xpoints,ypoints);
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(2.5);
  graph->Add(gr1);
  TGraph *gr2 = new TGraph(xinter,yinter);
  graph->Add(gr2);
  graph->Draw("ap");
  c1->cd();
  graph->GetXaxis()->SetTitle("x");
  graph->GetYaxis()->SetTitle("y");
  graph->Draw("ap");


  //turns off the program with mous clic
  theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
  //starts the canvas
  theApp.Run();

  return 0;
}
