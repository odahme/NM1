
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

  if(argc != 3)
  {
    cout<<"Usage: ./Ex1 <method(L,N)> <number of points>"<<endl;
    return -1;
  }

  const char* method = argv[1];
  unsigned int interpolationN = atoi(argv[2]);

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
  TVectorD xpoints(number);
  TVectorD ypoints(number);
  input.close();
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
  unsigned int N = interpolationN;
  TVectorD xinter(N+1);
  TVectorD yinter(N+1);
  double deltax = (xmax - xmin)/N;
  ofstream output;
  output.open("out.txt");
  unsigned int k = 0;
  double xint = xmin;
  if (strcmp("L", method) == 0){
    while (xint < xmax) {
      xint = xmin + deltax*k;
      double yint = 0;
      for (size_t i = 0; i < number; i++) {
        double val = 1;
        for (size_t j = 0; j < number; j++) {
          if (j==i) {
            /* code */
          }else{
            RooRealVar *xi = (RooRealVar*) pointList[i]->at(0);
            RooRealVar *xj = (RooRealVar*) pointList[j]->at(0);
            val *= (xint - xj->getVal())/(xi->getVal() - xj->getVal());
          }
        }
        RooRealVar *yi = (RooRealVar*) pointList[i]->at(1);
        yint += yi->getVal()*val;
      }
      output<< xint <<"\t"<< yint <<"\n";
      xinter[k] = xint;
      yinter[k] = yint;
      k++;
    }
  }

  if (strcmp("N", method) == 0){
    RooRealVar *var1 = (RooRealVar*) pointList[0]->at(0);
    RooRealVar *var2 = (RooRealVar*) pointList[1]->at(0);
    double h = abs(var1->getVal() - var2->getVal());
    for (size_t i = 1; i < number-1; i++) {
      RooRealVar *xi = (RooRealVar*) pointList[i]->at(0);
      RooRealVar *xj = (RooRealVar*) pointList[i+1]->at(0);
      double dif2 = abs(xi->getVal() - xj->getVal());
      if (h - dif2 > 0.01) {
        std::cout << "The grid is not uniform, plz use L for the interpolation" << std::endl;
        output.close();
        return 0;
      }
    }
    double delta[number];
    for (size_t n = 1; n < number; n++) {
      double deltacurr = 0;
      for (size_t s = 0; s <= n; s++) {
        RooRealVar *yval = (RooRealVar*) pointList[s]->at(1);
        deltacurr += evensign(s+n) * fak(n)/(fak(s)*fak(n-s)) * yval->getVal();
      }
      delta[n] = deltacurr/(fak(n)*pow(h,n));
    }
    while (xint < xmax) {
      xint = xmin + deltax*k;
      RooRealVar *yval = (RooRealVar*) pointList[0]->at(1);
      double yint = yval->getVal();
      for (size_t n = 1; n < number; n++) {
        double xdiff = 1;
        for (size_t s = 0; s < n; s++) {
          RooRealVar *xn = (RooRealVar*) pointList[s]->at(0);
          xdiff *= (xint - xn->getVal());
        }
        yint +=  (delta[n] * xdiff);
      }
      output<< xint <<"\t"<< yint <<"\n";
      xinter[k] = xint;
      yinter[k] = yint;
      k++;
    }
  }


  output.close();

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
