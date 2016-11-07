
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

  if(argc != 1)
  {
    cout<<"Usage: ./extraplation <>"<<endl;
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

  ofstream exstream;
  exstream.open("extrapolationInput.txt");
  pval = 1;
  for (size_t i = 0; i < 4; i++) {
    pval *= y[i]* pow((hk - x[i])/(x[i] - x[i-1]),i);
  }



  ifstream input("extrapolationInput.txt");
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

  nsigned int N = interpolationN;
  TVectorD xinter(N+1);
  TVectorD yinter(N+1);
  double deltax = (xmax - xmin)/N;
  ofstream output;
  output.open("out.txt");
  unsigned int k = 0;
  double xint = xmin;
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
