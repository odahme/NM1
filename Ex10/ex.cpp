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

double exactODE(double x) {
  double calc = 0;
  double e = 0.1;
  double w0 = 5;
  double c1 = 1.3;
  double c2 = -0.0748969;
  calc = c2*sin(sqrt(1-pow(e,2))*w0*x);
  calc += c1*cos(sqrt(1-pow(e,2))*w0*x);
  calc *= exp(-e*w0*x);
  return calc;
}

TMatrixD* add_b(TMatrixD *A, TVectorD *b){
  size_t rows = A->GetNrows();
  size_t cols = A->GetNcols();
  TMatrixD *M = new TMatrixD(rows,cols+1);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      (*M)[i][j] = (*A)[i][j];
    }
  }
  for (size_t i = 0; i < rows; i++) {
    (*M)[i][cols] = (*b)[i];
  }
  return M;
}

TVectorD* gaussmod(TMatrixD* A, TVectorD *b)
{
  TMatrixD *M = add_b(A,b);
  int n = M->GetNrows();
  for (int i=0;i<n;i++)
    {

      double max = abs((*M)[i][i]);
      int maxrow=i;
      for (int j=i+1;j<n;j++)
	{
	  if (abs((*M)[j][i])>max)
	    {
	      max = abs((*M)[j][i]);
	      maxrow = j;
	    }
	}


      for (int k=i;k<n+1;k++)
	{
	  double tmp = (*M)[maxrow][k];
	  (*M)[maxrow][k] = (*M)[i][k];
	  (*M)[i][k] = tmp;
	}

      for (int k=i+1; k<n;k++)
	{
	  double c = -(*M)[k][i]/(*M)[i][i];
	  for (int j=1;j<n+1;j++)
	    {
	      if(i==j)
		(*M)[k][j] = 0;
	      else (*M)[k][j] += c*(*M)[i][j];
	    }
	}
  if (i%100 == 0) {
    std::cout << "row "<<i<<" of "<<n<< " done" << std::endl;
  }
    }
  //print(A);
  TVectorD *m = new TVectorD(n);
  for (int i=n-1;i>=0;i--)
    {
      (*m)[i] = (*M)[i][n]/(*M)[i][i];
      for (int k=n-1;k>=0;k--)
	{
	  (*M)[k][n] -= (*M)[k][i]*(*m)[i];
	}
    }
  return m;
}

TVectorD* getxGS(TMatrixD *A, TVectorD *b,double under,double upper, double seed = 42) {

  if (A->Determinant() == 0) {
    std::cout << "can't solve this" << std::endl;
    TVectorD *sol = new TVectorD(1);
    return sol;
  }
  TRandom3 *rnd = new TRandom3(seed);
  size_t nrows = A->GetNrows();
  size_t ncols = A->GetNcols();
  TVectorD *sol = new TVectorD(nrows);
  for (size_t i = 0; i < nrows; i++) {
    (*sol)[i] = rnd->Uniform(under,upper);
  }

  double epsilon = 1e-5;
  for (size_t x = 0; x < 1000000; x++) {
    unsigned int delta = 0;
    TVectorD xk(nrows);
    xk = (*sol);
    for (size_t i = 0; i < nrows; i++) {
      double temp1 = 0;
      double temp2 = 0;

      for (size_t j = 0; j < ncols; j++) {
        if (j < i) {
          temp1 += (*A)[i][j]* (*sol)[j];
        }
        if (j>i) {
          temp2 += (*A)[i][j] * xk[j];
        }
      }
      (*sol)[i] = 1./(*A)[i][i] * ( (*b)[i] - temp1 - temp2);
    }
    for (size_t i = 0; i < nrows; i++) {
      if (abs( (*sol)[i] - xk[i]) < epsilon) {
        delta++;
      }
    }
    if (delta == nrows) {
      break;
    }
  }
  return sol;
}

int main(int argc, char **argv)
{
  if(argc != 3)
    {
      cout<<"Usage: ./ex <method> <N>"<<endl;
      return -1;
    }
  const char* method = argv[1];
  unsigned int N = atoi(argv[2]);

 TApplication theApp("demoApplication",&argc,argv);
 TCanvas c1("c1","c1",1,1,1024,768);

  TMatrixD *A = new TMatrixD(N,N);
  TVectorD *x = new TVectorD(N);
  TVectorD *ex = new TVectorD(N);
  TVectorD *b = new TVectorD(N);
  double upper = 2;
  double under = 0;
  double y0 = 1.3;
  double yn = -0.4;
  double h = (double) abs(upper-under)/N;
  string filename;
  //define lin equ system
  if (strcmp("f", method) == 0) {
    filename = "forward.png";
    for (size_t i = 1; i < N-1; i++) {
      (*A)[i][i-1] = pow(h,-2) - pow(h,-1) + 25;
      (*A)[i][i] = -2*pow(h,-2) + pow(h,-1);
      (*A)[i][i+1] = pow(h,-2);
    }
    //define starting conditions
    (*b)[0] = -y0*(pow(h,-2) - pow(h,-1) + 25);
    (*b)[N-1] = -yn*pow(h,-2);
    (*A)[0][0] = -2*pow(h,-2) + pow(h,-1);
    (*A)[0][1] = pow(h,-2);
    (*A)[N-1][N-1] = -2*pow(h,-2) + pow(h,-1);
    (*A)[N-1][N-2] = pow(h,-2) - pow(h,-1) + 25;
  }
  if (strcmp("c", method) == 0) {
    filename = "central.png";
    for (size_t i = 2; i < N-2; i++) {
      (*A)[i][i-2] = 1./(4.*pow(h,2));
      (*A)[i][i-1] = - 1./(2*h);
      (*A)[i][i] = 25. - 1./(2*pow(h,2));
      (*A)[i][i+1] = 1./(2*h);
      (*A)[i][i+2] = 1./(4.*pow(h,2));
    }
    (*b)[0] = - y0*(1./(4*pow(h,2)));
    (*A)[0][0] = - 1./(2*h);
    (*A)[0][1] = 25. - 1./(2*pow(h,2));
    (*A)[0][2] = 1./(2*h);
    (*A)[0][3] = 1./(4*pow(h,2));

    (*A)[1][1] = pow(h,-2) - pow(h,-1) + 25;
    (*A)[1][2] = -2*pow(h,-2) + pow(h,-1);
    (*A)[1][3] = pow(h,-2);

    (*b)[N-1] = - yn*(1./(4*pow(h,2)));
    (*A)[N-1][N-1] = 1./(2*h);
    (*A)[N-1][N-2] = 25. - 1./(2*pow(h,2));
    (*A)[N-1][N-3] = - 1./(2*h);
    (*A)[N-1][N-4] = 1./(4*pow(h,2));

    (*A)[N-2][N-2] = pow(h,-2);
    (*A)[N-2][N-3] = -2*pow(h,-2) + pow(h,-1);
    (*A)[N-2][N-4] = pow(h,-2) - pow(h,-1) + 25;
  }


  for (size_t i = 0; i < N; i++) {
    (*x)[i] = under +i*h;
  }




  //TVectorD *y = getxGS(A,b,under,upper,42);
  // A->Print();
  // b->Print();
  TVectorD *y = gaussmod(A,b);
  TVectorD *ey = new TVectorD(N);
  double maxerr = -1e32;
  for (size_t i = 0; i < N; i++) {
    double calc = abs((*y)[i] - exactODE((*x)[i]));
    (*ey)[i] = calc;
    if (calc > maxerr) {
      maxerr = calc;
    }
  }
  // TVectorD *y = getxGS(A,b,under,upper,42);
  // y->Print();

  // TGraph *gr = new TGraph(*x,*y);
  // c1.cd();
  // gr->Draw();

  TGraphErrors *grEr = new TGraphErrors(*x,*y,*ex,*ey);
  c1.cd();
  grEr->Draw();
  c1.SaveAs(filename.c_str());
  std::cout << "maxerr = "<< maxerr << std::endl;

  //turns off the program with mous clic
 theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
  //starts the canvas
 theApp.Run();

  return 1;
}
