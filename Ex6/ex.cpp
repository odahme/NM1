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

TVectorD* getx(TMatrixD *A, TVectorD *b) {
  size_t nrows = A->GetNrows();
  size_t ncols = A->GetNcols();
  TMatrixD L(nrows,ncols);
  TMatrixD U(nrows,ncols);
  TMatrixD D(nrows,ncols);
  TVectorD *sol = new TVectorD(nrows);
  double r = 0;

  for (size_t i = 0; i < nrows; i++) {
    for (size_t j = 0; j < ncols; j++) {
      if (i>j) {
        L[i][j] = (*A)[i][j];
      }
      if (i==j) {
        if ((*A)[i][j] == 0) {
          D[i][j] = 0;
        }else{
          D[i][j] = 1./(*A)[i][j];
        }
      }
      if (i<j) {
        U[i][j] = (*A)[i][j];
      }
    }
  }

  double epsilon = 1e-5;
  for (size_t x = 0; x < 100000; x++) {
    TVectorD xk(nrows);
    unsigned int delta = 0;
    for (size_t i = 0; i < nrows; i++) {
      double val = 0;
      for (size_t j = 0; j < ncols; j++) {
        val += -D[i][j]*(L[i][j] + U[i][j])* (*sol)[j] + D[i][j] * (*b)[j];
      }
      xk[i] = val;
    }
    for (size_t i = 0; i < nrows; i++) {
      if (abs( (*sol)[i] - xk[i]) < epsilon) {
        delta++;
      }
    }
    (*sol) = xk;
    if (delta == nrows) {
      break;
    }
  }
  return sol;
}

TVectorD* getxGS(TMatrixD *A, TVectorD *b) {
  size_t nrows = A->GetNrows();
  size_t ncols = A->GetNcols();
  TVectorD *sol = new TVectorD(nrows);

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
      std::cout << "/* message */" << std::endl;
    }
  }
  return sol;
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

  TMatrixD *A = new TMatrixD(10,10);
  TVectorD *b = new TVectorD(10);

  ifstream inputstream ("A_Matrix.txt",ios::in);
  size_t row = 0;
  if(inputstream.is_open()) {
    double val1;
    double val2;
    double val3;
    double val4;
    double val5;
    double val6;
    double val7;
    double val8;
    double val9;
    double val10;
    string buffer;
    getline(inputstream,buffer);
    while(inputstream >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7 >> val8 >> val9 >> val10 || row != 10 ) {
      (*A)[row][0] = val1;
      (*A)[row][1] = val2;
      (*A)[row][2] = val3;
      (*A)[row][3] = val4;
      (*A)[row][4] = val5;
      (*A)[row][5] = val6;
      (*A)[row][6] = val7;
      (*A)[row][7] = val8;
      (*A)[row][8] = val9;
      (*A)[row][9] = val10;
      row++;
    }
      inputstream.close();
    }
  else {
    cout << "Unable to open file" << endl;
  }
  ifstream bstream ("b_Vector.txt",ios::in);
  if(bstream.is_open()) {
    string buffer;
    double val1;
    size_t row = 0;
    getline(bstream,buffer);
    while( bstream >> val1 ) {
      (*b)[row] = val1;
      row++;
    }
      bstream.close();
    }
  else {
    cout << "Unable to open file" << endl;
  }

  A->Print();
  b->Print();
  //double detA = det(A);
  double detA = A->Determinant();
  std::cout << "detA = "<< detA << std::endl;
  //gauss(*A);
  //A->Print();
  //TMatrixD *M = add_b(A,b);
  //M->Print();
  TVectorD *x = getxGS(A,b);
  x->Print();






  return 1;
}
