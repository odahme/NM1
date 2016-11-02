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
  A->Print();
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
    std::cout << "A11 = "<< (*A)[1][1] << std::endl;
    while(inputstream >> val1 >> val2 >> val3 >> val4 >> val5 >> val6 >> val7 >> val8 >> val9 >> val10 || row != 10 ) {
      std::cout << "val1 = "<< val1 << std::endl;
      std::cout << "row = "<< row << std::endl;
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
  TVectorD *x = getx(A,b);
  x->Print();






  return 1;
}
