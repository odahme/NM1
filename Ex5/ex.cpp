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

int sub_i(TMatrixD &A, TVectorD& b,size_t x){

  size_t nrows = A.GetNrows();
  for (size_t i = 0; i < nrows; i++) {
    A[i][x] = b[i];
  }
  return 1;
}

double det2(TMatrixD *A) {
  double calc = 0;
  calc = (*A)[0][0] * (*A)[1][1] - (*A)[0][1] * (*A)[1][0];
  return calc;
}


double det(TMatrixD *A) {
  double calc = 0;
  size_t nrows = A->GetNrows();
  size_t ncols = A->GetNcols();
  if (nrows != ncols) {
    std::cout << "Matrix is not quadratic" << std::endl;
    return -1;
  }
  if (nrows == 2) {
    //A->Print();
    calc = det2(A);
    //std::cout << "calc = "<< calc << std::endl;
    //cout << "Press Enter to Continue";
    //cin.ignore();
    return calc;
  } else {

    for (size_t j = 0; j < ncols; j++) {
      for (size_t i = 0; i < nrows; i++) {
        TMatrixD *B = new TMatrixD(nrows-1,ncols-1);
        size_t ia = 0;
        size_t ja = 0;
        for (size_t ib = 0; ib < ncols-1; ib++) {
          if (ib == i) {
            ia++;
          }
          ja = 0;
          for (size_t jb = 0; jb < nrows-1; jb++) {
            if (jb == j) {
              ja++;
            }
            (*B)[ib][jb] = (*A)[ia][ja];
            ja++;
          }
          ia++;
        }
        //std::cout << "A("<<i<<","<<j<<") = "<<(*A)[i][j] << std::endl;
        //B->Print();
        calc += pow(-1,i+j)*(*A)[i][j]*det(B);
        delete B;
      }
    }
  }
  return calc;
}

int swapRow(TMatrixD &M,int k, int l){
  size_t cols = M.GetNcols();
  double temp[cols];
  for (size_t i = 0; i < cols; i++) {
    temp[i] = M[k][i];
    M[k][i] = M[l][i];
  }
  for (size_t i = 0; i < cols; i++) {
    M[l][i] = temp[i];
  }
  return 1;
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

struct colum {
  int i;
  double val;
};

int sortCol(TMatrixD &M, int k){
  std::vector<colum> v;
  size_t rows = M.GetNrows();
  size_t cols = M.GetNcols();
  if (k <= rows-2) {
    return 1;
  }
  v.reserve(rows);
  for (size_t i = 0; i < rows; i++) {
    colum temp;
    temp.i = i;
    temp.val = M[i][k];
    v.push_back(temp);
  }
  std::sort(v.begin(),v.end(), [](colum a, colum b) {
    if (abs(a.val) > abs(b.val)) {
      return true;
    } else {
      return false;
    }
  });
  TMatrixD tempM(M);
  for (size_t i = k; i < rows; i++) {
    for (size_t j = k; j < cols; j++) {
      M[i][j] = tempM[v[i].i][j];
    }
  }
  return 1;
}

int gauss(TMatrixD &M) {
  size_t rows = M.GetNrows();
  size_t cols = M.GetNcols();
  for (size_t k = 0; k < rows-1; k++) {
    sortCol(M,k);
    double a11 = M[k][k];
    for (size_t i = k; i < rows-1; i++) {
      double a22 = M[i+1][k];
      for (size_t j = k; j < cols; j++) {
        M[k][j] *= a22;
        M[i+1][j] *= a11;
      }
      for (size_t j = k; j < cols; j++) {
        M[i+1][j] -= M[k][j];
        M[k][j] /= a22;
        M[i+1][j] /= a11;
      }
    }
  }
  double maxval = -1e32;
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (M[i][j] > maxval) {
        maxval = M[i][j];
      }
    }
  }
  maxval = abs(maxval);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (abs(M[i][j]) < maxval*1e-10) {
        M[i][j] = 0;
      }
    }
  }
  return 1;
}



TVectorD* getx(TMatrixD *A, TVectorD *b) {
  TMatrixD *M = add_b(A,b);
  gauss(*M);
  M->Print();
  size_t rows = M->GetNrows();
  size_t cols = M->GetNcols();
  TVectorD *x = new TVectorD(rows);
  (*x)[rows-1] = (*M)[rows-1][cols-1]/(*M)[rows-1][cols-2];
  for (int i = rows-2; i >= 0; i--) {
    double temp = 0;
    for (size_t k = i+1; k < rows; k++) {
      temp += (*M)[i][k] * (*x)[k];
    }
    (*x)[i] = 1/(*M)[i][i] * ((*M)[i][cols-1] - temp);
  }
  return x;
}

TMatrixD* getLk(TMatrixD *A, size_t k) {
  size_t rows = A->GetNrows();
  size_t cols = A->GetNcols();
  TMatrixD *L = new TMatrixD(rows,cols);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (i == j) {
        (*L)[i][j] = 1.;
      }
      if (k == j && j > i) {
        (*L)[i][j] = - (*A)[i][j]/(*A)[j][j];
      }
    }
  }
  return L;
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

  TMatrixD *A3 = new TMatrixD(*A);

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

  TMatrixD *L = new TMatrixD(*A);
  std::vector<TMatrixD*> MLA;
  TMatrixD *unity = new TMatrixD(A->GetNrows(),A->GetNcols());
  unity->UnitMatrix();
  MLA.reserve(A->GetNrows());
  for (size_t i = 0; i < A->GetNrows(); i++) {
    std::cout << "L"<<i<<" = " << std::endl;
    TMatrixD *Ltemp = getLk(L,i);
    Ltemp->Print();
    *L = (*Ltemp)*(*L);
    MLA.push_back(Ltemp);
  }
  std::cout << "A10 = " << std::endl;
  L->Print();

  TMatrixD *U = new TMatrixD(*L);
  L->UnitMatrix();
  for (size_t i = 0; i < MLA.size(); i++) {
    TMatrixD* tempM = MLA[i];
    *L *= (-1.* (*tempM) +2. * (*unity));
  }
  L->Print();

  return 1;
}
