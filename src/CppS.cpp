#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double CppS(int n,NumericMatrix Kz, NumericVector Wz,
            NumericMatrix Kyz, NumericVector Wyz, List L){
  double T = 0;
  for(int p = 0; p<n; p++){
    NumericMatrix s(n,n);
    for(int i = 0; i<n; i++){
      NumericMatrix A = L(i);
      for(int j = 0; j<n; j++){
        double s1 = Kz(p,i)*Kz(p,j)/(Wz(p)*Wz(p))+Kyz(p,i)*Kyz(p,j)/(Wyz(p)*Wyz(p));
        double s2 = 0;
        for(int k = 0; k<n; k++){
          s2 = s2 + (Kz(p,k)/Wz(p)-Kyz(p,k)/Wyz(p))*A(k,j);
        }
        s(i,j) = s2*s2*s1;
      }
    }
    T = T + sum(s);
  }
  return(T/n);
}

//[[Rcpp::export]]
NumericMatrix  rowcolsampler(NumericMatrix A, NumericVector s){
  int n = A.ncol();
  NumericMatrix B(n,n);
  for(int k=0;k<n;k++){
    for(int l=0; l<n;l++){
      B(k,l) = A(s(k)-1,s(l)-1);
    }
  }
  return(B);
}

//[[Rcpp::export]]
NumericVector resample(int n,NumericMatrix Kz, NumericVector Wz,
                       NumericMatrix Kyz, NumericVector Wyz, List L,
                       NumericMatrix Pi){
  int R = Pi.ncol();
  NumericVector T(R);
  for(int i = 0; i<R; i++){
    NumericVector s = Pi(_,i);
    List L1(n);
    for(int j = 0; j<n;j++){
      NumericMatrix A = L[s(j)-1];
      L1[j] = rowcolsampler(A,s);
    }
    T[i] = CppS(n,Kz,Wz,Kyz,Wyz,L1);
  }
  return(T);
}
