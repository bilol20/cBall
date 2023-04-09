#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double CppS(int n, NumericMatrix Kz, NumericVector Wz,
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
      NumericMatrix A = L(s(j)-1);
      L1(j) = rowcolsampler(A,s);
    }
    T[i] = CppS(n,Kz,Wz,Kyz,Wyz,L1);
  }
  return(T);
}

//[[Rcpp::export]]
double kernel(double x,double h){
  if(-h <x && x<h){
    return(0.75*(1-x*x/(h*h)));
  }else{
    return(0);
  }
}

//[[Rcpp::export]]
int fun2(NumericVector x)
{
  int s = 0;
  int n = x.length();
  for(int i = 0; i<n; i++){
    if(x(i)!=0)
      s = s + 1;
  }
  return(s);
}

//[[Rcpp::export]]
double eu(NumericVector x, NumericVector y){
  return(std::sqrt(sum((x-y)*(x-y))));
}

//[[Rcpp::export]]
NumericMatrix dist_cpp(NumericMatrix X){
  int n = X.nrow();
  NumericMatrix D(n,n);
  for(int i = 0; i<n; i++){
    for(int j = 0; j<n; j++){
      D(i,j) = eu(X(i,_),X(j,_));
    }
  }
  return(D);
}

