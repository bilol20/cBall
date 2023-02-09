#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double CppS(int n,NumericMatrix Kz, NumericVector Wz, NumericMatrix Kyz, NumericVector Wyz, List L){
  double T = 0;
  for(int p = 0; p<n; p++){
    NumericMatrix s(n,n);
    for(int i = 0; i<n; i++){
      NumericMatrix A = L[i];
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
