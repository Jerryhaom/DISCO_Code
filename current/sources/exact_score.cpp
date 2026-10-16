#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector organ_score(NumericMatrix x, NumericVector mu, NumericVector sd,
                          NumericMatrix cc, NumericMatrix weight,
                          int nref, bool legacy=false, double coefficient=NA_REAL) {
  const int n=x.nrow(), p=x.ncol();
  const double nr=nref;
  const double a=Rcpp::NumericVector::is_na(coefficient) ? (legacy ? nr*nr/((nr+1)*(nr+1)*(nr-1)) : nr/((nr+1)*(nr-1))) : coefficient;
  NumericVector result(n);
  std::vector<double> z(p), inv(p);
  for (int i=0;i<n;i++) {
    for (int j=0;j<p;j++) { z[j]=(x(i,j)-mu[j])/sd[j]; inv[j]=1/std::sqrt(1+a*z[j]*z[j]); }
    double ds=0;
    for (int j=0;j<p;j++) for (int k=0;k<j;k++) {
      const double diff=cc(j,k)-(cc(j,k)+a*z[j]*z[k])*inv[j]*inv[k];
      ds+=diff*diff*(weight(j,k)+weight(k,j));
    }
    result[i]=std::log(ds*nr*nr);
  }
  return result;
}
