#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// Exact rank-one covariance update from sim_r27_rep50_nested.R. The
// sample covariance of rbind(ref, query) is ((n-1)/n) S + dd dd'/(n+1).
// [[Rcpp::export]]
NumericMatrix score_two_weights(NumericMatrix queries, NumericVector mu,
  NumericMatrix S, NumericMatrix Rref, NumericMatrix Wa, NumericMatrix Wu,
  int nref) {
  const int nr=queries.nrow(),p=queries.ncol();
  if (mu.size()!=p || S.nrow()!=p || S.ncol()!=p || Rref.nrow()!=p ||
      Rref.ncol()!=p || Wa.nrow()!=p || Wa.ncol()!=p || Wu.nrow()!=p || Wu.ncol()!=p)
    stop("Inconsistent score dimensions");
  const double c1=(nref-1.0)/nref,c2=1.0/(nref+1.0);
  NumericMatrix out(nr,2); std::vector<double> dd(p),ss(p);
  for(int i=0;i<nr;i++) {
    for(int j=0;j<p;j++) {
      dd[j]=queries(i,j)-mu[j];
      ss[j]=std::sqrt(c1*S(j,j)+c2*dd[j]*dd[j]);
      if(!std::isfinite(ss[j]) || ss[j]<=0) stop("Nonfinite or zero query variance");
    }
    double a=0.0,u=0.0;
    for(int j=0;j<p;j++) for(int k=0;k<j;k++) {
      const double r=(c1*S(j,k)+c2*dd[j]*dd[k])/(ss[j]*ss[k]);
      const double delta=Rref(j,k)-r, d=delta*delta;
      a += (Wa(j,k)+Wa(k,j))*d; u += (Wu(j,k)+Wu(k,j))*d;
    }
    out(i,0)=std::log(nref*double(nref)*a);
    out(i,1)=std::log(nref*double(nref)*u);
  }
  return out;
}
