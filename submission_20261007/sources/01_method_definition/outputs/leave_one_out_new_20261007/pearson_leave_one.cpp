#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// Exact augmented Pearson correlations. For each query, calculate every pair
// once and subtract the pairs containing the omitted marker. Re-normalise the
// retained fixed age weights, exactly as recomputing a reduced feature matrix.
// [[Rcpp::export]]
NumericMatrix pearson_all_omissions(NumericMatrix ref, NumericMatrix query,
                                   NumericMatrix weight) {
  int n=ref.nrow(),p=ref.ncol(),m=query.nrow();
  if(n<3 || p<3 || query.ncol()!=p || weight.nrow()!=p || weight.ncol()!=p)
    stop("Reference/query/weight dimensions are invalid");
  std::vector<double> mu(p,0),scatter(p*p,0),delta(p),variance(p),removed(p),w_removed(p,0);
  double total_weight=0;
  for(int j=0;j<p;++j) {
    for(int i=0;i<n;++i) {
      if(!R_finite(ref(i,j))) stop("Non-finite reference input");
      mu[j]+=ref(i,j);
    }
    mu[j]/=n;
    for(int k=0;k<p;++k) {
      if(!R_finite(weight(j,k)) || weight(j,k)<0) stop("Invalid weight");
      if(j==k && weight(j,k)!=0) stop("Diagonal weights must be zero");
      total_weight+=weight(j,k);
      if(j!=k) w_removed[j]+=weight(j,k)+weight(k,j);
    }
  }
  if(!(total_weight>0)) stop("Zero total weight");
  for(int j=0;j<p;++j) for(int k=j;k<p;++k) {
    double s=0;
    for(int i=0;i<n;++i) s+=(ref(i,j)-mu[j])*(ref(i,k)-mu[k]);
    scatter[j+p*k]=scatter[k+p*j]=s;
  }
  for(int j=0;j<p;++j) if(!(scatter[j+p*j]>0) || !(total_weight-w_removed[j]>0))
    stop("Constant reference feature or zero retained weight");
  NumericMatrix scores(m,p+1);
  const double fac=n/(n+1.0);
  for(int i=0;i<m;++i) {
    if(i%1000==0) checkUserInterrupt();
    for(int j=0;j<p;++j) {
      if(!R_finite(query(i,j))) stop("Non-finite query input");
      delta[j]=query(i,j)-mu[j];
      variance[j]=scatter[j+p*j]+fac*delta[j]*delta[j]; removed[j]=0;
    }
    double ds=0;
    for(int j=0;j<p;++j) for(int k=j+1;k<p;++k) {
      double rb=scatter[j+p*k]/std::sqrt(scatter[j+p*j]*scatter[k+p*k]);
      double ra=(scatter[j+p*k]+fac*delta[j]*delta[k])/std::sqrt(variance[j]*variance[k]);
      double d=ra-rb;
      double term=d*d*(weight(j,k)+weight(k,j));
      ds+=term;removed[j]+=term;removed[k]+=term;
    }
    if(!(ds>0)) stop("Non-positive DISCO perturbation; no pseudocount added");
    scores(i,0)=std::log(n*(double)n*ds/total_weight);
    for(int j=0;j<p;++j) {
      double d=ds-removed[j];
      if(!(d>0)) stop("Non-positive reduced DISCO perturbation");
      scores(i,j+1)=std::log(n*(double)n*d/(total_weight-w_removed[j]));
    }
  }
  return scores;
}
