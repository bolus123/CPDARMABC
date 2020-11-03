#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List optimCpp(const Rcpp::NumericVector& x, const Rcpp::NumericVector& order, const Rcpp::LogicalVector& include_mean) {
  Rcpp::Function optim_("arima");
  Rcpp::List out = optim_(Rcpp::Named("x", x),
                          Rcpp::Named("order", order),
                          Rcpp::Named("include.mean", include_mean));
  return(out);
}