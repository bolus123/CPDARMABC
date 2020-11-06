#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <omp.h>


//using namespace Rcpp;

// [[Rcpp::export()]]
void omp2 (int t = 1) {
  omp_set_num_threads(t) ;
# pragma omp parallel for
  for (int i = 0 ; i < 10 ; i++) {
    Rcpp::Rcout << " " << i << " " ;
  }
  Rcpp::Rcout << std::endl ;
}
