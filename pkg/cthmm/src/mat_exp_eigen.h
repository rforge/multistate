#ifndef MAT_EXP_EIGEN_H
#define MAT_EXP_EIGEN_H
#include <RcppArmadillo.h>
RcppExport SEXP mat_exp_eigen(SEXP eigen_decomp, SEXP time_int);
arma::mat mat_exp_eigen_cpp(Rcpp::List rate_eigen, double time_int);
#endif