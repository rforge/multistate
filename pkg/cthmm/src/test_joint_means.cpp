#include "auxmat_cpp.h"
#include "joint_mean_markov_rewards.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP test_joint_means(SEXP rate_eigen);
RcppExport SEXP test_joint_means(SEXP rate_eigen){
 Rcpp::List rrate_eigen(rate_eigen);
 arma::mat m1;
 m1=joint_mean_markov_rewards_cpp(1, .5, rrate_eigen);
 return(Rcpp::wrap(m1));
}
