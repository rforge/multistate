//requires auxmat_cpp.h
#ifndef JOINT_MEAN_MARKOV_JUMPS_CPP_H
#define JOINT_MEAN_MARKOV_JUMPS_CPP_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_mean_markov_jumps_cpp(Rcpp::List rate_eigen, arma::mat regist_matrix, double 
									  interval_len);
#endif

#ifndef JOINT_MEAN_MARKOV_JUMPS_H
#define JOINT_MEAN_MARKOV_JUMPS_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP joint_mean_markov_jumps(SEXP rate_eigen, SEXP regist_matrix, SEXP 
										interval_len);
#endif
