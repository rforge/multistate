#ifndef GET_MEAN_H
#define GET_MEAN_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP get_mean(SEXP single_joint_mean, SEXP obs_data, SEXP forward, SEXP backward, SEXP LL, SEXP emission_matrix);
#endif

#ifndef GET_MEAN_CPP_H
#define GET_MEAN_CPP_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
double get_mean_cpp(std::vector<arma::mat> joint_mean, arma::icolvec obs_data, arma::mat forward, arma::mat backward, double LL, arma::mat emission_matrix);
#endif