//requires auxmat_cpp.h
#ifndef JOINT_MEAN_MARKOV_REWARDS_CPP_H
#define JOINT_MEAN_MARKOV_REWARDS_CPP_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_mean_markov_rewards_cpp(int the_index, double interval_len, Rcpp::List rate_eigen);
#endif

#ifndef JOINT_MEAN_MARKOV_REWARDS_H
#define JOINT_MEAN_MARKOV_REWARDS_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP joint_mean_markov_rewards(SEXP i, SEXP interval_len, SEXP rate_eigen);
#endif
