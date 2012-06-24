#ifndef UNIFORMIZATION_MEAN
#define UNIFORMIZATION_MEAN_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat uniformization_mean(double t, arma::mat B, arma::mat rate);
#endif