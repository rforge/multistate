#ifndef AUXMAT_H
#define AUXMAT_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::cx_mat auxmat_cpp(arma::cx_colvec xx, double yy);
#endif

#ifndef AUXMAT2_H
#define AUXMAT2_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::cx_cube auxmat2_cpp(arma::cx_colvec v, double t);
#endif

#ifndef HOLBOLTH_FUN_H
#define HOLBOLTH_FUN_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::colvec hobolth_fun(Rcpp::List eigen_dec, double t, int a, int b, int c, int d, int e, int f);
#endif