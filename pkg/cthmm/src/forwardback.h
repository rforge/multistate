#ifndef FORWARDBACK_H
#define FORWARDBACK_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP forwardback(SEXP problist, SEXP emission, SEXP x, SEXP delta);
#endif 
