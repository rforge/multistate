#ifndef JOINT_DURATION_2MOMENT_H
#define JOINT_DURATION_2MOMENT_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_duration_2moment(int i, int j,Rcpp::List eigen_decomp, double interval_len);
#endif

#ifndef JOINT_DUR_TRANS_2MOMENT_H
#define JOINT_DUR_TRANS_2MOMENT_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma:: mat joint_dur_trans_2moment(int d1, int tr1, int tr2, Rcpp::List eigen_decomp, double interval_len);
#endif

#ifndef JOINT_TRANSITION_2MOMENT_H
#define JOINT_TRANSITION_2MOMENT_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_transition_2moment(int j1, int k1, int j2, int k2,Rcpp::List eigen_decomp, double interval_len);
#endif

#ifndef JOINT_DUR_DUR_EXACT_TIME_H
#define JOINT_DUR_DUR_EXACT_TIME_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_dur_dur_exact_time(int i,int j, double interval_len, int absorb_state,Rcpp::List eigen_decomp);
#endif

#ifndef JOINT_DUR_TRANS_EXACT_TIME_H
#define JOINT_DUR_TRANS_EXACT_TIME_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_dur_trans_exact_time(int d1, int tr1, int tr2, int absorb_state, Rcpp::List eigen_decomp, double interval_len);
#endif

#ifndef JOINT_TRANS_TRANS_EXACT_TIME_H
#define JOINT_TRANS_TRANS_EXACT_TIME_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
arma::mat joint_trans_trans_exact_time(int j1,int k1,int j2,int k2,int absorb_state,Rcpp::List eigen_decomp,double interval_len);
#endif

#ifndef JOINT_TRANS_TRANS_TIMES_H
#define JOINT_TRANS_TRANS_TIMES_H
#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP joint_trans_trans_times(SEXP time_diffs,SEXP j1_,SEXP k1_,SEXP j2_,SEXP k2_,SEXP absorb_state,SEXP eigen_decomp,SEXP exact_time_rank);
#endif

#ifndef JOINT_DUR_DUR_TIMES_H
#define JOINT_DUR_DUR_TIMES_H
#include <RcppArmadillo.h>
RcppExport SEXP joint_dur_dur_times(SEXP time_diffs,SEXP i_,SEXP j_,SEXP absorb_state,SEXP eigen_decomp,SEXP exact_time_rank);
#include <Rcpp.h>
#endif

#ifndef JOINT_TRANS_DUR_TIMES_H
#define JOINT_TRANS_DUR_TIMES_H
#include <RcppArmadillo.h>
RcppExport SEXP joint_trans_dur_times(SEXP time_diffs,SEXP d1_,SEXP tr1_,SEXP tr2_,SEXP absorb_state,SEXP eigen_decomp,SEXP exact_time_rank);
#include <Rcpp.h>
#endif