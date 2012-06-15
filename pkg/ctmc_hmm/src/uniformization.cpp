#include "uniformization_mean.h"

arma::mat uniformization_mean(double t, arma::mat B, arma::mat rate){

 arma::colvec the_diag=-1*rate.diag();
 double mu=the_diag.max();
 arma::mat A=arma::eye<arma::mat>(B.n_cols,B.n_cols)+rate/mu;

 int s=0;
 arma::mat out_sum=arma::zeros<arma::mat>(B.n_cols,B.n_cols);
 double error=mu*t*(1-exp(-mu*t));
 arma::mat A_s=arma::eye(B.n_cols,B.n_cols);
 arma::mat B_star=B/mu;
 arma::mat D_s=B_star;
 
 while(error>1e-7){
  double mu_t=mu*t;
  out_sum=out_sum+exp(-mu*t)*D_s*pow(mu_t,s+1)/(Rcpp::internal::factorial(s+1)); 
  A_s=A_s*A;
  D_s=D_s*A+A_s*B_star;
  s++;
  error=error-(mu*t)*exp(-mu*t)*pow(mu*t,s)/(Rcpp::internal::factorial(s));
}
 return(out_sum);
}

RcppExport SEXP uniformization_mean_R(SEXP time_,SEXP B_,SEXP rate_){
	double t=Rcpp::as<double>(time_);
	arma::mat B=Rcpp::as<arma::mat>(B_);
	arma::mat rate=Rcpp::as<arma::mat>(rate_);
	arma::mat out=uniformization_mean(t,B,rate);
	return(Rcpp::wrap(out));
}