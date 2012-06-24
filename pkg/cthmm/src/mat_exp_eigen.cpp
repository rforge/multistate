#include "mat_exp_eigen.h"
using namespace Rcpp;
using namespace RcppArmadillo;
arma::mat uniformization_prob2(double t, arma::mat rate);

RcppExport SEXP mat_exp_eigen(SEXP eigen_decomp, SEXP time_int){
  Rcpp::List rate_eigen(eigen_decomp);
  double t=Rcpp::as<double>(time_int);
  arma::mat out=mat_exp_eigen_cpp(rate_eigen,t); 
  return(Rcpp::wrap(out));
}

arma::mat mat_exp_eigen_cpp(Rcpp::List eigen_decomp, double time_int){
	arma::cx_mat vectors=Rcpp::as<arma::cx_mat>(eigen_decomp["vectors"]);
	arma::cx_mat invvectors=Rcpp::as<arma::cx_mat>(eigen_decomp["invvectors"]);
	arma::cx_colvec values=Rcpp::as<arma::cx_colvec>(eigen_decomp["values"]);
	arma::mat rate=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	
	arma::mat out;
	bool replicate=Rcpp::as<bool>(eigen_decomp["replicate"]);
	
	if(replicate==1){
		out=uniformization_prob2(time_int,rate);
		
	}else{
		out=real(vectors * diagmat(exp(values * time_int)) * invvectors);
	}
	return(out);
}
arma::mat uniformization_prob2(double t, arma::mat rate){
	
	arma::colvec the_diag=-1*rate.diag();
	double mu=the_diag.max();
	arma::mat A=arma::eye<arma::mat>(rate.n_cols,rate.n_cols)+rate/mu;
	
	int s=0;
	arma::mat out_sum=arma::zeros<arma::mat>(rate.n_cols,rate.n_cols);
	double error=(1-exp(-mu*t));
	arma::mat A_s=arma::eye<arma::mat>(rate.n_cols,rate.n_cols);
	
	while(error>1e-9){
		double mu_t=mu*t;
		out_sum=out_sum+exp(-mu*t)*A_s*pow(mu_t,s)/(Rcpp::internal::factorial(s)); 
		A_s=A_s*A;
		s++;
		error=error-exp(-mu*t)*pow(mu*t,s)/(Rcpp::internal::factorial(s));
	}
	return(out_sum);
}








