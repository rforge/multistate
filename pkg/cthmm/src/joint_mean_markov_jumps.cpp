#include "joint_mean_markov_jumps.h"
#include "auxmat_cpp.h"
#include "uniformization_mean.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>


arma::mat joint_mean_markov_jumps_cpp(Rcpp::List rate_eigen, arma::mat regist_matrix, double 
									  interval_len){
	int space_size=regist_matrix.n_cols;
	arma::cx_mat factorial_moments = arma::zeros<arma::cx_mat>(space_size,space_size);
	arma::cx_mat rate_reg=Rcpp::as<arma::cx_mat>(rate_eigen["rate"])%regist_matrix;
	arma::cx_mat vectors=Rcpp::as<arma::cx_mat>(rate_eigen["vectors"]);
	arma::cx_mat invvectors=Rcpp::as<arma::cx_mat>(rate_eigen["invvectors"]);
	arma::cx_colvec values=Rcpp::as<arma::cx_colvec>(rate_eigen["values"]);
	
	arma::mat rate=Rcpp::as<arma::mat>(rate_eigen["rate"]);
	arma::colvec realvec=arma::sort(real(values));
	bool replicate=Rcpp::as<bool>(rate_eigen["replicate"]);
	//bool replicate=0;
	arma::mat out;
	//for(int i=0; i<realvec.n_elem-1;i++){
	//	if((realvec(i+1)-realvec(i))<1e-5||(realvec(i)-realvec(i+1)<1e-5)){
	//		replicate=1;
	//		exit;
	//	}
	//}
	if(replicate==1){
		arma::mat real_reg=real(rate_reg);
		out=uniformization_mean(interval_len, real_reg, rate);  
	}else{
		arma::cx_mat int_matrix=auxmat_cpp(values,interval_len);
		factorial_moments=vectors * (int_matrix % (invvectors * rate_reg* vectors)) *invvectors;
		out=arma::real(factorial_moments);
	}
	return(out);
}


RcppExport SEXP joint_mean_markov_jumps(SEXP rate_eigen, SEXP regist_matrix, SEXP 
										interval_len){
	double iinterval_len = Rcpp::as<double>(interval_len);
	arma::mat rregist_matrix=Rcpp::as<arma::mat>(regist_matrix);
	Rcpp::List rrate_eigen(rate_eigen);
	arma::mat out = joint_mean_markov_jumps_cpp(rrate_eigen, rregist_matrix, iinterval_len);

	return(Rcpp::wrap(out));
}

