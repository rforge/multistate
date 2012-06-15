#include "auxmat_cpp.h"
#include "joint_mean_markov_rewards.h"
#include "uniformization_mean.h"


RcppExport SEXP joint_mean_markov_rewards(SEXP i, SEXP interval_len, SEXP rate_eigen){
 double iinterval_len = Rcpp::as<double>(interval_len);
 int the_index=Rcpp::as<int>(i);
 Rcpp::List rrate_eigen(rate_eigen);
 arma::mat mean_rewards=joint_mean_markov_rewards_cpp(the_index-1, iinterval_len, rrate_eigen);
 return(Rcpp::wrap(mean_rewards));
 }


arma::mat joint_mean_markov_rewards_cpp(int di, double interval_len, Rcpp::List rate_eigen){
	arma::cx_colvec values=Rcpp::as<arma::cx_colvec>(rate_eigen["values"]);
	arma::cx_mat vectors=Rcpp::as<arma::cx_mat>(rate_eigen["vectors"]);
	arma::cx_mat invvectors=Rcpp::as<arma::cx_mat>(rate_eigen["invvectors"]);
	arma::cx_mat int_matrix=auxmat_cpp(values,interval_len);
	arma::colvec reward_vector=arma::zeros(values.n_elem);
	reward_vector(di)=1;
	arma::mat reward_matrix=arma::diagmat(reward_vector);
	
	
	arma::mat rate=Rcpp::as<arma::mat>(rate_eigen["rate"]);
	arma::colvec realvec=arma::sort(real(values));
	arma::mat out;
	bool replicate=Rcpp::as<bool>(rate_eigen["replicate"]);
	//bool replicate=0;

	//for(int i=0; i<realvec.n_elem-1;i++){
//		if((realvec(i+1)-realvec(i))<1e-5||(realvec(i)-realvec(i+1)<1e-5)){
//			replicate=1;
//			exit;
//		}
	//}
	if(replicate==1){
		out=uniformization_mean(interval_len, reward_matrix, rate);  
	}else{
	arma::cx_mat mean_rewards=vectors * (int_matrix % (invvectors *reward_matrix*vectors))*invvectors;
	out=arma::real(mean_rewards);
	}
	return(out);
}





