#include <RcppArmadillo.h>
#include <Rcpp.h>

RcppExport SEXP E_step_emission_init_all(SEXP likelihood_forward_backward_list,SEXP obs_data_list,SEXP num_states_,SEXP num_obs_states_,SEXP num_subjects_);
Rcpp::List get_emission_init_states(Rcpp::List llfb, arma::icolvec obs, int num_states, int num_obs_states);

RcppExport SEXP E_step_emission_init_all(SEXP likelihood_forward_backward_list,SEXP obs_data_list,SEXP num_states_,SEXP num_obs_states_,SEXP num_subjects_){

 Rcpp::List llfb_list(likelihood_forward_backward_list);
 Rcpp::List obsdata_list(obs_data_list);
 int num_subjects=Rcpp::as<int>(num_subjects_);

 int num_states=Rcpp::as<int>(num_states_);
 int num_obs_states=Rcpp::as<int>(num_obs_states_);
 arma::cube emission_array=arma::zeros<arma::cube>(num_states,num_obs_states,num_subjects);
 arma::mat init_array=arma::zeros<arma::mat>(num_states,num_subjects);
  for(int i=0; i<num_subjects;i++){
    arma::icolvec obs=Rcpp::as<arma::icolvec>(obsdata_list[i]);
    Rcpp::List out=get_emission_init_states(llfb_list[i], obs, num_states, num_obs_states);
    emission_array.slice(i)=Rcpp::as<arma::mat>(out["emission_states"]);
		init_array.col(i)=arma::strans(Rcpp::as<arma::mat>(out["init_states"]));
    
 }
return(Rcpp::List::create(Rcpp::Named("emission_array")=emission_array,Rcpp::Named("init_array")=init_array));

}

Rcpp::List get_emission_init_states(Rcpp::List llfb, arma::icolvec obs, int num_states, int num_obs_states){
 
 arma::mat alpha=Rcpp::as<arma::mat>(llfb["logalpha"]);
 arma::mat beta=Rcpp::as<arma::mat>(llfb["logbeta"]);
 double LL=Rcpp::as<double>(llfb["LL"]);
 
 arma::mat probs=exp(alpha+beta-LL);
 arma::mat init=probs.row(0);
 arma::mat emission_out=arma::zeros<arma::mat>(num_states,num_obs_states);
 
 for(int k=1;k<=num_obs_states;k++){
   arma::mat probs_temp=probs;
   for(int j=0;j<obs.n_elem;j++){
     
     if(obs(j)!=k){
        probs_temp.row(j)=arma::zeros<arma::rowvec>(num_states);
     }
     
   }
   emission_out.col(k-1)=strans(sum(probs_temp,0));
 
  
}
 return(Rcpp::List::create(Rcpp::Named("emission_states")=emission_out,Rcpp::Named("init_states")=init));
}
