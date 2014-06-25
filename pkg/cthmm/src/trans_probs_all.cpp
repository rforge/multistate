#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat uniformization_prob(double t, arma::mat rate);
arma::mat mat_exp_eigen_cpp2(Rcpp::List eigen_decomp, double time_int);
Rcpp::List trans_prob(Rcpp::List eigen_decomp, Rcpp::NumericVector time_int,int exact_time_ranks);
RcppExport SEXP trans_probs_all(SEXP eigen_decomp_list_, SEXP time_int_list_, SEXP exact_times_list_);
Rcpp::List trans_prob_h(Rcpp::List eigen_decomp, Rcpp::NumericVector time_int,int exact_time_ranks, Rcpp::NumericVector indivH, arma::mat indivQ);
RcppExport SEXP trans_probs_all_h(SEXP eigen_decomp_list_, SEXP time_int_list_, SEXP exact_times_list_, SEXP h_list_, SEXP Q_list_);

arma::mat uniformization_prob(double t, arma::mat rate){

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

arma::mat mat_exp_eigen_cpp2(Rcpp::List eigen_decomp, double time_int){
    arma::cx_mat vectors=Rcpp::as<arma::cx_mat>(eigen_decomp["vectors"]);
	arma::cx_mat invvectors=Rcpp::as<arma::cx_mat>(eigen_decomp["invvectors"]);
	arma::cx_colvec values=Rcpp::as<arma::cx_colvec>(eigen_decomp["values"]);
	arma::mat rate=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
  
  arma::mat out;
  bool replicate=Rcpp::as<bool>(eigen_decomp["replicate"]);

  if(replicate==1){
  out=uniformization_prob(time_int,rate);

  }else{
	out=real(vectors * diagmat(exp(values * time_int)) * invvectors);
  }
  return(out);
}

Rcpp::List trans_prob(Rcpp::List eigen_decomp, Rcpp::NumericVector time_int,int exact_time_ranks){
 	Rcpp::List out(time_int.size());
	for(int m=0; m<time_int.size();m++){
		out[m]=mat_exp_eigen_cpp2(eigen_decomp,time_int[m]);
	} 
  
	
      if(exact_time_ranks>1){
      arma::mat rate=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
      rate.diag()=arma::zeros<arma::colvec>(rate.n_rows);
      int exact_time_index=exact_time_ranks-2; 
      arma::mat out_temp=Rcpp::as<arma::mat>(out[exact_time_index]);
      out[exact_time_index]=out_temp*rate;
     }
	
	return(wrap(out));
}

Rcpp::List trans_prob_h(Rcpp::List eigen_decomp, Rcpp::NumericVector time_int,int exact_time_ranks, Rcpp::NumericVector indivH, arma::mat indivQ){
 	Rcpp::List out(time_int.size());
	
	for(int m=0; m<time_int.size();m++){
		arma::mat rate=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
		rate.diag()=arma::zeros<arma::colvec>(rate.n_rows);

		out[m]=mat_exp_eigen_cpp2(eigen_decomp,time_int[m]);
			 
		if(indivH[m+1]==1){
			arma::mat out_temp=Rcpp::as<arma::mat>(out[m]);
			out[m]=out_temp*abs(indivQ);	
		}
	/*	if(indivH[m+1]==3){ //3 codes known time of X(t) into an absorbing state
			arma::mat out_temp=Rcpp::as<arma::mat>(out[m]);
			out[m]=out_temp*rate;
		} */
		
  	} 
	
	
	if(exact_time_ranks>1){
		arma::mat rate=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
		rate.diag()=arma::zeros<arma::colvec>(rate.n_rows);
		int exact_time_index=exact_time_ranks-2; 
		arma::mat out_temp=Rcpp::as<arma::mat>(out[exact_time_index]);
		out[exact_time_index]=out_temp*rate;
	}
	
	return(wrap(out));
}



RcppExport SEXP trans_probs_all(SEXP eigen_decomp_list_, SEXP time_int_list_, SEXP exact_times_list_){

 Rcpp::List eigen_list=Rcpp::as<Rcpp::List>(eigen_decomp_list_);
 Rcpp::List time_list=Rcpp::as<Rcpp::List>(time_int_list_);
 Rcpp::List exact_list=Rcpp::as<Rcpp::List>(exact_times_list_);
 int length=time_list.size();
 Rcpp::List out(length);

for(int i=0;i<length;i++){
   Rcpp::NumericVector indivTime=Rcpp::as<Rcpp::NumericVector>(time_list[i]);
		
  int exact_time=0;
   if(Rf_isNull(exact_list[i])==0){
   exact_time=Rcpp::as<int>(exact_list[i]);
  }
  out[i]=trans_prob(eigen_list[i], indivTime,exact_time);

}
return(Rcpp::wrap(out));

}



RcppExport SEXP trans_probs_all_h(SEXP eigen_decomp_list_, SEXP time_int_list_, SEXP exact_times_list_, SEXP h_list_, SEXP Q_list_){
	
	Rcpp::List eigen_list=Rcpp::as<Rcpp::List>(eigen_decomp_list_);
	Rcpp::List time_list=Rcpp::as<Rcpp::List>(time_int_list_);
	Rcpp::List exact_list=Rcpp::as<Rcpp::List>(exact_times_list_);
	Rcpp::List h_list=Rcpp::as<Rcpp::List>(h_list_);
	Rcpp::List Q_list=Rcpp::as<Rcpp::List>(Q_list_);
	

	int length=time_list.size();
	Rcpp::List out(length);
	
	for(int i=0;i<length;i++){
		Rcpp::NumericVector indivTime=Rcpp::as<Rcpp::NumericVector>(time_list[i]);
		Rcpp::NumericVector indivH=Rcpp::as<Rcpp::NumericVector>(h_list[i]);
		arma::mat indivQ=Rcpp::as<arma::mat>(Q_list[i]);

		int exact_time=0;
		if(Rf_isNull(exact_list[i])==0){
			exact_time=Rcpp::as<int>(exact_list[i]);
		}
		
		out[i]=trans_prob_h(eigen_list[i], indivTime,exact_time,indivH,indivQ);
		
	}
	return(Rcpp::wrap(out));
	
}
