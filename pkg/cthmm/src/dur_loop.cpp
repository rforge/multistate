#include "joint_mean_markov_rewards.h"
#include "auxmat_cpp.h"
#include "get_mean.h"

RcppExport SEXP dur_loop(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, 
  					 SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP time_dep_emission);
RcppExport SEXP dur_loop_h(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP obs_data_list, 
						SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP h_list_, SEXP Q_list_,SEXP time_dep_emission);

RcppExport SEXP dur_times(SEXP the_index, SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP the_state_size, SEXP exact_time_ranks_list,SEXP absorb_state);
arma::mat exact_dur(std::vector<arma::mat> joint_means_dur, arma::mat rate_matrix, arma::ivec absorb_states, int state, int exact_time_index);
arma::mat exact_dur2(arma::cube joint_means_dur, arma::mat rate_matrix, arma::ivec absorb_states, int state, int exact_time_index);

RcppExport SEXP dur_loop_h(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, 
						   SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state, SEXP h_list_, SEXP Q_list_,SEXP time_dep_emission){
	
	Rcpp::List llfb(likelihood_forward_backward);
	Rcpp::List time_diff(time_diffs_list);
	Rcpp::List eigen_list(eigen_decomp_list);
	Rcpp::List obsdata_list(obs_data_list);
	Rcpp::List emission(emission_list);
	Rcpp::List exact_time_ranks(exact_time_ranks_list);
	Rcpp::List h_list=Rcpp::as<Rcpp::List>(h_list_);
	Rcpp::List Q_list=Rcpp::as<Rcpp::List>(Q_list_);
	bool time_dep=Rcpp::as<bool>(time_dep_emission);
	int state_size=Rcpp::as<int>(the_state_size);
	int num_subjects=eigen_list.size();
	arma::mat expected_dur_array(state_size,num_subjects);
	
	for(int i=1; i<=state_size; i++){
		for(int k=0; k<num_subjects;k++){
			Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[k]);
			arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[k]);
			arma::icolvec obsdata=Rcpp::as<arma::icolvec>(obsdata_list[k]);
			Rcpp::List indivLLFB=Rcpp::as<Rcpp::List>(llfb[k]);
			arma::mat logalpha=Rcpp::as<arma::mat>(indivLLFB["logalpha"]);
			arma::mat logbeta=Rcpp::as<arma::mat>(indivLLFB["logbeta"]);
			arma::mat rate_matrix=Rcpp::as<arma::mat>(indivEigen["rate"]); 
			Rcpp::NumericVector indivH=Rcpp::as<Rcpp::NumericVector>(h_list[k]);
			arma::mat indivQ=Rcpp::as<arma::mat>(Q_list[k]);
			
			double LL=Rcpp::as<double>(indivLLFB["LL"]);
			
			int size=indivTime.n_elem;
			std::vector<arma::mat> joint_means_dur(size);
			for(int l=0; l<size;l++){
				joint_means_dur[l]=joint_mean_markov_rewards_cpp(i-1, indivTime(l), indivEigen);
				if(indivH[l+1]==1){
					joint_means_dur[l]=joint_means_dur[l]*abs(indivQ);
				}
			}
			
			
			if(Rf_isNull(exact_time_ranks[k])==0){
				int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
				arma::ivec the_absorb_states=Rcpp::as<arma::ivec>(absorb_state);
				joint_means_dur[exact_time_index]=exact_dur(joint_means_dur,rate_matrix,the_absorb_states,i,exact_time_index);
			}
      
			if(time_dep==false){
				arma::mat indivEmission=Rcpp::as<arma::mat>(emission[k]);
				expected_dur_array(i-1,k)=get_mean_cpp(joint_means_dur,obsdata, logalpha, logbeta, LL, indivEmission);
			}else{
				Rcpp::List indivEmissionList=Rcpp::as<Rcpp::List>(emission[k]);  
				expected_dur_array(i-1,k)=get_mean_cpp_time_dep_emission(joint_means_dur,obsdata, logalpha, logbeta, LL, indivEmissionList);
			}
					
		}
	}
	return(Rcpp::wrap(expected_dur_array));
}


RcppExport SEXP dur_loop(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, 
						 SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP time_dep_emission){
	
	Rcpp::List llfb(likelihood_forward_backward);
	Rcpp::List time_diff(time_diffs_list);
	Rcpp::List eigen_list(eigen_decomp_list);
	Rcpp::List obsdata_list(obs_data_list);
	Rcpp::List emission(emission_list);
	Rcpp::List exact_time_ranks(exact_time_ranks_list);
	bool time_dep=Rcpp::as<bool>(time_dep_emission);
	
	int state_size=Rcpp::as<int>(the_state_size);
	int num_subjects=eigen_list.size();
	
	arma::mat expected_dur_array(state_size,num_subjects);
	
	for(int i=1; i<=state_size; i++){
		for(int k=0; k<num_subjects;k++){
			
			
			Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[k]);
			arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[k]);
			arma::icolvec obsdata=Rcpp::as<arma::icolvec>(obsdata_list[k]);
			Rcpp::List indivLLFB=Rcpp::as<Rcpp::List>(llfb[k]);
			arma::mat logalpha=Rcpp::as<arma::mat>(indivLLFB["logalpha"]);
			arma::mat logbeta=Rcpp::as<arma::mat>(indivLLFB["logbeta"]);
			arma::mat rate_matrix=Rcpp::as<arma::mat>(indivEigen["rate"]); 
			
			double LL=Rcpp::as<double>(indivLLFB["LL"]);
			
			int size=indivTime.n_elem;
			std::vector<arma::mat> joint_means_dur(size);
			for(int l=0; l<size;l++){
				joint_means_dur[l]=joint_mean_markov_rewards_cpp(i-1, indivTime(l), indivEigen);
			}
			
			if(Rf_isNull(exact_time_ranks[k])==0){
				int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
				arma::ivec the_absorb_states=Rcpp::as<arma::ivec>(absorb_state);
				joint_means_dur[exact_time_index]=exact_dur(joint_means_dur,rate_matrix,the_absorb_states,i,exact_time_index);
			}
			if(time_dep==false){
				arma::mat indivEmission=Rcpp::as<arma::mat>(emission[k]);
				expected_dur_array(i-1,k)=get_mean_cpp(joint_means_dur,obsdata, logalpha, logbeta, LL, indivEmission);
			}else{
				Rcpp::List indivEmissionList=Rcpp::as<Rcpp::List>(emission[k]);  
				expected_dur_array(i-1,k)=get_mean_cpp_time_dep_emission(joint_means_dur,obsdata, logalpha, logbeta, LL, indivEmissionList);
			}
			
		}
	}
	return(Rcpp::wrap(expected_dur_array));
}


RcppExport SEXP dur_times(SEXP the_index, SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP the_state_size, SEXP exact_time_ranks_list,SEXP absorb_state){
	
	int test=5;
	Rcpp::List time_diff(time_diffs_list);
	Rcpp::List eigen_list(eigen_decomp_list);
	Rcpp::List exact_time_ranks(exact_time_ranks_list);
	
	int state_size=Rcpp::as<int>(the_state_size);
	int index=Rcpp::as<int>(the_index);
	int num_subjects=eigen_list.size();

	Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[0]);
	arma::mat rate_matrix=Rcpp::as<arma::mat>(indivEigen["rate"]); 
	
	arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[0]);
	int size=indivTime.n_elem;
	arma::cube joint_means_dur(state_size, state_size, size); 
	int k=0;
	int i=index;
	
	for(int l=0; l<size;l++){   
		joint_means_dur.slice(l)=joint_mean_markov_rewards_cpp(i-1, indivTime(l), indivEigen);
	}
	
	if(Rf_isNull(exact_time_ranks[k])==0){
		int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
		arma::ivec the_absorb_states=Rcpp::as<arma::ivec>(absorb_state);
		joint_means_dur.slice(exact_time_index)=exact_dur2(joint_means_dur,rate_matrix,the_absorb_states,i,exact_time_index);
		//std::cout<<joint_means_dur[exact_time_index]; 
	}
	
	
	return(Rcpp::wrap(joint_means_dur));

	
}


arma::mat exact_dur(std::vector<arma::mat> joint_means_dur, arma::mat rate_matrix, arma::ivec absorb_states, int state, int exact_time_index){
	arma::mat out=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	arma::mat temp=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	state=state-1;
	int k=0;
	int i=0;
	bool in_list=0;
	
	while(in_list==0 && i<absorb_states.size()){
		int test=absorb_states[i]-1;
		if(test==state){
			in_list=1;
			//std::cout<<in_list;
		}
		i++;
	}
	
	if(in_list==0){
		for(k=0;k<absorb_states.size();k++){
			int absorb_state=absorb_states[k]-1;
			temp.col(absorb_state)=rate_matrix.col(absorb_state);
		}
		
		out=joint_means_dur[exact_time_index]*temp;
	}
	return(out);
	
}

arma::mat exact_dur2(arma::cube joint_means_dur, arma::mat rate_matrix, arma::ivec absorb_states, int state, int exact_time_index){
	arma::mat out=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	arma::mat temp=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	state=state-1;
	int k=0;
	int i=0;
	bool in_list=0;
	
	while(in_list==0 && i<absorb_states.size()){
		int test=absorb_states[i]-1;
		if(test==state){
			in_list=1;
			//std::cout<<in_list;
		}
		i++;
	}
	
	if(in_list==0){
		for(k=0;k<absorb_states.size();k++){
			int absorb_state=absorb_states[k]-1;
			temp.col(absorb_state)=rate_matrix.col(absorb_state);
		}
		
		out=joint_means_dur.slice(exact_time_index)*temp;
	}
	return(out);
}
