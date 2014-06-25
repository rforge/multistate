#include "joint_mean_markov_jumps.h"
#include "auxmat_cpp.h"
#include "get_mean.h"
#include "mat_exp_eigen.h"


RcppExport SEXP trans_loop(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP ij_indices,SEXP time_dep_emission);
RcppExport SEXP trans_loop_h(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state, SEXP ij_indices,SEXP h_list_, SEXP Q_list_,SEXP time_dep_emission);

RcppExport SEXP trans_times(SEXP ni,SEXP nj, SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP the_state_size, SEXP exact_time_ranks_list,SEXP absorb_state);
arma::mat exact_trans(std::vector<arma::mat> joint_means_trans, Rcpp::List eigen_decomp, double time_int, arma::ivec absorb_states, int start_state, int end_state, int exact_time_index);
arma::mat exact_trans2(arma::cube joint_means_trans, Rcpp::List eigen_decomp, double time_int, arma::ivec absorb_states, int start_state, int end_state, int exact_time_index);

RcppExport SEXP trans_loop_h(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP obs_data_list, SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP ij_indices,SEXP h_list_, SEXP Q_list_,SEXP time_dep_emission){
	arma::imat state_pairs=Rcpp::as<arma::imat>(ij_indices);
	int num_state_pairs=state_pairs.n_rows;
	
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
	arma::cube expected_trans_array=arma::zeros<arma::cube>(state_size,state_size,num_subjects);
	
	for(int i=0; i<num_state_pairs; i++){
		for(int k=0; k<num_subjects;k++){
			Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[k]);
			arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[k]);
			arma::icolvec obsdata=Rcpp::as<arma::icolvec>(obsdata_list[k]);
			Rcpp::List indivLLFB=Rcpp::as<Rcpp::List>(llfb[k]);
			arma::mat logalpha=Rcpp::as<arma::mat>(indivLLFB["logalpha"]);
			arma::mat logbeta=Rcpp::as<arma::mat>(indivLLFB["logbeta"]);
			Rcpp::NumericVector indivH=Rcpp::as<Rcpp::NumericVector>(h_list[k]);
			arma::mat indivQ=Rcpp::as<arma::mat>(Q_list[k]);
			
			double LL=Rcpp::as<double>(indivLLFB["LL"]);
			
			int size=indivTime.n_elem;
			std::vector<arma::mat> joint_means_trans(size);
			arma::mat regist_matrix=arma::zeros<arma::mat>(state_size,state_size);
			regist_matrix(state_pairs(i,0)-1,state_pairs(i,1)-1)=1;
			for(int l=0; l<size;l++){
				
				
				joint_means_trans[l]=joint_mean_markov_jumps_cpp(indivEigen,regist_matrix, indivTime(l));
				if(indivH[l+1]==1){
					joint_means_trans[l]=joint_means_trans[l]*abs(indivQ);
				}
			//std::cout<<joint_means_trans[l];
			}
			
			if(Rf_isNull(exact_time_ranks[k])==0){
				int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
				arma::ivec the_absorb_state=Rcpp::as<arma::ivec>(absorb_state);
				int start_state=state_pairs(i,0);
				int end_state=state_pairs(i,1);
				joint_means_trans[exact_time_index]=exact_trans(joint_means_trans,indivEigen,indivTime(exact_time_index),the_absorb_state, start_state,end_state,exact_time_index);
			}
			
			if(time_dep==false){
				arma::mat indivEmission=Rcpp::as<arma::mat>(emission[k]);
				expected_trans_array(state_pairs(i,0)-1,state_pairs(i,1)-1,k)=get_mean_cpp(joint_means_trans,obsdata, logalpha, logbeta, LL, indivEmission);
			}else{
				Rcpp::List indivEmissionList=Rcpp::as<Rcpp::List>(emission[k]);  
				expected_trans_array(state_pairs(i,0)-1,state_pairs(i,1)-1,k)=get_mean_cpp_time_dep_emission(joint_means_trans,obsdata, logalpha, logbeta, LL, indivEmissionList);
			}
		}
	}
	return(Rcpp::wrap(expected_trans_array));
}


RcppExport SEXP trans_loop(SEXP likelihood_forward_backward,SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP obs_data_list, 
						   SEXP emission_list, SEXP exact_time_ranks_list,SEXP the_state_size,SEXP absorb_state,SEXP ij_indices,SEXP time_dep_emission){
	arma::imat state_pairs=Rcpp::as<arma::imat>(ij_indices);
	int num_state_pairs=state_pairs.n_rows;
	
	Rcpp::List llfb(likelihood_forward_backward);
	Rcpp::List time_diff(time_diffs_list);
	Rcpp::List eigen_list(eigen_decomp_list);
	Rcpp::List obsdata_list(obs_data_list);
	Rcpp::List emission(emission_list);
	Rcpp::List exact_time_ranks(exact_time_ranks_list);
	
	bool time_dep=Rcpp::as<bool>(time_dep_emission);
	
	int state_size=Rcpp::as<int>(the_state_size);
	int num_subjects=eigen_list.size();
	arma::cube expected_trans_array=arma::zeros<arma::cube>(state_size,state_size,num_subjects);
	
	for(int i=0; i<num_state_pairs; i++){
		for(int k=0; k<num_subjects;k++){
			Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[k]);
			arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[k]);
			arma::icolvec obsdata=Rcpp::as<arma::icolvec>(obsdata_list[k]);
			Rcpp::List indivLLFB=Rcpp::as<Rcpp::List>(llfb[k]);
			arma::mat logalpha=Rcpp::as<arma::mat>(indivLLFB["logalpha"]);
			arma::mat logbeta=Rcpp::as<arma::mat>(indivLLFB["logbeta"]);
			
			double LL=Rcpp::as<double>(indivLLFB["LL"]);
			
			int size=indivTime.n_elem;
			std::vector<arma::mat> joint_means_trans(size);
			arma::mat regist_matrix=arma::zeros<arma::mat>(state_size,state_size);
			regist_matrix(state_pairs(i,0)-1,state_pairs(i,1)-1)=1;
			for(int l=0; l<size;l++){
				joint_means_trans[l]=joint_mean_markov_jumps_cpp(indivEigen,regist_matrix, indivTime(l));
			}
			
			if(Rf_isNull(exact_time_ranks[k])==0){
				int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
				arma::ivec the_absorb_state=Rcpp::as<arma::ivec>(absorb_state);
				int start_state=state_pairs(i,0);
				int end_state=state_pairs(i,1);
				joint_means_trans[exact_time_index]=exact_trans(joint_means_trans,indivEigen,indivTime(exact_time_index),the_absorb_state, start_state,end_state,exact_time_index);
			}
			
			
			if(time_dep==false){
				arma::mat indivEmission=Rcpp::as<arma::mat>(emission[k]);
				expected_trans_array(state_pairs(i,0)-1,state_pairs(i,1)-1,k)=get_mean_cpp(joint_means_trans,obsdata, logalpha, logbeta, LL, indivEmission);
			}else{
				Rcpp::List indivEmissionList=Rcpp::as<Rcpp::List>(emission[k]);  
				expected_trans_array(state_pairs(i,0)-1,state_pairs(i,1)-1,k)=get_mean_cpp_time_dep_emission(joint_means_trans,obsdata, logalpha, logbeta, LL, indivEmissionList);
			}
			
		}
	}
	return(Rcpp::wrap(expected_trans_array));
}



RcppExport SEXP trans_times(SEXP ni,SEXP nj, SEXP time_diffs_list, SEXP eigen_decomp_list, SEXP the_state_size, SEXP exact_time_ranks_list,SEXP absorb_state){

    int test=3;
	Rcpp::List time_diff(time_diffs_list);
	Rcpp::List eigen_list(eigen_decomp_list);
	Rcpp::List exact_time_ranks(exact_time_ranks_list);
	
	int state_size=Rcpp::as<int>(the_state_size);
	int num_subjects=eigen_list.size();

	Rcpp::List indivEigen=Rcpp::as<Rcpp::List>(eigen_list[0]);
	arma::mat rate_matrix=Rcpp::as<arma::mat>(indivEigen["rate"]); 
	
	arma::vec indivTime=Rcpp::as<arma::vec>(time_diff[0]);
	int size=indivTime.n_elem;
    arma::cube joint_means_trans(state_size, state_size, size); 

	int k=0;
	int i=Rcpp::as<int>(ni);
	int j=Rcpp::as<int>(nj);
	
	arma::mat regist_matrix=arma::zeros<arma::mat>(state_size,state_size);
	
	regist_matrix(i-1,j-1)=1;
	for(int l=0; l<size;l++){
		joint_means_trans.slice(l)=joint_mean_markov_jumps_cpp(indivEigen,regist_matrix, indivTime(l));
	}
	
	if(Rf_isNull(exact_time_ranks[k])==0){
		int exact_time_index=Rcpp::as<int>(exact_time_ranks[k])-2;  
		arma::ivec the_absorb_state=Rcpp::as<arma::ivec>(absorb_state);
		int start_state=i;
		int end_state=j;
		joint_means_trans.slice(exact_time_index)=exact_trans2(joint_means_trans,indivEigen,indivTime(exact_time_index),the_absorb_state,   start_state,end_state,exact_time_index);
    }
	
return(Rcpp::wrap(joint_means_trans));

//	return(Rcpp::wrap(test));
}


arma::mat exact_trans(std::vector<arma::mat> joint_means_trans, Rcpp::List eigen_decomp, double time_int, arma::ivec absorb_states, int start_state, int end_state, int exact_time_index){
	arma::mat rate_matrix=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	arma::mat out=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	arma::mat temp=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	int i=start_state-1;
	int j=end_state-1;
	int k=0;  
	
	bool i_in_A=0;
	bool j_in_A=0;
	
	//std::cout<<absorb_states;
	
	while(i_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k]-1;
		if(test==i){
			i_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(j_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k]-1;
		if(test==j){
			j_in_A=1;
		}
		k++;
	}
	
	int in_either=i_in_A+j_in_A;
	
	if(in_either==0){
		for(int l=0;l<absorb_states.size();l++){
			int absorb_state=absorb_states[l]-1;
			temp.col(absorb_state)=rate_matrix.col(absorb_state);
		}
		out=joint_means_trans[exact_time_index]*temp;
	}
	if(i_in_A==0 && j_in_A==1){
		arma::mat prob_mat=mat_exp_eigen_cpp(eigen_decomp,time_int);
		out.col(j)=prob_mat.col(i)*rate_matrix(i,j);
	}
	
	
	return(out);
}

arma::mat exact_trans2(arma::cube joint_means_trans, Rcpp::List eigen_decomp, double time_int, arma::ivec absorb_states, int start_state, int end_state, int exact_time_index){
	arma::mat rate_matrix=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	arma::mat out=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	arma::mat temp=arma::zeros<arma::mat>(rate_matrix.n_rows,rate_matrix.n_rows);
	
	int i=start_state-1;
	int j=end_state-1;
	int k=0;  
	
	bool i_in_A=0;
	bool j_in_A=0;
	
	//std::cout<<absorb_states;
	
	while(i_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k]-1;
		if(test==i){
			i_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(j_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k]-1;
		if(test==j){
			j_in_A=1;
		}
		k++;
	}
	
	int in_either=i_in_A+j_in_A;

	
	if(in_either==0){
		for(int l=0;l<absorb_states.size();l++){
			int absorb_state=absorb_states[l]-1;
			temp.col(absorb_state)=rate_matrix.col(absorb_state);
		}
		out=joint_means_trans.slice(exact_time_index)*temp;
	}
	if(i_in_A==0 && j_in_A==1){
		arma::mat prob_mat=mat_exp_eigen_cpp(eigen_decomp,time_int);
		out.col(j)=prob_mat.col(i)*rate_matrix(i,j);
	}
	
	
	return(out);
}	

