#include "eigen_moments2.h"
#include "joint_mean_markov_jumps.h"
#include "joint_mean_markov_rewards.h"
#include "auxmat_cpp.h"
#include "mat_exp_eigen.h"


arma::mat joint_duration_2moment(int i, int j,Rcpp::List eigen_decomp, double interval_len){
	
	arma::cx_colvec v=Rcpp::as<arma::cx_colvec>(eigen_decomp["values"]);
	double t=interval_len;
	int size=v.n_elem;
	arma::cx_mat out(size,size);
	
	for(int a = 0; a < size; a++){
		for(int b = 0; b < size; b++){
			out(a,b)=hobolth_fun(eigen_decomp, t, a, i, i, j, j, b)(0)+hobolth_fun(eigen_decomp, t, a, j, j, i, i, b)(0);
		}
	}
	arma::mat real_out=arma::real(out);
	return(real_out);
}
arma:: mat joint_dur_trans_2moment(int d1, int tr1, int tr2, Rcpp::List eigen_decomp, double interval_len){
	
	double t=interval_len;
	
	arma::cx_colvec v=Rcpp::as<arma::cx_colvec>(eigen_decomp["values"]);
	arma::mat Q=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	
	int size=v.n_elem;
	arma::cx_mat out(size,size);
	
	for(int a = 0; a < size; a++){
		for(int b = 0; b < size; b++){
			out(a,b)=Q(tr1,tr2)*(hobolth_fun(eigen_decomp, t, a,tr1,tr2,d1,d1,b)(0)+hobolth_fun(eigen_decomp, t, a,d1,d1,tr1,tr2,b)(0));
		}
	}
	arma::mat real_out=arma::real(out);
	return(real_out);
}

arma::mat joint_transition_2moment(int j1, int k1, int j2, int k2,Rcpp::List eigen_decomp, double interval_len){
	double t=interval_len;
	arma::cx_colvec v=Rcpp::as<arma::cx_colvec>(eigen_decomp["values"]);
	arma::mat Q=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	int size=v.n_elem;
	arma::mat mean_mat=arma::zeros<arma::mat>(size,size);
	arma::cx_mat out(size,size);
	if(j1==j2&k1==k2){
		arma::mat regist_matrix=arma::zeros<arma::mat>(size,size);
		regist_matrix(j1,k1)=1;
		mean_mat=joint_mean_markov_jumps_cpp(eigen_decomp,regist_matrix, t);
		
	}
	
	for(int a = 0; a < size; a++){
		for(int b = 0; b < size; b++){
			out(a,b)=mean_mat(a,b)+Q(j1,k1)*Q(j2,k2)*(hobolth_fun(eigen_decomp, t, a,j1,k1,j2,k2,b)(0)+hobolth_fun(eigen_decomp, t, a,j2,k2,j1,k1,b)(0));
		}
	}
	arma::mat real_out=arma::real(out);
	return(real_out);
}

arma::mat joint_dur_dur_exact_time(int i,int j, double interval_len, arma::ivec absorb_states,Rcpp::List eigen_decomp){
	arma::cx_mat Q=Rcpp::as<arma::cx_mat>(eigen_decomp["rate"]);
	double t=interval_len;
	arma::mat joint_dur=joint_duration_2moment(i, j,eigen_decomp,t);
	arma::cx_mat temp=arma::zeros<arma::cx_mat>(Q.n_rows,Q.n_rows);
	
	bool i_in_A=0;
	bool j_in_A=0;
	int k=0;	
	//std::cout<<absorb_states;
	
	while(i_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==i){
			i_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(j_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==j){
			j_in_A=1;
		}
		k++;
	}
	
	int in_either=i_in_A+j_in_A;
	
	if(i_in_A+j_in_A==0){
		
		for(int l=0;l<absorb_states.size();l++){
			int absorb_state=absorb_states[l];
			temp.col(absorb_state)=Q.col(absorb_state);
		}
		
		temp=joint_dur*temp;
		
		
	}else{
		temp=arma::zeros<arma::cx_mat>(Q.n_rows,Q.n_rows);
	}
	arma::mat out=arma::real(temp);   
	return(out);
	
}


arma::mat joint_dur_trans_exact_time(int d1, int tr1, int tr2, arma::ivec absorb_state, Rcpp::List eigen_decomp, double interval_len){
    
	arma::ivec absorb_states=absorb_state;
	double t=interval_len;
	
	arma::cx_mat Q=Rcpp::as<arma::cx_mat>(eigen_decomp["rate"]);
	arma::mat H=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
	arma::cx_mat temp=arma::zeros<arma::cx_mat>(Q.n_rows,Q.n_rows);
	
	bool d1_in_A=0;
	bool tr1_in_A=0;
	bool tr2_in_A=0;
	int k=0;  
	
	while(d1_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==d1){
			d1_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(tr1_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==tr1){
			tr1_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(tr2_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==tr2){
			tr2_in_A=1;
		}
		k++;
	}
	
	if(d1_in_A==0 && tr1_in_A==0 && tr2_in_A==0){
		for(int l=0;l<absorb_states.size();l++){
			int absorb_state=absorb_states[l];
			temp.col(absorb_state)=Q.col(absorb_state);
		}
		arma::mat joint_dur_trans=joint_dur_trans_2moment(d1, tr1, tr2, eigen_decomp, t);
		temp=joint_dur_trans*temp;
	}else if((d1_in_A==0 && tr1_in_A==0) && tr2_in_A==1){
		H.col(tr1)=joint_mean_markov_rewards_cpp(d1, t, eigen_decomp).col(tr1);
		temp.col(tr2)=Q.col(tr2);
		temp=H*temp;
	}
	
	arma::mat out=arma::real(temp);
	return(out);
}

arma::mat joint_trans_trans_exact_time(int j1,int k1,int j2,int k2,arma::ivec absorb_state,Rcpp::List eigen_decomp,double interval_len){
	double t=interval_len;
	arma::ivec absorb_states=absorb_state;
	
	arma::mat Q=Rcpp::as<arma::mat>(eigen_decomp["rate"]);
	arma::mat temp=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
	
	arma::mat M1=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
	
	bool j1_in_A=0;
	bool j2_in_A=0;
	bool k1_in_A=0;
	bool k2_in_A=0;
	
	int k=0;  
	
	while(j1_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==j1){
			j1_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(j2_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==j2){
			j2_in_A=1;
		}
		k++;
	}
	
	k=0;
	while(k1_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==k1){
			k1_in_A=1;
		}
		k++;
	}
	k=0;
	while(k2_in_A==0 && k<absorb_states.size()){
		int test=absorb_states[k];
		if(test==k2){
			k2_in_A=1;
		}
		k++;
	} 
	if(j1_in_A==0 && j2_in_A==0 && k1_in_A==0 && k2_in_A==0){
		arma::mat M2=joint_transition_2moment(j1, k1, j2, k2,eigen_decomp, t);
		for(int l=0;l<absorb_states.size();l++){
			int absorb_state=absorb_states[l];
			temp.col(absorb_state)=Q.col(absorb_state);
		}
		temp=M2*temp;
	}else if(j1_in_A==0 && j2_in_A==0 && k1_in_A==1 && k2_in_A==0){
		arma::mat regist_matrix=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
		regist_matrix(j2,k2)=1;
		M1.col(j1)=joint_mean_markov_jumps_cpp(eigen_decomp,regist_matrix, t).col(j1);
		temp.col(k1)=Q.col(k1);
		temp=M1*temp;
	}else if(j1_in_A==0 && j2_in_A==0 && k1_in_A==0 && k2_in_A==1){
		arma::mat regist_matrix=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
		regist_matrix(j1,k1)=1;
		M1.col(j2)=joint_mean_markov_jumps_cpp(eigen_decomp,regist_matrix, t).col(j2);
		temp.col(k2)=Q.col(k2);
		temp=M1*temp;
	}else if(j1_in_A==0 &&j2_in_A==0 &&j1==j2 && k1_in_A==1 && k2_in_A==1 && k1==k2){
		arma::mat P=arma::zeros<arma::mat>(Q.n_rows,Q.n_rows);
		P.col(j1)=mat_exp_eigen_cpp(eigen_decomp,t).col(j1);
		temp.col(k2)=Q.col(k2);	
		temp=P*temp;
    }
	
	return(temp);
}


RcppExport SEXP joint_trans_trans_times(SEXP time_diffs,SEXP j1_,SEXP k1_,SEXP j2_,SEXP k2_,SEXP absorb_state,SEXP eigen_decomp,SEXP exact_time_rank){
	
	int j1=Rcpp::as<int>(j1_);
	int k1=Rcpp::as<int>(k1_);
	int j2=Rcpp::as<int>(j2_);
	int k2=Rcpp::as<int>(k2_);
	Rcpp::List eigen_dec=Rcpp::as<Rcpp::List>(eigen_decomp);
	
	arma::colvec times=Rcpp::as<arma::colvec>(time_diffs);
	arma::mat Q=Rcpp::as<arma::mat>(eigen_dec["rate"]);
	arma::cube out=arma::zeros<arma::cube>(Q.n_rows,Q.n_rows,times.n_elem);
	
	for(int l=0; l<times.n_elem;l++){   
		if(Rf_isNull(exact_time_rank)==0){
			arma::ivec absorb=Rcpp::as<arma::ivec>(absorb_state);
			int exact_time_index=Rcpp::as<int>(exact_time_rank)-2;   
			if(l==exact_time_index){
				out.slice(l)=joint_trans_trans_exact_time(j1,k1,j2,k2,absorb,eigen_decomp,times(l));
			}else{
				out.slice(l)=joint_transition_2moment(j1, k1, j2, k2,eigen_decomp, times(l));
			}
		}else{
			out.slice(l)=joint_transition_2moment(j1, k1, j2, k2,eigen_decomp, times(l));
		} 
	}
	
	return(Rcpp::wrap(out));
}

RcppExport SEXP joint_dur_dur_times(SEXP time_diffs,SEXP i_,SEXP j_,SEXP absorb_states,SEXP eigen_decomp,SEXP exact_time_rank){ 
	int i=Rcpp::as<int>(j_);
	int j=Rcpp::as<int>(i_);
	Rcpp::List eigen_dec=Rcpp::as<Rcpp::List>(eigen_decomp);
	
	arma::colvec times=Rcpp::as<arma::colvec>(time_diffs);
	arma::mat Q=Rcpp::as<arma::mat>(eigen_dec["rate"]);
	arma::cube out=arma::zeros<arma::cube>(Q.n_rows,Q.n_rows,times.n_elem);
	
	for(int l=0; l<times.n_elem;l++){   
		if(Rf_isNull(exact_time_rank)==0){
			arma::ivec absorb=Rcpp::as<arma::ivec>(absorb_states);
			int exact_time_index=Rcpp::as<int>(exact_time_rank)-2;   
			if(l==exact_time_index){
				out.slice(l)=joint_dur_dur_exact_time(i,j, times(l), absorb,eigen_decomp);
			}else{
				out.slice(l)=joint_duration_2moment(i, j,eigen_decomp, times(l));
			}
		}else{
			out.slice(l)=joint_duration_2moment(i, j,eigen_decomp, times(l));
		} 
	}
	
	return(Rcpp::wrap(out));
}


RcppExport SEXP joint_trans_dur_times(SEXP time_diffs,SEXP d1_,SEXP tr1_,SEXP tr2_,SEXP absorb_state,SEXP eigen_decomp,SEXP exact_time_rank){
	int d1=Rcpp::as<int>(d1_);
	int tr1=Rcpp::as<int>(tr1_);
	int tr2=Rcpp::as<int>(tr2_);
	Rcpp::List eigen_dec=Rcpp::as<Rcpp::List>(eigen_decomp);
	
	arma::colvec times=Rcpp::as<arma::colvec>(time_diffs);
	arma::mat Q=Rcpp::as<arma::mat>(eigen_dec["rate"]);
	arma::cube out=arma::zeros<arma::cube>(Q.n_rows,Q.n_rows,times.n_elem);
	
	for(int l=0; l<times.n_elem;l++){   
		if(Rf_isNull(exact_time_rank)==0){
			arma::ivec absorb=Rcpp::as<arma::ivec>(absorb_state);
			int exact_time_index=Rcpp::as<int>(exact_time_rank)-2;   
			if(l==exact_time_index){
				out.slice(l)= joint_dur_trans_exact_time(d1, tr1, tr2, absorb, eigen_decomp, times(l));
			}else{
				out.slice(l)=joint_dur_trans_2moment(d1, tr1, tr2, eigen_decomp, times(l));
			}
		}else{
			out.slice(l)=joint_dur_trans_2moment(d1, tr1, tr2, eigen_decomp, times(l));
		} 
	}
	
	return(Rcpp::wrap(out));
}


