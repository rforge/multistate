trans_loop_R<-function(likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
     exact_time_ranks_list, the_state_size,absorb_state, ij_indices,h_list=NULL,Q=NULL,time_dep_emission=F){
	
    if(!is.null(h_list)){
		out<-.Call( "trans_loop_h",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
				   exact_time_ranks_list, the_state_size,absorb_state, ij_indices,h_list,Q,time_dep_emission)
	}else{
		out<-.Call( "trans_loop",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, 
				 obs_data_list, emission_list, exact_time_ranks_list, the_state_size,absorb_state, ij_indices,time_dep_emission)
    }
		return(out)
}

dur_loop_R<-function(likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
exact_time_ranks_list,the_state_size,absorb_state=NULL,h_list=NULL,Q=NULL,time_dep_emission=F){
	if(!is.null(h_list)){
	  out<-.Call( "dur_loop_h",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
			   exact_time_ranks_list,the_state_size,absorb_state,h_list,Q,time_dep_emission)
	}else{
	  out<-.Call( "dur_loop",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
				   exact_time_ranks_list,the_state_size,absorb_state,time_dep_emission)
    }
	return(out)
}

#These functions do not allow for informative sampling times or time dependent emission matrices
trans_times_R<-function(i,j,time_diffs,eigen_decomp,num_states,exact_time_ranks,absorb_state){
	
	time_diffs_list=list(time_diffs)
	eigen_decomp_list=list(eigen_decomp)
	exact_time_ranks_list=list(exact_time_ranks)
	
	out<-.Call( "trans_times",i,j,time_diffs_list, eigen_decomp_list, num_states,exact_time_ranks_list,absorb_state)
	return(out)
}

dur_times_R<-function(i,time_diffs,eigen_decomp,num_states,exact_time_ranks,absorb_state){
	
	time_diffs_list=list(time_diffs)
	eigen_decomp_list=list(eigen_decomp)
	exact_time_ranks_list=list(exact_time_ranks)
	
	out<-.Call( "dur_times",i,time_diffs_list, eigen_decomp_list, num_states,exact_time_ranks_list,absorb_state)
	return(out)
}

