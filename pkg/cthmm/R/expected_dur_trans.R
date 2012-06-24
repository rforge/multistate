
trans_loop_R<-function(likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, exact_time_ranks_list, the_state_size,absorb_state, ij_indices){
      out<-.Call( "trans_loop",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, exact_time_ranks_list, the_state_size,absorb_state, ij_indices)
     return(out)
}


trans_times_R<-function(i,j,time_diffs,eigen_decomp,num_states,exact_time_ranks,absorb_state){
	
	time_diffs_list=list(time_diffs)
	eigen_decomp_list=list(eigen_decomp)
	exact_time_ranks_list=list(exact_time_ranks)
	
	out<-.Call( "trans_times",i,j,time_diffs_list, eigen_decomp_list, num_states,exact_time_ranks_list,absorb_state)
	return(out)
}

dur_loop_R<-function(likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
exact_time_ranks_list,the_state_size,absorb_state=NULL){
	out<-.Call( "dur_loop",likelihood_forward_backward,time_diffs_list, eigen_decomp_list, obs_data_list, emission_list, 
			   exact_time_ranks_list,the_state_size,absorb_state)
	return(out)
}


dur_times_R<-function(i,time_diffs,eigen_decomp,num_states,exact_time_ranks,absorb_state){
	
	time_diffs_list=list(time_diffs)
	eigen_decomp_list=list(eigen_decomp)
	exact_time_ranks_list=list(exact_time_ranks)
	
	out<-.Call( "dur_times",i,time_diffs_list, eigen_decomp_list, num_states,exact_time_ranks_list,absorb_state)
	return(out)
}

