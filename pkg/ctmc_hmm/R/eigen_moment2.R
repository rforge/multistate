joint_trans_trans_times <- function(time_diffs, j1,k1,j2,k2,absorb_state,eigen_decomp, exact_time_rank){
        j1_=j1-1
        j2_=j2-1
        k1_=k1-1
        k2_=k2-1
        absorb_state=absorb_state-1
        out<-.Call( "joint_trans_trans_times", time_diffs, j1_,k1_,j2_,k2_,absorb_state,eigen_decomp, exact_time_rank)
        return(out)
}

joint_dur_dur_times <- function(time_diffs, di,dj,absorb_state,eigen_decomp, exact_time_rank){
       i_=di-1
       j_=dj-1      
       absorb_state=absorb_state-1
       out<-.Call( "joint_dur_dur_times", time_diffs, i_,j_,absorb_state,eigen_decomp, exact_time_rank)
       return(out)
}

joint_dur_trans_times <- function(time_diffs, di, tr1, tr2,absorb_state,eigen_decomp, exact_time_rank){
       d1_=di-1
       tr1_=tr1-1
       tr2_=tr2-1
       absorb_state=absorb_state-1
       out<-.Call( "joint_trans_dur_times", time_diffs, d1_,tr1_,tr2_,absorb_state,eigen_decomp, exact_time_rank)
       return(out)
}


