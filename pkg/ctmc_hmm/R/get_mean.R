get_mean<-function(single_joint_mean, obs_data, forward, backward, LL, emission_matrix){
      out<-.Call( "get_mean", single_joint_mean, obs_data, forward, backward, LL, emission_matrix)
     return(out)
}
      
