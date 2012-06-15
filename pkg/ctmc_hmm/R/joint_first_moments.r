joint_mean_markov_jumps <- function(rate.eigen, regist.matrix, interval.len){
        out<-.Call( "joint_mean_markov_jumps", rate.eigen, regist.matrix, interval.len)
        return(out)
}

joint_mean_markov_rewards <- function(i,interval.len, rate.eigen){
	out<-.Call( "joint_mean_markov_rewards", i, interval.len, rate.eigen) 
	return(out)
} 


uniformization_mean <- function(B,interval.len, rate){
	out<-.Call( "uniformization_mean_R", interval.len, B, rate) 
	return(out)
} 


