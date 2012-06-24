#' @include mat_exp_eigen_R.R
NULL
#' Hazard of reaching absorbing state 
#'
#' Get h(t|inital dist) in phase type model based on given initial distribution and given time. 
#' May need to be modified for multiple absorbing states!
#'
#' @param time evaluation time
#' @param initial.dist initial state occupancy distribution
#' @param rate.matrix rate matrix
#' @param absorb.state index of absorbing state
#' @return h(t|initial dist) 
#' @export
#' @author Jane Lange
phase.hazard<-function(time, initial.dist,rate.matrix,absorb.state){
############################################################## 
#Author: JL, 11/21/2011
#function the hazard of attaining the absorbing state based on the initial dist for a given time (WRT duration in the state)
#INPUTS:time, initial.dist,rate.matrix,absorb.state
#OUTPUTS:hazard evaluated at time
##############################################

transient.matrix=rate.matrix[-absorb.state,-absorb.state]
exit.vector=rate.matrix[-absorb.state,absorb.state]

matexp=mat.exp(mat=transient.matrix,t=time)
e=rep(1,times=dim(transient.matrix)[1])

out=(initial.dist %*%matexp%*%exit.vector)/(initial.dist %*%matexp%*%e)
return(out)
}

#' @include mat_exp_eigen_R.R
NULL
#' Survivial in phase type model
#'
#' Get S(t|inital dist) in phase type model based on given initial distribution and given time. 
#' May need to be modified for multiple absorbing states!
#'
#' @inheritParams phase.hazard
#' @return S(t|initial dist)
#' @export
#' @author Jane Lange
phase.surv<-function(time, initial.dist,rate.matrix,absorb.state){
############################################################## 
#Author: JL, 11/21/2011
#function the hazard of attaining the absorbing state based on the initial dist for a given time (WRT duration in the state)
#INPUTS:time, initial.dist,rate.matrix,absorb.state
#OUTPUTS:hazard evaluated at time
##############################################
transient.matrix=rate.matrix[-absorb.state,-absorb.state]
exit.vector=rate.matrix[-absorb.state,absorb.state]

matexp=mat.exp(mat=transient.matrix,t=time)
e=rep(1,times=dim(transient.matrix)[1])

out=(initial.dist %*%matexp%*%e)
return(out)
}
