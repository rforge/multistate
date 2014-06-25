NULL
#' @title Latent CTMC hazard functions
#'
#' @description Obtain hazard function based on specified starting states and destination state.
#' 
#'
#' @param start start time for the hazard estimates
#' @param end end time for the hazard estimates
#' @param state_of_interest destination state
#' @param at_risk_states potential states of origin
#' @param alpha starting distribution (vector)
#' @param rate latent CTMC rate matrix to estimate the hazard
#' @param length.out number of points between start and end times to estimate the hazard
#' @return 
#' \item{hazard}{estimated hazard at different times}
#' \item{times}{times at which hazard is estimated}.
#' @examples
#' \dontrun{
#' library(cthmm)
#' data(panel_demodata)
#' 
#' rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
#' rate.setup=rates.setup, do.list=T)
#' rate.firstpassage.AB_0=rates.list[[200]]
#' rate.firstpassage.AB_0[3:4,]=0
#' hazardA.B_0=hazard_times(start=.01,end=10,state_of_interest=c(3),
#' at_risk_states=c(1,2),alpha=c(1,0,0,0),rate.firstpassage.AB_0)
#' 
#' plot(hazardA.B_0$times,hazardA.B_0$haz,type="l",ylab="hazard",xlab="time since entrance in state A",
#'      main="Hazard function for first passage time from state A to B, for covariate X=0")
#' }
#' @author Jane Lange
#' @export
hazard_times<-function(start=0,end,state_of_interest,at_risk_states,alpha,rate,length.out=1000){
	times=seq(start,end,length.out=length.out)
	out=unlist(lapply(times,FUN=get_hazard,state_of_interest=state_of_interest,at_risk_states=at_risk_states,alpha=alpha,rate=rate))
	return(list(hazard=out,times=times))
}

NULL
#' @title Latent CTMC subdistribution functions
#' @description Obtain subdistribution functions or transition probabilities for specified starting and destination states. 
#' 
#'
#' @param start start time for the hazard estimates
#' @param end end time for the hazard estimates
#' @param states destination state(s) 
#' @param alpha starting distribution (vector)
#' @param rate latent CTMC rate matrix to estimate the hazard
#' @param length.out number of points between start and end times to estimate the subdistribution function.
#' @return 
#' \item{dist}{estimated subdistribution at different times}
#' \item{times}{times at which subdistribution function is esitmated}.
#' @examples
#'\dontrun{
#' library(cthmm)
#' data(panel_demodata)
#' 
#' rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
#' rate.setup=rates.setup, do.list=T)
#' rate.firstpassage.AB_0=rates.list[[200]]
#' rate.firstpassage.AB_0[3:4,]=0
#' cdfA.B_0=sub_dist_times(start=.01,end=10,states=3,
#' alpha=c(1,0,0,0),rate.firstpassage.AB_0)
#' plot(cdfA.B_0$times,cdfA.B_0$dist,type="l",ylab="CDF",xlab="time since entrance in state A",
#'      main="Cumulative distribution function for first passage time from state A to B, for covariate X=0")
#' }
#' @author Jane Lange
#' @export
sub_dist_times<-function(start=0,end,states,alpha,rate,length.out=1000){
	times=seq(start,end,length.out=length.out)
	out=unlist(lapply(times,FUN=get_sub_distribution,states=states,alpha=alpha,rate=rate))
	return(list(dist=out,times=times))
}


#############################################################################
NULL
#' @title Subdistribution Standard Errors 
#'
#' @description Get standard errors for the subdistribution function estimates
#' 
#'
#' @param start start time for the hazard estimates
#' @param end end time for the hazard estimates
#' @param covar ovariance matrix for the rate estimates only
#' @param rate latent CTMC rate matrix to estimate the hazard
#' @param states destination state(s)
#' @param alpha starting distribution (vector)
#' @param param.deriv derivative array for the individual
#' @param num.transitions number of transtions between states in the latent CTMC model with non-zero rates
#' @param num.params number of rate parameters in the latent CTMC rate model
#' @param transitions transition codes for the latent CTMC model
#' @param num.states number of states in the latent CTMC model
#' @param length.out number of points between start and end times to estimate the standard error
#' @return standard errors for the subdistribution function 
#' @examples
#' \dontrun{
#' library(cthmm)
#' data(panel_demodata)
#' 
#' rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
#'            rate.setup=rates.setup, do.list=T)
#'
#' rate.firstpassage.AB_0=rates.list[[200]]
#' rate.firstpassage.AB_0[3:4,]=0
#'
#' covar=covariance_4state$covariance[1:8,1:8]
#'
#' num.transitions=dim(rates.setup$transition.codes)[1]
#'
#' num.params=8
#'
#' transitions=rates.setup$transition.codes
#'
#' num.states=4
#'
#' param.deriv=rates.setup$deriv.array[,,1]
#'
#' se_AB_0=se.dist.times(start=.01,end=10,
#'           covar[1:8,1:8],
#'           rate.firstpassage.AB_0,
#'           states=3,alpha=c(1,0,0,0), 
#'           param.deriv,num.transitions,
#'           num.params,transitions,
#'            num.states,length.out=1000)
#' }
#' @author Jane Lange
#' @export
#get the SE of the distribution fun at different times
se.dist.times<-function(start,end,covar,rate,states,alpha, param.deriv,num.transitions,num.params,transitions,num.states,length.out=100){
times=seq(start,end,length.out=length.out)
deriv.list=lapply(times,FUN="dist.fun.deriv",rate,states,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
}

#############################################################################
NULL
#' @title Hazard Standard Errors 
#'
#' @description Get standard errors for the hazard function estimates
#' 
#'
#' @param start start time for the hazard estimates
#' @param end end time for the hazard estimates
#' @param covar covariance matrix for the rate estimates only
#' @param rate latent CTMC rate matrix to estimate the hazard
#' @param at_risk_states potential states of origin
#' @param state_of_interest destination state
#' @param alpha starting distribution (vector)
#' @param param.deriv derivative array for the individual
#' @param num.transitions number of transtions between states in the latent CTMC model with non-zero rates
#' @param num.params number of rate parameters in the latent CTMC rate model
#' @param transitions transition codes for the latent CTMC model
#' @param num.states number of states in the latent CTMC model
#' @param length.out number of points between start and end times to estimate the standard error
#' @return standard errors for the subdistribution function 
#' @examples
#' \dontrun{
#' library(cthmm)
#' data(DDO_data)
#'
#' rates.list=get.rate.matrix.list(current.params=DDO_EM$params[1:7],rate.setup=rates.setup,do.list=T)
#'
#' se=se.haz.times(start=0, end=10,
#' covar=covariance_DDO$covariance[1:7,1:7],
#' rate=rates.list[[1]],
#' at_risk_states=c(1,2),
#' state_of_interest=3,
#' alpha=c(1,0,0,0),
#' param.deriv=rates.setup$deriv.array[,,1],
#' num.transitions=dim(rates.setup$transition.codes)[1],
#' num.params=7,
#' transitions=rates.setup$transition.codes,
#' num.states=4,
#' length.out = 100)
#' }
#' @author Jane Lange
#' @export
se.haz.times<-function(start,end,covar,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states,length.out=100){
times=seq(start,end,length.out=length.out)
deriv.list=lapply(times,FUN="haz.fun.deriv",rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
se=sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
return(se)
}




#density
#parameters: state_of_interest, at_risk_states, alpha, rate,t_
get_density<-function(t_,state_of_interest,at_risk_states,alpha,rate){
  rate_eigen=eigen_decomp_list_R(list(rate))
  P_t=mat_exp_eigen_R(rate_eigen[[1]],t_)
  return(t(alpha)%*%P_t[,at_risk_states]%*%(rate[at_risk_states,state_of_interest]))
}

get_sub_distribution<-function(t_,states,alpha,rate){
   rate_eigen=eigen_decomp_list_R(list(rate))
  P_t=mat_exp_eigen_R(rate_eigen[[1]],t_)
  sum(t(alpha)%*%P_t[,states])
}

get_hazard<-function(t_,state_of_interest,at_risk_states,alpha,rate){
  dens=get_density(t_,state_of_interest,at_risk_states,alpha,rate)
  dist=get_sub_distribution(t_,at_risk_states,alpha,rate)
  return(dens/dist)
}
  
###########################################################
#a list of matrices corresponding to the derivatives of each of the rate parameters
rate.deriv<-function(rate,param.deriv,num.transitions,num.params,transitions,num.states){
 rate.deriv.list<-list()
 for(l in 1:num.params){
  rate.deriv.list[[l]]=matrix(0,nrow=num.states,ncol=num.states)
  for(k in 1:num.transitions){
    rate.deriv.list[[l]][transitions[k,"ni"],transitions[k,"nj"]]=param.deriv[l,k]
     }
  rate.deriv.list[[l]]=rate.deriv.list[[l]]*rate
   diag(rate.deriv.list[[l]])=-1*apply(rate.deriv.list[[l]],1,"sum")
}
  return(rate.deriv.list)
}

#get the derivative of the matrix exponential in the direction of the matrix deriv.mat
mat.exp.deriv<-function(deriv.mat,rate,t_){
  return(uniformization_mean(deriv.mat,t_, rate))
}
  
dist.fun<-function(P,alpha,states){
  sum(t(alpha)%*%P[,states])
}
  
#get a vector of the derivatives of the dist function WRT to each of the B parameters
dist.fun.deriv<-function(t_,rate,states,alpha, param.deriv,num.transitions,num.params,transitions,num.states){
  rate.derivs.list=rate.deriv(rate,param.deriv,num.transitions,num.params,transitions,num.states)
  mat.deriv=lapply(rate.derivs.list,FUN="mat.exp.deriv",rate=rate,t_=t_)
  out=unlist(lapply(mat.deriv,FUN="dist.fun",alpha=alpha,states=states))
  return(out)
}
#the derivative of the density function WRT to a single B parameter
dens.fun.deriv.single=function(rate.deriv,mat.deriv,rate,P_t,at_risk_states,state_of_interest,alpha){
out=t(alpha)%*%(P_t[,at_risk_states]%*%rate.deriv[at_risk_states,state_of_interest]+mat.deriv[,at_risk_states]%*%rate[at_risk_states,state_of_interest])
return(as.numeric(out))
}
#get a vector of the derivatives of the density function WRT to each of the B parameters
dens.fun.deriv<-function(t_,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states){
  rate.derivs.list=rate.deriv(rate,param.deriv,num.transitions,num.params,transitions,num.states)
  mat.derivs.list=lapply(rate.derivs.list,FUN="mat.exp.deriv",rate=rate,t_=t_)
  rate_eigen=eigen_decomp_list_R(list(rate))
  P_t=mat_exp_eigen_R(rate_eigen[[1]],t_)
  out=mapply(rate.derivs.list,mat.derivs.list,FUN="dens.fun.deriv.single",MoreArgs = list(rate=rate,P_t=P_t,at_risk_states=at_risk_states,state_of_interest=state_of_interest,alpha=alpha))
  return(out)
}
haz.fun.deriv.single<-function(dist.deriv,dens.deriv,dist,dens){
 (dist*dens.deriv-dens*dist.deriv)/(dist^2)
}
haz.fun.deriv<-function(t_,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states){
  dens.derivs=dens.fun.deriv(t_,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
  dist.derivs=dist.fun.deriv(t_,rate,states=at_risk_states,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
  dist=get_sub_distribution(t_,at_risk_states,alpha,rate)
  dens=get_density(t_,state_of_interest,at_risk_states,alpha,rate)
  out=mapply(dist.derivs,dens.derivs,FUN="haz.fun.deriv.single",MoreArgs=list(dist=dist,dens=dens))
  return(out)
}

delta.method.mat=function(x,covar){
  return(x%*%covar%*%x)
}



se.dens.times<-function(start,end,covar,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states){
 times=seq(start,end,length.out=100)
 deriv.list=lapply(times,FUN="dens.fun.deriv",rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
 se=sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
 return(se)
}
  
##############################################################################
