#' @include covariates_max.R
NULL
#' Performs Newton Raphson to get MLE for rate parameters.
#'
#' Important: currently only runs NR max.it times (not until convergence).
#' Thus this function can be used to obtain the rate matrices 
#' based on the provided values of current.params, if transitions=NULL, durations=NULL, and max.it=1.
#'
#' @inheritParams update.emission.params
#' @param rate.setup object: must have covariate array, param.types
#' @param transitions matrix of dim(num.states,num.states,num.subjects) with all i,j transtions for all subjects
#' @param durations matrix of dim(num,states, num.subjects) with the duration spent in each state
#'               
#' @return list with 
#'  \item{params}{an array of dim number of parameters, max.it} 
#'  \item{rates}{array with each column representing the vector of calculated rates for a single subject, based on the last value of the parameters}
#'  \item{transition.codes}{matrix columns "ni" and "nj" i,j entries in the rate matrix}
#' @export
#' @author Jane Lange 
update.rate.params<-function(current.params, rate.setup, transitions=NULL,durations=NULL,max.it){
##################################################
#Author: JL, 11/15/2011
#gets updates of rate parameters based on NR algorithm (run max.it, not to convergence!)
#INPUTS:    
#           current.params=initial values for parameters
#           rate.setup=rate object with:
#              rate.setup$covariate.array=array of dim (num.betas,num.thetas, num.subjects) that contains the deriv. of theta WRT beta
#              rate.setup$param.types = vector for each parameter coding NR (0), simple update (1), or fixed (2)
#              rate.setup$transition.codes= matrix with ni, and nj columns encoding the transitions 
#                             corresponding to each rate
#           transitions=matrix of dim(num.states,num.states,num.subjects) with all i,j transtions for all subjects
#           durations= matrix of dim(num,states, num.subjects) with the duration spent in each state
#           max.it=maximum number of NR updates
#OUTPUTS: a list with params: NR updates of MLEs and rates: array with rate matrices calculated for each individualx
###################################################
covariate.array=rate.setup$covariate.array
param.types=rate.setup$param.types
transition.codes=rate.setup$transition.codes
init.params=current.params
no.NR=0
if(is.null(transitions)){
no.NR=1
  }

#number of parameters for NR update
 num.NR.params=length(param.types[param.types==0])
 #deriv array for NR update

 #CHANGE
 # deriv.array=array(covariate.array[param.types==0,,],dim=c(sum(param.types==0),dim(covariate.array)[2],dim(covariate.array)[3]))
 deriv.array=rate.setup$deriv.array
 score=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params[,1]=init.params[param.types==0]

 
 #updates for rate parameters that do not require NR
 params=array(0,dim=c(length(init.params),(max.it+1)))
 params[,1]=init.params
 
 for(i in 1:(max.it)){
  rates=get.rates(param.values=params[,i],covariate.array=covariate.array)
  if(!no.NR){
  score[,i]=rate_score(rates,durations, deriv.array, transitions, transition.codes)
  hessian=rate_hessian(rates,durations, deriv.array, transition.codes)
  NR.params[,i+1]=NR_update(current.params=NR.params[,i],score[,i], hessian)
  params[,i+1]=init.params
  params[param.types==0,i+1]<-NR.params[,i+1]
  } 
  }
return(list(params=params,rates=rates,transition.codes=transition.codes))
}

