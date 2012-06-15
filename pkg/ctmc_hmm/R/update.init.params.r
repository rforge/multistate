#' @include covariates_max.R
NULL
#' Performs Newton Raphson to get MLE for initial dist parameters.
#'
#' Important: currently only runs NR max.it times (not until convergence).
#' Thus this function can be used to obtain the initial distributions based
#'  on the provided values of current.params, if max.it=1 and init.counts=NULL.
#'
#' @inheritParams update.emission.params
#' @param init.setup initial distributions setup object: must have covariate array, param.types
#' @param init.counts matrix (dim(num.hidden.states,num.subjects)) with the initial dist count data               
#' @return list with 
#'  \item{params}{an array of dim number of parameters, max.it} 
#'  \item{prob.matrix}{array with each column representing the vector of initial probabilities for a single subject, based on the last value of the parameters}
#'  \item{states}{hidden states corresponding to rows of prob.matrix: will not include ref states or states with zero prob!}
#'  \item{init.counts}{array with (dim=num.states,num.subjects) of the initial distribution for each individual.}
#' @export
#' @author Jane Lange
update.init.params<-function(current.params,init.setup,init.counts=NULL,max.it=1){
  no.NR=0
  init=init.setup
  if(!is.null(init.counts)){
  counts=matrix(init.counts[init$states[init$states!=init$ref],],ncol=dim(init.counts)[2])
 
  }else{
  counts=NULL
  no.NR=1
  }
  N_vector=rep(1,times=dim(init$covariate.array)[3])
  update_init=update.multinom.params(current.params=current.params,
                                     covariate.array=init$covariate.array,
                                     deriv.array=init$deriv.array,
                                     param.types=init$param.types,
                                      counts=counts,
                                     N_vector=N_vector,
                                     max.it=max.it,
                                     no.NR=no.NR)
 update_init$states=init$states[init$states!=init$ref]
  return(update_init)

}
