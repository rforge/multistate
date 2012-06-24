#' @include covariates_max.R
NULL
#' Performs Newton Raphson to get MLE for emission parameters.
#'
#' Important: currently only runs NR max.it times (not until convergence). 
#' Thus this function can be used to obtain the emission probabilties
#' based on the provided values of current.params, if max.it=1
#' and emission.counts=NULL.
#'
#' @param current.params initial values of the parameters
#' @param emission.setup object with covariate array, param.types
#' @param emission.counts matrix (dim(num.hidden.states,num.obs.states,num.subjects)) with the emission count data
#' @param max.it maximum number of Newton Raphson iterations. If max.it=1, no updates are performed.
#               
#' @return list with 
#'  \item{params}{an array of dim number of parameters, max.it} 
#'  \item{prob.matrix}{array with each column representing the vector of calculation emission probabilities for a single subject, based on the last value of the parameters}
#'  \item{states}{matrix columns "i" and "j" representing hidden(i) and observed(j) states corresponding to the rows of prob.matrix}
#' @export
#' @author Jane Lange 
update.emission.params<-function(current.params,emission.setup, emission.counts=NULL,max.it=1){
emission=emission.setup
  no.NR=0
  if(is.null(emission.counts)){
  no.NR=1
  }
  #number of parameters for NR update
 num.NR.params=length(emission$param.types[emission$param.types==0])
 #deriv array for NR update

#CHANGE THIS
# deriv.array=array(emission$covariate.array[emission$param.types==0,,],dim=c(sum(emission$param.types==0),dim(emission$covariate.array)[2],dim(emission$covariate.array)[3]))
 deriv.array=emission.setup$deriv.array
score=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params[,1]=current.params[emission$param.types==0]

 score=array(0,dim=c(length(current.params),(max.it+1)))
 params=array(0,dim=c(length(current.params),(max.it+1)))
 params[,1]=current.params
 
for(i in 1:(max.it)){
   
  out=emission.score.hessian(emission.counts,params[,i],emission,deriv.array,no.NR=no.NR)
  if(!no.NR){
  score[,i]=out$score
  hessian=out$hessian
  NR.params[,i+1]=NR_update(current.params=NR.params[,i],score[,i], hessian)
  params[,i+1]=current.params
  params[emission$param.types==0,i+1]<-NR.params[,i+1]
  }
}
return(list(params=params,states=out$states,prob.matrix=out$prob.array))
}
