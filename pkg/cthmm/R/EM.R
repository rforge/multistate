#' @include EM.R
NULL
#' @title EM
#'
#' @description EM algorithm to get MLEs for latent CTMC parameters 
#'
#' @details Run the expectation-maximization algorithm to get the maximum likelihood
#' TEST
#' estimates for latent CTMC model parameters for multistate disease models. Data are discretely observed disease
#' processes from multiple independent individuals. Times may be non-informative
#' scheduled times or disease-dependent (informative) times. Uses the SQUAREM package for EM acceleration. 
#'
#' @param rates.setup list with rate setup information
#' @param init.setup list with initial distribution setup information
#' @param emission.setup list with emission distribution setup information
#' @param the.data list with the observed data, with one entry per individual 
#' @param num.subjects number of individuals in the study
#' @param num.states number of latent states in the CTMC
#' @param num.obs.states number of observed states
#' @param tol convergence tolerance for the EM
#' @param absorb.state vector of one or more absorbing latent states
#' @param maxiter=500 maximum number of iterations of the EM
#' @param DDO.setup setup object for disease driven observation model
#' @param do.DDO  Indicator (T/F) if there are disease driven observation times in the model?
#' @param time.dep.emission (Not supported) Indicator if emission distribution has time dependent covariates.
#' @return 
#' \item{params}{MLEs for model parameters: rates/emission/initial/DDO, parameters listed in order of model specification}
#' \item{details}{Details of the EM optimization from the \code{SQUAREM} package output}.
#' \item{updates}{The iterative oupdates for the parmaeter esitmates at each EM step}
#' \item{LL.updates}{Updates of observed data likelihood at each step}
#' \item{time}{Run time}
#' \item{num.evals}{Number of iterations of the EM} 
#' \item{LL}{Final complete data log likelihood}
#' @export
#' @examples
#'\dontrun{
#' library(cthmm)
#'#load the model setup and example data for a competing risks model
#'#observed at scheduled and informative observation times.
#'data(DDO_data)
#' #run the EM on the example data
#'DDO_EM=EM(rates.setup=rates.setup,
#'init.setup=init.setup,
#'emission.setup=emission.setup,
#'the.data=DDO_data, 
#'num.subjects=500,
#'num.states=4,
#'num.obs.states=3, 
#'tol = 1e-07, 
#'absorb.state=c(3,4), 
#'maxiter = 500, 
#'DDO.setup = DDO.setup,
#'do.DDO = T)}
#' @author Jane Lange
#################################################
EM<-function(rates.setup,init.setup,emission.setup,the.data,num.subjects,num.states,
             num.obs.states,tol=1e-7,absorb.state,maxiter=500,DDO.setup=NULL,do.DDO=F,time.dep.emission=F){
  #################################################
   #################################################

  Q=NULL
  h.list=NULL
  rates<-rates.setup
  init<-init.setup 
  emission<-emission.setup
  DDO<-DDO.setup
  ij.indices=matrix(unlist(rates$transition.codes),ncol=2,byrow=F)
  colnames(ij.indices)=c("i","j")
  ij.indices<-ij.indices
  obs.data.list<-lapply(the.data, "[[",c("obs.data"))
  exact.time.ranks.list<-lapply(the.data,"[[",c("exact.times.ranks"))
  time.diffs.list<-lapply(lapply(the.data,"[[","obs.times"),FUN="diff")

  if(do.DDO){
    h.list=lapply(the.data,"[[",c("h"))
  }
  rep_item<-function(i,item){
    return(item)
  }

#par=c(rates$param.values)
  par=c(0)

  if(!is.null(rates$fixed.rates)){
    rates.list=lapply(seq(1:num.subjects), FUN="rep_item",item=rates$fixed.rates)
    num.rate.params <-0
    do.rates=F
  }else{
    par=c(par,rates$param.values)
    num.rate.params <-length(rates$param.values)
    do.rates=T
  }
 if(!is.null(emission$fixed.dist)){
    emission.list=lapply(seq(1:num.subjects),FUN="rep_item",item=emission$fixed.dist)
    num.emission.params<-0
    do.emission=F
 }else{
    par=c(par,emission$param.values)
    num.emission.params=length(emission$param.values)
    do.emission=T
 }
 if(!is.null(init$fixed.dist)){
    delta.list=lapply(seq(1:num.subjects),FUN="rep_item",item=init$fixed.dist)
    num.init.params=0
    do.init=F
 }else{
    par=c(par,init$param.values)
    num.init.params<-length(init$param.values)
    do.init=T
 }
 if(do.DDO){
    if(!is.null(DDO$fixed.rates)){
      Q=lapply(seq(1:num.subjects),FUN="rep_item",item=DDO$fixed.rates)
      fixed.DDO=T
      num.DDO.params=0
    }else{
      par=c(par,DDO$param.values)
      num.DDO.params<-length(DDO$param.values)
      fixed.DDO=F
    }
  }

  par=par[-1]
  param.updates=par
  LL.updates=0
  start=proc.time()

EM.update<-function(par){

 if(do.rates){
   rates.list=get.rate.matrix.list(current.params=par[1:num.rate.params],rate.setup=rates,do.list=T)
 }
 if(do.emission){
   emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 1):(num.rate.params +
   num.emission.params)], emission.setup = emission, 
   num.states = num.states, num.obs.states = num.obs.states, 
   num.subjects = num.subjects,time.dep.emission=time.dep.emission)      
}
if(do.init){
   delta.list = get.init.matrix.list(current.params = par[(num.rate.params +num.emission.params + 1):(num.rate.params + num.emission.params + 
										num.init.params)], init.setup = init)
 }
#Change rates to Lambda-Q if do.DDO=T
 if(do.DDO){
   if(!fixed.DDO){
    Q=get.rate.matrix.list(current.params=par[(num.rate.params + 
										num.emission.params + num.init.params +1):(num.rate.params 
										+ num.emission.params + num.init.params+num.DDO.params)],
						                rate.setup=DDO, do.list = T)    
   }
    rates.list=mapply("+",rates.list,Q,SIMPLIFY=F)
    
 }

 eigen.decomp.list<-eigen_decomp_list_R(rates.list) 

 transition.probabilities.list=trans_probs_all(eigen.decomp.list,time.diffs.list,exact.time.ranks.list,h.list,Q)
 likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilities.list,delta.list,emission.list,FUN="forwardback_R",MoreArgs=list(time.dep=time.dep.emission),SIMPLIFY=F)

 expected.trans<-trans_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list, exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state,ij.indices,h.list,Q,time_dep_emission=time.dep.emission)
 expected.dur<-dur_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list,exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state,h.list,Q,time_dep_emission=time.dep.emission)
 expected.trans[is.na(expected.trans)]=0
 expected.dur[is.na(expected.dur)]=0
 expected.trans[is.nan(expected.trans)]=0
 expected.dur[is.nan(expected.dur)]=0
 
 if (do.emission|do.init) {
   if(time.dep.emission==F){
     E.emission.init = E_step_emission_init_all(likelihood.forward.backward.list, 
     obs.data.list, num.states, num.obs.states, num.subjects)
   }else{
     E.emission.init=E_step_emission_init_all_time_dep(likelihood.forward.backward.list, 
     obs.data.list, num.states, num.obs.states, num.subjects,
     total.obs.times=dim(emission$deriv.array)[3],limits=emission$limits)
}
  expected.emission = E.emission.init$emission_array
  expected.init = E.emission.init$init_array
  expected.emission[is.nan(expected.emission)]=0
  expected.init[is.nan(expected.init)]=0
} 

if (do.DDO) {
   expected.DDO =E_step_DDO(likelihood.forward.backward.list, h.list,num.states,num.subjects)
   expected.DDO[is.nan(expected.DDO)]=0
}

params=c(0)
  if(do.rates){
    if(sum(rates$param.types>0)==length(rates$param.values)){
      new=rates$param.values
    }else{
      update_rates = update.rate.params(rate.setup = rates, 
                                    transitions = expected.trans, durations = expected.dur, 
                                    current.params = par[1:num.rate.params], max.it = 10)
      new = update_rates$param[, dim(update_rates$param)[2]]
    }
	names(new) = rep("r", times = length(new))
	if(min(new)<(-16)){
    new[new<(-16)]=-16
    } 
   params = c(params,new)
  
}
 if(do.emission){
   if(sum(emission$param.types>0)==length(emission$param.values)){
     new=emission$param.values
   }else{
 
     update_emission = update.emission.params(current.params = par[(num.rate.params + 
													1):(num.rate.params + num.emission.params)], 
                                           emission.setup = emission, emission.counts = expected.emission, max.it = 10)
     new = update_emission$param[, dim(update_emission$param)[2]]
  }
  names(new) = rep("e", times = length(new))
  params = c(params, new)
 }
 if (do.init) {
  
   if(sum(init$param.types>0)==length(init$param.values)){
    new=init$param.values
   }else{
    update_init = update.init.params(current.params = par[(num.rate.params + 
							num.emission.params + 1):(num.rate.params + num.emission.params + num.init.params)], 
							init.setup = init, init.counts = expected.init, max.it = 2)
    new = update_init$param[, dim(update_init$param)[2]]
   }
  
  
   names(new) = rep("i", times = length(new))
   params = c(params, new)
 }
 if(do.DDO){
   if(is.null(DDO$fixed.rates)){
    if(sum(DDO$param.types>0)==length(DDO$param.values)){
      new=DDO$param.values
    }else{  
      update_DDO = update.rate.params(rate.setup = DDO, 
                                  transitions = expected.DDO, durations = expected.dur, 
                                  current.params = par[(num.rate.params +num.emission.params 
                                                        + num.init.params +1):(num.rate.params 
								+ num.emission.params +num.init.params+num.DDO.params)], max.it = 10)
  
	  new = update_DDO$param[, dim(update_DDO$param)[2]]
    }
  names(new) = rep("d", times = length(new))
  if(min(new)<(-50)){
    new[new<(-50)]=-50
  }

  params = c(params,new)
  }
}

params=params[-1]
LL.vec=(sapply(likelihood.forward.backward.list,"[[","LL"))
LL=sum(LL.vec[LL.vec!=-Inf&LL.vec!="NaN"])
LL.updates<<-c(LL.updates,LL)
print(LL)

print(params)
param.updates<<-rbind(param.updates,params)

return(params)

}

  
objfn<-function(par){
    
    rep_item<-function(i,item){
      return(item)
    }  
    if(!is.null(rates.setup$fixed.rates)){
      rates.list=lapply(seq(1:num.subjects), FUN="rep_item",item=rates.setup$fixed.rates)
    }else{
      rates.list=get.rate.matrix.list(current.params=par[1:num.rate.params],rate.setup=rates)
    }
    if(!is.null(emission.setup$fixed.dist)){
      emission.list=lapply(seq(1:num.subjects),FUN="rep_item",item=emission.setup$fixed.dist)
    }else{
		emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 1):(num.rate.params +
		num.emission.params)], emission.setup = emission, 
		num.states = num.states, num.obs.states = num.obs.states, 
		num.subjects = num.subjects,time.dep.emission=time.dep.emission)      
    }
    if(!is.null(init.setup$fixed.dist)){
      delta.list=lapply(seq(1:num.subjects),FUN="rep_item",item=init.setup$fixed.dist)
    }else{
      delta.list = get.init.matrix.list(current.params = par[(num.rate.params + 
                                                                num.emission.params + 1):(num.rate.params + num.emission.params + 
                                                                                            num.init.params)], init.setup = init)
    }
    if(do.DDO){
	
	 if(!is.null(DDO$fixed.rates)){
	     Q=lapply(seq(1:num.subjects),FUN="rep_item",item=DDO$fixed.rates)
	 }else{		
        Q=get.rate.matrix.list(current.params=par[(num.rate.params + 
						num.emission.params + num.init.params +1):(num.rate.params + num.emission.params + 
                                                                                                  num.init.params+num.DDO.params)], rate.setup=DDO, do.list = T)    
     } 
     rates.list=mapply("+",rates.list,Q,SIMPLIFY=F) 
    }
    
    eigen.decomp.list<-eigen_decomp_list_R(rates.list)
    transition.probabilities.list=trans_probs_all(eigen.decomp.list,time.diffs.list,exact.time.ranks.list,h.list,Q)
	likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilities.list,delta.list,emission.list,FUN="forwardback_R",MoreArgs=list(time.dep=time.dep.emission),SIMPLIFY=F)
    LL=sum(sapply(likelihood.forward.backward.list,"[[","LL"))
    return(LL)
  }
  
  #out <- fpiter(par=par, fixptfn=EM.update, control=list(tol=tol), objfn=objfn)
  #out <- squarem(par, fixptfn=EM.update,objfn=objfn,control=list(tol=tol,step.min0=1,maxiter=2500))
  out <- squarem(par, fixptfn=EM.update,control=list(tol=tol,step.min0=1,maxiter=maxiter))
  
  par=out$par
  rate=get.rate.matrix.list(current.params=par[1:num.rate.params],rate.setup=rates)[[1]]
  
  the.time=proc.time()-start
  
  return(list(params=par,rate=rate,details=out,updates=param.updates,LL.updates=LL.updates,time=the.time,num.evals=length(LL.updates),LL=max(LL.updates[-1])))
  
}