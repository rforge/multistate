#################################################
#This functions runs the EM for parameterized rates, emission, and initial distributions
# using the SQUAREM package for acceleration of updates
#Author: JL 11-24-2011
#Inputs: rates.setup, init.setup, emission.setup, the.date,num.subjects,num.states, num.obs.states,tol
#Outputs: List consisting of params=EM estimates for parameters and LL=final observed data log likelihood
################################################
#' @include EM.r
NULL
#' EM algorithm for HMM based on CTMC for parameterized rates, emission, and initial distributions
#' Using the SQUAREM package for accleration of convergence.
#' Get the MLEs for the parameters in the model (rate, initial distribution, and emission)
#' This implementation does not use the recursive filtering method, but rather uses forward and backward probabilities.
#'
#' @param rates.setup list rate setup information
#' @param init.setup list with initial distribution setup information
#' @param emission.setup list with emission distribution setup information
#' @param the.data list with the data 
#' @param num.subjects number of subjects in the study
#' @param num.states number of states in the CTMC
#' @param num.obs.states number of observed states
#' @param tol convergence tolerance
#' @param absorb.state vector of one or more absorbing states
#' @param maxiter=500
#' @return a list of estimates from the EM:
#' \item{params}{MLEs for model parameters}
#' \item{param.updates}{The iterative oupdates for the parmaeter esitmates at each EM step}
#' \item{LL.updates}{Updates of observed data likelihood at each step}
#' \item{time}{Run time}
#' \item{LL}{Final complete data log likelihood}
#' @export
#' @author Jane Lange


#################################################
EM<-function(rates.setup,init.setup,emission.setup,the.data,num.subjects,num.states,
num.obs.states,tol=1e-7,absorb.state,maxiter=500){
#################################################
#This functions runs the EM for parameterized rates, emission, and initial distributions
# using the SQUAREM package for acceleration of updates
#Author: JL 11-24-2011
#Inputs: rates.setup, init.setup, emission.setup, the.date,num.subjects,num.states, num.obs.states,tol
#Outputs: List consisting of params=EM estimates for parameters and LL=final observed data log likelihood
#################################################
print("new package")
rates<-rates.setup
init<-init.setup
emission<-emission.setup
num.states<-num.states
num.obs.states<-num.obs.states

ij.indices=matrix(unlist(rates$transition.codes),ncol=2,byrow=F)
colnames(ij.indices)=c("i","j")
ij.indices<-ij.indices
obs.data.list<-lapply(the.data, "[[",c("obs.data"))
exact.time.ranks.list<-lapply(the.data,"[[",c("exact.times.ranks"))
time.diffs.list<-lapply(lapply(the.data,"[[","obs.times"),FUN="diff")

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

par=par[-1]
param.updates=par
LL.updates=0
start=proc.time()

EM.update<-function(par){
  
   if(do.rates){
     rates.list=get.rate.matrix.list(current.params=par[1:num.rate.params],rate.setup=rates)
   }
   if(do.emission){
           emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 
            1):(num.rate.params + num.emission.params)], emission.setup = emission, 
            num.states = num.states, num.obs.states = num.obs.states, 
            num.subjects = num.subjects)      
  }
  if(do.init){
    delta.list = get.init.matrix.list(current.params = par[(num.rate.params + 
          num.emission.params + 1):(num.rate.params + num.emission.params + 
          num.init.params)], init.setup = init)
  }
  
  eigen.decomp.list<-eigen_decomp_list_R(rates.list) 
  transition.probabilities.list=trans_probs_all(eigen.decomp.list,time.diffs.list,exact.time.ranks.list)
  likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilities.list,delta.list,emission.list,FUN="forwardback_R",SIMPLIFY=F)

  expected.trans<-trans_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list, exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state,ij.indices)
  expected.dur<-dur_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list,exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state)
  expected.trans[is.na(expected.trans)]=0
  expected.dur[is.na(expected.dur)]=0
   
  if (do.emission|do.init) {
            E.emission.init = E_step_emission_init_all(likelihood.forward.backward.list, 
                obs.data.list, num.states, num.obs.states, num.subjects)
            expected.emission = E.emission.init$emission_array
            expected.init = E.emission.init$init_array
  } 

   params=c(0)
 if(do.rates){
#  print(par)
   update_rates = update.rate.params(rate.setup = rates, 
                transitions = expected.trans, durations = expected.dur, 
                current.params = par[1:num.rate.params], max.it = 2)
                new = update_rates$param[, dim(update_rates$param)[2]]
            names(new) = rep("r", times = length(new))
             if(min(new)<(-50)){
                new[new<(-50)]=-50
              }
 
            params = c(params,new)
            
 }
if(do.emission){
              update_emission = update.emission.params(current.params = par[(num.rate.params + 
                1):(num.rate.params + num.emission.params)], 
                emission.setup = emission, emission.counts = expected.emission, 
                max.it = 2)
            new = update_emission$param[, dim(update_emission$param)[2]]
            names(new) = rep("e", times = length(new))
            params = c(params, new)
}
if (do.init) {
    
            update_init = update.init.params(current.params = par[(num.rate.params + 
                num.emission.params + 1):(num.rate.params + num.emission.params + 
                num.init.params)], init.setup = init, init.counts = expected.init, 
                max.it = 2)
            new = update_init$param[, dim(update_init$param)[2]]
            names(new) = rep("i", times = length(new))
            params = c(params, new)
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

# objfn<-function(par){
#   rates.list=get.rate.matrix.list(current.params=par,rate.setup=rates)
#   eigen.decomp.list<-eigen_decomp_list_R(rates.list)
#   transition.probabilities.list<-mapply(time.diffs.list,eigen.decomp.list, exact.time.ranks.list,FUN=transition_prob_list,SIMPLIFY=F)
#   likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilities.list,delta.list,emission.list,FUN="forwardback_R",SIMPLIFY=F)
#   LL.vec=(sapply(likelihood.forward.backward.list,"[[","LL"))
#   LL=sum(LL.vec[LL.vec!=-Inf])
#  
#   return(LL)
# }

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
      emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 
            1):(num.rate.params + num.emission.params)], emission.setup = emission, 
            num.states = num.states, num.obs.states = num.obs.states, 
            num.subjects = num.subjects)
  }
  if(!is.null(init.setup$fixed.dist)){
    delta.list=lapply(seq(1:num.subjects),FUN="rep_item",item=init.setup$fixed.dist)
  }else{
    delta.list = get.init.matrix.list(current.params = par[(num.rate.params + 
          num.emission.params + 1):(num.rate.params + num.emission.params + 
          num.init.params)], init.setup = init)
  }
 
  eigen.decomp.list<-eigen_decomp_list_R(rates.list)
  transition.probabilities.list<-mapply(time.diffs.list,eigen.decomp.list, exact.time.ranks.list,FUN=transition_prob_list,SIMPLIFY=F)
  likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilities.list,delta.list,emission.list,FUN="forwardback_R",SIMPLIFY=F)
  LL=sum(sapply(likelihood.forward.backward.list,"[[","LL"))
  return(LL)
}
  
#out <- fpiter(par=par, fixptfn=EM.update, control=list(tol=tol),objfn=objfn)
#out <- squarem(par, fixptfn=EM.update,objfn=objfn,control=list(tol=tol,step.min0=1,maxiter=2500))
out <- squarem(par, fixptfn=EM.update,control=list(tol=tol,step.min0=1,maxiter=maxiter))

par=out$par
rate=get.rate.matrix.list(current.params=par[1:num.rate.params],rate.setup=rates)[[1]]
 
the.time=proc.time()-start

return(list(params=par,rate=rate,details=out,updates=param.updates,LL.updates=LL.updates,time=the.time,num.evals=length(LL.updates),LL=max(LL.updates[-1])))
}

 



