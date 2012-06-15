##################################################################
#Functions for the parameterized m-step
#################################################################
get.rates<-function(param.values, covariate.array){
#############################################
#Author: JL, 11/15/2011
#gets the linear predictors corresponding to exp(theta)=exp(beta^TZ)
#INPUTS:  covariate.array=array of dim(num.beta, num.theta, num.subjects)
#                          that for each subject and each transition contains the values of the covariates
#                 param.values =vector with parameter values
#OUTPUTS:  matrix with rates for each subject dim(num.transitions.types, num.subjects)   
##############################################################       
linear.pred=apply(covariate.array*param.values,c(2,3),"sum")

return(exp(linear.pred))
}

rate_score<-function(rates,durations, deriv.array, transitions, transition.codes){
##############################################
#Author: JL, 11/15/2011
#Function to get the  rate matrix score function
#INPUTS:    rates=matrix with rates for each transition and each subject (dim(num.transitions, num.subjects))
#           transition.codes= matrix with ni, and nj columns encoding the transitions 
#                             corresponding to each rate
#           transition= array of dim(num.states,num.states,num.subjects) with the transitions for each subject
#           deriv.array = array of dim (num.betas,num.thetas, num.subjects) that contains the deriv. of theta WRT beta
#           duration= matrix of dim(num,states, num.subjects) with the duration spent in each state
##############################################
  suff.stat.score=array(NA,dim=c(dim(transition.codes)[1],dim(durations)[2]))
  score=vector(length=c(dim(deriv.array)[1]))
  for(l in 1:dim(transition.codes)[1]){
   suff.stat.score[l,]=transitions[transition.codes[l,"ni"],transition.codes[l,"nj"],]-rates[l,]*durations[transition.codes[l,"ni"],]
  }
  
  for(m in 1:dim(durations)[2]){
  score=score+deriv.array[,,m]%*%suff.stat.score[,m]
   }
  return(score) 
}
##############################################
rate_hessian<-function(rates,durations, deriv.array, transition.codes){
##############################################
#Author: JL, 11/15/2011
#INPUTS:    rates=matrix with rates for each transition and each subject (dim(num.transitions, num.subjects))
#           deriv.array = array of dim (num.betas,num.thetas, num.subjects) that contains the deriv. of theta WRT beta
#           transition.codes= matrix with ni, and nj columns encoding the transitions 
#                             corresponding to each rate
#           durations= matrix of dim(num,states, num.subjects) with the duration spent in each state
#OUTPUTS: hessian matrix for rate model parameters
##############################################
hessian=array(0,dim=c(dim(deriv.array)[1],dim(deriv.array)[1]))
if(dim(deriv.array)[1]==0){
  return(hessian)
}else{
for(m in 1:dim(durations)[2]){
  hessian=hessian+deriv.array[,,m]%*%(t(deriv.array[,,m])*((-1)*rates[,m]*durations[transition.codes[,"ni"],m]))
} 
return(hessian)
}
}

NR_update<-function(current.params,score, hessian){
##############################################
#Author: JL, 11/15/2011
#INPUTS:    score= score vector
#           hessian= hessian matrix
#           current.params =vector of current values of the parameters
#OUTPUTS: next newton raphson update for parameters (vector of length num.parameters)
##############################################
return(current.params-solve(hessian)%*%score)
}
 
get.prob.array<-function(covariate.array, param.values){
#############################################
#get array of multinomial probabilities (not including the reference state) for all subects
##############################################
#Author: JL, 11/15/2011
#INPUTS:    
#           covariate_array: array of covariates for all subjects of dim(num.params,num.states,num.subjects) 
#           param.values: vector of values for model parameters (betas)
#OUTPUTS:   array of  dim(num.states,num.subects)
   linear.pred=apply(covariate.array*param.values,c(2,3),"sum")
   exp.linear.pred=exp(linear.pred)
   return(t(t(exp.linear.pred)*1/(1+apply(exp.linear.pred,2,"sum"))))
}

multinom.score<-function(prob.matrix,deriv.array,counts,N_vector){
#############################################
#multinomial matrix score
##############################################
#Author: JL, 11/15/2011
#INPUTS:    prob.matrix: matrix of dim(num.states, num.subjects) with probabilities calculated
#                           for each state/subject
#           deriv.array = array of dim (num.betas,num.thetas, num.subjects) that contains the deriv. of theta WRT beta
#           counts: array of category counts dim(num.states, num.sub) (excluding reference category)
#           N_vector: vector (length =num.subjects) of multinomial sample size n for each subject
#OUTPUTS:   score vector for multinomial parameters
       num.subjects=dim(counts)[2]
       suff.stat.score=array(0,dim=c(dim(counts)[1],dim(counts)[2]))  
       param.score=rep(0,times=dim(deriv.array)[1])
       suff.stat.score=(counts-prob.matrix*N_vector)
       for(m in 1:num.subjects){
        param.score=param.score+matrix(deriv.array[,,m])%*%suff.stat.score[,m]
         }
    return(param.score)
}
multinom.covariance.array<-function(prob.matrix,N_vector){
##############################################################
#multinomial covariance array
############################################################
#Author: JL, 11/15/2011
#INPUTS:     prob.matrix: matrix of dim(num.states, num.subjects) with probabilities calculated
#                           for each state/subject
#           N_vector: vector (length =num.subjects) of multinomial sample size n for each subject
#OUTPUT     array (nums.state,num.states,num.subjects) of multinomial variance-covariance matrices
out=array(0,dim=c(dim(prob.matrix)[1],dim(prob.matrix)[1],dim(prob.matrix)[2]))
for(m in 1:dim(prob.matrix)[2]){
  
  if(dim(out)[1]==1){
  out[,,m]=prob.matrix[,m]*(1-prob.matrix[,m])
  }else{
  out[,,m]=(-1)*prob.matrix[,m]%*%t(prob.matrix[,m])
  diag(out[,,m])=prob.matrix[,m]*(1-prob.matrix[,m])
  }
  out[,,m]=out[,,m]*N_vector[m]
 
}
 return(out)
}
multinom.hessian<-function(prob.matrix, deriv.array,N_vector){
#############################################
#multinomial matrix hessian
##############################################
#Author: JL, 11/15/2011
#INPUTS:    prob.matrix: matrix of dim(num.states, num.subjects) with probabilities calculated
#                           for each subject
#           deriv.array = array of dim (num.betas,num.thetas, num.subjects) 
#                         that contains the deriv. of theta WRT beta
#           N_vector: vector (length =num.subjects) of multinomial sample size n for each subject
#OUTPUTS:   array of hessian matrix for multinomial parameters (dim num.param, num.params, num.subjects) 
#covariance_array=array(0,dim=c(dim(prob.matrix)[1],dim(prob.matrix)[1],dim(prob.matrix)[2]))
#covariance_array=multinom.covariance.array(prob.matrix,N_vector)
#hessian=array(0,dim=c(dim(deriv.array)[1],dim(deriv.array)[1]))
#for(m in 1:dim(deriv.array)[3]){
#  hessian=hessian+deriv.array[,,m]%*%covariance_array[,,m]%*%t(deriv.array[,,m])
#}
#return(-1*hessian)
myDims=dim(deriv.array)
out=.Call("multinom_hessian",prob.matrix,deriv.array,myDims,N_vector)
return(out)

}

emission.score.hessian<-function(emission_updates,param.values,emission,deriv.array,no.NR=0){
#######################################################################
#emission score and hessian updates
##############################################
#Author: JL, 11/20/2011
#INPUTS: 
#        emission_updates:  matrix (dim(num.hidden.states,num.obs.states,num.subjects)) with the emission count data
#        param.values: vector of values for model parameters
#        emission: emission object with covariate array, param.types
#           deriv.array = array of dim (num.betas,num.thetas, num.subjects) 
#                         that contains the deriv. of theta WRT beta
#        no.NR (indicator =1 if not doing NR update )
#OUTPUTS: a list with two elements: score (vector of length(num.params)) and hessian (matrix of dim(num.params,num.params))
##########################################
if(!no.NR){
N_array=apply(emission_updates,c(1,3),"sum")
}
score=rep(0,times=dim(deriv.array)[1])
hessian=array(0,dim=c(dim(deriv.array)[1],dim(deriv.array)[1]))
hidden_states=unique(emission$emission.states["i"])
all.states=c(0,0)
prob.array=rep(0,length=dim(emission$covariate.array)[3])
#prob.array=array(0,dim=c(dim(emission$emission.states)[1],dim(emission$covariate.array)[3]))
 
for(k in 1:dim(hidden_states)[1]){
#  states=emission$emission.states[unlist(emission$emission.states["i"])==hidden_states[k,"i"],]
  states=emission$emission.states[emission$emission.states$i==hidden_states[k,"i"],]
  cov.array.states=array(emission$covariate.array[,as.numeric(rownames(states)),],dim=c(dim(emission$covariate.array)[1],dim(states)[1],dim(emission$covariate.array)[3]))
  deriv.array.states=array(deriv.array[,as.numeric(rownames(states)),],dim=c(dim(deriv.array)[1],dim(states)[1],dim(deriv.array)[3]))
  prob.matrix=get.prob.array(cov.array.states, param.values)
  
  if(!no.NR){
  N_vector=N_array[hidden_states[k,"i"],]
 # counts=matrix(emission_updates[unlist(states["i"]),unlist(states["j"]),],ncol=dim(emission_updates)[3])
  counts=matrix(emission_updates[states$i,states$j,],ncol=dim(emission_updates)[3])
  score=score+multinom.score(prob.matrix,deriv.array.states,counts,N_vector)
  hessian=hessian+multinom.hessian(prob.matrix,deriv.array.states,N_vector) 
  }
  all.states=rbind(all.states,states)
  prob.array=rbind(prob.array,prob.matrix)

}
rownames(all.states)=seq(0,(dim(all.states)[1]-1))
return(list(score=score,hessian=hessian,states=all.states[-1,],prob.array=prob.array[-1,]))
}

update.multinom.params<-function(current.params,covariate.array,deriv.array,param.types,counts=NULL,N_vector=NULL,max.it,no.NR=0){
#############################################
#update the multinomial parameters
##############################################
#Author: JL, 11/15/2011
#INPUTS:    current.params: initial values of the parameters
#           covariate_array: array of covariates for all subjects of dim(num.states,num.states,num.subjects) 
#           deriv_array: array of covariates for all subjects corresponding to non-fixed parameters
#           param.types:  vector for each parameter coding NR (0), simple update (1), or fixed (2)
#           counts: array of category counts dim(num.states, num.sub) (excluding reference category)
#           N_vector: vector (length =num.subjects) of multinomial sample size n for each subject
#           max.it=maximum number of NR updates
#        no.NR (indicator =1 if not doing NR update )
#OUTPUTS: updated values of the parameters
#######################################################
  #number of parameters for NR updaet
 num.NR.params=dim(deriv.array)[1]

 #deriv array for NR update
 ########CHANGE
 #deriv.array=array(covariate.array[param.types==0,,],dim=c(sum(param.types==0),dim(covariate.array)[2],dim(covariate.array)[3]))
 hessian=NULL
 score=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params=array(0,dim=c(num.NR.params,(max.it+1)))
 NR.params[,1]=current.params[param.types==0]

score=array(0,dim=c(length(current.params),(max.it+1)))
params=array(0,dim=c(length(current.params),(max.it+1)))
params[,1]=current.params
 
for(i in 1:(max.it)){
  prob.matrix=get.prob.array(covariate.array,params[,i])
 if(!no.NR){
  score[,i]=multinom.score(prob.matrix,deriv.array,counts,N_vector)
  hessian=multinom.hessian(prob.matrix,deriv.array,N_vector)
  NR.params[,i+1]=NR_update(current.params=NR.params[,i],score[,i], hessian)
  params[,i+1]=current.params
  params[param.types==0,i+1]<-NR.params[,i+1]
 }
}
return(list(params=params,prob.matrix=prob.matrix,hessian=hessian))
}
  

 