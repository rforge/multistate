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
  
hazard_times<-function(start=0,end,state_of_interest,at_risk_states,alpha,rate,length.out=1000){
  times=seq(start,end,length.out=length.out)
  out=unlist(lapply(times,FUN=get_hazard,state_of_interest=state_of_interest,at_risk_states=at_risk_states,alpha=alpha,rate=rate))
  return(list(hazard=out,times=times))
}
sub_dist_times<-function(start=0,end,states,alpha,rate,length.out=1000){
  times=seq(start,end,length.out=length.out)
  out=unlist(lapply(times,FUN=get_sub_distribution,states=states,alpha=alpha,rate=rate))
  return(list(dist=out,times=times))
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
#############################################################################
#get the SE of the distribution fun at different times
se.dist.times<-function(start,end,covar,rate,states,alpha, param.deriv,num.transitions,num.params,transitions,num.states,length.out=100){
 times=seq(start,end,length.out=length.out)
 deriv.list=lapply(times,FUN="dist.fun.deriv",rate,states,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
 sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
}

delta.method.mat=function(x,covar){
  return(x%*%covar%*%x)
}

se.haz.times<-function(start,end,covar,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states,length.out=100){
  times=seq(start,end,length.out=length.out)
  deriv.list=lapply(times,FUN="haz.fun.deriv",rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
  se=sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
  return(se)
}

se.dens.times<-function(start,end,covar,rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states){
 times=seq(start,end,length.out=100)
 deriv.list=lapply(times,FUN="dens.fun.deriv",rate,at_risk_states,state_of_interest,alpha, param.deriv,num.transitions,num.params,transitions,num.states)
 se=sqrt(unlist(lapply(deriv.list,FUN=delta.method.mat,covar=covar)))
 return(se)
}
  
##############################################################################
