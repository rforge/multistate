forward_modified<-function(obs.data, transition.probabilities, delta, emission.matrix){
#################################################################################
#Author: JL, 9/19/2011
#This function computes modfied forward probabilities P(O_t, X_t=j|o_{1:t-1})
#INPUTS: the.data= list with all subjects' data
#        transition.probabilities.list=list of transition probabilities from t1-t2, t2-t3,...,t_n-1-t_n
#        delta.list = vector correspondding to initial distribution of hidden states
#        emmision.list=list of matrices with the emission probablities for all subjects
#        the ith row corresponds to the hidden value X(t)=i, and the kth column to O(t)=k|X(t)=i
#        thus the rows sum to 1, and k columns correspond to the k possible observed states
#OUTPUTS: a_matrix, the s x n (num.states x n time points) matrix of modified forward probabiliites. 
####################################################################################
#modified forward algorithm
Pi=transition.probabilities
emission.matrix=emission.matrix
e=emission.matrix
o<-obs.data

num.states=dim(Pi[[1]])[1]
num.times=length(o)
a<-matrix(0, nrow=num.states, ncol=num.times)
a[,1]=e[,o[1]]*delta

for(k in 2:num.times){
for(j in 1:num.states){
  a[j,k]=e[j,o[k]]*sum(a[,k-1])^(-1)*sum(a[,k-1]*Pi[[k-1]][,j])
  }
}
return(a)
}

cappe.first.moment.tau<-function(transition.probabilities,emission.matrix,delta,obs.data,s,t1){
#Author: JL, 9/19/2011
#This function computes the recursive functions tau for complete data sufficient statistics based on Cappe 2011
#        tau_k(x_k)=E[I(X_k=x_k)t(x_1:k)|o_1:k] for subject k
#INPUTS:
#        transition.probabilities=list of transition probabilities from t1-t2, t2-t3,...,t_n-1-t_n
#        emmision.matrix=the  matrix with the emission probablities 
#        delta = initial distribution 
#        obs.data= the subject's observed data
#         s=4 dimensional array with dim=c(state.size,state.size,num.times,num.stats)), corresponding to 
#           s_k(x_k, x_{k+1})
#         t1=matrix (num.states x num.stats) with initial value t(x1) for all sufficient statistics
#OUTPUTS: tau, the recursive updates for online first moments for each of the sufficient statistics
#         tau is a 3-d array with dim=c(num.state,num.times,num.stats)
####################################################################################

 Pi=transition.probabilities
 delta=delta
 e=emission.matrix
 o<-obs.data

 num.states=dim(Pi[[1]])[1]
 num.times=length(o)
 num.stats= dim(s)[4]
 a=forward_modified(obs.data,transition.probabilities,delta,emission.matrix)
 filtering.prob=t(t(a)/apply(a,2,sum))
 conditional.likelihood=apply(a,2,sum)

 tau=array(0,c(num.states,num.times,num.stats))
  
 for(m in 1:num.stats){
  for(j in 1:num.states){
   tau[j,1,m]=e[j,o[1]]*delta[j]
  }
   tau[,1,m]=t1[,m]*tau[,1,m]/sum(tau[,1,m])
 }
 
 
for(m in 1:num.stats){
 for(l in 1:(num.times-1)){
  for(j in 1:num.states){
  tau[j,l+1,m]=conditional.likelihood[l+1]^(-1)*e[j,o[l+1]]*sum((tau[,l,m]*Pi[[l]][,j])+filtering.prob[,l]*s[,j,l,m]*Pi[[l]][,j])
 }
}
}
return(tau)
}
  
cappe.second.moment.tau<-function(transition.probabilities,emission.matrix,delta,obs.data,s1,s2,t2,tau1){
#Author: JL, 9/26/2011
#This function computes the recursive functions tau for complete data sufficient statistics based on Cappe 2011
#        tau_k(x_k)=E[I(X_k=x_k)t(x_1:k)|o_1:k] for subject k
#         where t corresponds to the matrix of second moments of sufficient statistics.
#INPUTS: 
#        transition.probabilities=list of transition probabilities from t1-t2, t2-t3,...,t_n-1-t_n
#        emmision.matrix=list of matrices with the emission probablities for all subjects
#        the ith row corresponds to the hidden value X(t)=i, and the kth column to O(t)=k|X(t)=i
#        thus the rows sum to 1, and k columns correspond to the k possible observed states
#        delta.list = vector correspondding to initial distribution of hidden states
#        the.data= list with all subjects' data
#         s1=4 dimensional array with dim=c(state.size,state.size,num.times,num.stats)), corresponding to 
#           s1_k(x_k, x_{k+1})
#         s2=5 dimensional array with dim=c(state.size,state.size,num.times,num.stats,num.stats)), corresponding to 
#           s2_k(x_k, x_{k+1})
#         t2=matrix (num.states x num.stats x num.stats) with initial value t(x1) for all sufficient statistics
#         tau1=precalculate matrix (num.states,num.times) with recursive functions for first moments.
#OUTPUTS: tau2, the recursive updates for online second moments for each of the sufficient statistics
#         tau2 is a 4-d array with dim=c(num.state,num.times,num.stats,num.stats)
# NOTE:  To get the expected value of second moments conditional on all observed data use:
#        apply(tau2[,num.times,,],c(2,3),"sum")
####################################################################################

Pi=transition.probabilities
 delta=delta
 e=emission.matrix
 o<-obs.data

 num.states=dim(Pi[[1]])[1]
 num.times=length(o)
 num.stats= dim(s1)[4]
 a=forward_modified(obs.data,transition.probabilities,delta,emission.matrix)
 filtering.prob=t(t(a)/apply(a,2,sum))
 conditional.likelihood=apply(a,2,sum)


 tau2=array(0,c(num.states,num.times,num.stats,num.stats))
 #initialize tau2
for(m in 1:num.stats){
       for(n in 1:num.stats){
          for(j in 1:num.states){
             tau2[j,1,m,n]=e[j,o[1]]*delta[j]
          }
          tau2[,1,m,n]=t2[,m,n]*tau2[,1,m,n]/sum(tau2[,1,m,n])
        }
}

#update tau2
for(l in 1:(num.times-1)){
  for(j in 1:num.states){
    temp=array(0,dim=c(num.stats,num.stats))
    for(i in 1:num.states){
      temp=temp+Pi[[l]][i,j]*(tau2[i,l,,]+tau1[i,l,]%*%t(s1[i,j,l,])+s1[i,j,l,]%*%t(tau1[i,l,])+s2[i,j,l,,]*filtering.prob[i,l])
    }
    tau2[j,l+1,,]=temp*conditional.likelihood[l+1]^-1*e[j,o[l+1]]
  }
}
return(tau2)
}





