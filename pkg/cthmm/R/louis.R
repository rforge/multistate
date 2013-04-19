##############################
#get the information and variance/covariance matrix!!
##############################
#' @include louis.R
NULL
#' Get the information matrix for the parameter esitmates from the CTMC HMM model
#' Note: does not actually provide the variance: you need to get the inverse (or pseudo-inverse) to get the variance. 
#' @param the.data list with the data 
#' @param num.subjects number of subjects in the study
#' @param num.states number of states in the CTMC
#' @param num.obs.states number of observed states
#' @param rate.param.values MLEs for the rate params (if not estimated, set to NULL)
#' @param init.param.values MLEs for the initial distribution params (if not estimated, set to NULL)
#' @param emission.param.values MLEs for the emission params (if not estimated, set to NULL)
#' @param rates.setup list rate setup information
#' @param init.setup list with initial distribution setup information
#' @param emission.setup list with emission distribution setup information
#'
#' @return information a list with the following components
#' \item{out}{Information matrix}
#' \item{complete_info} the complete data information}
#' \item{score_sq}{Conditional expectation of the square of complete data information}
#' \item{cross_term}{Cross term in the conditional expecation of the complete data information}
#' @export
#' @author Jane Lange
get_variance=function(the.data,num.subjects,num.states,num.obs.states,rate.param.values=NULL,
emission.param.values=NULL,
init.param.values=NULL,absorb.state=NULL,rates.setup,emission.setup,init.setup){
	print(rates.setup[[1]])
	ij.indices=matrix(unlist(rates.setup$transition.codes),ncol=2,byrow=F)
	colnames(ij.indices)=c("i","j")
	
	obs.data.list<-lapply(the.data, "[[",c("obs.data"))
	exact.time.ranks.list<-lapply(the.data,"[[",c("exact.times.ranks"))
	time.diffs.list<-lapply(lapply(the.data,"[[",c("obs.times")),FUN="diff")
	
	rep_item<-function(i,item){
		return(item)
	}
	
	if(!is.null(rates.setup$fixed.rates)){
		rates.list=lapply(seq(1:num.subjects), FUN="rep_item",item=rates.setup$fixed.rates)
	}else{
		rates.list=get.rate.matrix.list(rate.param.values,rate.setup=rates.setup)
	}
	if(!is.null(emission.setup$fixed.dist)){
		emission.list=lapply(seq(1:num.subjects),FUN="rep_item",item=emission.setup$fixed.dist)
	}else{
		emission.list=get.emission.matrix.list(emission.param.values,emission.setup=emission.setup,num.states=num.states,num.obs.states=num.obs.states,num.subjects=num.subjects)
	}
	if(!is.null(init.setup$fixed.dist)){
		delta.list=lapply(seq(1:num.subjects),FUN="rep_item",item=init.setup$fixed.dist)
	}else{
		delta.list=get.init.matrix.list(current.params=init.param.values,init.setup=init.setup)
	}
	
	eigen.decomp.list<-eigen_decomp_list_R(rates.list) 
	transition.probabilties.list<-mapply(time.diffs.list,eigen.decomp.list, exact.time.ranks.list,FUN=transition_prob_list,SIMPLIFY=F)
	likelihood.forward.backward.list<-mapply(obs.data.list,transition.probabilties.list,delta.list,emission.list,FUN="forwardback_R",SIMPLIFY=F)
	durations<-dur_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list,exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state)
	transitions<-trans_loop_R(likelihood.forward.backward.list,time.diffs.list,eigen.decomp.list,obs.data.list,emission.list, exact.time.ranks.list,the_state_size=num.states,absorb_state=absorb.state,ij.indices)
	
	E.emission.init=E_step_emission_init_all(likelihood.forward.backward.list, obs.data.list,num.states,num.obs.states,num.subjects)
	emission.counts=E.emission.init$emission_array
	init.counts=E.emission.init$init_array
	
#get M1 and M2 lists
	
	M=get_M1_M2_list(rates.setup,init.setup,emission.setup,delta.list,rates.list,emission.list, obs.data.list,eigen.decomp.list,time.diffs.list,
					 exact.time.ranks.list,transition.probabilties.list,absorb.state=absorb.state)
	M1.list=M$M1_list
	M2.list=M$M2_list
	
	information=louis_information(rates.setup,init.setup,emission.setup,delta.list,emission.list,rates.list,M1.list,M2.list,
								  rate.param.values=rate.param.values,init.param.values=init.param.values,emission.param.values=emission.param.values,durations, emission.counts,init.counts,num.subjects,transitions)
	return(list(information=information))
}

NULL
#' Get the information matrix using the Louis decomposition method.
#' @param rates.setup list rate setup object
#' @param init.setup init setup object
#' @param emission.setup emission setup object
#' @param init.list list with individual initial distributions
#' @param emission.list with individual emission distributions
#' @param rates.list with individual rate matrices
#' @param M1.list list with conditional means of sufficient statistics for each subject
#' @param M2.list list with 2nd moments of sufficient statistics for each subject
#' @param rate.param.values values of rate parameters
#' @param init.param.values values of initial distribution parameters
#' @param emission.param.values of emission parameters
#' @param durations array (dim num.states x num.subjects) with the expected durations 
#' @param emission.counts array expected nunber of o_ij 's dim (num states x num.obs.states x num.subjects)
#' @param num.subjects number of subjects in the study
#'
#' @return list with the following components
#' \item{out}{Information matrix}
#' \item{complete_info} {the complete data information}
#' \item{score_sq}{Conditional expectation of the square of complete data information}
#' \item{cross_term}{Cross term in the conditional expecation of the complete data information}
#' @export
#' @author Jane Lange

louis_information<-function(rates.setup,init.setup,emission.setup,delta.list,emission.list,rates.list,M1.list,M2.list,
                  rate.param.values=NULL,init.param.values=NULL,emission.param.values=NULL,durations, emission.counts,init.counts,num.subjects,transitions){
 #hessian component

   if(is.null(rates.setup$fixed.rates)){
    rate_hessian1=rate_hessian_louis(rates.setup,rate.param.values,durations)
    rate_score=rate_score_louis(rates.setup,rate.param.values,durations,transitions)
   }else{
    rate_hessian1=NA
    rate_score=matrix(NA,nrow=1,ncol=num.subjects)
   }
   if(is.null(init.setup$fixed.dist)){
    init_hessian1=init_hessian_louis(init.setup,init.param.values)
    init_score=init_score_louis(init.setup,init.param.values,init.counts)
   }else{
    init_hessian1=NA
    init_score=matrix(NA,nrow=1,ncol=num.subjects)
   }
   if(is.null(emission.setup$fixed.dist)){
    emission_hessian1=emission_hessian_louis(emission.setup,emission.param.values,emission.counts)
    emission_score=emission_score_louis(emission.setup,emission.param.values,emission.counts)  
   }else{
    emission_hessian1=NA
    emission_score=matrix(NA,nrow=1,ncol=num.subjects)
   }
 
   hessian_component=adiag(rate_hessian1,init_hessian1,emission_hessian1)
   
   hessian_component=hessian_component[!is.na(apply(hessian_component,1,"sum")),!is.na(apply(hessian_component,1,"sum"))] 
   print(hessian_component)
  #"cross term" component 
 score=rbind(rate_score,init_score,emission_score)
 score=score[!is.na(apply(score,1,"sum")),]
 score_out=array(0,dim=c(dim(hessian_component)[1],dim(hessian_component)[2]))  
 for(i in 1:num.subjects){
     score_out=score_out+(score[,i])%*%t(score[,i])
  }
  
 #E[SS^T|o] component
  diff_theta_beta=diff_theta_beta_array(rates.setup,emission.setup,init.setup,num.subjects)
  score_sq=diff_theta_beta[,,1]%*%louis_decomp_score_sq(rates.setup,init.setup,emission.setup,delta.list[[1]],emission.list[[1]],rates.list[[1]],M1.list[[1]],M2.list[[1]])%*%t(diff_theta_beta[,,1])

 for(i in 2:num.subjects){
  #for(i in 2:2){
  score_sq_new=louis_decomp_score_sq(rates.setup,init.setup,emission.setup,delta.list[[i]],emission.list[[i]],rates.list[[i]],M1.list[[i]],M2.list[[i]])
  score_sq=score_sq+diff_theta_beta[,,i]%*%score_sq_new%*%t(diff_theta_beta[,,i])
 }
 
  print(score_sq)

  #observed information
  out=-hessian_component-score_sq+score_out
  return(list(out=out, complete_info=hessian_component, score_sq=score_sq, cross_term=score_out))
}
NULL
#' Outerproduct
#' 
#' Get (c1*D)(c2*E)^T.
#' Suppose Y is the vector of complete data sufficient statistics.
#'
#' If l_1[k]=0 (l_2[k]=0) this signfies a value of 1 rather than a suffcient statistic in Y. 
#' Hence in this case the desired product is a function of the first moment vector M1 rather than M2.
#' @param l_1 vector of indices of Y in D, i.e. D=Y[l_a] 
#' @param l_2 vector of indices of Y in E, i.e. E=Y[l_2]
#' @param M2 matrix of E[YY^T|o]
#' @param M1 vector of E[Y|o]
#' @param c1 vector of constants to be multiplied by elements in A
#' @param c2 vector of constants to be multipled by elements in B
#' @return matrix c1*D(c2*E)^T
#' @export
#' @author Jane Lange
outerproduct<-function(l_1,l_2,c1,c2,M1,M2){
outmatrix<-array(0,dim=c(length(l_1),length(l_2)))
for(i in 1:length(l_1)){
for(j in 1:length(l_2)){
  if(l_1[i]>0&&l_2[j]>0){
    if(is.null(dim(M2))){
     outmatrix[i,j]=M2[l_1[i]]*c1[i]*c2[j] 
     }else{
     outmatrix[i,j]=M2[l_1[i],l_2[j]]*c1[i]*c2[j]
     }  
  }else if(l_2[j]==0&&l_1[i]>0){
     outmatrix[i,j]=M1[l_1[i]]*c1[i]*c2[j]
  #    print("that")
  }else if(l_1[i]==0&&l_2[j]>0){
    outmatrix[i,j]=M1[l_2[j]]*c1[i]*c2[j]
   # outmatrix[i,j]=0
   #  print("this")
    print(c(i,j))
  }else if(l_1[i]==0&&l_2[j]==0){
       outmatrix[i,j]=c1[i]*c2[j]
  }
 }
}

return(outmatrix)
}

NULL
#' get.suff.stat.table
#'
#' Get the information matrix using the Louis decomposition method.
#' @param rates.setup list rate setup object
#' @param init.setup init setup object
#' @param emission.setup emission setup object
#'
#'@return a data frame with num rows=number of sufficient stats, and columns ni nj di zi Ni oi oj. If the statistic is n_12, then the corresponding column is (1,2,0,0,0,0,0)
get.suff.stat.table<-function(rates.setup,init.setup,emission.setup){
 
  nij=matrix(nrow=1,ncol=2)
  di=matrix()
  zi=matrix()
  Ni=matrix()
  oij=matrix(nrow=1,ncol=2)  

  if(is.null(rates.setup$fixed.rates)){
     nij=matrix(unlist(rates.setup$transition.codes),byrow=F,ncol=2,dimnames=list(NULL, c("ni","nj")))
     di=matrix(unique(rates.setup$transition.codes[,"ni"]),ncol=1,dimnames=list(NULL,"di"))
  }
  if(is.null(init.setup$fixed.dist)){
     zi=matrix(init.setup$states[init.setup$states!=init.setup$ref],ncol=1,dimnames=list(NULL,"zi"))
  }
  if(is.null(emission.setup$fixed.dist)){
      Ni=matrix(unique(emission.setup$emission.states[,"i"]),ncol=1,dimnames=list(NULL,"Ni"))
      oij=matrix(unlist(emission.setup$emission.states),byrow=F,ncol=2,dimnames=list(NULL,c("oi","oj")))
  }

  suff.stats=adiag(nij,di,zi,Ni,oij)
  suff.stats=suff.stats[!is.na(apply(suff.stats,1,"sum")),]
  colnames(suff.stats)=c("ni","nj","di","zi","Ni","oi","oj")

  if(is.null(rates.setup$fixed.rates)){
  if(sum(rates.setup$param.types!=0)==length(rates.setup$param.types)){
   suff.stats<-suff.stats[suff.stats$ni==0,]
   suff.stats<-suff.stats[suff.stats$di==0,]
  }
  }
  if(is.null(init.setup$fixed.dist)){
  if(sum(init.setup$param.types!=0)==length(init.setup$param.types)){
   suff.stats<-suff.stats[suff.stats$zi==0,]
  }
  }
  if(is.null(emission.setup$fixed.dist)){
  if(sum(emission.setup$param.types!=0)==length(emission.setup$param.types)){
   suff.stats<-suff.stats[suff.stats$oi==0&suff.stats$Ni==0,]
  }
  }
  rownames(suff.stats)=seq(1,dim(suff.stats)[1])
  return(data.frame(suff.stats))
}

#' Get the score element indices
#' 
#' The complete data score equations consist of elements such as nij-di*lambda ij, zi-pi, and o_ij-Ni*e_ij.
#' Hence the score can be written as A-cB, where A is a vector of nij,zi,oij and B, di,1,Ni, and c, lambda ij,pi,eij
#' @param suff.stats table with coding all of the sufficient statistics
#' @param rate.matrix the rate matrix
#' @param emission.matrix the emission matrix
#' @param init.dist vector with the initial distribution
#' @return list with 
#'  \item{l_a}{row index of suff.stats corresponding to A}
#'  \item{l_b}{row index of suff.stats corresponding to B--for 1, this index is set to 0} 
#'  \item{c}{vector of constants}
#' @export
#' @author Jane Lange
element.indices=function(suff.stats,rate.matrix,emission.matrix,init.dist){
  Atable=suff.stats[suff.stats[,"ni"]>0|suff.stats[,"zi"]>0|suff.stats[,"oi"]>0,]
  l_a=as.numeric(rownames(Atable))
  l_b=rep(0,times=length(l_a))
  c=rep(1,times=length(l_a))
  for(k in 1:dim(Atable)[1]){
    if(Atable[k,"ni"]>0){
      l_b[k]=as.numeric(rownames(suff.stats[Atable[k,"ni"]==suff.stats[,"di"],]))
      c[k]=rate.matrix[Atable[k,"ni"],Atable[k,"nj"]]
    }
    if(Atable[k,"oi"]>0){
      l_b[k]=as.numeric(rownames(suff.stats[Atable[k,"oi"]==suff.stats[,"Ni"],]))
      c[k]=emission.matrix[Atable[k,"oi"],Atable[k,"oj"]]
    }
    if(Atable[k,"zi"]>0){
      l_b[k]=0
      c[k]=init.dist[Atable[k,"zi"]]
    }
  }
 return(list(l_a=l_a,l_b=l_b,c_b=c))
}

NULL
#' Louis decomposition score component
#' 
#'  Get E[SS^T|o] for a single subject: i.e. the louis decomposition score component
#'  @param rates.setup rate setup object
#'  @param init.setup initial distribution setup object
#'  @param emission.setup emission distribution setup object
#'  @param rate.matrix rate matrix for individual
#'  @param emission.matrix emission matrix for individual
#'  @param init.dist vector of initial dist for individual
#'  @param M1 vector of first moments of complete data sufficient statistics
#'  @param M2 matrix of second moments of complete data sufficient statistics
#'  @return matrix with second moments of the score function. 
louis_decomp_score_sq<-function(rates.setup,init.setup,emission.setup,init.dist,emission.matrix,rate.matrix,M1,M2){

     
    suff.stats=get.suff.stat.table(rates.setup,init.setup,emission.setup)
    
    
    
    the.indices=element.indices(suff.stats,rate.matrix,emission.matrix,init.dist)
    l_a=the.indices$l_a
    l_b=the.indices$l_b
    c_b=the.indices$c_b
    c_a=rep(1,times=length(l_a))
    out=outerproduct(l_a,l_a,c_a,c_a,M1,M2)-outerproduct(l_a,l_b,c_a,c_b,M1,M2) -outerproduct(l_b,l_a,c_b,c_a,M1,M2)+outerproduct(l_b,l_b,c_b,c_b,M1,M2)

         
   return(out)
    
}    

NULL
#'diff_theta_beta_array
#'get an array that corresponds to diff_theta_beta for each individual in the dataset
#'The "theta" params are the natural parameters encoding the rate, initial distribution, and emission matrices. 
#'The "beta"params are paramters in the linear predictors that relate to the natural parameters. 
#' @param rates.setup rate setup object
#' @param emission.setup emission setup object
#' @param init.setup init setup object
#' @param num.subjects number of subjects
#' @return an array of dimension (num theta param)x(num beta params)xnumber of subjects.
diff_theta_beta_array<-function(rates.setup,emission.setup,init.setup,num.subjects){

 if(is.null(rates.setup$fixed.rates)){
   num.rate.params=dim(rates.setup$deriv.array)[1]
   num.rate.stats=I(num.rate.params>0)*dim(rates.setup$deriv.array)[2]
 }else{
   num.rate.params=0
   num.emission.stats=0
 }
 
 if(is.null(emission.setup$fixed.dist)){
 num.emission.params=dim(emission.setup$deriv.array)[1]
 num.emission.stats=I(num.emission.params>0)*dim(emission.setup$deriv.array)[2]
 }else{
   num.emission.params=0
   num.emission.stats=0
 }
 
 if(is.null(init.setup$fixed.dist)){
 num.init.params=dim(init.setup$deriv.array)[1]
 num.init.stats=I(num.init.params>0)*dim(init.setup$deriv.array)[2]
 }else{
  num.init.params=0
  num.init.stats=0
 }
 #numbers of sufficient statistics
 out.array=array(0,dim=c(num.rate.params+num.emission.params+num.init.params,num.rate.stats+num.emission.stats+num.init.stats,num.subjects))
 
 if(num.rate.params>0){
 out.array[1:num.rate.params, 1:num.rate.stats,]<-rates.setup$deriv.array
 }
 if(num.init.stats>0){
 out.array[(num.rate.params+1):(num.rate.params+num.init.params),(num.rate.stats+1):(num.init.stats+num.rate.stats),]<-init.setup$deriv.array
 }
 if(num.emission.stats>0){
   out.array[(num.rate.params+num.init.params+1):(num.rate.params+num.init.params+num.emission.params),
           (num.init.stats+num.rate.stats+1):(num.init.stats+num.rate.stats+num.emission.stats),]<-emission.setup$deriv.array
 }        
 return(out.array)
}

NULL
#'get_M1_M2
#' get the list with first and second moments of complete data sufficient statistics for an individual, using the recursive smoothing method. 
#' @param stat.table data frame with num rows=number of sufficient stats, and columns ni nj di zi Ni oi oj. If the statistic is n_12, then the corresponding column is (1,2,0,0,0,0,0)
#' @param time.diffs observation time intervals (single subject)
#' @param obs.data observed data (single subject)
#' @param eigen.decomp eigen.decomp object for rate matrix (single subject)
#' @param exact.time.ranks number encoding the absorption time index
#' @param absorb.state one or more absorbing states in the model
#' @param transition.probs transition probabilities for each time interval
#' @param emission.matrix emission matrix for subject
#' @param init.dist initial distribution for sujbect
#' @return list consisting of 
#' \item {M1} {vector of first moments of complete data sufficient statistics}
#' \item {M2} {matrix of second moments of complete data sufficient statistics}
get_M1_M2<-function(stat.table,time.diffs,obs.data,eigen.decomp,exact.time.ranks,absorb.state,transition.probs,emission.matrix,init.dist){
     
 num.states=dim(eigen.decomp$rate)[1]
 suff.stats=stat.table

 t1=t1_1(suff.stats,obs_data=obs.data,num_states=num.states) 
 s1=s_1(suff.stats,num.states,obs.data,time.diffs,eigen.decomp, exact_time_ranks=exact.time.ranks, absorb_state=absorb.state,trans_probs_list=transition.probs)
   
 t2=t1_2(stat_table=suff.stats,obs_data=obs.data,num_states=num.states,t1_1=t1)
 s2=s_2(stat_table=suff.stats,num_states=num.states,time_diffs=time.diffs,eigen_decomp=eigen.decomp, exact_time_rank=exact.time.ranks, absorb_state=absorb.state,s1,transition.probs)
 tau1=cappe.first.moment.tau(transition.probs,emission.matrix,init.dist,obs.data,s1,t1)
 tau2=cappe.second.moment.tau(transition.probabilities=transition.probs,emission.matrix,init.dist,obs.data,s1,s2,t2,tau1)

 M1=apply(tau1, c(2,3),"sum")[dim(tau1)[2],]
 M2=apply(tau2, c(2,3,4),"sum")[dim(tau2)[2],,]

 return(list(M1=M1,M2=M2))
}

NULL
#' get_M1_M2_list 
#' Get a list of first and second moments for each subject
#' @param rates.setup rate setup object
#' @param init.setup init setup object
#' @param emission.setup emission setup object
#' @param delta.list list with intial distribution for each subject
#' @param rates.list list with rates for each subject
#' @param emission.list list with emission matrices for each subject
#' @param obs.data.list list with obs.data for each subject
#' @param eigen.decomp.list list with eigen.decomp objects for each subject
#' @param time.diffs.list list with interobservation intervals for each subject
#' @param exact.time.ranks.list list encoding the observation number of times of absorption for each subject
#' @param transition.probabilities.list list encoding the transition probabilities for each subject
#' @param absorb.state vector with one or more absorbing states
#' @return a list consisting of 
#' \item {M1_list} {list with vectors of first moments for complete data sufficient statistics}
#' \item M2_list {list with matrices of second moments for complete data sufficient statistics}
get_M1_M2_list<-function(rates.setup,init.setup,emission.setup,delta.list,rates.list,emission.list, obs.data.list,eigen.decomp.list,
                         time.diffs.list,exact.time.ranks.list,transition.probabilties.list,absorb.state){

  num.subjects=length(rates.list) 
 stat.table=get.suff.stat.table(rates.setup,init.setup,emission.setup)
 M1.list<-vector("list", num.subjects)
 M2.list<-vector("list",num.subjects)
 for(i in 1:num.subjects){
     out=get_M1_M2(stat.table,time.diffs.list[[i]],obs.data.list[[i]],eigen.decomp.list[[i]],exact.time.ranks.list[[i]],absorb.state,transition.probabilties.list[[i]],emission.list[[i]],delta.list[[i]])
     M1.list[[i]]=out$M1
     M2.list[[i]]=out$M2
     print(i)
  }
 return(list(M1_list=M1.list, M2_list=M2.list))
}
NULL
#'init_score_louis
#'Get the inital dist scores
#' @param init.setup list init setup object
#' @param init.param.values values of init dist parameters
#' @param  init.counts matrix (num.states x num.subjects) with multinomial vector encoding the initial dist.
#' @return score array (dim num params x number of subjects)
init_score_louis<-function(init.setup,param.values,init.counts){
     
# if(sum(init.setup$param.types==0)==0){return(c(0))}else{
	
	init=init.setup
	if(!is.null(init.counts)){
		counts=matrix(init.counts[init$states[init$states!=init$ref],],ncol=dim(init.counts)[2])
	}
	N_vector=rep(1,times=dim(init$covariate.array)[3])
	
	deriv.array=init.setup$deriv.array
	prob.matrix=get.prob.array(init.setup$covariate.array,param.values)
	num.subjects=dim(counts)[2]
	suff.stat.score=array(0,dim=c(dim(counts)[1],dim(counts)[2]))  
	param.score=array(0,dim=c(dim(deriv.array)[1],dim(counts)[2]))
	suff.stat.score=(counts-prob.matrix*N_vector)
	for(m in 1:num.subjects){
        param.score[,m]=matrix(deriv.array[,,m])%*%suff.stat.score[,m]
	}
	return(param.score)
# hessian=multinom.hessian(prob.matrix,deriv.array,N_vector)
	return(param.score)
# }
}

NULL
#'Get the complete data hessian matrix for initial dist parameters
#' @param init.setup list init setup object
#' @param param.values values of init parameters
#'
#' @return complete data hessian matrix for init parameters     
init_hessian_louis<-function(init.setup,param.values){
  # if(sum(init.setup$param.types==0)==0){return(c(0))}else{
  deriv.array=init.setup$deriv.array
  prob.matrix=get.prob.array(init.setup$covariate.array,param.values)
  N_vector=rep(1,times=dim(init.setup$covariate.array)[3])
  hessian=multinom.hessian(prob.matrix,deriv.array,N_vector)
  return(hessian)
  # }
}

NULL
#'emission_score_louis
#'Get the score matrix rate parameters
#' @param emission.setup list emission setup object
#' @param emission.param.values values of emission parameters
#' @param emission.counts array (dim num states x num obs states x num subjects) with the (expected) counts of the o_ij  
#' @return score array (dim num params x number of subjects)
emission_score_louis<-function(emission.setup,emission.param.values,emission.counts){
	out=emission.score.hessian(emission.counts,emission.param.values,emission.setup,emission.setup$deriv.array)
	state_order=emission.setup$emission.states
	states=out$states
	states=data.frame(states)
	state_order=data.frame(emission.setup$emission.states)
	reorder=as.numeric(paste(row.names(states[state_order[,"i"]==states[,"i"]&state_order[,"j"]==states[,"j"],])))
	
	N_vector=apply(emission.counts,c(1,3),"sum")
	
#get the score for each individual
	suff.stat.score=array(0,dim=c(dim(emission.setup$deriv.array)[2],dim(emission.setup$deriv.array)[3]))
	
	for(l in 1:dim(state_order)[1]){
		suff.stat.score[l,]=emission.counts[state_order[l,"i"],state_order[l,"j"],]-out$prob.array[reorder[l],]*N_vector[state_order[l,"i"],]
	}
	out_score=array(0,dim=c(dim(emission.setup$deriv.array)[1],dim(emission.setup$deriv.array)[3]))
	for(m in 1:dim(emission.setup$deriv.array)[3]){
		out_score[,m]=emission.setup$deriv.array[,,m]%*%suff.stat.score[,m]
	}
	return(out_score)
}

NULL
#'emission_hessian_louis
#'Get the complete data hessian matrix for emission parameters
#' @param emission.setup list emission setup object
#' @param param.values values of emission parameters
#' @param emission.counts array (dim num states x num obs states x num subjects) with the (expected) counts of the o_ij  
#'
#' @return complete data hessian matrix for emission parameters
emission_hessian_louis<-function(emission.setup, param.values,emission.counts){
  num.subjects=dim(emission.setup$covariate.array)[3]
  deriv.array=emission.setup$deriv.array
  hessian=emission.score.hessian(emission_updates=emission.counts,param.values=param.values,emission=emission.setup,deriv.array=deriv.array,no.NR=0)$hessian
  return(hessian)
}
NULL
#'rate_hessian_louis
#'Get the complete data hessian matrix for rate parameters
#' @param rates.setup list rate setup object
#' @param rate.param.values values of rate parameters
#' @param durations array (dim num.states x num.subjects) with the expected durations 
#'
#' @return complete data hessian matrix for rate parameters

rate_hessian_louis<-function(rates.setup,rate.param.values,durations){
  deriv.array=rates.setup$deriv.array
  rates=get.rates(param.values=rate.param.values,covariate.array=rates.setup$covariate.array)
  hessian=rate_hessian(rates,durations, deriv.array, rates.setup$transition.codes)
  return(hessian) 
}
NULL
#'rate_score_louis
#'Get the score matrix rate parameters
#' @param rates.setup list rate setup object
#' @param rate.param.values values of rate parameters
#' @param durations array (dim num.states x num.subjects) with the expected durations 
#' @param transitions array (dum num.states x num.states x num.subjects) with the expected transition counts 
#' @return score array (dim num params x number of subjects)
rate_score_louis<-function(rates.setup,rate.param.values,durations,transitions){
  
  deriv.array=rates.setup$deriv.array
  transition.codes=rates.setup$transition.codes
  num.trans=dim(transition.codes)[1]
  rates=get.rates(param.values=rate.param.values,covariate.array=rates.setup$covariate.array)
  
  suff.stat.score=array(NA,dim=c(dim(transition.codes)[1],dim(durations)[2]))
  score=array(0,dim=c(dim(deriv.array)[1],dim(transitions)[3]))
  
  
  for(l in 1:dim(transition.codes)[1]){
    suff.stat.score[l,]=transitions[transition.codes[l,"ni"],transition.codes[l,"nj"],]-rates[l,]*durations[transition.codes[l,"ni"],]
  }
  for(m in 1:dim(durations)[2]){
    score[,m]=deriv.array[,,m]%*%suff.stat.score[,m]
  }
  
  return(score)
  
  
}



