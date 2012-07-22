
NULL
#' Obtain emission distribution matrices for each subject based on given parameter values
#'
#' Can be output as an array or list. The initial distribution will include all hidden states. 
#'
#' @inheritParams update.emission.params
#' @param num.states number of states in the model
#' @param num.obs.states number of observed states in the model
#' @param num.subjects number of individuals in the data set
#' @param do.list T if output as a list, F if array
#' @return an array of (dim num.states,num.subjects) or list with the initial dist  for each subject based on covariates and current parameter values.
#' @export
#' @author Jane Lange
get.emission.matrix.list<-function(current.params,emission.setup,num.states,num.obs.states,num.subjects, do.list=T){
emission=emission.setup

#initialize emission
update_emission=update.emission.params(current.params=current.params,
emission.setup=emission,
max.it=1)
emission.list=populate_emission_matrix(update_emission,emission,num.states,num.obs.states,num.subjects,do.list=do.list)
return(emission.list)
}



NULL
#' Obtain rate matrices for each subject based on given parameter values
#'
#' Can be output as an array or list.  
#'
#' @inheritParams update.rate.params
#' @inheritParams get.emission.matrix.list
#' @return an array of (dim num.states,num.states, num.subjects) or list with the rate matrix for each subject based on covariates and current parameter values.
#' @export
#' @author Jane Lange
get.rate.matrix.list<-function(current.params,rate.setup,do.list=T){
 rates=rate.setup
 num.states=rates$num.states
  num.subjects=dim(rates$covariate.array)[3]
  update_rates=update.rate.params(rate.setup=rates,
                                 current.params=current.params,
                                 max.it=1)
  rates.list=populate_rate_matrix(update_rates,rates,num.states,num.subjects,do.list=do.list)
  return(rates.list)
}


NULL
#' Obtain initial distribution matrices for each subject based on given parameter values
#'
#' Can be output as an array or list. The initial distribution will include all hidden states. 
#'
#' @inheritParams update.init.params
#' @inheritParams get.emission.matrix.list
#' @return an array of (dim num.states,num.subjects) or list with the initial dist  for each subject based on covariates and current parameter values.
#' @export
#' @author Jane Lange
get.init.matrix.list<-function(current.params,init.setup,do.list=T){
init=init.setup
num.subjects=dim(init$covariate.array)[3]
num.states=init$num.states
#initialize initial dist
update_init=update.init.params(current.params=current.params,
init.setup=init,
max.it=1)
delta.list=populate_initial_dist(update_init,init,num.states,num.subjects,do.list=do.list)
return(delta.list)
}



#########################################################
#Functions to get individual level rate, emission, and transition matrices
###########################################################
#Documented functions: get.rate.matrix.list, get.init.matrix.list, get.emission.matrix.list
#######################################################################################

populate_emission_matrix=function(update_emission,emission,num.states,num.obs.states,num.subjects,do.list=F){ 
############################################################## 
#Author: JL, 11/21/2011
#function to repopulate the emission matrices for each subject
#INPUTS:update_emission,emission,num.states,num.obs.states,num.subjects,do.list
#OUTPUTS:an array (dim num.state,num.obs.state,num.subjects) with the populated emission matrices for each subject
##############################################
	emission_array=array(0,dim=c(num.states,num.obs.states,num.subjects))
	for(k in 1:dim(update_emission$states)[1]){
		emission_array[update_emission$states[k,"i"],update_emission$states[k,"j"],]=update_emission$prob.matrix[k,]
	}
	for(k in 1:dim(emission$ref)[1]){
		emission_array[emission$ref[k,"i"],emission$ref[k,"j"],]=1-apply(emission_array[emission$ref[k,"i"],,],2,"sum")
	}
	if(!is.null(emission$exact.states)){
	 for(k in 1:dim(emission$exact.states)[1]){
		emission_array[emission$exact.states[k,"i"],emission$exact.states[k,"j"],]=1
	 }
	}
	if(do.list){
		emission_array=apply(emission_array,c(3),"list")
		emission_array=lapply(emission_array,c(1),FUN="[[")
	}
	return(emission_array)
}
populate_initial_dist=function(update_init,init,num.states,num.subjects,do.list=F){ 
############################################################## 
#Author: JL, 11/21/2011
#function to repopulate the initial distributions for each subject
#INPUTS:update_init,init,num.states,num.subjects,do.list=F
#OUTPUTS:an array (dim num.state,num.subjects) with the initial distribution for each subject based on covariates and current parameterized values
##############################################
#function to repopulate the rate matrices for each subject
	
	init_array=array(0,dim=c(num.states,num.subjects))
	init_array[init$states[!init$states%in%init$ref],]=update_init$prob.matrix
	init_array[init$ref,]=1-apply(init_array,2,"sum")
	if(do.list){
		init_array=apply(init_array,c(2),"list")
		init_array=lapply(init_array,c(1),FUN="[[")
	}
	return(init_array)
	
}
populate_rate_matrix=function(update_rates,rates,num.states,num.subjects,do.list=F){ 
############################################################## 
#Author: JL, 11/21/2011
#function to repopulate the emission matrices for each subject
#INPUTS:update_rates,rates,num.states,num.subjects,do.list=F
#OUTPUTS:an array (dim num.state,num.states,num.subjects) with the rate matrix for each subject based on covariates and current parameterized values
##############################################
	rate_array=array(0,dim=c(num.states,num.states,num.subjects))
	for(k in 1:dim(rates$transition.codes)[1]){
		rate_array[rates$transition.codes[k,"ni"],rates$transition.codes[k,"nj"],]=update_rates$rates[k,]
	}
	row.totals=apply(rate_array,c(1,3),"sum")
	for(l in 1:num.states){
		rate_array[l,l,]=(-1)*row.totals[l,]
	}
	
	if(do.list){
		rate_array=apply(rate_array,c(3),"list")
		rate_array=lapply(rate_array,c(1),FUN="[[")
	}
	return(rate_array)
	
}


