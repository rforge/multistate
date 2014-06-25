NULL
#' @title Individual Emission Matrices
#'
#' @description Obtain initiatl distribution matrices for each individual based on given parameter values and individual-specific 
#' covariates.
#'
#' @details Can be output as an array or list.  
#'
#' @param current.params values of emission parameters
#' @param emission.setup list with emission distribution setup information
#' @param num.states number of latent states in the CTMC
#' @param num.obs.states number of observed states
#' @param num.subjects number of individuals in the study
#' @param do.list If T, return as list, otherwise return as array
#' @param time.dep.emission (not supported) Indicator if emission distribution has time dependent covariates.
#' @return an array of (dim num.states,num.subjects) or list with the fitted emission distribution for each subject.
#' @examples
#'\dontrun{
#' library(cthmm)
#' data(DDO_data)
#' emission.list=get.emission.matrix.list(current.params=DDO_EM$params[8:10],
#'             emission.setup=emission.setup,num.states=4, 
#'            num.obs.states=3, num.subjects=500, do.list=T)
#' }
#' @author Jane Lange
#' @export
get.emission.matrix.list<-function(current.params, emission.setup, num.states, num.obs.states, num.subjects, 
do.list = T,time.dep.emission=F){
	emission = emission.setup
	update_emission = update.emission.params(current.params = current.params, 
											 emission.setup = emission, max.it = 1)
	
	if(!time.dep.emission){
		emission.list = populate_emission_matrix(update_emission, 
												 emission, num.states, num.obs.states, 
												 num.subjects, do.list = do.list)
	}else{
		emission.list = populate_emission_matrix(update_emission, 
												 emission, num.states, num.obs.states, 
												 num.subjects=dim(emission$deriv.array)[3], do.list = F)
		emission.list=lapply(seq(1,dim(emission.setup$limits)[1]),FUN="extract_array",
							 limits=emission.setup$limits, the_array=emission.list)
	}
	
	return(emission.list)
}


NULL
#' @title Individual Rate Matrices
#'
#' @description Obtain rate matrices for each individual based on given parameter values and individual-specific 
#' covariates.
#'
#' @details Can be output as an array or list.  
#'
#' @param current.params values of rate parameters
#' @param rate.setup list with rate matrix setup information
#' @param do.list If T, return as list, otherwise return as array
#' @return an array of (dim num.states,num.subjects) or list with the fitted rate matrix for each subject.
#' @examples
#'\dontrun{
#' library(cthmm)
#' data(DDO_data)
#' rate.list=get.rate.matrix.list(current.params=DDO_EM$params[1:7],rate.setup=rates.setup,do.list=T)
#' }
#' @author Jane Lange
#' @export
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
#' @title Individual Initial Distributions
#'
#' @description Obtain initial distributions for each individual based on given parameter values and individual-specific 
#' covariates.
#'
#' @details Can be output as an array or list.  
#'
#' @param current.params values of emission parameters
#' @param init.setup list with initial distribution setup information
#' @param do.list If T, return as list, otherwise return as array
#' @return an array of (dim num.states,num.subjects) or list with the fitted initial distribution for each subject.
#' @examples
#'\dontrun{
#' library(cthmm)
#' data(panel_demodata)
#' init.list=get.init.matrix.list(current.params=fit.4state$params[12:13],init.setup=init.setup,do.list=T)
#' }
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

#############################################
#functions to extract subarrays from an array
#############################################
extract_array<-function(index,limits,the_array){
temp=(the_array[,,seq(limits[index,1],limits[index,2])])
if(limits[index,1]==limits[index,2]){
out=list(temp)
}else{
out=lapply(1:dim(temp)[3], FUN="extract_mat",the_array=temp)
}
return(out)
}

extract_mat<-function(i,the_array){
return(the_array[,,i])
}




