

E_step_DDO<-function(likelihood.forward.backward.list, h.list,num.states,num.subjects){
	out=.Call("E_step_DDO_all",likelihood.forward.backward.list,h.list,num.states,num.subjects)
	return(out)
}

E_step_emission_init_all<-function(likelihood.forward.backward.list, obs.data.list,num.states,num.obs.states,num.subjects){

		out=.Call("E_step_emission_init_all",likelihood.forward.backward.list,obs.data.list,num.states,num.obs.states,num.subjects)
        return(out)
}


E_step_emission_init_all_time_dep<-function(likelihood.forward.backward.list, obs.data.list,num.states,num.obs.states,num.subjects,total.obs.times,limits){
  out=.Call("E_step_emission_init_all_time_dep",likelihood.forward.backward.list,obs.data.list,num.states,num.obs.states,num.subjects,total.obs.times,limits)
  return(out)
}



###########################################
#Functions to do the E-step for emission and inital distributions
#Used for documented functions: E.step.emission.init.all
################################################
#' @include E_step_init_emission.r
NULL
#' E step for emission and initial distributions
#'
#'Get the expected emission counts E[(O_ij)|o] and initial distribution 
#' counts E[(Z_i)|o] given observed data for all individuals.
#' This implementation does not use the recursive filtering method, but rather uses forward and backward probabilities.
#' @param likelihood.forward.backward.list list with LL, forward, and backward probs for each inividual
#' @param the.data list with suject level obs.data
#' @param num.states number of hidden states in the model
#' @param num.obs.state number of observed states in the model
#' @param num.subjects number of individuals
#' @return a list with two elements:
#' \item{emission.array}{(dim (num.states, num.obs.states, num.subjects))}
#' \item{init.array}{(dim (num.states,num.subjects))}
#' @export
#' @author Jane Lange
################################################
E.step.emission.init.all<-function(likelihood.forward.backward.list, the.data,num.states,num.obs.states,num.subjects){
############################################################## 
#Author: JL, 11/21/2011
#E-step for initial and emission distributions for all subjects
#INPUTS:num.states,num.obs.states,num.subjects,likelihood.forward.backward.list, the.data
#OUTPUTS:a list with two elements: emission.array(dim (num.states, num.obs.states, num.subjects)) and init.array (dim (num.states,num.subjects))
###############################################################
	emission_array=array(0,dim=c(num.states,num.obs.states,num.subjects))
	init_array=array(0,dim=c(num.states,num.subjects))
	
	for(i in 1:num.subjects){
		out=get.emission.init.states(likelihood.forward.backward.list[[i]],the.data[[i]]$obs.data,num.obs.states=num.obs.states,num.states=num.states)
		emission_array[,,i]=out$emission.states
		init_array[,i]=out$init.states
	}
	return(list(emission.array=emission_array, init.array=init_array))
}

################################################
get.emission.init.states=function(likelihood.forward.backward,obs.data,num.obs.states,num.states){
############################################################## 
#Author: JL, 11/21/2011
#E-step for initial and emission distributions for an inidividual
#INPUTS:likelihood.forward.backward,obs.data,num.obs.states, num.states
#OUTPUTS:a list with two elements: emission.states(dim (num.states, num.obs.states)) and init.states (dim (num.states))
###############################################################
  lfb=likelihood.forward.backward
  alpha=exp(lfb$logalpha)
  beta=exp(lfb$logbeta)
  likelihood=exp(lfb$LL)

 probs=alpha*beta/likelihood
 emission.out=array(0,dim=c(num.states,num.obs.states))
 init.out=probs[1,]

 #Make indicator vectors for each of the observed states
 for(k in 1:num.obs.states){
 emission.out[,k]=apply(probs*I(obs.data==k),2,"sum")
 }
 return(list(init.states=init.out,emission.states=emission.out))
}

