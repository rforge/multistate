########################################################
#Helper functions in getting the rate.setup, init.setup, and emission.setup
############################################################
#Documentation: get_covariate_array, construct.init.setup, construct. emission.setup, construct.rates.setup
###############################################################
NULL
#' Get the covariate arrays
#'
#' Creates matrices with covariate values for all subjects for use in parameterized rate, initial dist., and emission matrices. 
#'
#' @param design.matrix a matrix that codes the constants and variables for the parameterization
#'                       with rows corresponding to theta and columns corresponding to betas
#' @param cov.data a data frame with baseline covariate data for each individual
#' @param variable.entries matrix with rows for each parameter corresponding to variable covariates, with the names for
#'                      each.
#'  @return an array (dim=c(num.betas,num.thetas, num.subjects)) with covariates values for all subjects
#' @export
#' @author Jane Lange
get.covariate.array<-function(design.matrix,cov.data,variable.entries=NULL){
#####################################################################################
#Author: JL, 11/15/2011
# This function creates an array with matrix entries of covariate values for all subjects
#INPUTS: design.matrix: matrix that codes the constants and variables for the parameterization
#                       with rows corresponding to theta and columns corresponding to betas
#        variable.entries=matrix with rows for each parameter corresponding to variable covariates, with the names for
#                      each.  
#OUTPUTS: an array (dim=c(num.betas,num.thetas, num.subjects)) with covariates values for all subjects
#########################################################################
	covariate.array=array(0,dim=c(dim(design.matrix)[2],dim(design.matrix)[1],dim(cov.data)[1]))
	if(!is.null(variable.entries)){
		for(i in 1:dim(cov.data)[1]){
			covariate.array[,,i]=fill.variables(i,cov.data,design.matrix,variable.entries)
		}
	}
	for(i in 1:dim(design.matrix)[1]){
		for(j in 1:dim(design.matrix)[2]){
			if(is.numeric(design.matrix[i,j])&!is.na(design.matrix[i,j])){
				print(c(i,j))
				covariate.array[j,i,]=design.matrix[i,j]
			}  
		}
	}
	return(covariate.array)
}

NULL
#'Get deriv arrays based on the non-fixed parameters
#'@param param.types (not fixed if=0)
#'@param covariate.array
get.deriv.array<-function(covariate.array, param.types){
deriv.array=array(covariate.array[param.types==0,,],dim=c(sum(param.types==0),dim(covariate.array)[2],dim(covariate.array)[3]))
return(deriv.array)
}




NULL
#' Get an initial distribution design setup object
#' 
#' @param num.states number of hidden states
#' @param states vector of index of states with non-zero initial probabilities, not including ref state
#' @param ref index of reference state
#' @param fixed.dist option fixed initial distribution (same for all subjects)
#' @param design.matrix array of dim(num beta parameters, initial states) that encode the variables/constants used in linear predictors of log(emission_ij).  Variables are named, and constants are entered as constants.  The rows correspond to the transitions in the states entry.  
#' @param covariate.array (dim=c(num.betas,num.thetas, num.subjects)) with covariates values for all subjects
#' @param deriv.array
#' @param param.values inital or fixed values of beta parameters
#' @param param.types a vector of length(num.beta.parameters), with entries=0 if parameter is free, 1 if fixed at param.value.
#' @param variable.predictors matrix with columns "i" and "j" that encode the variable predictors in the design.matrix
#' @return init.setup object
#' @export
#' @author Jane Lange
construct.init.setup<-function(num.states=NULL,states=NULL,ref=1, fixed.dist=NULL, design.matrix=NULL, covariate.array=NULL,deriv.array=NULL, param.values=NULL,param.types=NULL,variable.predictors=NULL){
	init.setup=list(num.states=NULL,states=c(1,3),ref=1, fixed.dist=NULL, design.matrix=NULL, covariate.array=NULL, param.values=NULL,param.types=NULL,variable.predictors=NULL)
	class(init.setup)<-"initSetup" 
}
NULL
#' Get an emission distribution design setup object
#' 
#' The design matrix needs to all i,j entries in the emission matrix that are non-zero, except for the reference state for row i. 
#' @param emission.states number of hidden states
#' @param ref.states vector of index of states with non-zero initial probabilities
#' @param exact.states index of reference state
#' @param fixed.dist option emission matrix set at fixed value for all individuals. 
#' @param design.matrix array of dim(num beta parameters, num.emissions) that encodes the variables/constants used in linear predictors of log(emission_ij).  Variables are named, and constants are entered as constants.  The rows correspond to the transitions in the transition.codes entry.  
#' @param covariate.array (dim=c(num.betas,num.thetas, num.subjects)) with covariates values for all subjects
#' @param deriv.array 
#' @param param.values inital or fixed values of beta parameters
#' @param param.types a vector of length(num.beta.parameters), with entries=0 if parameter is free, 1 if fixed at param.value.
#' @param variable.predictors matrix with columns "i" and "j" that encode the variable predictors in the design.matrix
#' @return emission.setup object (a list)
#' @author Jane Lange
construct.emission.setup<-function(emission.states=NULL,ref.states=NULL, exact.states=NULL, fixed.dist=NULL, design.matrix=NULL, covariate.array=NULL, deriv.array=NULL, param.values=NULL,param.types=NULL,variable.predictors=NULL){
	emission.setup=list(emission.states=NULL,ref.states=NULL, exact.states=NULL, fixed.dist=NULL, design.matrix=NULL, covariate.array=NULL, param.values=NULL,param.types=NULL,variable.predictors=NULL)
	class(emission.setup)<-"emissionSetup" 
}
NULL  
#' Get a rate matrix design setup object
#' 
#' Rates are defined as log(theta_{ij})=W^T Beta, where W are predictors for rate entry i,j
#' @param num.states number of hidden states
#' @param transition.codes matrix with columns "ni" and "nj" for each non-zero i,j transition intensity
#' @param fixed.rates optional rate matrix with fixed entries (same for all individuals)
#' @param design.matrix array of dim(num beta parameters, num.non.zero.transitions) that encodes the variables/constants used in linear predictors of log(rates).  Variables are named, and constants are entered as constants.  The rows correspond to the transitions in the transition.codes entry.  
#' @param covariate.array (dim=c(num.betas,num.thetas, num.subjects)) with covariates values for all subjects
#' @param deriv.array
#' @param param.values inital or fixed values of beta parameters
#' @param param.types a vector of length(num.beta.parameters), with entries=0 if parameter is free, 1 if fixed at param.value.
#' @param variable.predictors matrix with columns "i" and "j" that encode the variable predictors in the design.matrix
#' @return rate.setup object
#' @export
#' @author Jane Lange
construct.rates.setup<-function(num.states=NULL, transition.codes=NULL, fixed.rates=NULL, design.matrix=NULL, covariate.array=NULL, deriv.array=NULL,param.values=NULL, param.types=NULL,variable.predictors=NULL){
	rates.setup=list(num.states=NULL, transition.codes=NULL, fixed.rates=NULL, design.matrix=NULL, covariate.array=NULL, param.values=NULL, param.types=NULL,variable.predictors=NULL)
	class(rates.setup)<-"ratesSetup" 
}

fill.variables<-function(row,cov.data,design.matrix,variable.entries){
###################################################################
#Author: JL, 11/15/2011
#This function fills the covariate.matrix with a single subject's variable values
#INPUTS: row=row number for subject in cov.data
#        cov.data=data frame with baseline covariate data by subject
#        variable.entries=matrix with rows for each parameter corresponding to variable covariates, with the names for
#                      each.  
#OUTPUTS: 
# matrix with covariate entries: rows correspond to beta parameters and columns correspond to the theta parameters
####################################################################################
  cov.data.row=cov.data[row,]
  matrix.out=matrix(0,nrow=dim(design.matrix)[1],ncol=dim(design.matrix)[2])  
  for(k in 1:dim(variable.entries)[1]){
   matrix.out[variable.entries[k,"i"],variable.entries[k,"j"]]=as.numeric(cov.data.row[paste(design.matrix[variable.entries[k,"i"],variable.entries[k,"j"]])])
    }
      
   return(t(matrix.out))
}

