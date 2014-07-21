#' @include Numeric_hessian.R
NULL
#' @title Covariance estimate for latent CTMC parameters
#'
#' @description Obtain the covariance matrix for parameter estimates from the latent CTMC model. 
#'
#' @details Run after the EM algorithm that provides maximum likelihood estimates. Evaluates the hessian of the observed data likelihood, using numeric differentiation, 
#' 
#' @param par the parameter estimates (rates, emission, initial distribution, DDO (in that order)
#' @param the.data list with the observed data, with one entry per individual 
#' @param num.subjects number of individuals in the study
#' @param num.states number of latent states in the CTMC
#' @param num.obs.states number of observed states
#' @param rates.setup list with rate setup information
#' @param emission.setup list with emission distribution setup information
#' @param init.setup list with initial distribution setup information
#' @param DDO.setup setup object for disease driven observation model
#' @param do.DDO  indicator (T/F) if there are disease driven observation times in the model
#' @param time.dep.emission (not supported) Indicator if emission distribution has time dependent covariates.
#' @return 
#' \item{covariance}{Estimated covariance of parameter estiamtes}
#' \item{information}{Observed Fisher information of parameter estimates}.
#' @examples
#'\dontrun{
#' library(cthmm)
#' data(DDO_data)
#' #run the EM on the example data
#' covariance_DDO=get_covariance(par=DDO_EM$param,
#' the.data=DDO_data,
#' num.subjects=500,
#' num.states=4,
#' num.obs.states=3,
#' rates.setup=rates.setup,
#' emission.setup=emission.setup,
#' init.setup=init.setup,
#' DDO.setup=DDO.setup,
#' do.DDO=T)}
#' @author Jane Lange
#' @export
#################################################
#Function to evaluate the likelihood
get_covariance<-function(par,
                               the.data,
                               num.subjects,
                               num.states,
                               num.obs.states,
                               rates.setup,
                               emission.setup=NULL,
                               init.setup=NULL,
                               DDO.setup=NULL,
                               do.DDO=F,
                               time.dep.emission=F){
  the.data<-the.data
  num.subjects<-num.subjects
  num.states<-num.states
  num.obs.states<-num.obs.states
  rates.setup<-rates.setup
  emission.setup<-emission.setup
  init.setup<-init.setup
  DDO.setup<-DDO.setup
  do.DDO<-do.DDO
  time.dep.emission<-time.dep.emission
  
  get_LL_reduced_version<-function(par){  
    Q = NULL
    h.list = NULL
    rates <- rates.setup
    init <- init.setup
    emission <- emission.setup
    DDO <- DDO.setup
    ij.indices = matrix(unlist(rates$transition.codes), ncol = 2, 
                        byrow = F)
    colnames(ij.indices) = c("i", "j")
    ij.indices <- ij.indices
    obs.data.list <- lapply(the.data, "[[", c("obs.data"))
    exact.time.ranks.list <- lapply(the.data, "[[", c("exact.times.ranks"))
    time.diffs.list <- lapply(lapply(the.data, "[[", "obs.times"), 
                              FUN = "diff")
    if (do.DDO) {
      h.list = lapply(the.data, "[[", c("h"))
    }
    rep_item <- function(i, item) {
      return(item)
    }
   # par = c(0)
    if (!is.null(rates$fixed.rates)) {
      rates.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                          item = rates$fixed.rates)
      num.rate.params <- 0
      do.rates = F
    }else {
    #  par = c(par, rates$param.values)
      num.rate.params <- length(rates$param.values)
      do.rates = T
    }
    if (!is.null(emission$fixed.dist)) {
      emission.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                             item = emission$fixed.dist)
      num.emission.params <- 0
      do.emission = F
    }else {
     # par = c(par, emission$param.values)
      num.emission.params = length(emission$param.values)
      do.emission = T
    }
    if (!is.null(init$fixed.dist)) {
      delta.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                          item = init$fixed.dist)
      num.init.params = 0
      do.init = F
    }else {
      #par = c(par, init$param.values)
      num.init.params <- length(init$param.values)
      do.init = T
    }
    if (do.DDO) {
      if (!is.null(DDO$fixed.rates)) {
        Q = lapply(seq(1:num.subjects), FUN = "rep_item", 
                   item = DDO$fixed.rates)
        fixed.DDO = T
        num.DDO.params = 0
      }else {
       # par = c(par, DDO$param.values)
        num.DDO.params <- length(DDO$param.values)
        fixed.DDO = F
      }
    }
#    par = par[-1]
    param.updates = par
    LL.updates = 0
    start = proc.time()
    
    if (do.rates) {
      rates.list = get.rate.matrix.list(current.params = par[1:num.rate.params], 
                                        rate.setup = rates, do.list = T)
    }
    if (do.emission) {
      emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 
                                                                       1):(num.rate.params + num.emission.params)], 
                                               emission.setup = emission, num.states = num.states, 
                                               num.obs.states = num.obs.states, num.subjects = num.subjects, 
                                               time.dep.emission = time.dep.emission)
    }
    if (do.init) {
      delta.list = get.init.matrix.list(current.params = par[(num.rate.params + 
                                                                num.emission.params + 1):(num.rate.params + num.emission.params + 
                                                                                            num.init.params)], init.setup = init)
    }
    if (do.DDO) {
      if (!fixed.DDO) {
        Q = get.rate.matrix.list(current.params = par[(num.rate.params + 
                                                         num.emission.params + num.init.params + 1):(num.rate.params + 
                                                                                                       num.emission.params + num.init.params + num.DDO.params)], 
                                 rate.setup = DDO, do.list = T)
      }
      rates.list = mapply("+", rates.list, Q, SIMPLIFY = F)
    }
    eigen.decomp.list <- eigen_decomp_list_R(rates.list)
    transition.probabilities.list = trans_probs_all(eigen.decomp.list, 
                                                    time.diffs.list, exact.time.ranks.list, h.list, Q)
    likelihood.forward.backward.list <- mapply(obs.data.list, 
                                               transition.probabilities.list, delta.list, emission.list, 
                                               FUN = "forwardback_R", MoreArgs = list(time.dep = time.dep.emission), 
                                               SIMPLIFY = F)
    LL.vec=(sapply(likelihood.forward.backward.list,"[[","LL"))
    LL=sum(LL.vec[LL.vec!=-Inf&LL.vec!="NaN"]) 
    LL.updates=sum(sapply(likelihood.forward.backward.list,"[[","LL"))
    return(LL)
    
  }
  information=-1*hessian(func=get_LL_reduced_version,x=par)
  covariance=pseudoinverse(information)
								   
  return(list(covariance=covariance, information=information))
  
}
                               
                               
get_LL<-function(par,
                 the.data,
                 num.subjects,
                 num.states,
                 num.obs.states,
                 rates.setup=NULL,
                 emission.setup=NULL,
                 init.setup=NULL,
                 DDO.setup=NULL,
                 do.DDO=F,
                 time.dep.emission=NULL){  
  Q = NULL
  h.list = NULL
  rates <- rates.setup
  init <- init.setup
  emission <- emission.setup
  DDO <- DDO.setup
  ij.indices = matrix(unlist(rates$transition.codes), ncol = 2, 
                      byrow = F)
  colnames(ij.indices) = c("i", "j")
  ij.indices <- ij.indices
  obs.data.list <- lapply(the.data, "[[", c("obs.data"))
  exact.time.ranks.list <- lapply(the.data, "[[", c("exact.times.ranks"))
  time.diffs.list <- lapply(lapply(the.data, "[[", "obs.times"), 
                            FUN = "diff")
  if (do.DDO) {
    h.list = lapply(the.data, "[[", c("h"))
  }
  rep_item <- function(i, item) {
    return(item)
  }
  #par = c(0)
  if (!is.null(rates$fixed.rates)) {
    rates.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                        item = rates$fixed.rates)
    num.rate.params <- 0
    do.rates = F
  }else {
  #  par = c(par, rates$param.values)
    num.rate.params <- length(rates$param.values)
    do.rates = T
  }
  if (!is.null(emission$fixed.dist)) {
    emission.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                           item = emission$fixed.dist)
    num.emission.params <- 0
    do.emission = F
  }else {
   # par = c(par, emission$param.values)
    num.emission.params = length(emission$param.values)
    do.emission = T
  }
  if (!is.null(init$fixed.dist)) {
    delta.list = lapply(seq(1:num.subjects), FUN = "rep_item", 
                        item = init$fixed.dist)
    num.init.params = 0
    do.init = F
  }else {
  #  par = c(par, init$param.values)
    num.init.params <- length(init$param.values)
    do.init = T
  }
  if (do.DDO) {
    if (!is.null(DDO$fixed.rates)) {
      Q = lapply(seq(1:num.subjects), FUN = "rep_item", 
                 item = DDO$fixed.rates)
      fixed.DDO = T
      num.DDO.params = 0
    }else {
   #   par = c(par, DDO$param.values)
      num.DDO.params <- length(DDO$param.values)
      fixed.DDO = F
    }
  }
  #par = par[-1]
  param.updates = par
  LL.updates = 0
  start = proc.time()
  
  if (do.rates) {
    rates.list = get.rate.matrix.list(current.params = par[1:num.rate.params], 
                                      rate.setup = rates, do.list = T)
  }
  if (do.emission) {
    emission.list = get.emission.matrix.list(current.params = par[(num.rate.params + 
                                                                     1):(num.rate.params + num.emission.params)], 
                                             emission.setup = emission, num.states = num.states, 
                                             num.obs.states = num.obs.states, num.subjects = num.subjects, 
                                             time.dep.emission = time.dep.emission)
  }
  if (do.init) {
    delta.list = get.init.matrix.list(current.params = par[(num.rate.params + 
                                                              num.emission.params + 1):(num.rate.params + num.emission.params + 
                                                                                          num.init.params)], init.setup = init)
  }
  if (do.DDO) {
    if (!fixed.DDO) {
      Q = get.rate.matrix.list(current.params = par[(num.rate.params + 
                                                       num.emission.params + num.init.params + 1):(num.rate.params + 
                                                                                                     num.emission.params + num.init.params + num.DDO.params)], 
                               rate.setup = DDO, do.list = T)
    }
    rates.list = mapply("+", rates.list, Q, SIMPLIFY = F)
  }
  eigen.decomp.list <- eigen_decomp_list_R(rates.list)
  transition.probabilities.list = trans_probs_all(eigen.decomp.list, 
                                                  time.diffs.list, exact.time.ranks.list, h.list, Q)
  likelihood.forward.backward.list <- mapply(obs.data.list, 
                                             transition.probabilities.list, delta.list, emission.list, 
                                             FUN = "forwardback_R", MoreArgs = list(time.dep = time.dep.emission), 
                                             SIMPLIFY = F)
   LL.vec=(sapply(likelihood.forward.backward.list,"[[","LL"))
  LL=sum(LL.vec[LL.vec!=-Inf&LL.vec!="NaN"]) 
  LL.updates=sum(sapply(likelihood.forward.backward.list,"[[","LL"))
  return(LL)
}
