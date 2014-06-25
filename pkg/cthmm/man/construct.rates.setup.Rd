\name{construct.rates.setup}
\alias{construct.rates.setup}
\title{Get a rate matrix design setup object}
\usage{
  construct.rates.setup(num.states = NULL,
    transition.codes = NULL, fixed.rates = NULL,
    design.matrix = NULL, covariate.array = NULL,
    deriv.array = NULL, param.values = NULL,
    param.types = NULL, variable.predictors = NULL)
}
\arguments{
  \item{num.states}{number of hidden states}

  \item{transition.codes}{matrix with columns "ni" and "nj"
  for each non-zero i,j transition intensity}

  \item{fixed.rates}{optional rate matrix with fixed
  entries (same for all individuals)}

  \item{design.matrix}{array of dim(num beta parameters,
  num.non.zero.transitions) that encodes the
  variables/constants used in linear predictors of
  log(rates).  Variables are named, and constants are
  entered as constants.  The rows correspond to the
  transitions in the transition.codes entry.}

  \item{covariate.array}{(dim=c(num.betas,num.thetas,
  num.subjects)) with covariates values for all subjects}

  \item{deriv.array}{}

  \item{param.values}{inital or fixed values of beta
  parameters}

  \item{param.types}{a vector of
  length(num.beta.parameters), with entries=0 if parameter
  is free, 1 if fixed at param.value.}

  \item{variable.predictors}{matrix with columns "i" and
  "j" that encode the variable predictors in the
  design.matrix}
}
\value{
  rate.setup object
}
\description{
  Rates are defined as log(theta_{ij})=W^T Beta, where W
  are predictors for rate entry i,j
}
\author{
  Jane Lange
}

