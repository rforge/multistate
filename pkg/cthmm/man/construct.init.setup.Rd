\name{construct.init.setup}
\alias{construct.init.setup}
\title{Get an initial distribution design setup object}
\usage{
  construct.init.setup(num.states = NULL, states = NULL,
    ref = 1, fixed.dist = NULL, design.matrix = NULL,
    covariate.array = NULL, deriv.array = NULL,
    param.values = NULL, param.types = NULL,
    variable.predictors = NULL)
}
\arguments{
  \item{num.states}{number of hidden states}

  \item{states}{vector of index of states with non-zero
  initial probabilities, not including ref state}

  \item{ref}{index of reference state}

  \item{fixed.dist}{option fixed initial distribution (same
  for all subjects)}

  \item{design.matrix}{array of dim(num beta parameters,
  initial states) that encode the variables/constants used
  in linear predictors of log(emission_ij).  Variables are
  named, and constants are entered as constants.  The rows
  correspond to the transitions in the states entry.}

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
  init.setup object
}
\description{
  Get an initial distribution design setup object
}
\author{
  Jane Lange
}

