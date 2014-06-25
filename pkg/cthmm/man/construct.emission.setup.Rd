\name{construct.emission.setup}
\alias{construct.emission.setup}
\title{Get an emission distribution design setup object}
\usage{
  construct.emission.setup(emission.states = NULL,
    ref.states = NULL, exact.states = NULL,
    fixed.dist = NULL, design.matrix = NULL,
    covariate.array = NULL, deriv.array = NULL,
    param.values = NULL, param.types = NULL,
    variable.predictors = NULL)
}
\arguments{
  \item{emission.states}{number of hidden states}

  \item{ref.states}{vector of index of states with non-zero
  initial probabilities}

  \item{exact.states}{index of reference state}

  \item{fixed.dist}{option emission matrix set at fixed
  value for all individuals.}

  \item{design.matrix}{array of dim(num beta parameters,
  num.emissions) that encodes the variables/constants used
  in linear predictors of log(emission_ij).  Variables are
  named, and constants are entered as constants.  The rows
  correspond to the transitions in the transition.codes
  entry.}

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
  emission.setup object (a list)
}
\description{
  The design matrix needs to all i,j entries in the
  emission matrix that are non-zero, except for the
  reference state for row i.
}
\author{
  Jane Lange
}

