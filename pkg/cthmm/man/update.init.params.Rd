\name{update.init.params}
\alias{update.init.params}
\title{Performs Newton Raphson to get MLE for initial dist parameters.}
\usage{
  update.init.params(current.params, init.setup,
    init.counts = NULL, max.it = 1)
}
\arguments{
  \item{init.setup}{initial distributions setup object:
  must have covariate array, param.types}

  \item{init.counts}{matrix
  (dim(num.hidden.states,num.subjects)) with the initial
  dist count data}

  \item{current.params}{initial values of the parameters}

  \item{max.it}{maximum number of Newton Raphson
  iterations. If max.it=1, no updates are performed.}
}
\value{
  list with \item{params}{an array of dim number of
  parameters, max.it} \item{prob.matrix}{array with each
  column representing the vector of initial probabilities
  for a single subject, based on the last value of the
  parameters} \item{states}{hidden states corresponding to
  rows of prob.matrix: will not include ref states or
  states with zero prob!} \item{init.counts}{array with
  (dim=num.states,num.subjects) of the initial distribution
  for each individual.}
}
\description{
  Important: currently only runs NR max.it times (not until
  convergence). Thus this function can be used to obtain
  the initial distributions based on the provided values of
  current.params, if max.it=1 and init.counts=NULL.
}
\author{
  Jane Lange
}

