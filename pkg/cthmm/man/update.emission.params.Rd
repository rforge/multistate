\name{update.emission.params}
\alias{update.emission.params}
\title{Performs Newton Raphson to get MLE for emission parameters.}
\usage{
  update.emission.params(current.params, emission.setup,
    emission.counts = NULL, max.it = 1)
}
\arguments{
  \item{current.params}{initial values of the parameters}

  \item{emission.setup}{object with covariate array,
  param.types}

  \item{emission.counts}{matrix
  (dim(num.hidden.states,num.obs.states,num.subjects)) with
  the emission count data}

  \item{max.it}{maximum number of Newton Raphson
  iterations. If max.it=1, no updates are performed.}
}
\value{
  list with \item{params}{an array of dim number of
  parameters, max.it} \item{prob.matrix}{array with each
  column representing the vector of calculation emission
  probabilities for a single subject, based on the last
  value of the parameters} \item{states}{matrix columns "i"
  and "j" representing hidden(i) and observed(j) states
  corresponding to the rows of prob.matrix}
}
\description{
  Important: currently only runs NR max.it times (not until
  convergence). Thus this function can be used to obtain
  the emission probabilties based on the provided values of
  current.params, if max.it=1 and emission.counts=NULL.
}
\author{
  Jane Lange
}

