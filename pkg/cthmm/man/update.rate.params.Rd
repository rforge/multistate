\name{update.rate.params}
\alias{update.rate.params}
\title{Performs Newton Raphson to get MLE for rate parameters.}
\usage{
  update.rate.params(current.params, rate.setup,
    transitions = NULL, durations = NULL, max.it)
}
\arguments{
  \item{rate.setup}{object: must have covariate array,
  param.types}

  \item{transitions}{matrix of
  dim(num.states,num.states,num.subjects) with all i,j
  transtions for all subjects}

  \item{durations}{matrix of dim(num,states, num.subjects)
  with the duration spent in each state}

  \item{current.params}{initial values of the parameters}

  \item{max.it}{maximum number of Newton Raphson
  iterations. If max.it=1, no updates are performed.}
}
\value{
  list with \item{params}{an array of dim number of
  parameters, max.it} \item{rates}{array with each column
  representing the vector of calculated rates for a single
  subject, based on the last value of the parameters}
  \item{transition.codes}{matrix columns "ni" and "nj" i,j
  entries in the rate matrix}
}
\description{
  Important: currently only runs NR max.it times (not until
  convergence). Thus this function can be used to obtain
  the rate matrices based on the provided values of
  current.params, if transitions=NULL, durations=NULL, and
  max.it=1.
}
\author{
  Jane Lange
}

