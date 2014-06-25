\name{get.rate.matrix.list}
\alias{get.rate.matrix.list}
\title{Individual Rate Matrices}
\usage{
  get.rate.matrix.list(current.params, rate.setup,
    do.list = T)
}
\arguments{
  \item{current.params}{values of rate parameters}

  \item{rate.setup}{list with rate matrix setup
  information}

  \item{do.list}{If T, return as list, otherwise return as
  array}
}
\value{
  an array of (dim num.states,num.subjects) or list with
  the fitted rate matrix for each subject.
}
\description{
  Obtain rate matrices for each individual based on given
  parameter values and individual-specific covariates.
}
\details{
  Can be output as an array or list.
}
\examples{
\dontrun{
library(cthmm)
data(DDO_data)
rate.list=get.rate.matrix.list(current.params=DDO_EM$params[1:7],rate.setup=rates.setup,do.list=T)
}
}
\author{
  Jane Lange
}

