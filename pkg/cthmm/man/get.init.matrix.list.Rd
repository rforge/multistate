\name{get.init.matrix.list}
\alias{get.init.matrix.list}
\title{Individual Initial Distributions}
\usage{
  get.init.matrix.list(current.params, init.setup,
    do.list = T)
}
\arguments{
  \item{current.params}{values of emission parameters}

  \item{init.setup}{list with initial distribution setup
  information}

  \item{do.list}{If T, return as list, otherwise return as
  array}
}
\value{
  an array of (dim num.states,num.subjects) or list with
  the fitted initial distribution for each subject.
}
\description{
  Obtain initial distributions for each individual based on
  given parameter values and individual-specific
  covariates.
}
\details{
  Can be output as an array or list.
}
\examples{
\dontrun{
library(cthmm)
data(panel_demodata)
init.list=get.init.matrix.list(current.params=fit.4state$params[12:13],init.setup=init.setup,do.list=T)
}
}
\author{
  Jane Lange
}

