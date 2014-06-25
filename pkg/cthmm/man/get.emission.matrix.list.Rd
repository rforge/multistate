\name{get.emission.matrix.list}
\alias{get.emission.matrix.list}
\title{Individual Emission Matrices}
\usage{
  get.emission.matrix.list(current.params, emission.setup,
    num.states, num.obs.states, num.subjects, do.list = T,
    time.dep.emission = F)
}
\arguments{
  \item{current.params}{values of emission parameters}

  \item{emission.setup}{list with emission distribution
  setup information}

  \item{num.states}{number of latent states in the CTMC}

  \item{num.obs.states}{number of observed states}

  \item{num.subjects}{number of individuals in the study}

  \item{do.list}{If T, return as list, otherwise return as
  array}

  \item{time.dep.emission}{(not supported) Indicator if
  emission distribution has time dependent covariates.}
}
\value{
  an array of (dim num.states,num.subjects) or list with
  the fitted emission distribution for each subject.
}
\description{
  Obtain initiatl distribution matrices for each individual
  based on given parameter values and individual-specific
  covariates.
}
\details{
  Can be output as an array or list.
}
\examples{
\dontrun{
library(cthmm)
data(DDO_data)
emission.list=get.emission.matrix.list(current.params=DDO_EM$params[8:10],
            emission.setup=emission.setup,num.states=4,
           num.obs.states=3, num.subjects=500, do.list=T)
}
}
\author{
  Jane Lange
}

