\name{EM}
\alias{EM}
\title{EM}
\usage{
  EM(rates.setup, init.setup, emission.setup, the.data,
    num.subjects, num.states, num.obs.states, tol = 1e-07,
    absorb.state, maxiter = 500, DDO.setup = NULL,
    do.DDO = F, time.dep.emission = F)
}
\arguments{
  \item{rates.setup}{list with rate setup information}

  \item{init.setup}{list with initial distribution setup
  information}

  \item{emission.setup}{list with emission distribution
  setup information}

  \item{the.data}{list with the observed data, with one
  entry per individual}

  \item{num.subjects}{number of individuals in the study}

  \item{num.states}{number of latent states in the CTMC}

  \item{num.obs.states}{number of observed states}

  \item{tol}{convergence tolerance for the EM}

  \item{absorb.state}{vector of one or more absorbing
  latent states}

  \item{maxiter=500}{maximum number of iterations of the
  EM}

  \item{DDO.setup}{setup object for disease driven
  observation model}

  \item{do.DDO}{Indicator (T/F) if there are disease driven
  observation times in the model?}

  \item{time.dep.emission}{(Not supported) Indicator if
  emission distribution has time dependent covariates.}
}
\value{
  \item{params}{MLEs for model parameters:
  rates/emission/initial/DDO, parameters listed in order of
  model specification} \item{details}{Details of the EM
  optimization from the \code{SQUAREM} package output}.
  \item{updates}{The iterative oupdates for the parmaeter
  esitmates at each EM step} \item{LL.updates}{Updates of
  observed data likelihood at each step} \item{time}{Run
  time} \item{num.evals}{Number of iterations of the EM}
  \item{LL}{Final complete data log likelihood}
}
\description{
  EM algorithm to get MLEs for latent CTMC parameters
}
\details{
  Run the expectation-maximization algorithm to get the
  maximum likelihood estimates for latent CTMC model
  parameters for multistate disease models. Data are
  discretely observed disease processes from multiple
  independent individuals. Times may be non-informative
  scheduled times or disease-dependent (informative) times.
  Uses the SQUAREM package for EM acceleration.
}
\examples{
\dontrun{
library(cthmm)
#load the model setup and example data for a competing risks model
#observed at scheduled and informative observation times.
data(DDO_data)
#run the EM on the example data
DDO_EM=EM(rates.setup=rates.setup,
init.setup=init.setup,
emission.setup=emission.setup,
the.data=DDO_data,
num.subjects=500,
num.states=4,
num.obs.states=3,
tol = 1e-07,
absorb.state=c(3,4),
maxiter = 500,
DDO.setup = DDO.setup,
do.DDO = T)}
}
\author{
  Jane Lange
}

