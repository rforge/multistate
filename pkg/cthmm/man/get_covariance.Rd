\name{get_covariance}
\alias{get_covariance}
\title{Covariance estimate for latent CTMC parameters}
\usage{
  get_covariance(par, the.data, num.subjects, num.states,
    num.obs.states, rates.setup, emission.setup = NULL,
    init.setup = NULL, DDO.setup = NULL, do.DDO = F,
    time.dep.emission = F)
}
\arguments{
  \item{par}{the parameter estimates (rates, emission,
  initial distribution, DDO (in that order)}

  \item{the.data}{list with the observed data, with one
  entry per individual}

  \item{num.subjects}{number of individuals in the study}

  \item{num.states}{number of latent states in the CTMC}

  \item{num.obs.states}{number of observed states}

  \item{rates.setup}{list with rate setup information}

  \item{emission.setup}{list with emission distribution
  setup information}

  \item{init.setup}{list with initial distribution setup
  information}

  \item{DDO.setup}{setup object for disease driven
  observation model}

  \item{do.DDO}{indicator (T/F) if there are disease driven
  observation times in the model}

  \item{time.dep.emission}{(not supported) Indicator if
  emission distribution has time dependent covariates.}
}
\value{
  \item{covariance}{Estimated covariance of parameter
  estiamtes} \item{information}{Observed Fisher information
  of parameter estimates}.
}
\description{
  Obtain the covariance matrix for parameter estimates from
  the latent CTMC model.
}
\details{
  Run after the EM algorithm that provides maximum
  likelihood estimates. Evaluates the hessian of the
  observed data likelihood, using numeric differentiation,
}
\examples{
\dontrun{
library(cthmm)
data(DDO_data)
#run the EM on the example data
covariance_DDO=get_covariance(par=DDO_EM$param,
the.data=DDO_data,
num.subjects=500,
num.states=4,
num.obs.states=3,
rates.setup=rates.setup,
emission.setup=emission.setup,
init.setup=init.setup,
DDO.setup=DDO.setup,
do.DDO=T)}
}
\author{
  Jane Lange
}

