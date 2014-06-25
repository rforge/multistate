\name{se.haz.times}
\alias{se.haz.times}
\title{Hazard Standard Errors}
\usage{
  se.haz.times(start, end, covar, rate, at_risk_states,
    state_of_interest, alpha, param.deriv, num.transitions,
    num.params, transitions, num.states, length.out = 100)
}
\arguments{
  \item{start}{start time for the hazard estimates}

  \item{end}{end time for the hazard estimates}

  \item{covar}{covariance matrix for the rate estimates
  only}

  \item{rate}{latent CTMC rate matrix to estimate the
  hazard}

  \item{at_risk_states}{potential states of origin}

  \item{state_of_interest}{destination state}

  \item{alpha}{starting distribution (vector)}

  \item{param.deriv}{derivative array for the individual}

  \item{num.transitions}{number of transtions between
  states in the latent CTMC model with non-zero rates}

  \item{num.params}{number of rate parameters in the latent
  CTMC rate model}

  \item{transitions}{transition codes for the latent CTMC
  model}

  \item{num.states}{number of states in the latent CTMC
  model}

  \item{length.out}{number of points between start and end
  times to estimate the standard error}
}
\value{
  standard errors for the subdistribution function
}
\description{
  Get standard errors for the hazard function estimates
}
\examples{
\dontrun{
library(cthmm)
data(DDO_data)

rate.list=get.rate.matrix.list(current.params=DDO_EM$params[1:7],rate.setup=rates.setup,do.list=T)

se=se.haz.times(start=0, end=10,
covar=covariance_DDO$covariance[1:7,1:7],
rate=rates.list[[1]],
at_risk_states=c(1,2),
state_of_interest=3,
alpha=c(1,0,0,0),
param.deriv=rates.setup$deriv.array[,,1],
num.transitions=dim(rates.setup$transition.codes)[1],
num.params=7,
transitions=rates.setup$transition.codes,
num.states=4,
length.out = 1000)
}
}
\author{
  Jane Lange
}

