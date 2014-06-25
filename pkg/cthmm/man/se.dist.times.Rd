\name{se.dist.times}
\alias{se.dist.times}
\title{Subdistribution Standard Errors}
\usage{
  se.dist.times(start, end, covar, rate, states, alpha,
    param.deriv, num.transitions, num.params, transitions,
    num.states, length.out = 100)
}
\arguments{
  \item{start}{start time for the hazard estimates}

  \item{end}{end time for the hazard estimates}

  \item{covar}{ovariance matrix for the rate estimates
  only}

  \item{rate}{latent CTMC rate matrix to estimate the
  hazard}

  \item{states}{destination state(s)}

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
  Get standard errors for the subdistribution function
  estimates
}
\examples{
\dontrun{
library(cthmm)
data(panel_demodata)

rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
           rate.setup=rates.setup, do.list=T)

rate.firstpassage.AB_0=rates.list[[200]]
rate.firstpassage.AB_0[3:4,]=0

covar=covariance_4state$covariance[1:8,1:8]

num.transitions=dim(rates.setup$transition.codes)[1]

num.params=8

transitions=rates.setup$transition.codes

num.states=4

param.deriv=rates.setup$deriv.array[,,1]

se_AB_0=se.dist.times(start=.01,end=10,
          covar[1:8,1:8],
          rate.firstpassage.AB_0,
          states=3,alpha=c(1,0,0,0),
          param.deriv,num.transitions,
          num.params,transitions,
           num.states,length.out=1000)
}
}
\author{
  Jane Lange
}

