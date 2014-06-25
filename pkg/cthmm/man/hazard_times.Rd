\name{hazard_times}
\alias{hazard_times}
\title{Latent CTMC hazard functions}
\usage{
  hazard_times(start = 0, end, state_of_interest,
    at_risk_states, alpha, rate, length.out = 1000)
}
\arguments{
  \item{start}{start time for the hazard estimates}

  \item{end}{end time for the hazard estimates}

  \item{state_of_interest}{destination state}

  \item{at_risk_states}{potential states of origin}

  \item{alpha}{starting distribution (vector)}

  \item{rate}{latent CTMC rate matrix to estimate the
  hazard}

  \item{length.out}{number of points between start and end
  times to estimate the hazard}
}
\value{
  \item{hazard}{estimated hazard at different times}
  \item{times}{times at which hazard is estimated}.
}
\description{
  Obtain hazard function based on specified starting states
  and destination state.
}
\examples{
\dontrun{
library(cthmm)
data(panel_demodata)

rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
rate.setup=rates.setup, do.list=T)
rate.firstpassage.AB_0=rates.list[[200]]
rate.firstpassage.AB_0[3:4,]=0
hazardA.B_0=hazard_times(start=.01,end=10,state_of_interest=c(3),
at_risk_states=c(1,2),alpha=c(1,0,0,0),rate.firstpassage.AB_0)

plot(hazardA.B_0$times,hazardA.B_0$haz,type="l",ylab="hazard",xlab="time since entrance in state A",
     main="Hazard function for first passage time from state A to B, for covariate X=0")
}
}
\author{
  Jane Lange
}

