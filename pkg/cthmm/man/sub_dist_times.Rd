\name{sub_dist_times}
\alias{sub_dist_times}
\title{Latent CTMC subdistribution functions}
\usage{
  sub_dist_times(start = 0, end, states, alpha, rate,
    length.out = 1000)
}
\arguments{
  \item{start}{start time for the hazard estimates}

  \item{end}{end time for the hazard estimates}

  \item{states}{destination state(s)}

  \item{alpha}{starting distribution (vector)}

  \item{rate}{latent CTMC rate matrix to estimate the
  hazard}

  \item{length.out}{number of points between start and end
  times to estimate the subdistribution function.}
}
\value{
  \item{dist}{estimated subdistribution at different times}
  \item{times}{times at which subdistribution function is
  esitmated}.
}
\description{
  Obtain subdistribution functions or transition
  probabilities for specified starting and destination
  states.
}
\examples{
\dontrun{
library(cthmm)
data(panel_demodata)

rates.list=get.rate.matrix.list(current.params=fit.4state$param[1:8],
rate.setup=rates.setup, do.list=T)
rate.firstpassage.AB_0=rates.list[[200]]
rate.firstpassage.AB_0[3:4,]=0
cdfA.B_0=sub_dist_times(start=.01,end=10,states=3,
alpha=c(1,0,0,0),rate.firstpassage.AB_0)
plot(cdfA.B_0$times,cdfA.B_0$dist,type="l",ylab="CDF",xlab="time since entrance in state A",
     main="Cumulative distribution function for first passage time from state A to B, for covariate X=0")
}
}
\author{
  Jane Lange
}

