\name{phase.surv}
\alias{phase.surv}
\title{Survivial in phase type model}
\usage{
  phase.surv(time, initial.dist, rate.matrix, absorb.state)
}
\arguments{
  \item{time}{evaluation time}

  \item{initial.dist}{initial state occupancy distribution}

  \item{rate.matrix}{rate matrix}

  \item{absorb.state}{index of absorbing state}
}
\value{
  S(t|initial dist)
}
\description{
  Get S(t|inital dist) in phase type model based on given
  initial distribution and given time. May need to be
  modified for multiple absorbing states!
}
\author{
  Jane Lange
}

