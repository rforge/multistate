\name{phase.hazard}
\alias{phase.hazard}
\title{Hazard of reaching absorbing state}
\usage{
  phase.hazard(time, initial.dist, rate.matrix,
    absorb.state)
}
\arguments{
  \item{time}{evaluation time}

  \item{initial.dist}{initial state occupancy distribution}

  \item{rate.matrix}{rate matrix}

  \item{absorb.state}{index of absorbing state}
}
\value{
  h(t|initial dist)
}
\description{
  Get h(t|inital dist) in phase type model based on given
  initial distribution and given time. May need to be
  modified for multiple absorbing states!
}
\author{
  Jane Lange
}

