\name{rate_score_louis}
\alias{rate_score_louis}
\title{rate_score_louis
Get the score matrix rate parameters}
\usage{
  rate_score_louis(rates.setup, rate.param.values,
    durations, transitions)

  rate_score_louis(rates.setup, rate.param.values,
    durations, transitions)
}
\arguments{
  \item{rates.setup}{list rate setup object}

  \item{rate.param.values}{values of rate parameters}

  \item{durations}{array (dim num.states x num.subjects)
  with the expected durations}

  \item{transitions}{array (dum num.states x num.states x
  num.subjects) with the expected transition counts}

  \item{rates.setup}{list rate setup object}

  \item{rate.param.values}{values of rate parameters}

  \item{durations}{array (dim num.states x num.subjects)
  with the expected durations}

  \item{transitions}{array (dum num.states x num.states x
  num.subjects) with the expected transition counts}
}
\value{
  score array (dim num params x number of subjects)

  score array (dim num params x number of subjects)
}
\description{
  rate_score_louis Get the score matrix rate parameters

  rate_score_louis Get the score matrix rate parameters
}

