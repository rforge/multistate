\name{get_M1_M2}
\alias{get_M1_M2}
\title{get_M1_M2
get the list with first and second moments of complete data sufficient statistics for an individual, using the recursive smoothing method.}
\usage{
  get_M1_M2(stat.table, time.diffs, obs.data, eigen.decomp,
    exact.time.ranks, absorb.state, transition.probs,
    emission.matrix, init.dist)
}
\arguments{
  \item{stat.table}{data frame with num rows=number of
  sufficient stats, and columns ni nj di zi Ni oi oj. If
  the statistic is n_12, then the corresponding column is
  (1,2,0,0,0,0,0)}

  \item{time.diffs}{observation time intervals (single
  subject)}

  \item{obs.data}{observed data (single subject)}

  \item{eigen.decomp}{eigen.decomp object for rate matrix
  (single subject)}

  \item{exact.time.ranks}{number encoding the absorption
  time index}

  \item{absorb.state}{one or more absorbing states in the
  model}

  \item{transition.probs}{transition probabilities for each
  time interval}

  \item{emission.matrix}{emission matrix for subject}

  \item{init.dist}{initial distribution for sujbect}
}
\value{
  list consisting of \item{M1}{vector of first moments of
  complete data sufficient statistics} \item{M2}{matrix of
  second moments of complete data sufficient statistics}
}
\description{
  get_M1_M2 get the list with first and second moments of
  complete data sufficient statistics for an individual,
  using the recursive smoothing method.
}

