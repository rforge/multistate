\name{get_M1_M2_list}
\alias{get_M1_M2_list}
\title{get_M1_M2_list
Get a list of first and second moments for each subject}
\usage{
  get_M1_M2_list(rates.setup, init.setup, emission.setup,
    delta.list, rates.list, emission.list, obs.data.list,
    eigen.decomp.list, time.diffs.list,
    exact.time.ranks.list, transition.probabilties.list,
    absorb.state)
}
\arguments{
  \item{rates.setup}{rate setup object}

  \item{init.setup}{init setup object}

  \item{emission.setup}{emission setup object}

  \item{delta.list}{list with intial distribution for each
  subject}

  \item{rates.list}{list with rates for each subject}

  \item{emission.list}{list with emission matrices for each
  subject}

  \item{obs.data.list}{list with obs.data for each subject}

  \item{eigen.decomp.list}{list with eigen.decomp objects
  for each subject}

  \item{time.diffs.list}{list with interobservation
  intervals for each subject}

  \item{exact.time.ranks.list}{list encoding the
  observation number of times of absorption for each
  subject}

  \item{transition.probabilities.list}{list encoding the
  transition probabilities for each subject}

  \item{absorb.state}{vector with one or more absorbing
  states}
}
\value{
  a list consisting of \item{M1_list}{list with vectors of
  first moments for complete data sufficient statistics}
  \item{M2_list}{list with matrices of second moments for
  complete data sufficient statistics}
}
\description{
  get_M1_M2_list Get a list of first and second moments for
  each subject
}

