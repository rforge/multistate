\name{get_variance}
\alias{get_variance}
\title{Get the information matrix for the parameter esitmates from the CTMC HMM model
Note: does not actually provide the variance: you need to get the inverse (or pseudo-inverse) to get the variance.}
\usage{
  get_variance(the.data, num.subjects, num.states,
    num.obs.states, rate.param.values = NULL,
    emission.param.values = NULL, init.param.values = NULL,
    absorb.state = NULL, rates.setup, emission.setup,
    init.setup)
}
\arguments{
  \item{the.data}{list with the data}

  \item{num.subjects}{number of subjects in the study}

  \item{num.states}{number of states in the CTMC}

  \item{num.obs.states}{number of observed states}

  \item{rate.param.values}{MLEs for the rate params (if not
  estimated, set to NULL)}

  \item{init.param.values}{MLEs for the initial
  distribution params (if not estimated, set to NULL)}

  \item{emission.param.values}{MLEs for the emission params
  (if not estimated, set to NULL)}

  \item{rates.setup}{list rate setup information}

  \item{init.setup}{list with initial distribution setup
  information}

  \item{emission.setup}{list with emission distribution
  setup information}
}
\value{
  information a list with the following components
  \item{out}{Information matrix} \item{complete_info}{the
  complete data information} \item{score_sq}{Conditional
  expectation of the square of complete data information}
  \item{cross_term}{Cross term in the conditional
  expecation of the complete data information}
}
\description{
  Get the information matrix for the parameter esitmates
  from the CTMC HMM model Note: does not actually provide
  the variance: you need to get the inverse (or
  pseudo-inverse) to get the variance.
}
\author{
  Jane Lange
}

