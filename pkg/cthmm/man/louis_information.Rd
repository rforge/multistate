\name{louis_information}
\alias{louis_information}
\title{Get the information matrix using the Louis decomposition method.}
\usage{
  louis_information(rates.setup, init.setup,
    emission.setup, delta.list, emission.list, rates.list,
    M1.list, M2.list, rate.param.values = NULL,
    init.param.values = NULL, emission.param.values = NULL,
    durations, emission.counts, init.counts, num.subjects,
    transitions)
}
\arguments{
  \item{rates.setup}{list rate setup object}

  \item{init.setup}{init setup object}

  \item{emission.setup}{emission setup object}

  \item{init.list}{list with individual initial
  distributions}

  \item{emission.list}{with individual emission
  distributions}

  \item{rates.list}{with individual rate matrices}

  \item{M1.list}{list with conditional means of sufficient
  statistics for each subject}

  \item{M2.list}{list with 2nd moments of sufficient
  statistics for each subject}

  \item{rate.param.values}{values of rate parameters}

  \item{init.param.values}{values of initial distribution
  parameters}

  \item{emission.param.values}{of emission parameters}

  \item{durations}{array (dim num.states x num.subjects)
  with the expected durations}

  \item{emission.counts}{array expected nunber of o_ij 's
  dim (num states x num.obs.states x num.subjects)}

  \item{num.subjects}{number of subjects in the study}
}
\value{
  list with the following components \item{out}{Information
  matrix} \item{complete_info}{the complete data
  information} \item{score_sq}{Conditional expectation of
  the square of complete data information}
  \item{cross_term}{Cross term in the conditional
  expecation of the complete data information}
}
\description{
  Get the information matrix using the Louis decomposition
  method.
}
\author{
  Jane Lange
}

