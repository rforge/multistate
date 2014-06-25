\name{emission_score_louis}
\alias{emission_score_louis}
\title{emission_score_louis
Get the score matrix rate parameters}
\usage{
  emission_score_louis(emission.setup,
    emission.param.values, emission.counts)
}
\arguments{
  \item{emission.setup}{list emission setup object}

  \item{emission.param.values}{values of emission
  parameters}

  \item{emission.counts}{array (dim num states x num obs
  states x num subjects) with the (expected) counts of the
  o_ij}
}
\value{
  score array (dim num params x number of subjects)
}
\description{
  emission_score_louis Get the score matrix rate parameters
}

