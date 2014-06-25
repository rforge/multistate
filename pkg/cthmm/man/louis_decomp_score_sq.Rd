\name{louis_decomp_score_sq}
\alias{louis_decomp_score_sq}
\title{Louis decomposition score component}
\usage{
  louis_decomp_score_sq(rates.setup, init.setup,
    emission.setup, init.dist, emission.matrix,
    rate.matrix, M1, M2)
}
\arguments{
  \item{rates.setup}{rate setup object}

  \item{init.setup}{initial distribution setup object}

  \item{emission.setup}{emission distribution setup object}

  \item{rate.matrix}{rate matrix for individual}

  \item{emission.matrix}{emission matrix for individual}

  \item{init.dist}{vector of initial dist for individual}

  \item{M1}{vector of first moments of complete data
  sufficient statistics}

  \item{M2}{matrix of second moments of complete data
  sufficient statistics}
}
\value{
  matrix with second moments of the score function.
}
\description{
  Get E[SS^T|o] for a single subject: i.e. the louis
  decomposition score component
}

