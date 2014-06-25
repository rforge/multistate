\name{get.suff.stat.table}
\alias{get.suff.stat.table}
\title{get.suff.stat.table}
\usage{
  get.suff.stat.table(rates.setup, init.setup,
    emission.setup)
}
\arguments{
  \item{rates.setup}{list rate setup object}

  \item{init.setup}{init setup object}

  \item{emission.setup}{emission setup object}
}
\value{
  a data frame with num rows=number of sufficient stats,
  and columns ni nj di zi Ni oi oj. If the statistic is
  n_12, then the corresponding column is (1,2,0,0,0,0,0)
}
\description{
  Get the information matrix using the Louis decomposition
  method.
}

