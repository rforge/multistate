\name{element.indices}
\alias{element.indices}
\title{Get the score element indices}
\usage{
  element.indices(suff.stats, rate.matrix, emission.matrix,
    init.dist)
}
\arguments{
  \item{suff.stats}{table with coding all of the sufficient
  statistics}

  \item{rate.matrix}{the rate matrix}

  \item{emission.matrix}{the emission matrix}

  \item{init.dist}{vector with the initial distribution}
}
\value{
  list with \item{l_a}{row index of suff.stats
  corresponding to A} \item{l_b}{row index of suff.stats
  corresponding to B--for 1, this index is set to 0}
  \item{c}{vector of constants}
}
\description{
  The complete data score equations consist of elements
  such as nij-di*lambda ij, zi-pi, and o_ij-Ni*e_ij. Hence
  the score can be written as A-cB, where A is a vector of
  nij,zi,oij and B, di,1,Ni, and c, lambda ij,pi,eij
}
\author{
  Jane Lange
}

