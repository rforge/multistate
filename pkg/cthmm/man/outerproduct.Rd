\name{outerproduct}
\alias{outerproduct}
\title{Outerproduct}
\usage{
  outerproduct(l_1, l_2, c1, c2, M1, M2)
}
\arguments{
  \item{l_1}{vector of indices of Y in D, i.e. D=Y[l_a]}

  \item{l_2}{vector of indices of Y in E, i.e. E=Y[l_2]}

  \item{M2}{matrix of E[YY^T|o]}

  \item{M1}{vector of E[Y|o]}

  \item{c1}{vector of constants to be multiplied by
  elements in A}

  \item{c2}{vector of constants to be multipled by elements
  in B}
}
\value{
  matrix c1*D(c2*E)^T
}
\description{
  Get (c1*D)(c2*E)^T. Suppose Y is the vector of complete
  data sufficient statistics.
}
\details{
  If l_1[k]=0 (l_2[k]=0) this signfies a value of 1 rather
  than a suffcient statistic in Y. Hence in this case the
  desired product is a function of the first moment vector
  M1 rather than M2.
}
\author{
  Jane Lange
}

