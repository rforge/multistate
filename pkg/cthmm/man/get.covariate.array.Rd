\name{get.covariate.array}
\alias{get.covariate.array}
\title{Get the covariate arrays}
\usage{
  get.covariate.array(design.matrix, cov.data,
    variable.entries = NULL)
}
\arguments{
  \item{design.matrix}{a matrix that codes the constants
  and variables for the parameterization with rows
  corresponding to theta and columns corresponding to
  betas}

  \item{cov.data}{a data frame with baseline covariate data
  for each individual}

  \item{variable.entries}{matrix with rows for each
  parameter corresponding to variable covariates, with the
  names for each.}
}
\value{
  an array (dim=c(num.betas,num.thetas, num.subjects)) with
  covariates values for all subjects
}
\description{
  Creates matrices with covariate values for all subjects
  for use in parameterized rate, initial dist., and
  emission matrices.
}
\author{
  Jane Lange
}

