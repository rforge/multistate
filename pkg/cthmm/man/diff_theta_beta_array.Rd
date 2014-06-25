\name{diff_theta_beta_array}
\alias{diff_theta_beta_array}
\title{diff_theta_beta_array
get an array that corresponds to diff_theta_beta for each individual in the dataset
The "theta" params are the natural parameters encoding the rate, initial distribution, and emission matrices.
The "beta"params are paramters in the linear predictors that relate to the natural parameters.}
\usage{
  diff_theta_beta_array(rates.setup, emission.setup,
    init.setup, num.subjects)
}
\arguments{
  \item{rates.setup}{rate setup object}

  \item{emission.setup}{emission setup object}

  \item{init.setup}{init setup object}

  \item{num.subjects}{number of subjects}
}
\value{
  an array of dimension (num theta param)x(num beta
  params)xnumber of subjects.
}
\description{
  diff_theta_beta_array get an array that corresponds to
  diff_theta_beta for each individual in the dataset The
  "theta" params are the natural parameters encoding the
  rate, initial distribution, and emission matrices. The
  "beta"params are paramters in the linear predictors that
  relate to the natural parameters.
}

