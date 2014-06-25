\name{E.step.emission.init.all}
\alias{E.step.emission.init.all}
\title{E step for emission and initial distributions}
\usage{
  E.step.emission.init.all(likelihood.forward.backward.list,
    the.data, num.states, num.obs.states, num.subjects)
}
\arguments{
  \item{likelihood.forward.backward.list}{list with LL,
  forward, and backward probs for each inividual}

  \item{the.data}{list with suject level obs.data}

  \item{num.states}{number of hidden states in the model}

  \item{num.obs.state}{number of observed states in the
  model}

  \item{num.subjects}{number of individuals}
}
\value{
  a list with two elements: \item{emission.array}{(dim
  (num.states, num.obs.states, num.subjects))}
  \item{init.array}{(dim (num.states,num.subjects))}
}
\description{
  Get the expected emission counts E[(O_ij)|o] and initial
  distribution counts E[(Z_i)|o] given observed data for
  all individuals. This implementation does not use the
  recursive filtering method, but rather uses forward and
  backward probabilities.
}
\author{
  Jane Lange
}

