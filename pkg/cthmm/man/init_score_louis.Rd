\name{init_score_louis}
\alias{init_score_louis}
\title{init_score_louis
Get the inital dist scores}
\usage{
  init_score_louis(init.setup, param.values, init.counts)
}
\arguments{
  \item{init.setup}{list init setup object}

  \item{init.param.values}{values of init dist parameters}

  \item{init.counts}{matrix (num.states x num.subjects)
  with multinomial vector encoding the initial dist.}
}
\value{
  score array (dim num params x number of subjects)
}
\description{
  init_score_louis Get the inital dist scores
}

