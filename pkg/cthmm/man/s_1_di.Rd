\name{s_1_di}
\alias{s_1_di}
\title{s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is d_i}
\usage{
  s_1_di(di, num_states, time_diffs, rate_eigen,
    exact_time_ranks, absorb_state, trans_probs_list)
}
\arguments{
  \item{di}{state d_i}

  \item{num_states}{}

  \item{time_diffs}{}

  \item{rate_eigen}{}

  \item{exact_time_ranks}{}

  \item{absorb_state}{}

  \item{trans_probs_list}{list with the transition
  probabilities}
}
\value{
  s array of dim(num_state,num_states,num.times,1)
}
\description{
  s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if
  the type of stat [l] is d_i
}

