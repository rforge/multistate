\name{s_1_nij}
\alias{s_1_nij}
\title{s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is n_ij}
\usage{
  s_1_nij(ni, nj, num_states, time_diffs, rate_eigen,
    exact_time_ranks, absorb_state, trans_probs_list)
}
\arguments{
  \item{time_diffs}{vector of time intervals}

  \item{ni}{}

  \item{nj}{}

  \item{num_states}{number of hidden states}

  \item{exact_time_ranks}{}

  \item{absorb_state}{}

  \item{rate_eigen}{}

  \item{trans_probs_list}{list with the transition
  probabilities}
}
\value{
  s array of dim(num_state,num_states,num.times,1)
}
\description{
  s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if
  the type of stat [l] is n_ij
}

