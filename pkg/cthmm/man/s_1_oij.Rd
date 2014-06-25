\name{s_1_oij}
\alias{s_1_oij}
\title{s function for transitions s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is o_ij}
\usage{
  s_1_oij(oi, oj, num_states, num_times, obs_data)
}
\arguments{
  \item{obs_data}{}

  \item{num_times}{}

  \item{oi}{hidden state oi}

  \item{oj}{observed state oj}
}
\value{
  s array of dim(num_state,num_states,num.times,1)
}
\description{
  s function for transitions
  s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of
  stat [l] is o_ij
}

