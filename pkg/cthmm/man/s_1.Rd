\name{s_1}
\alias{s_1}
\title{Get s matrix for first moments of complete data sufficient statistics for an individual}
\usage{
  s_1(stat_table, num_states, obs_data, time_diffs,
    rate_eigen, exact_time_ranks, absorb_state,
    trans_probs_list)
}
\arguments{
  \item{stat_table}{table with sufficient statistic
  information (need to run suff.stat.table)}

  \item{obs_data}{obseved data}

  \item{time_diffs}{t2-t1, t3-t2,...,tn-t_{n-1}}

  \item{num_states}{number of hidden states}

  \item{trans_probs_list}{list with the transition
  probabilities for the individual}
}
\value{
  s array( dim(num_states, num_states,num.times,num.stats)
}
\description{
  This functions gets the s_k(x_k,x_{k+1}) function for the
  recursive smoothing algorithm.
}

