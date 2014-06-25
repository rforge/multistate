\name{t1_1}
\alias{t1_1}
\title{Get t1 for first moments of complete data sufficient statistics for an inidvidual}
\usage{
  t1_1(stat_table, obs_data, num_states)
}
\arguments{
  \item{stat_table}{table with sufficient statistic
  information (need to run suff.stat.table)}

  \item{obs_data}{obseved data}

  \item{num_states}{number of hidden states}
}
\value{
  t array(num.states,num.stats)
}
\description{
  Get t1 for first moments of complete data sufficient
  statistics for an inidvidual
}

