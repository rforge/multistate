NULL
#'Get t1 for first moments of complete data sufficient statistics for an inidvidual
#'
#' @param stat_table table with sufficient statistic information (need to run suff.stat.table)
#' @param obs_data obseved data
#' @param num_states number of hidden states
#' @return t array(num.states,num.stats)
t1_1<-function(stat_table,obs_data,num_states){
  num_stats=dim(stat_table)[1]
  stat_array=stat_table
  t1=array(0,dim=c(num_states,num_stats))
  
  for(i in 1:num_stats){
  
	if(!is.null(stat_array[i,"ni"])){
	if(stat_array[i,"ni"]>0){
            t1[,i]=t1_1_nij(num_states)
		}
	}
	if(!is.null(stat_array[i,"di"])){
	  if(stat_array[i,"di"]>0){
       t1[,i]=t1_1_di(num_states)
      }
    }
    if(!is.null(stat_array[i,"zi"])){
		if(stat_array[i,"zi"]>0){
       t1[,i]=t1_1_zi(stat_array[i,"zi"],num_states)
	}
	}
	  if(!is.null(stat_array[i,"oi"])){
		  if(stat_array[i,"oi"]>0){
       t1[,i]=t1_1_oij(stat_array[i,"oi"],stat_array[i,"oj"],num_states,obs_data)
    }
	}
	if(!is.null(stat_array[i,"Ni"])){
	   if(stat_array[i,"Ni"]>0){
       t1[,i]=t1_1_Ni(stat_array[i,"Ni"],num_states)
	}
    }
  
  }
  return(t1)
}
  
NULL
#'Get first moment t1 for nij
#' @param num_states number of hidden states
#' @return t1 array (num.states, 1)
t1_1_nij<-function(num_states){
  out=array(0, dim=c(num_states,1))
  return(out)
}

NULL
#'Get first moment t1 for di
#' @param num_states number of hidden states
#' @return t1 array (num.states, 1)
t1_1_di<-function(num_states){
   out=array(0, dim=c(num_states,1))
  return(out)
}

NULL
#'Get first moment t1 for zi
#' @param num_states
#' @param zi 
#' @return t1 array dim=(num_states, 1)
t1_1_zi<-function(zi,num_states){
  out=array(0, dim=c(num_states,1))
  out[zi,1]=1
  return(out)
}

NULL
#'Get first moment t1 for Ni
#' @param Ni
#' @param num_states
#' @return t1 array (num_states, 1)
t1_1_Ni<-function(Ni,num_states){
  out=array(0, dim=c(num_states,1))
  out[Ni,1]=1
  return(out)
}

NULL
#'Get first moment t1 for oij
#' @param oi
#' @param oj
#' @param num_states
#' @param obs_data
#' @return t1 array (num.states, 1)
t1_1_oij<-function(oi,oj,num_states, obs_data){
  out=array(0, dim=c(num_states,1))
  out[oi,1]=1*I(obs_data[1]==oj)
  return(out)
}


NULL
#'Get s matrix for first moments of complete data sufficient statistics for an individual
#'
#' This functions gets the s_k(x_k,x_{k+1}) function for the recursive smoothing algorithm. 
#' @param stat_table table with sufficient statistic information (need to run suff.stat.table)
#' @param obs_data obseved data
#' @param time_diffs t2-t1, t3-t2,...,tn-t_{n-1}
#' @param num_states number of hidden states
#' @param trans_probs_list list with the transition probabilities for the individual
#' @return s array( dim(num_states, num_states,num.times,num.stats)
#'
s_1<-function(stat_table,num_states,obs_data,time_diffs,rate_eigen, exact_time_ranks, absorb_state,trans_probs_list){
  num_times=length(obs_data)-1
  obs_data=obs_data[-1]
  num_stats=dim(stat_table)[1]
  stat_array=stat_table
  s_1=array(0,dim=c(num_states,num_states,num_times,num_stats))
  
  for(i in 1:num_stats){
	  if(!is.null(stat_array[i,"ni"])){
	if(stat_array[i,"ni"]>0){
      s_1[,,,i]= s_1_nij(stat_array[i,"ni"], stat_array[i,"nj"], num_states, time_diffs, rate_eigen, exact_time_ranks, absorb_state,trans_probs_list)
     }
	}
	if(!is.null(stat_array[i,"di"])){
    if(stat_array[i,"di"]>0){
       s_1[,,,i]=s_1_di(stat_array[i,"di"], num_states, time_diffs, rate_eigen, exact_time_ranks, absorb_state,trans_probs_list)
    }
	}
	if(!is.null(stat_array[i,"zi"])){
    if(stat_array[i,"zi"]>0){
       s_1[,,,i]=s_1_zi(stat_array[i,"zi"],num_states,num_times)
    }
	}
	if(!is.null(stat_array[i,"oi"])){
	if(stat_array[i,"oi"]>0){
       s_1[,,,i]=s_1_oij(stat_array[i,"oi"],stat_array[i,"oj"],num_states,num_times,obs_data)
    }
	   }
	   if(!is.null(stat_array[i,"Ni"])){				
    if(stat_array[i,"Ni"]>0){
       s_1[,,,i]=s_1_Ni(stat_array[i,"Ni"],num_states,num_times)
    }
	}
  
  }
  return(s_1)
}


NULL
#' s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is n_ij
#'
#' @param time_diffs vector of time intervals
#' @param ni 
#' @param nj
#' @param num_states number of hidden states
#' @param exact_time_ranks 
#' @param absorb_state
#' @param rate_eigen 
#' @param trans_probs_list list with the transition probabilities
#' @return s array of dim(num_state,num_states,num.times,1)
s_1_nij<-function(ni,nj,num_states,time_diffs,rate_eigen,exact_time_ranks,absorb_state,trans_probs_list){
  num_times=length(time_diffs)
  out=array(0,dim=c(num_states,num_states,num_times,1))
  regist_matrix=matrix(0,nrow=num_states,ncol=num_states)
  regist_matrix[ni,nj]=1
  out[,,,1]=trans_times_R(i=ni,j=nj,time_diffs,rate_eigen, num_states=num_states,exact_time_ranks,absorb_state)
  for(i in 1:num_times){out[,,i,1]=out[,,i,1]/trans_probs_list[[i]]}
    out[is.na(out)]=0

  return(out)

}



NULL
#' s function for s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is d_i
#'
#' @param di state d_i
#' @param num_states
#' @param time_diffs
#' @param rate_eigen 
#' @param exact_time_ranks
#' @param absorb_state
#' @param trans_probs_list list with the transition probabilities
#' @return s array of dim(num_state,num_states,num.times,1)
s_1_di<-function(di,num_states, time_diffs,rate_eigen,exact_time_ranks,absorb_state,trans_probs_list){
   num_times=length(time_diffs)
   out=array(0,dim=c(num_states,num_states,num_times,1)) 
   out[,,,1]= dur_times_R(di, time_diffs, rate_eigen, num_states, exact_time_ranks, absorb_state)
   for(i in 1:num_times){out[,,i,1]=out[,,i,1]/trans_probs_list[[i]]}
    out[is.na(out)]=0
  return(out)
}


NULL
#' s function for transitions s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is Z_i
#'
#' @param zi state i
#' @param num_states
#' @param num_times
#' @return s array of dim(num_states,num_states,num.times,1)
s_1_zi<-function(zi,num_states,num_times){
  out=array(0,dim=c(num_states,num_states,num_times,1))
  return(out)
}


NULL
#' s function for transitions s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is o_ij
#'
#' @param obs_data
#' @param num_times
#' @param oi hidden state oi
#' @param oj observed state oj
#' @return s array of dim(num_state,num_states,num.times,1)
s_1_oij<-function(oi,oj,num_states, num_times, obs_data){
  out=array(0,dim=c(num_states,num_states,num_times,1))
  out[,oi,,1]=1
  out[,oi,,1]=t(t(out[,oi,,1])*I(obs_data==oj))
  return(out)
}

NULL
#' s function for transitions s_k(x_k,x_{k+1})[1:num.times,stat[l]], if the type of stat [l] is Ni
#' @param Ni hidden state i
#' @param num_states
#' @param num_times
#' @return s array of dim(num_states,num_states,num.times,1)
s_1_Ni<-function(Ni,num_states,num_times){
   out=array(0,dim=c(num_states,num_states,num_times,1))
  out[,Ni,,1]=1
  return(out)
}
