s_2<-function(stat_table,num_states,time_diffs,eigen_decomp, exact_time_rank, absorb_state,s1,trans_probs){
    num_times=length(time_diffs)
    num_stats=dim(stat_table)[1]
    s2=array(0,dim=c(num_states,num_states,num_times,num_stats,num_stats))
    s1
      stat_array=stat_table
     for(k in 1: num_stats){
       for(m in 1:k){
    if(stat_array[k,"ni"]>0&stat_array[m,"ni"]>0){   
       j1=stat_array[k,"ni"]
       k1=stat_array[k,"nj"]
       j2=stat_array[m,"ni"]
       k2=stat_array[m,"nj"]
       s2[,,,m,k]=joint_trans_trans_times(time_diffs, j1,k1,j2,k2,absorb_state,eigen_decomp, exact_time_rank) 
           for(l in 1:num_times){
         s2[,,l,m,k]=s2[,,l,m,k]/trans_probs[[l]]
       }
    }else if(stat_array[k,"di"]>0&stat_array[m,"di"]>0){  
       di1=stat_array[k,"di"]
       di2=stat_array[m,"di"]
       s2[,,,m,k]=joint_dur_dur_times(time_diffs, di1,di2,absorb_state,eigen_decomp, exact_time_rank)
          for(l in 1:num_times){
         s2[,,l,m,k]=s2[,,l,m,k]/trans_probs[[l]]
       }
    }else if(stat_array[k,"ni"]>0&stat_array[m,"di"]>0){
       tr1=stat_array[k,"ni"]
       tr2=stat_array[k,"nj"]
       di=stat_array[m,"di"]
       s2[,,,m,k]=joint_dur_trans_times(time_diffs, di, tr1, tr2,absorb_state,eigen_decomp, exact_time_rank)        
        for(l in 1:num_times){
        s2[,,l,m,k]=s2[,,l,m,k]/trans_probs[[l]]
       }
    }else if(stat_array[m,"ni"]>0&stat_array[k,"di"]>0){
       tr1=stat_array[m,"ni"]
       tr2=stat_array[m,"nj"]
       di=stat_array[k,"di"]
       s2[,,,m,k]=joint_dur_trans_times(time_diffs, di, tr1, tr2,absorb_state,eigen_decomp, exact_time_rank)   
         for(l in 1:num_times){
         s2[,,l,m,k]=s2[,,l,m,k]/trans_probs[[l]]
       }
    }else{
      s2[,,,m,k]=s1[,,,m]*s1[,,,k]
     }  
      s2[,,,k,m]=s2[,,,m,k]
    }
    }
  s2[is.na(s2)]=0
return(s2)
}

NULL
t1_2<-function(stat_table,obs_data,num_states,t1_1){
#'Get t1 for second moments of complete data sufficient statistics for an inidvidual
#'
#' @param stat_table table with sufficient statistic information (need to run suff.stat.table)
#' @param obs_data obseved data
#' @param num_states number of hidden states
#' @param t1_1 initial matrix for the parameter
#' @return t array(num.states,num.stats, num.stats)
num_stats=dim(stat_table)[1]
stat_array=stat_table
t2=array(0, dim=c(num_states,num_stats,num_stats))
   for(k in 1: num_stats){
    for(m in 1:k){
        t2[,k,m]=t1_1[,k]*t1_1[,m]
        t2[,m,k]=t2[,k,m]
      }
    }
  
return(t2)
}


  