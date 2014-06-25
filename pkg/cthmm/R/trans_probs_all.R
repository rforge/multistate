 trans_probs_all <- function (eigen_decomp_list_, time_int_list_, exact_times_list_,h_list_=NULL,Q_list_=NULL){
   if(is.null(h_list_)){
     out=.Call("trans_probs_all", eigen_decomp_list_, time_int_list_, exact_times_list_)
   }else{
     out=.Call("trans_probs_all_h", eigen_decomp_list_, time_int_list_, exact_times_list_,h_list_,Q_list_)
   }
   return(out)
} 

