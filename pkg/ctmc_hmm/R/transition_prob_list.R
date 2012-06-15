trans_probs_all <- function (eigen_decomp_list_, time_int_list_, exact_times_list_) 
.Call("trans_probs_all", eigen_decomp_list_, time_int_list_, exact_times_list_, 
PACKAGE = "test4")



transition_prob_list<-function(time.intervals, eigen.decomp, exact.times.ranks=NULL){
#####################################################################
#author JL 2/7/2011
#
#This function gets a list of transition probability matrices for several
#time intervals, assuming a time-homogeneous CTMC
#i.e., P(t2-t1), P(t3-t2), ...P(tn-t_{n-1})
#If exact observation times is indicated at time l>1, then set
# equal toP[[l-1]] times rate matrix(diag=0)
#
#INPUTS: time.intervals t2-t1, t3-t2, ...,t_n-t_n-1
#        rate.matrix=rate matrix
#        exact times.ranks=vector of indices corresponding to exact observation times (all indices>1!)
#OUTPUTS: rate.matrix = transition intensity matrix
###########################################################

probs.list<-lapply(time.intervals, FUN="mat_exp", ev=eigen.decomp)

#probs.list=transition_probs(eigen.decomp,time.intervals)

	
if(!is.null(exact.times.ranks)){
 rate.temp<-eigen.decomp$rate
 diag(rate.temp)=0
 for(i in 1:length(exact.times.ranks)){
  probs.list[[exact.times.ranks[i]-1]]<-probs.list[[exact.times.ranks[i]-1]]%*%rate.temp
 }
}
return(probs.list)
}
#####################################################################################



#####################################################################################
mat_exp<-function (t = 1,ev=NULL)
{
#####################################################################################
#
#
#This function is from the msm library (orginally called MatrixExp) that computes the matrix expnential,
# using the eigen decomposition method if there are no duplicated eigen values (=>matrix is diagonalizable)
# or otherwise the series approximaton method
#to get transition probability matrices
#INPUTS: t=time interval, mat=rate matrix, n=number of items in series approximation, k=Underflow correction factor, ev=eigen decomposition of rate matrix
#OUTPUTS: the transition probability matrix
#Note: need to have msm loaded
#
#
#####################################################################################
    mat=ev$rate
    nr <- nrow(mat)
    if(is.null(ev)){
    ev <- eigen(mat)
    }
    if (any(duplicated(ev$values))) {             
                resi <- .C("MatrixExpPadeR", res = double(length(mat)),
                as.double(mat), as.integer(nr), as.double(t))$res
                resi <- matrix(resi, nrow = nrow(mat))
            
            res <- resi
    } else {
             resi <-mat_exp_eigen_R(ev,t)
                 res <- resi
     
    }
    
   return(res);
}


