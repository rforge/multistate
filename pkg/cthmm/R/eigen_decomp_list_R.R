eigen_decomp_list_R<-function(rate_matrix_list){
      out<-.Call( "eigen_decomp_list",rate_matrix_list)
     return(out)
}
