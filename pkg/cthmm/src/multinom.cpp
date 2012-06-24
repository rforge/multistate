#include <RcppArmadillo.h>
#include <Rcpp.h>
RcppExport SEXP multinom_hessian(SEXP prob_matrix,SEXP deriv_,SEXP myDims, SEXP N_vector);
arma::cube multinom_cov(arma::mat prob,arma::colvec N_vector);


RcppExport SEXP multinom_hessian(SEXP prob_matrix,SEXP deriv_,SEXP myDims, SEXP N_vector_){
 
 arma::mat prob=Rcpp::as<arma::mat>(prob_matrix);
 arma::colvec N_vector=Rcpp::as<arma::colvec>(N_vector_);

 Rcpp::NumericVector deriv_array(deriv_);
 Rcpp::IntegerVector arrayDims(myDims);
 arma::cube deriv(deriv_array.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
 
 arma::cube cov=multinom_cov(prob,N_vector);
 arma::mat hessian=arma::zeros<arma::mat>(deriv.n_rows,deriv.n_rows);

 for(int m=0;m<deriv.n_slices;m++){
  hessian=hessian+deriv.slice(m)*cov.slice(m)*arma::strans(deriv.slice(m));
 }
 hessian=-1*hessian;
 return(Rcpp::wrap(hessian));
}



arma::cube multinom_cov(arma::mat prob,arma::colvec N_vector){
 arma::cube out=arma::zeros<arma::cube>(prob.n_rows,prob.n_rows,prob.n_cols);

for(int m=0;m<prob.n_cols;m++){
if(prob.n_rows==1){
out.slice(m)=prob.col(m)%(1-prob.col(m));
}else{
  out.slice(m)=(-1)*prob.col(m)*strans(prob.col(m));
  out.slice(m).diag()=prob.col(m)%(1-prob.col(m));
}
 out.slice(m)=out.slice(m)*N_vector(m);
}
 return(out);
}
