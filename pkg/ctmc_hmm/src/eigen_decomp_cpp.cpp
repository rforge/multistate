#include<RcppArmadillo.h>
#include<Rcpp.h>
Rcpp::List eigen_decomp_cpp(SEXP the_rate);
RcppExport SEXP eigen_decomp_list(SEXP rate_matrix_list);


Rcpp::List eigen_decomp_cpp(SEXP the_rate){	 
	arma::mat rate=Rcpp::as<arma::mat>(the_rate);
	arma::cx_vec eigval;
	arma::cx_mat eigvec;
	arma::cx_mat invvectors;
	
	arma::eig_gen(eigval, eigvec, rate);
	invvectors=inv(eigvec);
	Rcpp::NumericVector values=Rcpp::wrap(eigval);
	values.attr( "dim" ) = R_NilValue ;

	Rcpp::List out=Rcpp::List::create(Rcpp::Named("rate") =rate,Rcpp::Named("values") =values, Rcpp::Named("vectors")=eigvec,Rcpp::Named("invvectors")=invvectors);
     return(out);
}

RcppExport SEXP eigen_decomp_list_old(SEXP rate_matrix_list){
	Rcpp::List rates_list(rate_matrix_list);
    int size=rates_list.size();
    Rcpp::List out_list(size);
	out_list = Rcpp::lapply(rates_list, eigen_decomp_cpp);
	return(Rcpp::wrap(out_list));
}
	
RcppExport SEXP eigen_decomp_list(SEXP rate_matrix_list){
    Rcpp::List rates_list(rate_matrix_list);
    int size=rates_list.size();
    Rcpp::List out_list(size);
	
	for(int i=0;i<size;i++){
		arma::mat rate=Rcpp::as<arma::mat>(rates_list[i]);
		arma::cx_colvec eigval;
		arma::cx_mat eigvec;
		arma::eig_gen(eigval, eigvec, rate);
		arma::cx_mat invvectors=arma::inv(eigvec);
		
		bool replicate=0;
		arma::colvec realvec=arma::real(eigval);
		for(int j=0; j<realvec.n_elem-1;j++){
			if(((realvec(j+1)-realvec(j))<1e-5)||((realvec(j)-realvec(j+1))<1e-5)){
				replicate=1;
				exit;
			}
		}
		//just adding this for now
		replicate=0;
		out_list[i]=Rcpp::List::create(Rcpp::Named("values") =eigval, Rcpp::Named("invvectors") =invvectors, Rcpp::Named("vectors") =eigvec,
									   Rcpp::Named("rate") =real(rate),Rcpp::Named("replicate")=replicate);
		
	}
	
	return(wrap(out_list));
}

