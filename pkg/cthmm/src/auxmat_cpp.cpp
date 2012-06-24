#include "auxmat_cpp.h"
arma::cx_mat auxmat_cpp(arma::cx_colvec xx, double yy){
	int n_xx = xx.n_elem;
	arma::cx_mat ans(n_xx,n_xx);
	
	for(int i = 0; i < n_xx; i++) {
		for(int j = 0; j < n_xx; j++){
			if (xx(i)==xx(j)){
				ans(i,j) = exp(xx(j)*yy)*yy;
			}else{
				ans(i,j) = (exp(xx(i)*yy) - exp(xx(j)*yy))/(xx(i) - xx(j));
			}
		}
	}
	return (ans);
}

arma::cx_cube auxmat2_cpp(arma::cx_colvec v, double t){
	
	int n_v=v.n_elem;
	arma::cx_cube ans(n_v,n_v,n_v);
	
	for(int i = 0; i < n_v; i++) {
		for(int j = 0; j < n_v; j++){
			for(int k= 0; k < n_v; k++){
				if (v(i)==v(j) & v(j)==v(k)){
					ans(i,j,k)=exp(v(i)*t)*(t*t/2);
				}else if((v(i)==v(j)) &(v(j)!=v(k))){
					ans(i,j,k)=t*exp(v(i)*t)/(v(j)-v(k))-(exp(v(j)*t)-exp(v(k)*t))/((v(j)-v(k))*(v(j)-v(k)));
				}else if((v(i)!=v(j)) &(v(i)==v(k))){
					ans(i,j,k)=(t*exp(v(i)*t)-(exp(v(j)*t)-exp(v(k)*t))/(v(j)-v(k)))/(v(i)-v(j));
				}else if((v(i)!=v(j)) &(v(j)==v(k))){
					ans(i,j,k)=((exp(v(i)*t)-exp(v(k)*t))/(v(i)-v(k))-t*exp(v(j)*t))/(v(i)-v(j));
				}else{
					ans(i,j,k)=(exp(v(i)*t)-exp(v(k)*t))/((v(i)-v(j))*(v(i)-v(k)))-(exp(v(j)*t)-exp(v(k)*t))/((v(i)-v(j))*(v(j)-v(k)));
				}
			}
		}
	}
	
	
	return((ans));
}


arma::colvec hobolth_fun(Rcpp::List eigen_dec, double t, int a, int b, int c, int d, int e, int f){
	arma::cx_colvec v=Rcpp::as<arma::cx_colvec>(eigen_dec["values"]);
	arma::cx_mat U=Rcpp::as<arma::cx_mat>(eigen_dec["vectors"]);
	arma::cx_mat U_inv=Rcpp::as<arma::cx_mat>(eigen_dec["invvectors"]);
	int size=v.n_elem;
	
	arma::cx_cube int_array=auxmat2_cpp(v,t);
	arma::cx_colvec temp_sum(1);
	temp_sum(0)=0;
	
	
	
	for(int i=0; i<size;i++){
		for(int j=0; j<size;j++){
			for(int k=0; k<size;k++){
				temp_sum(0)=temp_sum(0)+U(a,i)*U_inv(i,b)*U(c,j)*U_inv(j,d)*U(e,k)*U_inv(k,f)*int_array(i,j,k);     
			}    
		}
	}
	arma::mat out=arma::real(temp_sum);
	return(out);
	
}

