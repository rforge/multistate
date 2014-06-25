#include "forwardback.h"

RcppExport SEXP forwardbacktimedep(SEXP problist, SEXP emission, SEXP x, SEXP delta);
RcppExport SEXP forwardback(SEXP problist, SEXP emission, SEXP x, SEXP delta){

  Rcpp::List Pi(problist);
  arma::mat eemission=Rcpp::as<arma::mat>(emission);
  arma::icolvec xx	= Rcpp::as<arma::icolvec>(x);
  arma::rowvec phi	= Rcpp::as<arma::rowvec>(delta);

  int m = phi.n_elem;
  int n= xx.n_elem;
  arma::mat logalpha(n,m);
  arma::mat logbeta(n,m);

  double lscale=0;
  double sumphi;
  double LL1;
  int i;
  for(int i=0; i<n; i++){

    if(i>0){
      phi=phi*Rcpp::as<arma::mat>(Pi[i-1]);
     }
    phi=phi*(arma::diagmat(eemission.col(xx(i)-1)));
    sumphi=arma::sum(phi);
    phi=phi/sumphi;
    lscale =lscale + log(sumphi);
    logalpha.row(i) = log(phi) + lscale;
  }
   LL1=lscale;

   arma::colvec phi2= arma::ones<arma::vec>(m,1)/m;
   logbeta.row(n-1)=arma::zeros<arma::rowvec>(1,m);
   lscale=log(m);

    for(int k=(n-2); k>=0; k--){
      phi2= Rcpp::as<arma::mat>(Pi[k])*arma::diagmat(eemission.col(xx(k+1)-1))*phi2;
      logbeta.row(k) = arma::trans(log(phi2)) + lscale;
      sumphi=arma::sum(phi2);
      phi2 = phi2/sumphi;
      lscale = lscale + log(sumphi);
     }

    return Rcpp::List::create(logalpha, logbeta, LL1) ;
}

RcppExport SEXP forwardbacktimedep(SEXP problist, SEXP emission, SEXP x, SEXP delta){
	Rcpp::List Pi(problist);
	Rcpp::List emission_list(emission); 
	arma::icolvec xx	= Rcpp::as<arma::icolvec>(x);
	arma::rowvec phi	= Rcpp::as<arma::rowvec>(delta);
	
	int m = phi.n_elem;
	int n= xx.n_elem;
	
	arma::mat logalpha(n,m);
	arma::mat logbeta(n,m);
	
	double lscale=0;
	double sumphi;
	double LL1;
	int i;
	for(int i=0; i<n; i++){
		
		arma::mat eemission=Rcpp::as<arma::mat>(emission_list[i]);
		
		//std::cout << i;
		//std::cout << eemission;
		
		if(i>0){
			phi=phi*Rcpp::as<arma::mat>(Pi[i-1]);
		}
		phi=phi*(arma::diagmat(eemission.col(xx(i)-1)));
		sumphi=arma::sum(phi);
		phi=phi/sumphi;
		lscale =lscale + log(sumphi);
		logalpha.row(i) = log(phi) + lscale;
	}
	LL1=lscale;
	
	arma::colvec phi2= arma::ones<arma::vec>(m,1)/m;
	logbeta.row(n-1)=arma::zeros<arma::rowvec>(1,m);
	lscale=log(m);
	
    for(int k=(n-2); k>=0; k--){
		arma::mat eemission=Rcpp::as<arma::mat>(emission_list[k+1]);
        //std::cout << k;
        //std::cout << eemission;
		phi2= Rcpp::as<arma::mat>(Pi[k])*arma::diagmat(eemission.col(xx(k+1)-1))*phi2;
		logbeta.row(k) = arma::trans(log(phi2)) + lscale;
		sumphi=arma::sum(phi2);
		phi2 = phi2/sumphi;
		lscale = lscale + log(sumphi);
	}
	
    return Rcpp::List::create(logalpha, logbeta, LL1) ;
}