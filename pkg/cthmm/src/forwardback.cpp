#include "forwardback.h"
RcppExport SEXP likelihood_list(SEXP trans_prob_list, SEXP emission_list_,SEXP delta_list_,SEXP obs_list_);
Rcpp::List forwardback_list(Rcpp::List Pi, arma::mat eemission, arma::icolvec xx, arma::rowvec phi);

RcppExport SEXP likelihood_list(SEXP trans_prob_list, SEXP emission_list_,SEXP delta_list_,SEXP obs_list_){
	Rcpp::List Pi_list(trans_prob_list);
	Rcpp::List emission_list(emission_list_);
	Rcpp::List delta_list(delta_list_);
	Rcpp::List obs_list(obs_list_);
	Rcpp::List out(obs_list.size());
	for(int m=0;m<obs_list.size();m++){
		Rcpp::List Pi=Rcpp::as<Rcpp::List>(Pi_list[m]);
		arma::mat eemission=Rcpp::as<arma::mat>(emission_list[m]);
		arma::icolvec xx=Rcpp::as<arma::icolvec>(obs_list[m]);
		arma::rowvec phi=Rcpp::as<arma::rowvec>(delta_list[m]);
		out[m]=forwardback_list(Pi, eemission, xx, phi);
	}
	return Rcpp::wrap(out);
}
Rcpp::List forwardback_list(Rcpp::List Pi, arma::mat eemission, arma::icolvec xx, arma::rowvec phi){
	
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
	
	return Rcpp::List::create(Rcpp::Named("logalpha") =logalpha, Rcpp::Named("logbeta") =logbeta, Rcpp::Named("LL") =LL1);
}


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

