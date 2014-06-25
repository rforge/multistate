#include "get_mean.h"
RcppExport SEXP get_mean(SEXP single_joint_mean, SEXP obs_data, SEXP forward, SEXP backward, SEXP LL, SEXP emission_matrix){

double temp_sum=0;
double LogLik= Rcpp::as<double>(LL);
Rcpp::List joint_mean(single_joint_mean);
arma::mat eemission_matrix=Rcpp::as<arma::mat>(emission_matrix);
arma::mat fforward=Rcpp::as<arma::mat>(forward);
arma::mat bbackward=Rcpp::as<arma::mat>(backward);
arma::icolvec obsdata=Rcpp::as<arma::icolvec>(obs_data);
int state_size = eemission_matrix.n_rows;
arma::mat joint=Rcpp::as<arma::mat>(joint_mean[0]);

for(int r=0; r<state_size; r++){
for(int s=0; s<state_size; s++){
for(int l=0; l<(obsdata.n_elem-1);l++){
joint=Rcpp::as<arma::mat>(joint_mean[l]);
//std::cout <<joint(r,s)*eemission_matrix(s,(obsdata(l+1)-1))*exp(fforward(l,r)+bbackward(l+1,s));
temp_sum=temp_sum+joint(r,s)*eemission_matrix(s,(obsdata(l+1)-1))*exp(fforward(l,r)+bbackward(l+1,s));
}
}
}
return(Rcpp::wrap(temp_sum/exp(LogLik)));
}


double get_mean_cpp(std::vector<arma::mat> joint_mean, arma::icolvec obs_data, arma::mat forward, arma::mat backward, double LL, arma::mat emission_matrix){
	double temp_sum=0;
	int state_size = emission_matrix.n_rows;
	arma::mat joint=joint_mean[0];
	
	for(int r=0; r<state_size; r++){
		for(int s=0; s<state_size; s++){
			for(int l=0; l<(obs_data.n_elem-1);l++){
				joint=joint_mean[l];
				//std::cout <<joint(r,s)*emission_matrix(s,(obs_data(l+1)-1))*exp(forward(l,r)+backward(l+1,s));
				temp_sum=temp_sum+joint(r,s)*emission_matrix(s,(obs_data(l+1)-1))*exp(forward(l,r)+backward(l+1,s));
			}
		}
	}
	return(temp_sum/exp(LL));
}

double get_mean_cpp_time_dep_emission(std::vector<arma::mat> joint_mean, arma::icolvec obs_data, arma::mat forward, 
									  arma::mat backward, double LL, Rcpp::List indivEmission){
	double temp_sum=0;
	int state_size = forward.n_cols;
	arma::mat joint=joint_mean[0];
	
	for(int r=0; r<state_size; r++){
		for(int s=0; s<state_size; s++){
			for(int l=0; l<(obs_data.n_elem-1);l++){
				joint=joint_mean[l];
				arma::mat emission_matrix=Rcpp::as<arma::mat>(indivEmission[l+1]);
				temp_sum=temp_sum+joint(r,s)*emission_matrix(s,(obs_data(l+1)-1))*exp(forward(l,r)+backward(l+1,s));
			}
		}
	}
	return(temp_sum/exp(LL));
}


