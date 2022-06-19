#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "rejection.h"
#include "noreg.h"


//' COM-Poisson sampling under mu and nu vectors
//' @param beta_mu location regression parameters
//' @param beta_nu dispersion regression parameter
//' @param X_mu model matrix for mu (first column=1: intercept)
//' @param X_nu model matrix for nu (first column=1: intercept)
// [[Rcpp::export]]
arma::vec rcompoisreg(arma::vec& beta_mu, arma::vec& beta_nu, arma::mat& X_mu, arma::mat& X_nu){

  arma::mat mu = exp( beta_mu.as_row()*X_mu.t());
  arma::mat nu = exp( beta_nu.as_row()*X_nu.t());
  
  int n =X_mu.n_rows;
  arma::vec y_sample(n);
  for(int i =0; i< n; i ++){
    y_sample[i] =  rcompois(1, mu[i], nu[i])[0];
  }
  return y_sample;
}


List logqregression(arma::vec& y, arma::vec& beta_mu, arma::vec& beta_nu, arma::mat& X_mu, arma::mat& X_nu){
  
  arma::mat mu = exp( beta_mu.as_row()*X_mu.t());
  arma::mat nu = exp( beta_nu.as_row()*X_nu.t());
  
  NumericVector mu_nv = as<NumericVector>(wrap(mu));
  NumericVector nu_nv = as<NumericVector>(wrap(nu));
  NumericVector y_nv = as<NumericVector>(wrap(y));
  
  double mu_term = sum((nu_nv*y_nv)*log(mu_nv));
  double nu_term = -sum(nu_nv*log(factorial(y_nv)));
    
  return List::create(_["mu"] = mu_term,
                    _["nu"] = nu_term);
}
  
//' Normal proposal
//' @param index starting from 0, the parameter to propose
//' @param param vector
//' @param sigma vector
//' @return sigma param vector with updated position 'index'
// [[Rcpp::export]]
arma::vec proposal_normal(int index, arma::vec& param, arma::vec& sigma){
   arma::vec prime = param;
   prime[index] = R::rnorm(param[index], sigma[index]); //mean, sd
   return prime;
}

//' Exchange move for location regression parameter
//' @param index index of parameter to update, from 0
//' @param beta_mu current beta_mu
//' @param beta_nu current beta_nu
//' @param y response vector
//' @param X_mu location model matrix
//' @param X_nu dispersion model matrix
//' @param sigma_mu proposal variances
//' @param hyperparams_mu named list with 'mean' and 'sd'
// [[Rcpp::export]]
List update_beta_mu(int index, arma::vec& beta_mu, arma::vec& beta_nu, arma::vec& y, arma::mat& X_mu, arma::mat& X_nu, arma::vec& sigma_mu, List hyperparams_mu){
  
  arma::vec beta_mu_prime = proposal_normal(index, beta_mu, sigma_mu);

  //Priors:
  double prior_mean = hyperparams_mu["mean"]; double prior_sd = hyperparams_mu["sd"];
  double log_prior_current = -pow( (beta_mu[index] -prior_mean)/prior_sd, 2)/(2*prior_sd);
  double log_prior_prime = -pow( (beta_mu_prime[index] - prior_mean)/prior_sd, 2)/(2*prior_sd);
  
  //Exchange draw:
  arma::vec y_prime = rcompoisreg(beta_mu_prime, beta_nu, X_mu, X_nu);
  
  double logq_current = logqregression(y, beta_mu, beta_nu, X_mu, X_nu)["mu"];
  double logq_prime = logqregression(y, beta_mu_prime, beta_nu, X_mu, X_nu)["mu"];
  
  double logq_current_aux = logqregression(y_prime, beta_mu, beta_nu, X_mu, X_nu)["mu"];
  double logq_prime_aux = logqregression(y_prime, beta_mu_prime, beta_nu, X_mu, X_nu)["mu"];
  
  double alpha_num = logq_prime + log_prior_prime + logq_current_aux;
  double alpha_den = logq_current + log_prior_current + logq_prime_aux;
  
  double alpha = alpha_num - alpha_den;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  arma::vec beta_mu_accepted;
  
  if( compare[0]==true ){
    beta_mu_accepted = beta_mu_prime;
    accept=1;
  } else {
    beta_mu_accepted = beta_mu;
    accept = 0;
  }
  
  return List::create(_["beta_mu"] = beta_mu_accepted, 
                      _["accepted"] = accept);
}

//' Exchange move for dispersion regression parameters
//' @param index index of parameter to update, from 0
//' @param beta_mu current beta_mu
//' @param beta_nu current beta_nu
//' @param y response vector
//' @param X_mu location model matrix
//' @param X_nu dispersion model matrix
//' @param sigma_nu proposal variances
//' @param hyperparams_nu named list with 'mean' and 'sd'
// [[Rcpp::export]]
List update_beta_nu(int index, arma::vec& beta_mu, arma::vec& beta_nu, arma::vec& y, arma::mat& X_mu, arma::mat& X_nu, arma::vec& sigma_nu, List hyperparams_nu){
  
  arma::vec beta_nu_prime = proposal_normal(index, beta_nu, sigma_nu);
  
  //Priors:
  double prior_mean = hyperparams_nu["mean"]; double prior_sd = hyperparams_nu["sd"];
  double log_prior_current = -pow( (beta_nu[index] -prior_mean)/prior_sd, 2)/(2*prior_sd);
  double log_prior_prime = -pow( (beta_nu_prime[index] - prior_mean)/prior_sd, 2)/(2*prior_sd);
  
  //Exchange draw:
  arma::vec y_prime = rcompoisreg(beta_mu, beta_nu_prime, X_mu, X_nu);
  
  double logq_current = logqregression(y, beta_mu, beta_nu, X_mu, X_nu)["nu"];
  double logq_prime = logqregression(y, beta_mu, beta_nu_prime, X_mu, X_nu)["nu"];
  
  double logq_current_aux = logqregression(y_prime, beta_mu, beta_nu, X_mu, X_nu)["nu"];
  double logq_prime_aux = logqregression(y_prime, beta_mu, beta_nu_prime, X_mu, X_nu)["nu"];
  
  double alpha_num = logq_prime + log_prior_prime + logq_current_aux;
  double alpha_den = logq_current + log_prior_current + logq_prime_aux;
  
  double alpha = alpha_num - alpha_den;
  
  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  arma::vec beta_nu_accepted;
  
  if( compare[0]==true ){
    beta_nu_accepted = beta_nu_prime;
    accept=1;
  } else {
    beta_nu_accepted = beta_nu;
    accept = 0;
  }
  
  return List::create(_["beta_nu"] = beta_nu_accepted, 
                      _["accepted"] = accept);
}


//' MCMC: regression case
//' @param beta_mu_init initial beta_mu
//' @param beta_nu_init initial beta_nu
//' @param y response vector
//' @param X_mu location model matrix
//' @param X_nu dispersion model matrix
//' @param burn_in integer
//' @param n_iter integer
//' @param sigma_mu proposal variances
//' @param sigma_nu proposal variances
//' @param hyperparams_mu named list with 'mean' and 'sd'
//' @param hyperparams_nu named list with 'mean' and 'sd'
// [[Rcpp::export]]
List exchange_reg(arma::vec& beta_mu_init, arma::vec& beta_nu_init, arma::vec& y, arma::mat& X_mu, arma::mat& X_nu, int burn_in, int n_iter, arma::vec& sigma_mu, arma::vec& sigma_nu, List hyperparams_mu,  List hyperparams_nu){
  
  int n_beta= X_mu.n_cols; int n_nu = X_nu.n_cols;
  
  //Initialize chains and acceptance rate
  arma::mat betamu_chain(n_iter-burn_in, n_beta);
  arma::mat betanu_chain(n_iter-burn_in, n_nu);
  
  arma::vec ac_counter_mu(n_beta); 
  ac_counter_mu.zeros();
  arma::vec ac_counter_nu(n_nu); 
  ac_counter_nu.zeros();
  
  arma::vec ac_rate_mu(n_beta); 
  arma::vec ac_rate_nu(n_nu); 
  ac_rate_mu.zeros();
  ac_rate_nu.zeros();
  
  arma::vec betamu_current = beta_mu_init;
  arma::vec betanu_current = beta_nu_init;
  
  int step =1;
  while(step <= n_iter){
  //Update beta mu elements:
  for(int j =0; j < n_beta; j++){
    List update1 =  update_beta_mu(j, betamu_current, betanu_current, y, X_mu, X_nu, sigma_mu, hyperparams_mu);
    betamu_current =  as<arma::vec>(update1["beta_mu"]);
    
    ac_counter_mu[j] += int(update1["accepted"]);
    ac_rate_mu = (1.0*ac_counter_mu)/(1.0*step);
    
    sigma_mu[j] = proposal_adjust(sigma_mu[j], ac_rate_mu[j], step);
  }
  
  
  //Update beta nu elements:
  for(int k =0; k < n_nu; k++){
    List update2 =  update_beta_nu(k, betamu_current, betanu_current, y, X_mu, X_nu, sigma_nu, hyperparams_nu);
    betanu_current =  as<arma::vec>(update2["beta_nu"]);
    
    ac_counter_nu[k] += int(update2["accepted"]);
    ac_rate_nu = (1.0*ac_counter_nu)/(1.0*step);
    
    sigma_nu[k] = proposal_adjust(sigma_nu[k], ac_rate_nu[k], step);
  }
  
  //Storing
  if(step > burn_in){
    betamu_chain.row(step-burn_in-1) = betamu_current.as_row();
    betanu_chain.row(step-burn_in-1) = betanu_current.as_row();
  }
  step +=  1;
  if((step%1000) == 0){ 
    Rprintf("Step: %i \n", step);
  }
  } // End interations

  List output = List::create(_["betamu"] = betamu_chain, 
                             _["betanu"] = betanu_chain,
                             _["ac_rates_mu"] = ac_rate_mu,
                             _["ac_rates_nu"] = ac_rate_nu);
  
  
  return(output);
              
}          
              
     