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

//' Normal proposal
//' @param index starting from 0, the parameter to propose
//' @param param vector
//' @param sigma vector
//' @return param vector with updated position 'index'
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
  double prior_current = -pow( (beta_mu[index] -prior_mean)/prior_sd, 2)/(2*prior_sd);
  double prior_prime = -pow( (beta_mu_prime[index] - prior_mean)/prior_sd, 2)/(2*prior_sd);
  
  //Exchange draw:
  arma::vec y_prime = rcompoisreg(beta_mu_prime, beta_nu, X_mu, X_nu);
  
  return List::create(_["prior_prime"] = prior_prime,_["beta_prime"] = beta_mu_prime,
                      _["prior_current"] = prior_current, _["y_draw"] = y_prime);
  
}

