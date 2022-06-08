#include <Rcpp.h>
using namespace Rcpp;
#include "rejection.h"

NumericVector proposal_ln(double param, double sigma){
  NumericVector param_prime = rlnorm(1, log(param), sigma);
  return(param_prime);
}

double logqcomp(NumericVector y, double mu, double nu){
  NumericVector log_q = nu*(y*log(mu) - log(factorial(y)));
  return sum(log_q);
}

//' Update positive parameters in the no-regression case
//' @param type 0 or 1: updates mu, or nu
//' @param params c(mu, nu) vector
//' @param y data
//' @param sigma c(sigma_mu, sigma_nu) vector of proposal parameters
//' @param shape c(shape_mu, shape_nu) Gamma prior shape parameters
//' @param rate c(rate_mu, rate_nu) Gamma prior rate parameters
// [[Rcpp::export]]
List update_positive(int type, NumericVector params, NumericVector y, NumericVector sigma, NumericVector shape, NumericVector rate){
  
  //Proposed draw
  NumericVector params_prime = clone(params);
  params_prime[type] = proposal_ln(params[type], sigma[type])[0];
 
  //Priors
  double log_prior_current = (shape[type]-1)*log(params[type]) - rate[type]*params[type];
  double log_prior_prime = (shape[type]-1)*log(params_prime[type]) - rate[type]*params_prime[type];
  
  //Auxiliary draws
  NumericVector y_prime = rcompois(y.size(), params_prime[0], params_prime[1]);
  
  //Unnormalised COM-Poisson likelihood
  double logq_current = logqcomp(y, params[0], params[1]);
  double logq_prime = logqcomp(y, params_prime[0], params_prime[1]);
  
  double logq_current_aux = logqcomp(y_prime, params[0], params[1]);
  double logq_prime_aux = logqcomp(y_prime, params_prime[0], params_prime[1]);
  
  double alpha_num = logq_prime + log_prior_prime + logq_current_aux;
  double alpha_den = logq_current + log_prior_current + logq_prime_aux;
  
  double alpha = alpha_num - alpha_den + log(params_prime[type]) - log(params[type]);

  LogicalVector compare = log(runif(1)) < alpha;
  int accept;
  NumericVector params_accepted;
  
  if( compare[0]==true ){
    params_accepted = params_prime;
    accept=1;
  } else {
    params_accepted = params;
    accept = 0;
  }
  
  return List::create(_["prime"] = params_prime, _["alpha"] = alpha, _["accepted"] = accept);
  
}

