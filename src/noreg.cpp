#include <Rcpp.h>
using namespace Rcpp;
#include "rejection.h"

List proposal_ln(double param, double sigma){
  NumericVector param_prime = rlnorm(1, log(param), sigma);
  double log_q_ratio = log(param_prime[0]) - log(param);
  List output =  List::create(_["proposed"] = param_prime, _["log_q_ratio"] = log_q_ratio);
  return(output);
}

//' Update positive parameters in the no-regression case
//' @param type 0 or 1: updates mu, or nu
//' @param params c(mu, nu) vector
//' @param y data
//' @param sigma c(sigma_mu, sigma_nu) vector of proposal parameters
// [[Rcpp::export]]
List update_positive(int type, NumericVector params, NumericVector y, NumericVector sigma){
  
  NumericVector params_prime = clone(params);
  List proposed = proposal_ln(params[type], sigma[type]);
  
  params_prime[type] = proposed["proposed"];
  List output = List::create(_["proposed"] = params_prime,
                        _["log_q_ratio"] = proposed["log_q_ratio"]);
  
  return output;
  
}

