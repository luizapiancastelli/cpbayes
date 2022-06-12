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
  
  double alpha = alpha_num - alpha_den + log(params[type]) - log(params_prime[type]) ;

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
  
  return List::create(_["params"] = params_accepted, 
                      _["accepted"] = accept);
  
}

double proposal_adjust(double current_var, double accept_rate, int nprops){
  double adj = (1.0/(1.0*nprops))*(accept_rate - 0.44);
  double new_sd = current_var*exp(adj);
  return new_sd;
}


//' MCMC sampling: no regression case
//' @param y data
//' @param params0 initial parameter values
//' @param sigma c(sigma_mu, sigma_nu) vector of proposal parameters
//' @param n_iter number of iterations to store
//' @param burn_in burn-in period
//' @param hyperparams named list containing 'shape' c(shape_mu, shape_nu) gamma prior shape parameters and
//' 'rate' c(rate_mu, rate_nu) Gamma prior rate parameters.
// [[Rcpp::export]]
List exchange_noreg(NumericVector y, NumericVector params0, NumericVector sigma, int n_iter, int burn_in, List hyperparams){

  NumericMatrix chain(n_iter, 2); 
  
  //Acceptance rates
  int ac_counter_mu=0; double ac_rate_mu;
  int ac_counter_nu=0; double ac_rate_nu;
  
  NumericVector params = clone(params0);
  int step =1;
  while(step <= (n_iter + burn_in) ){
    
    //mu move
    List mu_update = update_positive(0, params, y, sigma, hyperparams["shape"], hyperparams["rate"]);

    params= mu_update["params"];
    ac_counter_mu += int(mu_update["accepted"]);
    
    ac_rate_mu = (1.0*ac_counter_mu)/(1.0*step);
    sigma[0] = proposal_adjust(sigma[0], ac_rate_mu, step);
    
    //nu move
    List nu_update = update_positive(1, params, y, sigma, hyperparams["shape"], hyperparams["rate"]);
    
    params= nu_update["params"];
    ac_counter_nu += int(nu_update["accepted"]);
    
    ac_rate_nu = (1.0*ac_counter_nu)/(1.0*step);
    sigma[1] = proposal_adjust(sigma[1], ac_rate_nu, step);
    
    //Storing
    if(step > burn_in){
      chain( step-burn_in-1, _) = params;
    }
    step +=  1;
   if((step%1000) == 0){ 
        Rprintf("Step: %i \n", step);
    }
  }
  
  NumericVector ac_rates(2);
  ac_rates[0] = ac_rate_mu; ac_rates[1] = ac_rate_nu;
  List output = List::create(_["mu"] = chain(_,0), _["nu"] = chain(_,1),
                             _["ac_rates"] = ac_rates);
  
  
  return(output);
  
}
