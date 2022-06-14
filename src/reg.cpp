#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "rejection.h"
#include "noreg.h"


//' COM-Poisson sampling under mu and nu vectors
//' @param mu vector from linking function exp(beta_{mu,0} + beta_{mu,1}*x1 + ...)
//' @param nu vector from linking function exp(beta_{nu,0} + beta_{nu,1}*x1 + ...)
// [[Rcpp::export]]
NumericVector rcompoisreg(NumericVector mu, NumericVector nu){
  //Sample of COM-Poisson given mu and nu vectors
  int n = mu.length();
  NumericVector y_sample(n);
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

