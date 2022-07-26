#include <Rcpp.h>
using namespace Rcpp;

NumericVector log_bound(double mu, double nu, double p){
  
  NumericVector a; NumericVector b;
  NumericVector c; NumericVector B;
  NumericVector denom;
  
  if(nu < 1){
    a = pow((1-p),(1/nu));
    b= floor(mu/a);
    c = nu*b;
    denom = b*log(1-p) + nu*log(factorial(b));
    B = log(1/p) + c*log(mu) - denom;
    
  } else if(nu >=1){
    a = floor(mu);
    B = (nu-1)*( a*log(mu) - log(factorial(a)) );
  }
  return(B);
}

//' COM-Poisson distribution
//' 
//' Exact sampling, approximate density, distribution and quantile computation.
//' 
//' @param n Integer. Number of observations.
//' @param mu Numeric value. Location parameter, \eqn{\mu >0}.
//' @param nu Numeric value. Dispersion parameter, \eqn{\nu >0}.
//' @return Integer count vector of \code{n} independent draws.
//'
//' @details The COM-Poisson distribution has density 
//' \deqn{f(y|\mu, \nu) = (\mu^y/y!)^\nu 1/Z(\mu, \nu) } 
//' for \eqn{y = 0, 1, 2, ...}, where \eqn{Z(\mu, \nu) = \sum_{y=0}^\infty (\mu^y/y!)^\nu }. 
//' 
//' If \eqn{\nu <1} the distribution is overdispersed (\eqn{E(Y) < Var(Y)}), exhibits underdispersion if \eqn{\nu >1} or reduces to a Poisson (equidispersed) if \eqn{\nu=1}.
//' The Geometric and Bernoulli distributions can be obtained as limiting cases. \cr 
//' 
//' Exact sampling \eqn{f(y|\mu, \nu)} is done with the \href{10.1214/20-BA1230}{fast-rejection sampler} sampler of Benson and Friel (2021),
//' for which C++ implementation is provided in \code{rcompois}.
//' 
//' 
//' Density, distribution and quantiles are unavailable from intractability of \eqn{Z(\mu, \nu)}
//' and can be approximated using a truncated sum \eqn{Z(\mu, \nu) = \sum_{y=0}^T (\mu^y/y!)^\nu }.
//' \code{dcompois}, \code{qcompois} and \code{pcompois} apply a \eqn{Z(\mu, \nu)} approximation that
//' recursively increments \eqn{T} until difference between successive probabilities is below \eqn{10^{-5}}. \cr
//'    
//' Asymptotic moments of the distribution ((Shmueli et al. (2005))) are \eqn{E(Y) ~ \mu + 1/(2\nu) -1/2}, \eqn{Var(Y) ~ \mu/\nu}. 
//' 
//' For further details on the fast-rejection sampler, see Benson and Friel (2021).
//' 
//' 
//' @examples
//' y1 = rcompois(100, 1, 0.5) ## independent overdispersed counts
//' y2 = rcompois(100, 2, 1) ## Poisson random draws
//' y3 = rcompois(100, 0.5, 1.5) ## an underdispersed count vector
//' 
//' @references
//' Benson, A. and Friel, N. (2021) Bayesian Inference, Model Selection and Likelihood Estimation using Fast Rejection Sampling: The Conway-Maxwell-Poisson Distribution.
//' Bayesian Analysis (16) 905-931.
//' 
//' Shmueli, G., Minka, T., Kadane, J, Borle, S. and Boatwright, P. (2005) A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution.
//' Journal of the Royal Statistical Society (C) (54) 127-142.
//' @rdname compois
// [[Rcpp::export]]
NumericVector rcompois(int n, double mu, double nu){
  
  int n_accepts = 0;
  NumericVector x(n);
  NumericVector B;
  
  if(nu<1){
    
    double p = 2*nu/( 2*mu*nu +1+nu );
    B = log_bound(mu,nu,p);
    
    while(n_accepts < n) {
      NumericVector y = rgeom(1, p);
      
      //Calculate acceptance
      NumericVector numer =nu*( y*log(mu) - log(factorial(y)) );
      NumericVector denom = B + y*log(1-p) + log(p);
      LogicalVector comp = log(runif(1)) <= numer - denom;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      }
      
    }
    
  } else {
    
    B = log_bound(mu,nu,1);
    
    while(n_accepts < n) {
      
      NumericVector y= rpois(1, mu); 
      NumericVector alpha = (nu-1)*( y*log(mu) - log(factorial(y)) ) - B;
      LogicalVector comp = log(runif(1)) <= alpha;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      }
    }
  }
  return(x);
}

// [[Rcpp::export]]
List rcompois_internal(int n, double mu, double nu){
  
  int n_accepts = 0;
  int n_proposed = 0;
  NumericVector x(n);
  NumericVector B;
  double logZ_gamma;
  
  if(nu<1){
    
    double p = 2*nu/( 2*mu*nu +1+nu );
    B = log_bound(mu,nu,p);
    logZ_gamma = 0.0;
    
    while(n_accepts < n) {
      NumericVector y = rgeom(1, p);
      
      //Calculate acceptance
      NumericVector numer =nu*( y*log(mu) - log(factorial(y)) );
      NumericVector denom = B + y*log(1-p) + log(p);
      LogicalVector comp = log(runif(1)) <= numer - denom;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      } 
      n_proposed += 1;
      
    }
    
  } else {
    
    B = log_bound(mu,nu,1);
    logZ_gamma = mu;
    
    while(n_accepts < n) {
      
      NumericVector y= rpois(1, mu); 
      NumericVector alpha = (nu-1)*( y*log(mu) - log(factorial(y)) ) - B;
      LogicalVector comp = log(runif(1)) <= alpha;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      }
      n_proposed += 1;
    }
  }
  
  double Mhat = 1.0*n_proposed/(1.0*n);
  double invZ_est = log(Mhat)- B[0] - logZ_gamma;
  
  return List::create(_["sampled"] = x, _["invZ_est"] = invZ_est, 
                      _["n_proposed"] = n_proposed,
                      _["Mhat"] = Mhat);
}


//' Estimate COM-Poisson normalising constant with the fast-rejection sampler
//'
//' @param r number of fast-rejection sampler draws
//' @param mu location parameter
//' @param nu dispersion parameter
//' @return double 
// [[Rcpp::export]]
double Zhat(int r, double mu, double nu){
  
  List sampler_output = rcompois_internal(r, mu, nu);
  double Mhat = sampler_output["Mhat"];
  double p = 2*nu/( 2*mu*nu +1+nu );
  double log_B= log_bound( mu,  nu,  p)[0];
  
  double output;
  
  if(nu>= 1){
    output = log(Mhat)- log_B - mu;
  } else{
    output = log(Mhat) - log_B;
  }
  
  return 1/exp(output);
  
}

