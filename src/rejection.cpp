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

//' COM-Poisson rejection sampling
//'
//' @param n number of samples
//' @param mu location parameter
//' @param nu dispersion parameter
//' @return integer vector
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
                      _["n_proposed"] = n_proposed);
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

