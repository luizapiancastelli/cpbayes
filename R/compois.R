
#' @param y Count vector.
#' @param log Logical; if TRUE, returns the log-probabilities.
#' @inheritParams rcompois
#' @rdname compois
dcompois = function(y, mu, nu, log=FALSE){
  log_q = nu*(y*log(mu) - log(factorial(y)));
  z = Ztrunc(10^(-5), mu, nu);
  if(log){
    return(log_q - log(z))  
  } else{
    return(exp(log_q)/z)
  }

}

#' @param q Numeric value. Quantile of probability mass function.
#' @inheritParams dcompois
#' @param lower.tail Logical; if TRUE \eqn{P(X >x)}, otherwise \eqn{P(X \le x)}.
#' @param log.p Logical; if TRUE, returns the log-probabilities.
#' @rdname compois
pcompois = function(q, mu, nu, lower.tail = TRUE, log.p = FALSE){
  
  y = 0:q
  log_lower = dcompois(y, mu, nu, TRUE);
  p_upper = sum(exp(log_lower));
  
  if(lower.tail){ ## P(X <= x)
    
    if(log.p){
      return( log(p_upper) )
    }else{
      return(p_upper)
    }
  } else {
    
    p_lower= 1- p_upper
    
    if(log.p){
      return(log(p_lower) )
    } else {
      return( p_lower )
    }
    
  }
  
}

#' @param p Numeric value; probability.
#' @inheritParams rcompois
#' @rdname compois
qcompois = function(p, mu, nu){
  
  y = 0
  p_sum =  dcompois(y, mu, nu)
 
  while(p_sum < p){
    y = y + 1;
    p_sum = p_sum + dcompois(y, mu, nu)
  }
  return(y)
  
}

