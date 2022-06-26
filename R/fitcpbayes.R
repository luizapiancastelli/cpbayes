#' Run COM-Poisson GLM
#' @importFrom stats model.matrix model.frame model.response
#' 
#' @param formula_beta regression formula for location
#' @param formula_nu regression formula for dispersion
#' @param data data frame 
#' @param burnin integer. Number of iterations to discard as burn-in period.
#' @param niter integer. Number of iterations to store.
#' @param prior_mu optional named list with 'shape' and 'rate' if no regression on location, or 'mean' and 'sd' if running CP GLM. 
#' Defaults to ('shape' =2, 'rate' =2) or ('mean' =0, 'sd' =1).
#' @param prior_nu optional named list with 'shape' and 'rate' if no regression on dispersion, or 'mean' and 'sd' if running CP GLM.
#' Defaults to ('shape' =2, 'rate' =2) or ('mean' =0, 'sd' =1).
#' @export
fitcpbayes = function(formula_beta, formula_nu, data, burnin, niter, prior_mu, prior_nu){
  
  
  X_mu = model.matrix(formula_beta, data = data)
  X_nu = model.matrix(formula_nu, data = data)
  y = c(model.response(model.frame(formula_beta,data=data)))
  
  #### Regression case  ----------------------------------------------------------
  if(ncol(X_mu) + ncol(X_nu) >2){
  
  sigma_mu = rep(0.1, ncol(X_mu))
  sigma_nu = rep(0.2, ncol(X_nu))
  
  if(missing(prior_mu)){
    prior_mu = list('mean'=0, 'sd'=1)
  }
  if(missing(prior_nu)){
    prior_nu = list('mean'=0, 'sd'=1)
  }
  
  mcmc_raw = exchange_reg( rep(1, ncol(X_mu)), rep(1, ncol(X_nu)), y, X_mu, X_nu, burnin, niter, sigma_mu, sigma_nu, prior_mu, prior_nu)
  ac_matrix = cbind(mcmc_raw$ac_rates_mu, mcmc_raw$ac_rates_nu)
  colnames(ac_matrix) = c('mu', 'nu')
  mcmc = list('mu' = mcmc_raw$betamu, 'nu' = mcmc_raw$betanu, 
              'ac_rates' = ac_matrix)
  
  
  } else { #No regression case --------------------------------------------------
    
    
    if(missing(prior_mu)){
      prior_mu = list('shape'=2, 'rate'=2)
    }
    if(missing(prior_nu)){
      prior_nu = list('shape'=2, 'rate'=2)
    }
    
    hyperparams = list('shape' = c(prior_mu$shape, prior_nu$shape),'rate' =  c(prior_mu$rate, prior_nu$rate))
    sigma = c(0.1, 0.1)
    
    mcmc =exchange_noreg(y, c(1, 1), sigma, burnin, niter, hyperparams) 

  }
  
  class(mcmc) = 'cpbayes'
  return(mcmc)
  
  
}
