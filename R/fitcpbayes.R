#' Runs a single chain of COM-Poisson GLM
#' 
#' @param X_mu parsed regression matrix for location
#' @param X_nu parsed regression matrix for dispersion
#' @param y parsed response variable
#' @param burnin integer. Number of iterations to discard as burn-in period.
#' @param niter integer. Number of iterations to store.
#' @param prior_mu optional named list with 'shape' and 'rate' if no regression on location, or 'mean' and 'sd' if running CP GLM. 
#' Defaults to ('shape' =2, 'rate' =2) or ('mean' =0, 'sd' =1).
#' @param prior_nu optional named list with 'shape' and 'rate' if no regression on dispersion, or 'mean' and 'sd' if running CP GLM.
#' Defaults to ('shape' =2, 'rate' =2) or ('mean' =0, 'sd' =1).
fitcpbayes_single = function(X_mu, X_nu, y, burnin, niter, prior_mu, prior_nu){
  

  #### Regression case  ----------------------------------------------------------
  if(ncol(X_mu) + ncol(X_nu) >2){
  
  sigma_mu = rep(0.1, ncol(X_mu))
  sigma_nu = rep(0.2, ncol(X_nu))
  
  mcmc_raw = exchange_reg( c(log(mean(y)), rep(0, ncol(X_mu)-1)),  c(0.1, rep(0, ncol(X_nu)-1)), y, X_mu, X_nu, burnin, niter, sigma_mu, sigma_nu, prior_mu, prior_nu)
  ac_matrix = cbind(mcmc_raw$ac_rates_mu, mcmc_raw$ac_rates_nu)
  colnames(ac_matrix) = c('mu', 'nu')
  mcmc = list('mu' = mcmc_raw$betamu, 'nu' = mcmc_raw$betanu, 
              'ac_rates' = ac_matrix)
  
  
  } else { #No regression case --------------------------------------------------
    
    hyperparams = list('shape' = c(prior_mu$shape, prior_nu$shape),'rate' =  c(prior_mu$rate, prior_nu$rate))
    sigma = c(0.1, 0.1)
    
    mcmc =exchange_noreg(y, c(1, 1), sigma, burnin, niter, hyperparams) 

  }
  return(mcmc)
}

#' MCMC for COM-Poisson GLM
#' 
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
#' @param nchains integer, number of MCMC chains to run in parallel. Defaults to a single chain.
#' @param ncores integer, number of cores to use when running parallel computation. Defaults to 2.
#' @export
fitcpbayes = function(formula_beta, formula_nu, data, burnin, niter, prior_mu, prior_nu, nchains, ncores){
  
  X_mu = model.matrix(formula_beta, data = data)
  X_nu = model.matrix(formula_nu, data = data)
  y = c(model.response(model.frame(formula_beta,data=data)))
  
  if(ncol(X_mu) + ncol(X_nu) >2){ #Regression case -------------
    
    if(missing(prior_mu)){
      prior_mu = list('mean'=0, 'sd'=1)
    }
    if(missing(prior_nu)){
      prior_nu = list('mean'=0, 'sd'=1)
    }
    
  } else { #No regression ---------------------------------------
    
    if(missing(prior_mu)){
      prior_mu = list('shape'=2, 'rate'=2)
    }
    if(missing(prior_nu)){
      prior_nu = list('shape'=2, 'rate'=2)
    }
    
  }
  
  #Execution:
  if(missing(nchains)){
    mcmc = fitcpbayes_single(X_mu, X_nu, y, burnin, niter, prior_mu, prior_nu)
    mcmc = list(mcmc)
  } else{
    
    if(missing(ncores)){
      ncores =2
    }
    
    cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
    parallel::clusterEvalQ(cl, library("cpbayes"))
    
    mcmc = parallel::parSapply(cl, 1:nchains, 
                                 function(times,X_mu, X_nu, y, burnin, niter, prior_mu, prior_nu){
                                   fitcpbayes_single(X_mu, X_nu, y, burnin, niter, prior_mu, prior_nu)},
                                 X_mu, X_nu, y, burnin, niter, prior_mu, prior_nu,
                                 simplify = F)
    
    parallel::stopCluster(cl) 
  }
  
  class(mcmc) = 'cpbayes'
  return(mcmc)
  
}