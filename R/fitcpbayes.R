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
              'ac_rates' = ac_matrix, 'loglik' = mcmc_raw$loglik, 'y' = y)
  
  
  } else { #No regression case --------------------------------------------------
    
    hyperparams = list('shape' = c(prior_mu$shape, prior_nu$shape),'rate' =  c(prior_mu$rate, prior_nu$rate))
    sigma = c(0.1, 0.1)
    
    mcmc_raw =exchange_noreg(y, c(1, 1), sigma, burnin, niter, hyperparams) 
    mcmc = list("mu" = matrix(mcmc_raw$mu, ncol =1), 
                "nu" = matrix(mcmc_raw$nu, ncol =1), 
                'ac_rates' = mcmc_raw$ac_rates, 
                'loglik' = mcmc_raw$loglik, 'y' = y)
  }
  return(mcmc)
}

#' Bayesian fit of COM-Poisson GLMs
#' 
#' \code{fitcpbayes} fits generalised linear models to COM-Poisson distributed count data, supporting
#' regression specification in the location and dispersion parameters.
#' 
#' @importFrom stats model.matrix model.frame model.response
#' 
#' @param formula_beta Object of class \code{\link[stats]{formula}}. Regression structure for \eqn{\mu}, the location parameter.
#' @param formula_nu Object of class \code{\link[stats]{formula}}. Regression structure for \eqn{\nu}, the dispersion parameter.
#' @param data Data frame containing the variables used in \code{formula_beta} and \code{formula_nu}. 
#' @param burnin Integer. Number of iterations to discard as burn-in period.
#' @param niter Integer. Number of iterations to store (additional to \code{burnin}).
#' @param prior_mu optional named list with \code{shape} and \code{rate} if modelling only the intercept, 
#' or \code{mean} and \code{sd} if any covariates in \code{formula_beta}. Default values are \code{shape} =2, \code{rate} =2, and \code{mean} =0, \code{sd} =1. 
#' @param prior_nu  optional named list with \code{shape} and \code{rate} if modelling only the intercept, 
#' or \code{mean} and \code{sd} if any covariates in \code{formula_nu}. Default values are \code{shape} =2, \code{rate} =2, and \code{mean} =0, \code{sd} =1. 
#' See details for more information on prior specification.
#' @param nchains Integer. Number of chains to run in parallel, defaults to a single chain.
#' @param ncores Integer. Number of cores for parallel computation.
#' 
#' @return Object of class \code{cpbayes}. \code{nchains} lists, each containing the following.
#' \describe{
#' \item{\code{mu}}{(\code{niter} x \eqn{n_\mu}) matrix. Stored MCMC samples of parameters in \code{formula_beta}.}
#' \item{\code{nu}}{(\code{niter} x \eqn{n_\nu}) matrix. Stored MCMC samples of parameters in \code{formula_nu}.}
#' \item{\code{ac_rates}}{( \eqn{max( n_{\mu}, n_{\nu}) x 2}) matrix) of acceptance rates.}
#' \item{\code{loglik}}{(\code{niter} x \eqn{n}) or (\code{niter} x \eqn{1}) matrix.
#'  Log-likelihood values of response data \eqn{y_1, ... , y_n} at the MCMC iterations.   }
#' \item{\code{y}}{Response data.}
#' }
#' 
#' @seealso summary.cpbayes, plot.cpbayes, BIC.cpbayes, rcompois
#' 
#' @examples 
#' \dontrun{
#' burnin= 10000
#' niter = 20000
#'data("inventory") ## No regression in mu or nu, single MCMC chain
#'mcmc_noreg =fitcpbayes(sales~1, sales~1, inventory, burnin, niter, nchains=1)
#'
#' data("takeoverbids")
#' library(dplyr)
#' data=  takeoverbids %>% select(numbids, whtknght, size)
#' formula_beta = numbids ~ whtknght
#' formula_nu = numbids ~ size
#' #Parallel computation, regression on both parameters  
#' mcmc_reg = fitcpbayes(formula_beta, formula_nu, data, burnin, niter, nchains =3, ncores =3)
#'}
#'
#' @details
#' MCMC for the posterior model with COM-Poisson likelihood regression of \href{10.1214/20-BA1230}{Benson and Friel (2021)} where
#' \deqn{y_i ~ COM-Poisson(\mu_i, \nu_i)}  
#' \eqn{\mu_i = \log( \beta_\mu  X_\mu^T)} and \eqn{\nu_i = \log( \beta_nu  X_\nu^T)}.
#' 
#' Covariates in the model matrices \eqn{ X_\mu}, \eqn{X_\nu} are modelled alongside intercepts, so
#' \eqn{\beta_\mu = (\beta_{\mu,0},\beta_{\mu,1},..., \beta_{\mu,n_{\mu}} )}, and \eqn{\beta_nu = (\beta_{\nu,0},\beta_{\nu,1},..., \beta_{\nu,n_{\nu} } )}. \cr
#' If no regression is supplied, \eqn{(\beta_\mu = (\beta_{\mu,0} )} or \eqn{(\beta_\nu = (\beta_{\nu,0} )} the parameter is modelled directly.
#' 
#' The exchange algorithm ...
#' In this case the priors ...
#' 
#' @author Luiza Piancastelli \email{luiza.piancastelli@@ucdconnect.ie}
#' 
#' 
#' @references
#' Benson, A. and Friel, N. (2021) Bayesian Inference, Model Selection and Likelihood Estimation using Fast Rejection Sampling: The Conway-Maxwell-Poisson Distribution.
#' Bayesian Analysis (16) 905-931.
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

#' Summarise cpbayes MCMC output
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise group_by
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom coda mcmc as.mcmc.list effectiveSize gelman.diag
#' 
#' @param object cpbayes object
#' @param ... additional arguments
summary.cpbayes = function(object, ...){
  
  mcmc = object
  nchains = length(mcmc)
  niter = nrow(mcmc[[1]][[1]])
  
  mu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "mu")))
  nu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "nu")))
  
  if(ncol(mu_chain)>1){
    names(mu_chain)= paste0('betamu', 1:ncol(mu_chain))
  }else{
    names(mu_chain) = 'mu'
  }
  if(ncol(nu_chain)>1){
    names(nu_chain)= paste0('betanu', 1:ncol(nu_chain))
  } else{
    names(nu_chain)='nu'
  }
  
  params = cbind(mu_chain, nu_chain)
  params$iter = 1:niter
  params$chain = sort(rep(1:nchains, niter))
  
  params =pivot_longer(params, cols = -c(.data$iter,.data$chain), names_to = 'param')
  
  stats = params %>% group_by(.data$param)%>%  dplyr::summarise(mean = mean(.data$value), 
                                                    "sd" = sd(.data$value),
                                                  'Q5' = quantile(.data$value, 0.05),
                                                  'Q50' = quantile(.data$value, 0.5),
                                                  'Q95' = quantile(.data$value, 0.95))
  mcmc_coda = list();
  for(i in 1:nchains){
    mcmc_coda[[i]] = coda::mcmc(cbind(mcmc[[i]][[1]],mcmc[[i]][[2]]))
  }
  mcmc_coda = coda::as.mcmc.list(mcmc_coda)
  
  neff = coda::effectiveSize(mcmc_coda)
  if(nchains>1){
  rhat = coda::gelman.diag(mcmc_coda)
  stats$rhat = rhat$psrf[,1]
  }
  
  stats$neff = neff
  return(stats)
  
}


#' Plot cpbayes MCMC output
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise group_by
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_density theme_bw facet_wrap
#' 
#' @param x cpbayes object
#' @param type trace or density plot 
#' @param ... additional arguments
plot.cpbayes = function(x, type = "density", ...){
  
  mcmc = x
  nchains = length(mcmc)
  niter = nrow(mcmc[[1]]$mu)
  
  mu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "mu")))
  nu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "nu")))
  
  if(ncol(mu_chain)>1){
    names(mu_chain)= paste0('betamu', 1:ncol(mu_chain))
  }else{
    names(mu_chain) = 'mu'
  }
  if(ncol(nu_chain)>1){
    names(nu_chain)= paste0('betanu', 1:ncol(nu_chain))
  } else{
    names(nu_chain)='nu'
  }
  
  params = cbind(mu_chain, nu_chain)
  params$iter = 1:niter
  params$chain = as.factor(sort(rep(1:nchains, niter)))
  
  params =pivot_longer(params, cols = -c(.data$iter,.data$chain), names_to = 'param')
  
  if(type == 'density'){
  cpbayesplot = ggplot(params, aes(x = .data$value))+
    geom_density(adjust = 5)+facet_wrap(~.data$param, scales = 'free')+
    theme_bw()
  }
  if(type == "trace"){
    cpbayesplot = ggplot(params, aes(x = .data$iter, y = .data$value, col = .data$chain))+
      facet_wrap(~.data$param, scales = 'free')+theme_bw()+
      geom_line()
  }

  return(cpbayesplot)
  
}

#' BIC information criteria 
#' @param object cpbayes object
#' @param ... additional arguments
BIC.cpbayes = function(object, ...){
 
  logliks = do.call(rbind,lapply(object, "[[", "loglik"))
  k = ncol(object[[1]][[1]]) + ncol(object[[1]][[2]])
  
  if("matrix" %in% class(logliks)){
    logliks = rowSums(logliks)
  }
  
  bic = k*log(length(object[[1]]$y)) - 2*max(logliks)
  return(bic)
}
