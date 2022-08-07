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
  

  #### Regression case (in either)  ----------------------------------------------------------
  if(ncol(X_mu) + ncol(X_nu) >2){
  
  sigma_mu = rep(0.1, ncol(X_mu))
  sigma_nu = rep(0.2, ncol(X_nu))
  
  mcmc_raw = exchange_reg( c(log(mean(y)), rep(0, ncol(X_mu)-1)),  c(0.1, rep(0, ncol(X_nu)-1)), y, X_mu, X_nu, burnin, niter, sigma_mu, sigma_nu, prior_mu, prior_nu)
  
  ac_matrix = matrix(NA, nrow = max(c(ncol(X_mu), ncol(X_nu))), ncol=2)
  ac_matrix[1:ncol(X_mu),1] = mcmc_raw$ac_rates_mu
  ac_matrix[1:ncol(X_nu),2] = mcmc_raw$ac_rates_nu
  
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
#' regression in the location and dispersion parameters.
#' 
#' @importFrom stats model.matrix model.frame model.response
#' 
#' @param formula_beta Object of class \code{formula}. Regression structure for \eqn{\mu}, the location parameter.
#' @param formula_nu Object of class \code{formula}. Regression structure for \eqn{\nu}, the dispersion parameter.
#' @param data Data frame containing the variables used in \code{formula_beta} and \code{formula_nu}. 
#' @param burnin Integer. Number of iterations to discard as burn-in period.
#' @param niter Integer. Number of iterations to store (additional to \code{burnin}).
#' @param prior_mu Named list. Optional argument with \code{mean} and \code{sd} hyperparameters of Normal prior. Unless there is no regression structure, 
#' in which case \code{shape} =2, \code{rate} =2 is expected.
#' @param prior_nu Named list. Equivalent to \code{prior_mu} for the dispersion regression. See details section for more information on prior specification.
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
#' \bold{COM-Poisson regression}\cr
#' Assumes the COM-Poisson regression structure specified by \href{10.1214/20-BA1230}{Benson and Friel (2021)} where
#' \deqn{y_i ~ COM-Poisson(\mu_i, \nu_i)}  
#' and \eqn{\mu_i = \log( \beta_\mu  X_\mu^T)}, \eqn{\nu_i = \log( \beta_nu  X_\nu^T)}.
#' Covariate data is included in the model matrices \eqn{ X_\mu}, \eqn{X_\nu} parsed with \code{formula} and
#' \eqn{\beta_\mu = (\beta_{\mu,0},\beta_{\mu,1},..., \beta_{\mu,n_{\mu}} )}, and \eqn{\beta_nu = (\beta_{\nu,0},\beta_{\nu,1},..., \beta_{\nu,n_{\nu} } )}.
#' 
#' The likelihood of count observations \eqn{y_i} given \eqn{\mu_i}, \eqn{\nu_i} is then
#' \deqn{f(y_i|\mu_i, \nu_i) = (\mu_i^y_i/y_i!)^\nu_i 1/Z(\mu_i, \nu_i) } 
#' where \eqn{Z(\mu_i, \nu_i) = \sum_{y=0}^\infty (\mu_i^y_i/y_i!)^\nu_i } is an intractable normalisation constant. 
#' 
#' \bold{Markov Chain Monte Carlo}\cr
#' Sampling the posterior model parameters is done with the exchange algorithm (Murray, Ghahramani and Ghahramani (2006)),
#' an approach for doubly-intractable Bayesian inference that overcomes unavailability of the likelihood's normalising 
#' constant via introduction of auxiliary draws from the intractable distribution.
#'
#' The augmented posterior model is
#' \deqn{\pi(\beta_\mu, \beta_\nu, \beta_\mu', \beta_\nu', y'| y) \propto f(y|\beta_\mu, \beta_\nu) f(y'|\beta_\mu', \beta_\nu')\pi(\beta_\mu)\pi(\beta_\nu)h(\beta_mu'|\beta_\mu)h(\beta_nu'|\beta_nu) } 
#' where \eqn{\pi(\beta_\mu) = \prod_{i=1}^{n_\mu} \pi(\beta_{\mu, i}|} \code{mean}, \code{sd}) and \eqn{\pi(\beta_\nu) = \prod_{j=1}^{n_\nu} \pi(\beta_{\nu, j}|} \code{mean}, \code{sd}) are
#' the prior distributions of \eqn{\beta_\mu}, \eqn{\beta_\nu} elements. \cr
#' The auxiliary likelihood draws \eqn{y'} are sampled using \code{\link{rcompois}} with parameters \eqn{\beta_\mu'}, \eqn{\beta_\nu'}
#' proposed according to \eqn{h(\beta_{mu', i}|\beta_{mu, i}) = h(\beta_{nu', j}|\beta_{nu,j})}. Specifically, \eqn{h(\beta_{\mu', i}|\beta_{\mu, i} \equiv} 
#' Normal(\eqn{\beta_{\mu,i}, \sigma_{\mu,i}}) for \eqn{i=1, ..., n_\mu}, and similarly \eqn{h(\beta_{\nu, j}, \beta_{\nu', j})} ~
#' Normal(\eqn{\beta_{\nu,j}, \sigma_{\nu, j}}) for \eqn{j=1, ..., n_\nu}. \cr
#'
#' 
#' An automatic vanishing calibration of \eqn{\sigma_{\mu,i}} and \eqn{\sigma_{\nu, j}} proposal variabilities is implemented targetting 
#' the rule of thumb acceptance rate of 44\% to accept/reject proposed values appropriately. \cr
#' 
#' The prior distribution on location regressor effects are \eqn{\pi(\beta_{\mu, i}|} \code{mean}, \code{sd}) for all \eqn{i} with
#' the \code{mean}, \code{sd} values supplied in \code{prior_mu}. Default values are \code{mean} =0, \code{sd} =1. Equivalently
#' for the dispersion regression \eqn{\pi(\beta_{\nu, j}|} \code{mean}, \code{sd}) priors are Normally distributed with the \code{mean} and \code{sd}
#' provided in the \code{prior_nu} list, that also defaults to \code{mean} =0, \code{sd} =1. \cr
#' 
#' \bold{No regression case}\cr
#' The case of no regression, i.e, \eqn{n_\mu = n_\nu =0} is treated as a special case, avoiding the logarithm link and updating
#' \eqn{\mu>0} and \eqn{\nu >0} in place of \eqn{exp(\beta_{\mu,0})} and  \eqn{exp(\beta_{\nu,0})}. The same algorithmic strategy applies
#' with a change in the \emph{proposal and prior distributions}. Candidate \eqn{\mu}, \eqn{\nu} values are simulated from log-Normal distributions
#' with median given by the current parameter state and scale parameters \eqn{\sigma_\mu}, \eqn{\sigma_\nu}. That is, \eqn{h(\mu'|\mu)},\eqn{h(\nu'|\nu)}
#' are log-Normal(\eqn{\log(\mu), \sigma_\mu}), log-Normal(\eqn{\log(\nu), \sigma_\nu}). The change in prior specification acommodating \eqn{\mu>0} and \eqn{\nu >0} 
#' is to \eqn{\pi(\mu)} ~ Gamma(\code{shape}, \code{rate}) with \code{prior_mu} hyperparameters and also for \eqn{\pi(\nu)}, with \code{shape} and \code{rate} in \code{prior_nu}.
#' 
#' 
#' 
#' @author Luiza Piancastelli \email{luiza.piancastelli@@ucdconnect.ie}
#' 
#'  
#' @references
#' Benson, A. and Friel, N. (2021) Bayesian Inference, Model Selection and Likelihood Estimation using Fast Rejection Sampling: The Conway-Maxwell-Poisson Distribution.
#' Bayesian Analysis (16) 905-931.
#' 
#' Murray, I., Ghahramani, Z. and MacKay, D. (2006) MCMC for doubly-intractable distributions.
#' Proceedings of the 22nd Annual Conference on Uncertainty in Artificial Intelligence (UAI-06), AUAI Press.
#' 
#' 
#' 
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

#' Summarising COM-Poisson GLM model fits 
#' 
#' Plotting and summarising \code{methods} for \code{cpbayes} objects obtained with \code{\link{fitcpbayes}}.
#' 
#' 
#' @details 
#' 
#' The \bold{\code{summary.cpbayes}} method provides the posterior mean, standard deviation
#' and \%5, \%50 and \%95 quantiles of \code{cpbayes} MCMC draws for each model parameter.
#' The diagnostic tools 'number of efficient samples' \code{neff} and the 'potential scale reduction factor' \code{Rhat} (Gelman and Rubin (1992))
#' are also provided, computed with \code{coda::effectiveSize} and \code{coda::gelman.diag}. \cr
#' 
#' \code{neff} is the sample size adjusted for autocorrelation of the Markov Chain. Please refer to \code{coda} for further details on computation.
#' \code{Rhat} is included if \code{cpbayes} takes \code{nchains} >1, a convergence diagnostic tool 
#' where values close to 1 usually indicate the multiple chains to have converged to the same stationary distribution. 
#' 
#' \bold{\code{plot.cpbayes}} produces posterior density plots of \code{cpbayes} model parameters if \code{type} = "density".
#' If \code{nchains}>1 MCMC draws are combined for plotting. Traces (iteration versus parameter value) are available with \code{type} = "trace"
#' with lines colored by chain runs. 
#' 
#' 
#' @references 
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
#' 
#' 
#' 
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr summarise group_by
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom coda mcmc as.mcmc.list effectiveSize gelman.diag
#' 
#' @param object cpbayes object
#' @param ... additional arguments
#' @rdname cpbayesmethods
summary.cpbayes = function(object, ...){
  
  mcmc = object
  nchains = length(mcmc)
  niter = nrow(mcmc[[1]][[1]])
  
  mu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "mu")))
  nu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "nu")))
  
  if(ncol(mu_chain)==1 & ncol(nu_chain)==1){
    names(mu_chain) = 'mu'
    names(nu_chain)='nu'
  } else {
    names(mu_chain)= paste0('betamu', 1:ncol(mu_chain))
    names(nu_chain)= paste0('betanu', 1:ncol(nu_chain))
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



#' @importFrom magrittr %>%
#' @importFrom dplyr summarise group_by
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_density theme_bw facet_wrap
#' 
#' @inheritParams summary.cpbayes
#' @param x cpbayes object
#' @param type trace or density plot 
#' @param ... additional arguments
#' @rdname cpbayesmethods
plot.cpbayes = function(x, type = "density", ...){
  
  mcmc = x
  nchains = length(mcmc)
  niter = nrow(mcmc[[1]]$mu)
  
  mu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "mu")))
  nu_chain = data.frame(do.call(rbind, lapply(mcmc, "[[", "nu")))
  
  if(ncol(mu_chain)==1 & ncol(nu_chain)==1){
    names(mu_chain) = 'mu'
    names(nu_chain)='nu'
  } else {
    names(mu_chain)= paste0('betamu', 1:ncol(mu_chain))
    names(nu_chain)= paste0('betanu', 1:ncol(nu_chain))
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

#' BIC estimate
#' 
#' Estimates the Bayesian Information Criteria of a \code{cpbayes} model fit.
#' 
#' @details 
#' 
#' Exact BIC computation is not possible with a COM-Poisson likelihood model due to intractability
#' of the normalisation constant \eqn{Z(\mu, \nu)}. Benson and Friel (2021) propose a statistical estimator
#' for this quantity based on the acceptance rate of their rejection sampling algorithm (\code{rcompois}). 
#' 
#' Using \eqn{\widehat{M}^{(r)} =n_r/r}, the number of envelope draws \eqn{n_r} required for \eqn{r} acceptances,
#' an unbiased estimator of \eqn{1/Z(\mu, \nu)} can be obtained and \eqn{\widehat{f}^{(r)}(y|\mu, \nu)} is computed.
#' Please see Section 3 of the indicated reference for more details.
#' 
#' Obtaining an approximate BIC follows from recording the maximizer of \eqn{\widehat{f}^{(r)}(y|\mu_i, \nu_i)} for the MCMC
#' iterations \eqn{i=1, 2, ...}. Note that there is no additional cost in this procedure since \code{rcompois}
#' draws are part of the \code{fitcpbayes}'s posterior sampling via the exchange algorithm. 
#' 
#' The statistical estimator for BIC is then 
#' \deqn{\widehat{BIC}^{(r)} = k\log(n) - 2\max \widehat{f}^{(r)}(y|\mu, \nu).}
#' 
#' 
#' 
#' @references 
#' Benson, A. and Friel, N. (2021) Bayesian Inference, Model Selection and Likelihood Estimation using Fast Rejection Sampling: The Conway-Maxwell-Poisson Distribution.
#' Bayesian Analysis (16) 905-931.
#' 
#' @param object cpbayes object
#' @param ... additional arguments
BIC= function(object, ...){
 
  logliks = do.call(rbind,lapply(object, "[[", "loglik"))
  k = ncol(object[[1]][[1]]) + ncol(object[[1]][[2]])
  
  if("matrix" %in% class(logliks)){
    logliks = rowSums(logliks)
  }
  
  bic = k*log(length(object[[1]]$y)) - 2*max(logliks)
  return(bic)
}
