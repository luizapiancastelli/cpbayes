% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{exchange_reg}
\alias{exchange_reg}
\title{MCMC: regression case}
\usage{
exchange_reg(
  beta_mu_init,
  beta_nu_init,
  y,
  X_mu,
  X_nu,
  burn_in,
  n_iter,
  sigma_mu,
  sigma_nu,
  hyperparams_mu,
  hyperparams_nu
)
}
\arguments{
\item{beta_mu_init}{initial beta_mu}

\item{beta_nu_init}{initial beta_nu}

\item{y}{response vector}

\item{X_mu}{location model matrix}

\item{X_nu}{dispersion model matrix}

\item{burn_in}{integer}

\item{n_iter}{integer}

\item{sigma_mu}{proposal variances}

\item{sigma_nu}{proposal variances}

\item{hyperparams_mu}{named list with 'mean' and 'sd'}

\item{hyperparams_nu}{named list with 'mean' and 'sd'}
}
\description{
MCMC: regression case
}
