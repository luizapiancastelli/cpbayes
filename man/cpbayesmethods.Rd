% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitcpbayes.R
\name{summary.cpbayes}
\alias{summary.cpbayes}
\alias{plot.cpbayes}
\title{Summarising COM-Poisson GLM model fits}
\usage{
\method{summary}{cpbayes}(object, ...)

\method{plot}{cpbayes}(x, type = "density", ...)
}
\arguments{
\item{object}{cpbayes object}

\item{...}{additional arguments}

\item{x}{cpbayes object}

\item{type}{trace or density plot}
}
\description{
Plotting and summarising \code{methods} for \code{cpbayes} objects obtained with \code{\link{fitcpbayes}}.
}
\details{
The \bold{\code{summary.cpbayes}} method provides the posterior mean, standard deviation
and \%5, \%50 and \%95 quantiles of \code{cpbayes} MCMC draws for each model parameter.
The diagnostic tools 'number of efficient samples' \code{neff} and the 'potential scale reduction factor' \code{Rhat} (Gelman and Rubin (1992))
are also provided, computed with \code{coda::effectiveSize} and \code{coda::gelman.diag}. \cr

\code{neff} is the sample size adjusted for autocorrelation of the Markov Chain. Please refer to \code{coda} for further details on computation.
\code{Rhat} is included if \code{cpbayes} takes \code{nchains} >1, a convergence diagnostic tool 
where values close to 1 usually indicate the multiple chains to have converged to the same stationary distribution. 

\bold{\code{plot.cpbayes}} produces posterior density plots of \code{cpbayes} model parameters if \code{type} = "density".
If \code{nchains}>1 MCMC draws are combined for plotting. Traces (iteration versus parameter value) are available with \code{type} = "trace"
with lines colored by chain runs.
}
\references{
Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
}
