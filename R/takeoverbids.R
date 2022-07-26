#' Takeoverbids data
#' 
#' The \code{takeoverbids} data (Jaggia and Thosar (1993)) is a record of the number of bids received by U.S. firms 
#' in the period 1978-1985 that were taken over within 1 year of the initial offer alongside potential covariates. \cr
#' 
#' Additional information is as described in Saez-Castillo & Conde-Sanchez, A. (2013).
#' 
#' @format
#' \describe{
#' \item{ \code{numbids}}{Number of bids excluding the initial offer.}
#' \item{ \code{bidprem} }{Bid price divided by price 14 working days before bid.}
#' \item{\code{finrest}}{Did management propose changes in ownership structure? (1=YES,0=NO)}
#' \item{\code{insthold}}{Percentage of stock held by institutions.}
#' \item{\code{leglrest}}{Indicator variable for legal defence by lawsuit.}
#' \item{\code{rearest}}{Did management propose changes in asset structure? (1=YES,0=NO)}
#' \item{\code{regulatn}}{Did federal regulators intervene? (1=YES,0=NO)}
#' \item{\code{size}}{Total book value of assets. (USD billions)}
#' \item{\code{takeover}}{Indicates if the company was being taken over.}
#' \item{\code{weeks}}{Time in weeks between the initial and final offers.}
#' \item{\code{whtknght}}{Did management invite friendly third-party bid?  (1=YES,0=NO)}
#' \item{\code{sizesq}}{Book value squared.}
#' }
#' 
#' @references 
#' Jaggia, S. and S. Thosar (1993). Multiple bids as a consequence of target management resistance: 
#' A count data approach. Review of Quantitative Finance and Accounting 3(4), 447-457.
#' 
#' Saez-Castillo, Antonio & Conde-Sanchez, A. (2013). A hyper-Poisson regression model for overdispersed 
#' and underdispersed count data. Computational Statistics & Data Analysis. 61. 148-157.
"takeoverbids"
