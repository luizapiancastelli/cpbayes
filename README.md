
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cpbayes

Fast rejection sampling and Bayesian inference for COM-Poisson
regression models [Benson and Friel
(2021)](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-3/Bayesian-Inference-Model-Selection-and-Likelihood-Estimation-using-Fast-Rejection/10.1214/20-BA1230.full).
Regression structure on the location and dispersion parameters is
supported and MCMC sampling relies on the exchange algorithm.

Please ensure proper setup of a C++ compiler as our package links to
Rcpp and RcppArmadillo (you might test this by installing these two
first). Uncomment the following command line to install `cpbayes` from
this GitHub repository.

``` r
#devtools::install_github("luizapiancastelli/cpbayes")
library(cpbayes)

ls("package:cpbayes") #List of currently exported functions
#> [1] "fitcpbayes"   "inventory"    "rcompois"     "takeoverbids" "Zhat"
```

Usage of `rcompois` for set parameter values (the fast-rejection
sampler):

``` r
Y = rcompois(n =100, mu =1, nu = 0.5)
```

### COM-Poisson regression

Function `fitcpbayes` fits a COM-Poisson GLM and a quick check of its
arguments is rendered by `help("fitcpbayes")` (please note that
documentation is under development!). The code below provides and
example of its usage with the `takeoverbids` data set (model 4 in
(Benson and Friel (2021)) ) and some methods currently available.

``` r
#help("fitcpbayes")

data("takeoverbids")

mcmc = fitcpbayes(numbids ~ whtknght, numbids ~ size, takeoverbids, 10000, 10000, nchains =3)

summary(mcmc)
#> # A tibble: 4 × 8
#>   param     mean     sd     Q5    Q50    Q95  rhat  neff
#>   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 betamu1  0.333 0.101   0.159  0.336  0.491  1.01  499.
#> 2 betamu2  0.451 0.111   0.270  0.450  0.635  1.01  575.
#> 3 betanu1  0.691 0.265   0.253  0.692  1.13   1.00  742.
#> 4 betanu2 -0.204 0.0644 -0.322 -0.196 -0.116  1.00 1130.
plot(mcmc)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
