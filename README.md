
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
#> 
#> Attaching package: 'cpbayes'
#> The following object is masked from 'package:stats':
#> 
#>     BIC
ls("package:cpbayes") #List of currently exported functions
#> [1] "BIC"          "dcompois"     "fitcpbayes"   "inventory"    "pcompois"    
#> [6] "qcompois"     "rcompois"     "takeoverbids" "Zhat"
```

### COM-Poisson distribution

`cpbayes` provides a family of functions designed to work with the
COM-Poisson distribution, namely `rcompois`, `dcompois`, `pcompois` and
`qcompois`. A call to `help` renders documentation where further details
on their implementation can be found.

`rcompois` is the key component among these, an implementation of the
fast-rejection sampler of Benson and Friel (2021) that gives from a
COM-Poisson with location
![\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
and dispersion parameter
![\nu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu "\nu").
Currently `dcompois`, `pcompois` and `qcompois` involve approximation of
the normalisation constant via truncation, included so that `cpbayes`
gives a complete framework to work with the COM-Poisson.

The next code chunk exemplifies simulating over and underdispersed
counts with `rcompois` and computation of some COM-Poisson
probabilities.

``` r
Y1 = rcompois(n =200, mu =1, nu = 0.2) #nu < 1 renders overdispersed count data
Y2 = rcompois(n =200, mu =1, nu = 1.5) #nu > 1 underdispersion
c(mean(Y1), var(Y1))
#> [1] 2.805000 6.127613
c(mean(Y2), var(Y2))
#> [1] 0.7450000 0.5929397
#approximate mean and variance are E(Y) ~ mu + 1/2nu -1/2 and Var(Y) ~ mu/nu

#P(Y = 0)
dcompois(0, mu =1, nu = 0.2)
#> [1] 0.1903481
dcompois(0, mu =1, nu = 1.5)
#> [1] 0.4113677

#P(Y <= 1) for Y1, Y2
pcompois(1, mu = 1, nu= 0.2)
#> [1] 0.3806962
pcompois(1, mu = 1, nu=1.5)
#> [1] 0.8227354

#Median of distribution 
qcompois(0.5, 1, 0.2)
#> [1] 2
qcompois(0.5, 1, 1.5)
#> [1] 1
```

### COM-Poisson regression

In addition to the COM-Poisson distribution functions, `cpbayes` fits
generalised linear models to COM-Poisson distributed counts in the
Bayesian framework. Regression structure via the logarithm link is
supported in the location and dispersion parameters. This model assumes
counts
![y_1, \cdots, y_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_1%2C%20%5Ccdots%2C%20y_n "y_1, \cdots, y_n")
to be COM-Poisson distributed with likelihood

![f(y_i\|\mu_i, \nu_i) =  \left( \dfrac{\mu_i^{y_i}}{y_i!} \right)^{\nu_i} \frac{1}{Z(\mu_i, \nu_i)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;f%28y_i%7C%5Cmu_i%2C%20%5Cnu_i%29%20%3D%20%20%5Cleft%28%20%5Cdfrac%7B%5Cmu_i%5E%7By_i%7D%7D%7By_i%21%7D%20%5Cright%29%5E%7B%5Cnu_i%7D%20%5Cfrac%7B1%7D%7BZ%28%5Cmu_i%2C%20%5Cnu_i%29%7D "f(y_i|\mu_i, \nu_i) =  \left( \dfrac{\mu_i^{y_i}}{y_i!} \right)^{\nu_i} \frac{1}{Z(\mu_i, \nu_i)}")

where
![Z(\mu_i, \nu_i) = \sum\_{y=0}^\infty q(y_i\|\mu_i, \nu_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z%28%5Cmu_i%2C%20%5Cnu_i%29%20%3D%20%5Csum_%7By%3D0%7D%5E%5Cinfty%20q%28y_i%7C%5Cmu_i%2C%20%5Cnu_i%29 "Z(\mu_i, \nu_i) = \sum_{y=0}^\infty q(y_i|\mu_i, \nu_i)"),
![q(y_i\|\mu_i, \nu_i) = \left( \frac{\mu_i^{y_i}}{y_i!} \right)^{\nu_i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q%28y_i%7C%5Cmu_i%2C%20%5Cnu_i%29%20%3D%20%5Cleft%28%20%5Cfrac%7B%5Cmu_i%5E%7By_i%7D%7D%7By_i%21%7D%20%5Cright%29%5E%7B%5Cnu_i%7D "q(y_i|\mu_i, \nu_i) = \left( \frac{\mu_i^{y_i}}{y_i!} \right)^{\nu_i}")
and the regression
![\mu_i = \exp \left( \beta\_{\mu,0} + \sum\_{j=1}^b \beta\_{\mu, j} x^\mu\_{ij} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu_i%20%3D%20%5Cexp%20%5Cleft%28%20%5Cbeta_%7B%5Cmu%2C0%7D%20%2B%20%5Csum_%7Bj%3D1%7D%5Eb%20%5Cbeta_%7B%5Cmu%2C%20j%7D%20x%5E%5Cmu_%7Bij%7D%20%5Cright%29 "\mu_i = \exp \left( \beta_{\mu,0} + \sum_{j=1}^b \beta_{\mu, j} x^\mu_{ij} \right)")
and
![\nu_i = \exp \left( \beta\_{\nu,0} + \sum\_{k=1}^c \beta\_{\nu, k} x^\nu\_{ik} \right)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu_i%20%3D%20%5Cexp%20%5Cleft%28%20%5Cbeta_%7B%5Cnu%2C0%7D%20%2B%20%5Csum_%7Bk%3D1%7D%5Ec%20%5Cbeta_%7B%5Cnu%2C%20k%7D%20x%5E%5Cnu_%7Bik%7D%20%5Cright%29 "\nu_i = \exp \left( \beta_{\nu,0} + \sum_{k=1}^c \beta_{\nu, k} x^\nu_{ik} \right)").
![\boldsymbol{x}^\mu_1, \cdots, \boldsymbol{x}^\mu_b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bx%7D%5E%5Cmu_1%2C%20%5Ccdots%2C%20%5Cboldsymbol%7Bx%7D%5E%5Cmu_b "\boldsymbol{x}^\mu_1, \cdots, \boldsymbol{x}^\mu_b")
and
![\boldsymbol{x}^\nu_1, \cdots, \boldsymbol{x}^\nu_c](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7Bx%7D%5E%5Cnu_1%2C%20%5Ccdots%2C%20%5Cboldsymbol%7Bx%7D%5E%5Cnu_c "\boldsymbol{x}^\nu_1, \cdots, \boldsymbol{x}^\nu_c")
are observed covariates accounted in
![\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
and
![\nu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu "\nu")
for explaining
![\boldsymbol{y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7By%7D "\boldsymbol{y}").

The function `fitcpbayes` runs Markov Chain Monte Carlo for
![\boldsymbol{\beta\_\mu}, \boldsymbol{\beta\_\nu}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cbeta_%5Cmu%7D%2C%20%5Cboldsymbol%7B%5Cbeta_%5Cnu%7D "\boldsymbol{\beta_\mu}, \boldsymbol{\beta_\nu}")
Bayesian inference using the exchange algorithm (Murray (2006)). Further
details on the exchange framework for intractable likelihoods and
arguments to `fitcpbayes` are rendered by `help("fitcpbayes")`.

Code below provides and example of its usage with the `takeoverbids`
data, a data set on the number of bids `numbids` received by U.S. firms
in the period 1978-1985 that were taken over within 1 year of the
initial offer. This is modeled alongside the potential covariates
`whtknght` and `size` described in `help(takeoverbids)`.

`fitcpbayes`’s basic arguments include `R` formulas for
![\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
and
![\nu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cnu "\nu"),
the data set (a `data.frame`), number of iterations to discard as
`burnin` and to store (`niter`). `nchains` specifies how many MCMC
chains to run and parallel computation can be enabled by setting
`ncores`. `fitcpbayes` outputs an object of class `cpbayes` for which
generic `summary` and `plot` methods are available.

``` r
data("takeoverbids")
mcmc = fitcpbayes(numbids ~ whtknght, numbids ~ size, takeoverbids, 10000, 10000, nchains =3)

summary(mcmc) #neff is the number of efficient draws, and Rhat the potential scale reduction factor, estimated if nchains>1
#> # A tibble: 4 × 8
#>   param     mean     sd     Q5    Q50    Q95  rhat  neff
#>   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
#> 1 betamu1  0.336 0.0990  0.169  0.341  0.489  1.00  520.
#> 2 betamu2  0.450 0.108   0.275  0.448  0.628  1.00  597.
#> 3 betanu1  0.699 0.263   0.266  0.699  1.13   1.00  897.
#> 4 betanu2 -0.208 0.0654 -0.330 -0.200 -0.118  1.00 1084.
```

``` r
plot(mcmc) #density plots using combined posterior draws
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
plot(mcmc, type = 'trace') #trace plots colored by chain run
```

<img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

`cpbayes` also provides BIC estimation for comparing between models
fits. This is automatic from a `cpbayes` object, where likelihood
estimation occurs is linked to the acceptance rate of COM-Poisson
rejection sampling. For more details, please see the indicated reference
or `help("BIC")`. Let’s consider another regression structure for the
`takeoverbids` problem and use this tool for model selection.

``` r
mcmc2 = fitcpbayes(numbids ~ whtknght + bidprem, numbids ~ size, takeoverbids, 10000, 10000, nchains =3)

BIC(mcmc)
#> [1] 374.1528
BIC(mcmc2)
#> [1] 380.2662
```

It looks like the data further supports `numbids ~ whtknght`
![\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
regression in comparison to `numbids ~ whtknght + bidprem`.

#### No regression

A special case of `fitcpbayes` is implemented for when we are interested
in Bayesian inference for the
![\mu, \nu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%2C%20%5Cnu "\mu, \nu")
parameters with no regression structure. As an example, let’s revisit
`Y2` and see how the model fit supports the true values used to generate
the data, these were
![\mu =1, \nu = 1.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20%3D1%2C%20%5Cnu%20%3D%201.5 "\mu =1, \nu = 1.5").

``` r
y2 = data.frame('y' = Y2)
mcmc3 = fitcpbayes(y ~ 1, y ~ 1, y2, 10000, 10000, nchains =1)

summary(mcmc3)
#> # A tibble: 2 × 7
#>   param  mean    sd    Q5   Q50   Q95  neff
#>   <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 mu    0.859 0.117 0.649 0.869  1.03  138.
#> 2 nu    1.32  0.280 0.881 1.30   1.80  185.
```

Looks like the true values have high posterior probability as expected!
