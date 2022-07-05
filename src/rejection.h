#ifndef REJECTION_H
#define REJECTION_H

#include <Rcpp.h>

NumericVector rcompois(int n, double mu, double nu);

List rcompois_internal(int n, double mu, double nu);

NumericVector log_bound(double mu, double nu, double p);

#endif