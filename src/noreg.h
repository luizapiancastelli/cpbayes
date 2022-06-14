#ifndef NOREG_H
#define NOREG_H

#include <Rcpp.h>

double proposal_adjust(double current_var, double accept_rate, int nprops);
double logqcomp(NumericVector y, double mu, double nu);

#endif