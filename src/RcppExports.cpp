// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// update_positive
List update_positive(int type, NumericVector params, NumericVector y, NumericVector sigma, NumericVector shape, NumericVector rate);
RcppExport SEXP _cpbayes_update_positive(SEXP typeSEXP, SEXP paramsSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP shapeSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(update_positive(type, params, y, sigma, shape, rate));
    return rcpp_result_gen;
END_RCPP
}
// exchange_noreg
List exchange_noreg(NumericVector y, NumericVector params0, NumericVector sigma, int n_iter, int burn_in, List hyperparams);
RcppExport SEXP _cpbayes_exchange_noreg(SEXP ySEXP, SEXP params0SEXP, SEXP sigmaSEXP, SEXP n_iterSEXP, SEXP burn_inSEXP, SEXP hyperparamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type params0(params0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< List >::type hyperparams(hyperparamsSEXP);
    rcpp_result_gen = Rcpp::wrap(exchange_noreg(y, params0, sigma, n_iter, burn_in, hyperparams));
    return rcpp_result_gen;
END_RCPP
}
// rcompoisreg
arma::vec rcompoisreg(arma::vec& beta_mu, arma::vec& beta_nu, arma::mat& X_mu, arma::mat& X_nu);
RcppExport SEXP _cpbayes_rcompoisreg(SEXP beta_muSEXP, SEXP beta_nuSEXP, SEXP X_muSEXP, SEXP X_nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type beta_mu(beta_muSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_nu(beta_nuSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_mu(X_muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_nu(X_nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rcompoisreg(beta_mu, beta_nu, X_mu, X_nu));
    return rcpp_result_gen;
END_RCPP
}
// proposal_normal
arma::vec proposal_normal(int index, arma::vec& param, arma::vec& sigma);
RcppExport SEXP _cpbayes_proposal_normal(SEXP indexSEXP, SEXP paramSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(proposal_normal(index, param, sigma));
    return rcpp_result_gen;
END_RCPP
}
// update_beta_mu
List update_beta_mu(int index, arma::vec& beta_mu, arma::vec& beta_nu, arma::vec& y, arma::mat& X_mu, arma::mat& X_nu, arma::vec& sigma_mu, List hyperparams_mu);
RcppExport SEXP _cpbayes_update_beta_mu(SEXP indexSEXP, SEXP beta_muSEXP, SEXP beta_nuSEXP, SEXP ySEXP, SEXP X_muSEXP, SEXP X_nuSEXP, SEXP sigma_muSEXP, SEXP hyperparams_muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_mu(beta_muSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta_nu(beta_nuSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_mu(X_muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_nu(X_nuSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma_mu(sigma_muSEXP);
    Rcpp::traits::input_parameter< List >::type hyperparams_mu(hyperparams_muSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta_mu(index, beta_mu, beta_nu, y, X_mu, X_nu, sigma_mu, hyperparams_mu));
    return rcpp_result_gen;
END_RCPP
}
// rcompois
NumericVector rcompois(int n, double mu, double nu);
RcppExport SEXP _cpbayes_rcompois(SEXP nSEXP, SEXP muSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rcompois(n, mu, nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cpbayes_update_positive", (DL_FUNC) &_cpbayes_update_positive, 6},
    {"_cpbayes_exchange_noreg", (DL_FUNC) &_cpbayes_exchange_noreg, 6},
    {"_cpbayes_rcompoisreg", (DL_FUNC) &_cpbayes_rcompoisreg, 4},
    {"_cpbayes_proposal_normal", (DL_FUNC) &_cpbayes_proposal_normal, 3},
    {"_cpbayes_update_beta_mu", (DL_FUNC) &_cpbayes_update_beta_mu, 8},
    {"_cpbayes_rcompois", (DL_FUNC) &_cpbayes_rcompois, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cpbayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
