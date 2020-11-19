// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// OmegaMat
arma::mat OmegaMat(const int& n, const arma::colvec& order, const arma::colvec& phi, const arma::colvec& theta);
RcppExport SEXP _CPDARMABC_OmegaMat(SEXP nSEXP, SEXP orderSEXP, SEXP phiSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(OmegaMat(n, order, phi, theta));
    return rcpp_result_gen;
END_RCPP
}
// simARIMA
arma::colvec simARIMA(const double& mu, const double& sigma2, const Rcpp::String& innovDist, const arma::mat& OmegaMat, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta);
RcppExport SEXP _CPDARMABC_simARIMA(SEXP muSEXP, SEXP sigma2SEXP, SEXP innovDistSEXP, SEXP OmegaMatSEXP, SEXP BetaFlgSEXP, SEXP WSEXP, SEXP BetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type innovDist(innovDistSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type OmegaMat(OmegaMatSEXP);
    Rcpp::traits::input_parameter< const int& >::type BetaFlg(BetaFlgSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Beta(BetaSEXP);
    rcpp_result_gen = Rcpp::wrap(simARIMA(mu, sigma2, innovDist, OmegaMat, BetaFlg, W, Beta));
    return rcpp_result_gen;
END_RCPP
}
// YeoJohnson
arma::colvec YeoJohnson(const arma::colvec& Y, const double& lambda);
RcppExport SEXP _CPDARMABC_YeoJohnson(SEXP YSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(YeoJohnson(Y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// invYeoJohnson
arma::colvec invYeoJohnson(const arma::colvec& X, const double& lambda);
RcppExport SEXP _CPDARMABC_invYeoJohnson(SEXP XSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(invYeoJohnson(X, lambda));
    return rcpp_result_gen;
END_RCPP
}
// distPars
Rcpp::List distPars(const int& n, const Rcpp::List& model0, const double& lowerLambda, const double& upperLambda, const int& breakPoint, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, const int& maxIter, const Rcpp::String& innovDist, const int& nsim, const int& maxErr);
RcppExport SEXP _CPDARMABC_distPars(SEXP nSEXP, SEXP model0SEXP, SEXP lowerLambdaSEXP, SEXP upperLambdaSEXP, SEXP breakPointSEXP, SEXP BetaFlgSEXP, SEXP WSEXP, SEXP BetaSEXP, SEXP critSEXP, SEXP max_pSEXP, SEXP max_dSEXP, SEXP max_qSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP innovDistSEXP, SEXP nsimSEXP, SEXP maxErrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model0(model0SEXP);
    Rcpp::traits::input_parameter< const double& >::type lowerLambda(lowerLambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type upperLambda(upperLambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type breakPoint(breakPointSEXP);
    Rcpp::traits::input_parameter< const int& >::type BetaFlg(BetaFlgSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type crit(critSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_p(max_pSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_q(max_qSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type innovDist(innovDistSEXP);
    Rcpp::traits::input_parameter< const int& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxErr(maxErrSEXP);
    rcpp_result_gen = Rcpp::wrap(distPars(n, model0, lowerLambda, upperLambda, breakPoint, BetaFlg, W, Beta, crit, max_p, max_d, max_q, tol, maxIter, innovDist, nsim, maxErr));
    return rcpp_result_gen;
END_RCPP
}
// distLoglikRatio
arma::colvec distLoglikRatio(const int& n, const int& t, const Rcpp::List& distPars0, const double& lowerLambda, const double& upperLambda, const int& breakPoint, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, const int& maxIter, const Rcpp::String& innovDist, const int& nsim, const int& maxErr);
RcppExport SEXP _CPDARMABC_distLoglikRatio(SEXP nSEXP, SEXP tSEXP, SEXP distPars0SEXP, SEXP lowerLambdaSEXP, SEXP upperLambdaSEXP, SEXP breakPointSEXP, SEXP BetaFlgSEXP, SEXP WSEXP, SEXP BetaSEXP, SEXP critSEXP, SEXP max_pSEXP, SEXP max_dSEXP, SEXP max_qSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP innovDistSEXP, SEXP nsimSEXP, SEXP maxErrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type distPars0(distPars0SEXP);
    Rcpp::traits::input_parameter< const double& >::type lowerLambda(lowerLambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type upperLambda(upperLambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type breakPoint(breakPointSEXP);
    Rcpp::traits::input_parameter< const int& >::type BetaFlg(BetaFlgSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type crit(critSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_p(max_pSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_q(max_qSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type innovDist(innovDistSEXP);
    Rcpp::traits::input_parameter< const int& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxErr(maxErrSEXP);
    rcpp_result_gen = Rcpp::wrap(distLoglikRatio(n, t, distPars0, lowerLambda, upperLambda, breakPoint, BetaFlg, W, Beta, crit, max_p, max_d, max_q, tol, maxIter, innovDist, nsim, maxErr));
    return rcpp_result_gen;
END_RCPP
}
// binSeg
Rcpp::List binSeg(const arma::colvec& Y, const double& alpha, const int& GLRSApprox, const int& minSize, const double& lowerLambda, const double& upperLambda, const int& breakPoint, const int& BetaFlg, const arma::mat& W, const arma::colvec& Beta, const Rcpp::String& crit, const int& max_p, const int& max_d, const int& max_q, const double& tol, const int& maxIter, const double& maxLev, const Rcpp::String& innovDist, const int& nsim1, const int& nsim2, const int& maxErr, const int& minNonNAN);
RcppExport SEXP _CPDARMABC_binSeg(SEXP YSEXP, SEXP alphaSEXP, SEXP GLRSApproxSEXP, SEXP minSizeSEXP, SEXP lowerLambdaSEXP, SEXP upperLambdaSEXP, SEXP breakPointSEXP, SEXP BetaFlgSEXP, SEXP WSEXP, SEXP BetaSEXP, SEXP critSEXP, SEXP max_pSEXP, SEXP max_dSEXP, SEXP max_qSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP maxLevSEXP, SEXP innovDistSEXP, SEXP nsim1SEXP, SEXP nsim2SEXP, SEXP maxErrSEXP, SEXP minNonNANSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type GLRSApprox(GLRSApproxSEXP);
    Rcpp::traits::input_parameter< const int& >::type minSize(minSizeSEXP);
    Rcpp::traits::input_parameter< const double& >::type lowerLambda(lowerLambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type upperLambda(upperLambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type breakPoint(breakPointSEXP);
    Rcpp::traits::input_parameter< const int& >::type BetaFlg(BetaFlgSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type crit(critSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_p(max_pSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_q(max_qSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double& >::type maxLev(maxLevSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type innovDist(innovDistSEXP);
    Rcpp::traits::input_parameter< const int& >::type nsim1(nsim1SEXP);
    Rcpp::traits::input_parameter< const int& >::type nsim2(nsim2SEXP);
    Rcpp::traits::input_parameter< const int& >::type maxErr(maxErrSEXP);
    Rcpp::traits::input_parameter< const int& >::type minNonNAN(minNonNANSEXP);
    rcpp_result_gen = Rcpp::wrap(binSeg(Y, alpha, GLRSApprox, minSize, lowerLambda, upperLambda, breakPoint, BetaFlg, W, Beta, crit, max_p, max_d, max_q, tol, maxIter, maxLev, innovDist, nsim1, nsim2, maxErr, minNonNAN));
    return rcpp_result_gen;
END_RCPP
}
