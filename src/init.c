#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _CPDARMABC_binSeg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_distLoglikRatio(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_distPars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_invYeoJohnson(SEXP, SEXP);
extern SEXP _CPDARMABC_loglikFoward(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_OmegaMat(SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_OptLambdaCritBisec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_optModel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_simARIMA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_SldWin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_YeoJohnson(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_CPDARMABC_binSeg",             (DL_FUNC) &_CPDARMABC_binSeg,             22},
    {"_CPDARMABC_distLoglikRatio",    (DL_FUNC) &_CPDARMABC_distLoglikRatio,    18},
    {"_CPDARMABC_distPars",           (DL_FUNC) &_CPDARMABC_distPars,           17},
    {"_CPDARMABC_invYeoJohnson",      (DL_FUNC) &_CPDARMABC_invYeoJohnson,       2},
    {"_CPDARMABC_loglikFoward",       (DL_FUNC) &_CPDARMABC_loglikFoward,       11},
    {"_CPDARMABC_OmegaMat",           (DL_FUNC) &_CPDARMABC_OmegaMat,            4},
    {"_CPDARMABC_OptLambdaCritBisec", (DL_FUNC) &_CPDARMABC_OptLambdaCritBisec, 14},
    {"_CPDARMABC_optModel",           (DL_FUNC) &_CPDARMABC_optModel,           13},
    {"_CPDARMABC_simARIMA",           (DL_FUNC) &_CPDARMABC_simARIMA,            7},
    {"_CPDARMABC_SldWin",             (DL_FUNC) &_CPDARMABC_SldWin,             23},
    {"_CPDARMABC_YeoJohnson",         (DL_FUNC) &_CPDARMABC_YeoJohnson,          2},
    {NULL, NULL, 0}
};

void R_init_CPDARMABC(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}