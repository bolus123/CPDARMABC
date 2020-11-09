#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _CPDARMABC_BCBD(SEXP, SEXP);
extern SEXP _CPDARMABC_binSeg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_distLoglikRatio(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_distPars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_fastLm(SEXP, SEXP);
extern SEXP _CPDARMABC_invBCBD(SEXP, SEXP);
extern SEXP _CPDARMABC_OmegaMat(SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_simARIMA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _CPDARMABC_simCPDARIMABC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_CPDARMABC_BCBD",            (DL_FUNC) &_CPDARMABC_BCBD,             2},
    {"_CPDARMABC_binSeg",          (DL_FUNC) &_CPDARMABC_binSeg,          20},
    {"_CPDARMABC_distLoglikRatio", (DL_FUNC) &_CPDARMABC_distLoglikRatio, 18},
    {"_CPDARMABC_distPars",        (DL_FUNC) &_CPDARMABC_distPars,        17},
    {"_CPDARMABC_fastLm",          (DL_FUNC) &_CPDARMABC_fastLm,           2},
    {"_CPDARMABC_invBCBD",         (DL_FUNC) &_CPDARMABC_invBCBD,          2},
    {"_CPDARMABC_OmegaMat",        (DL_FUNC) &_CPDARMABC_OmegaMat,         4},
    {"_CPDARMABC_simARIMA",        (DL_FUNC) &_CPDARMABC_simARIMA,         7},
    {"_CPDARMABC_simCPDARIMABC",   (DL_FUNC) &_CPDARMABC_simCPDARIMABC,   20},
    {NULL, NULL, 0}
};

void R_init_CPDARMABC(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}