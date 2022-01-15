#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _testrcpp_G_thres_cpp(SEXP, SEXP);
extern SEXP _testrcpp_fun_dev(SEXP, SEXP);
extern SEXP _testrcpp_fun_mul(SEXP, SEXP);
extern SEXP _testrcpp_loglikelihood_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _testrcpp_post_e_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _testrcpp_post_tau_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _testrcpp_rcpparma_bothproducts(SEXP);
extern SEXP _testrcpp_rcpparma_hello_world();
extern SEXP _testrcpp_rcpparma_innerproduct(SEXP);
extern SEXP _testrcpp_rcpparma_outerproduct(SEXP);
extern SEXP _testrcpp_sample_c_l_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _testrcpp_sample_gibbs_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _testrcpp_sample_thres_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_testrcpp_G_thres_cpp",           (DL_FUNC) &_testrcpp_G_thres_cpp,            2},
    {"_testrcpp_fun_dev",               (DL_FUNC) &_testrcpp_fun_dev,                2},
    {"_testrcpp_fun_mul",               (DL_FUNC) &_testrcpp_fun_mul,                2},
    {"_testrcpp_loglikelihood_cpp",     (DL_FUNC) &_testrcpp_loglikelihood_cpp,     10},
    {"_testrcpp_post_e_cpp",            (DL_FUNC) &_testrcpp_post_e_cpp,            13},
    {"_testrcpp_post_tau_cpp",          (DL_FUNC) &_testrcpp_post_tau_cpp,          12},
    {"_testrcpp_rcpparma_bothproducts", (DL_FUNC) &_testrcpp_rcpparma_bothproducts,  1},
    {"_testrcpp_rcpparma_hello_world",  (DL_FUNC) &_testrcpp_rcpparma_hello_world,   0},
    {"_testrcpp_rcpparma_innerproduct", (DL_FUNC) &_testrcpp_rcpparma_innerproduct,  1},
    {"_testrcpp_rcpparma_outerproduct", (DL_FUNC) &_testrcpp_rcpparma_outerproduct,  1},
    {"_testrcpp_sample_c_l_cpp",        (DL_FUNC) &_testrcpp_sample_c_l_cpp,        14},
    {"_testrcpp_sample_gibbs_cpp",      (DL_FUNC) &_testrcpp_sample_gibbs_cpp,      20},
    {"_testrcpp_sample_thres_cpp",      (DL_FUNC) &_testrcpp_sample_thres_cpp,      12},
    {NULL, NULL, 0}
};

void R_init_testrcpp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
