#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "kzft.h"
#include "kz.h"
#include "config.h"

#include <R_ext/Rdynload.h>

 R_CallMethodDef callMethods[] = {
    {"kz", (DL_FUNC) &kz, 3},
    {"kza", (DL_FUNC) &kza, 6},
    {"kzp_energy", (DL_FUNC) &R_kzp_energy, 2},
    {"kztp_smooth", (DL_FUNC) &R_smooth_kztp, 5},
    {"kzp_smooth", (DL_FUNC) &R_smooth_kzp, 5},
    {"smooth_kzp", (DL_FUNC) &R_smooth_kzp, 5},
    {"smooth_kztp", (DL_FUNC) &R_smooth_kztp, 5},
    {"kzs", (DL_FUNC) &kzs, 8},
#ifdef HAVE_FFTW
	{"kzftwz", (DL_FUNC) &kzftwz, 4},
	{"kzp", (DL_FUNC) &R_kzp_fftw, 3},
	{"kzp_1k", (DL_FUNC) &R_kzp_fftw_1k, 3},
#endif
	{"fftw", (DL_FUNC) &check_fftw, 0},
   {NULL, NULL, 0}
};

void R_init_kza(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
