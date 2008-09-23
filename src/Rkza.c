#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "kzft.h"
#include "kz.h"

#include <R_ext/Rdynload.h>

 R_CallMethodDef callMethods[] = {
    {"kz", (DL_FUNC) &kz, 3},
    {"kza", (DL_FUNC) &kza, 6},
    {"kzp_energy", (DL_FUNC) &R_kzp_energy, 2},
    {"kztp_smooth", (DL_FUNC) &R_smooth_kztp, 5},
    {"kzp_smooth", (DL_FUNC) &R_smooth_kzp, 5},
    {"smooth_kzp", (DL_FUNC) &R_smooth_kzp, 5},
    {"smooth_kztp", (DL_FUNC) &R_smooth_kztp, 5},
   {NULL, NULL, 0}
};

void R_init_kza(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
