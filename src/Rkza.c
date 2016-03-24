#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "kz.h"

#include <R_ext/Rdynload.h>

 R_CallMethodDef callMethods[] = {
    {"kz", (DL_FUNC) &kz, 3},
    {"kza", (DL_FUNC) &kza, 6},
    {"kzs", (DL_FUNC) &kzs, 8},
   {NULL, NULL, 0}
};

void R_init_kza(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
