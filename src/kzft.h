
#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif
#ifndef ABS
#define ABS(x)  ((x<0)?(-x):(x))
#endif

#define HAVE_FORTRAN_DOUBLE_COMPLEX 1

void cmatprod(Rcomplex *, int, int,
                Rcomplex *, int, int, Rcomplex *);

void F77_NAME (sgemm) (char *, char *, int *, int *, int *, float *, float *,
		        int *, float *, int *, float *, float *, int *);

void F77_NAME(zgemm)(const char *transa, const char *transb, const int *m,
		        const int *n, const int *k, const Rcomplex *alpha,
		        const Rcomplex *a, const int *lda,
		        const Rcomplex *b, const int *ldb,
		        const Rcomplex *beta, Rcomplex *c, const int *ldc);

SEXP R_kzft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_kztp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_kztp_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP b_kzft(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP b_kztp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_smooth_kzp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_bsmooth_kzp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_smooth_kztp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_bsmooth_kztp(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_kzw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP R_coeff(SEXP, SEXP, SEXP);
SEXP R_kzft_fft(SEXP, SEXP, SEXP, SEXP);
SEXP R_kzp_energy(SEXP, SEXP);

SEXP bkzp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


/*
 k is the power
 m is the number of coefficients
*/
struct coeff {
    int m, k;
};

void kzft(double *, int, int, int, int, Rcomplex *, double *);

