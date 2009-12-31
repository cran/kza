
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

SEXP bc_kzft(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP kzftwz(SEXP, SEXP, SEXP, SEXP);
SEXP big_kzftw(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP R_kzp_fftw(SEXP, SEXP, SEXP);
SEXP R_kzp_fftw_1k(SEXP, SEXP, SEXP);
SEXP R_kzftwzf_1d(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP check_fftw(void);

double coeff(int, int, int);

/*
 k is the power
 m is the number of coefficients
*/
struct coeff {
    int m, k;
};

void kzft(double *, int, int, int, int, Rcomplex *, double *);

