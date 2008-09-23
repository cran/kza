#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include "kzft.h"


int kzft_fft(Rcomplex *z, int m, int inv)
{
    int maxf, maxp;
    double *work;
    int *iwork;
    
    fft_factor(m, &maxf, &maxp);
    if (maxf == 0) {
        error("fft factorization error");
        return NILSXP;
    }
    work = (double*)R_alloc(4 * maxf, sizeof(double));
    iwork = (int*)R_alloc(maxp, sizeof(int));
    //&(COMPLEX(z)[0].r)
    fft_work(&(z[0].r), &(z[0].i), 1, m, 1, inv, work, iwork);
    Free(work);
    Free(iwork);
         
    return 0;
}    


SEXP R_kzft_fft(SEXP z, SEXP r_m, SEXP r_k, SEXP r_inv)
{
    SEXP ans;
    int m = INTEGER_VALUE(r_m);
    int k = INTEGER_VALUE(r_k);
    
    switch (TYPEOF(z)) {
    case INTSXP:
    case LGLSXP:
    case REALSXP:
	z = coerceVector(z, CPLXSXP);
	break;
    case CPLXSXP:
	if (NAMED(z)) z = duplicate(z);
	break;
    default:
	error("non-numeric argument");
    }
    PROTECT(z);

    /* -2 for forward transform, complex values */
    /* +2 for backward transform, complex values */

    int inv = asLogical(r_inv);
    if (inv)
	    inv = -2;
    else
	    inv = 2;

    PROTECT( ans = allocMatrix(CPLXSXP, m, length(z)-m+1)) ;
	if (ans == NULL) {
	    error("Warning: Not enough memory.");
	    return NILSXP;
	}

    Rcomplex *p_ans=COMPLEX(ans);
    Rcomplex *p_z=COMPLEX(z);

    int maxf, maxp;
    double *work;
    int *iwork;
    
    fft_factor(m, &maxf, &maxp);
    if (maxf == 0) {
        error("fft factorization error");
        return NILSXP;
    }
    work = (double*)R_alloc(4 * maxf, sizeof(double));
    iwork = (int*)R_alloc(maxp, sizeof(int));
    
    int T=length(z)-(m-1);
    for (int i=0; i<T; i++) {
        Memcpy((p_ans+i*m),(p_z+i),m);
        fft_work(&(p_ans[i*m].r), &(p_ans[i*m].i), 1, m, 1, inv, work, iwork);
    }
   
    UNPROTECT(2);
    return ans;
}

//assume i is 1 based index (starts at 1 not zero)
double coeff(int i, int m, int k)
{
    double c;
    if (k == 1) return (double) 1.0/m;
    if (k == 2) {
        int j = (i > m) ? 2*m-i : i;
        return (double) j/(m*m);
    } else {
        double n=k-2+i;
        double r=k-1;
        c=choose(n,r)/pow(m,k);
    }
    return c;    
}

SEXP R_coeff(SEXP r_index, SEXP r_m, SEXP r_k)
{
    SEXP ans;
    int i=INTEGER_VALUE(r_index);
    int m=INTEGER_VALUE(r_m);
    int k=INTEGER_VALUE(r_k);
    struct coeff p;
    
    p.m=m; p.k=k;
    double a = coeff(i,m,k);
    
    PROTECT(ans = NEW_NUMERIC(1));
    NUMERIC_POINTER(ans)[0] = a;
    UNPROTECT(1);
    return (ans);
}


void test_kzft(double *x, int nrx, int ncx,
		     int nry, int ncy, Rcomplex *z, const int m, const int r, double *p)
{
    int i, j, k;
    double xij_r, xij_i, yjk_r, yjk_i;
    double sum_i, sum_r;

    for (i = 0; i < nrx; i++)
	for (k = 0; k < ncy; k++) {
	    z[i + k * nrx].r = NA_REAL;
	    z[i + k * nrx].i = NA_REAL;
	    sum_r = 0.0;
	    sum_i = 0.0;
	    for (j = 0; j < ncx; j++) {
		    xij_r = x[i + j * nrx];
		    xij_i = 0; //x[i + j * nrx].i;
		    double c1=p[j];
		    double c=coeff(j+1,m,r);
		    yjk_r = cos((2*PI*(k+1)/ncy) * (-j)) * coeff(j+1,m,r);
		    yjk_i = sin((2*PI*(k+1)/ncy) * (-j)) * coeff(j+1,m,r);
		    if (ISNAN(xij_r) || ISNAN(xij_i)
		        || ISNAN(yjk_r) || ISNAN(yjk_i))
		    goto next_ik;
		    sum_r += (xij_r * yjk_r - xij_i * yjk_i);
		    sum_i += (xij_r * yjk_i + xij_i * yjk_r);
	    }
	    z[i + k * nrx].r = sum_r;
	    z[i + k * nrx].i = sum_i;
	next_ik:
	    ;
	}
}

void kzft(double *x, int nrx, int ncx,
		     int nry, int ncy, Rcomplex *z, double *coeff)
{
    int i, j, k;
    double xij_r, xij_i, yjk_r, yjk_i;
    double sum_i, sum_r;
	
    for (i = 0; i < nrx; i++)
	for (k = 0; k < ncy; k++) {
	    z[i + k * nrx].r = NA_REAL;
	    z[i + k * nrx].i = NA_REAL;
	    sum_r = 0.0;
	    sum_i = 0.0;
	    for (j = 0; j < ncx; j++) {
		    xij_r = x[i + j * nrx];
		    xij_i = 0; //x[i + j * nrx].i;
		    yjk_r = cos((2*PI*(k+1)/ncy) * (-j)) * coeff[j];
		    yjk_i = sin((2*PI*(k+1)/ncy) * (-j)) * coeff[j];
		    if (ISNAN(xij_r) || ISNAN(xij_i)
		        || ISNAN(yjk_r) || ISNAN(yjk_i))
		    goto next_ik;
		    sum_r += (xij_r * yjk_r - xij_i * yjk_i);
		    sum_i += (xij_r * yjk_i + xij_i * yjk_r);
	    }
	    z[i + k * nrx].r = sum_r;
	    z[i + k * nrx].i = sum_i;
	next_ik:
	    ;
	}
}

SEXP R_kzft(SEXP x, SEXP r_coeff, SEXP r_m, SEXP r_k, SEXP n, SEXP p)
{
    SEXP ans;
    int N, M, L, T;
    int w;

    N = LENGTH(x);
    M = (INTEGER_VALUE(r_m)-1)*INTEGER_VALUE(r_k)+1;
    L = rint(round(M * REAL(p)[0]));
    T = floor((N-M)/L)+1;
    w = rint(REAL(n)[0]*INTEGER_VALUE(r_m));
 
    double *a = REAL(x);
    PROTECT( ans = allocMatrix(CPLXSXP, w, T)) ;
	if (ans == NULL) {
	    error("Warning: Not enough memory.");
	    return NILSXP;
	}

    double *coeff = REAL(r_coeff);
    Rcomplex *c = COMPLEX(ans);
    struct coeff poly;
    int m=INTEGER_VALUE(r_m);
    int k=INTEGER_VALUE(r_k);
    for (int i=0; i<T; i++) {
        kzft(a+(i*L), 1, M, M, w, c+(i*w), coeff);
//        test_kzft(a+(i*L), 1, M, M, w, c+(i*w), m, k, coeff);
    }

    UNPROTECT(1);
    return(ans);
}
