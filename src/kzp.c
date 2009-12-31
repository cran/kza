#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include "kzft.h"

static double sumdiff(double *x, int start, int end)
{
    double a=0;

    for (int k=start; k<end; k++) {
        a += pow(x[k+1]-x[k],2);
    }
    return a;
}

static double mean(double *x, int start, int end)
{
    double a=0;
    int d=0;
    
    for (int i=start; i<end; i++) {
        a += x[i];
        d++;
    }
    if (d>0) a = a/d;
    
    return a;
}

static int variation(double *pg, double *spg, int N, int iK, double c) 
{
    double total=c*sumdiff(pg,0,N-1);

    // row mode
    for (int i=0; i<N; i++) { //ith row
        int m=0;
        for (int j=1; j<=iK; j++) { //jth col
            int max=MAX(0,(i-j)); 
            int min=MIN(N,(i+j));
            m += (sumdiff(pg, max, min) <= total) ? 1 : 0;
        }
        int max = MAX(0,(i-m));
        int min = MIN(N,(i+m+1));
        spg[i] = mean(pg,max,min);
    }

    return 1;
}

static double abssum(double *x, int start, int end)
{
    double a=0;

    for (int i=start; i<end; i++) {
        double b=x[i+1]-2*x[i]+x[i-1];
        a += ABS(b);
    }
        
    return a;
}
    
static int nonlinear(double *pg, double *spg, int N, int iK, double c) 
{
    double total=c*abssum(pg,1,N-1);

    // row mode
    for (int i=0; i<N; i++) { //ith row
        int m=0;
        for (int j=1; j<iK; j++) { //jth col
            int max=MAX(1,(i-j)); 
            int min=MIN(N-1,(i+j+1));
            m += (abssum(pg, max, min) <= total) ? 1 : 0;
        }
        int max = MAX(0,(i-m));
        int min = MIN(N,(i+m+1));
        spg[i] = mean(pg,max,min);
    }

    return 1;
}

SEXP R_smooth_kzp(SEXP r_pg, SEXP R_c, SEXP K, SEXP r_spg, SEXP method)
{
    double c=REAL(R_c)[0];
    int N=LENGTH(r_pg);
    int iK=INTEGER_VALUE(K);
    double *spg = REAL(r_spg);
    double *p_pg = REAL(r_pg);

    if (strcmp(CHAR(STRING_ELT(method, 0)),"DZ")==0) {
        if (variation(p_pg, spg, N, iK, c) == 0) return NILSXP;
    } else {
        if (nonlinear(p_pg, spg, N, iK, c) == 0) return NILSXP;
    }
    return (r_spg);
}

SEXP R_kzp_energy(SEXP z, SEXP M)
{
    SEXP r_ans;

    Rcomplex *p_z=COMPLEX(z);
	int r = nrows(z);
	int c = ncols(z);

    PROTECT( r_ans = allocVector(REALSXP, r)) ;
	if (r_ans == NULL) {
	    error("Warning: Not enough memory.");
	    return NILSXP;
	}
	double *ans = REAL(r_ans);
    for (int i=0; i<r; i++) {
        ans[i] = 0;
    }
    
    for (int i=0; i<c; i++) {
        for (int j=0; j<r; j++) {
            ans[j] += pow(hypot(p_z[j+(i*r)].r, p_z[j+(i*r)].i),2)*INTEGER_VALUE(M);
        }
    }
    for (int i=0; i<r; i++) {
        ans[i] /= c;
    }
	UNPROTECT(1);
	return(r_ans);
}


