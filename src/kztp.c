#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include "kzft.h"

static double variance_squared(double *pg, int rows, int cols, int row_start, int row_stop, int col_start, int col_stop)
{
    double sum1=0;
    double sum2=0;

    // col mode
    for (int j=col_start; j<col_stop; j++) {
        for (int i=row_start; i<row_stop-1; i++) {
            sum1 += pow(pg[(i+1)+j*rows]-pg[i+j*rows],2);
        }
    }

    // row mode
    for (int i=row_start; i<row_stop; i++) {
        for (int j=col_start; j<col_stop-1; j++) {
            sum2 += pow(pg[i+(j+1)*rows]-pg[i+j*rows],2);
        }
    }

    return sum1+sum2;    
}

static double mean(double *x, int N, int start_row, int end_row, int start_col, int end_col)
{
    double sum=0;
    int cnt=0;
    
    for (int i=start_row; i<end_row; i++) {
        for (int j=start_col; j<end_col; j++) {
            sum += x[i+j*N];
            cnt++;
        }
    }

    return sum/cnt;
}

static void variation(double *pg, int N, int K, double c, double *spg)
{
    double total = c*variance_squared(pg,N,N,0,N,0,N);
    double sq;
    int cnt;
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            cnt = 0;
            for (int k=1; k<K; k++) {
                //pg[(max(1,(i-k+1)):min(N,(i+k-1))),(max(1,(j-k+1)):min(N,(j+k-1)))]
                int max1=MAX(0,(i-k));
                int min1=MIN(N,(i+k+1));
                int max2=MAX(0,(j-k));
                int min2=MIN(N,(j+k+1));
                if (max1<=min1 && max2<=min2) {
                    sq = variance_squared(pg,N,N,max1,min1,max2,min2);
                    if (sq <= total) cnt++;
                }
            }
            int max1=MAX(0,(i-cnt));
            int min1=MIN(N,(i+cnt+1));
            int max2=MAX(0,(j-cnt));
            int min2=MIN(N,(j+cnt+1));
            spg[i+j*N] = mean(pg,N,max1,min1,max2,min2);
        }
    }
    return;
}

//#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
SEXP R_smooth_kztp(SEXP R_pg, SEXP R_c, SEXP R_K, SEXP R_spg, SEXP method)
{
    int nrow=INTEGER(GET_DIM(R_pg))[0];
    int ncol=INTEGER(GET_DIM(R_pg))[1];
    int K=INTEGER_VALUE(R_K);
    double c=REAL(R_c)[0];
    double *spg=REAL(R_spg);
    int N=nrow;
    double *pg=REAL(R_pg);

    variation(pg, N, K, c, spg);
    
    return R_spg;    
}

static Rcomplex cmult(Rcomplex x1, Rcomplex x2)
{
    Rcomplex ans;

	ans.r = x1.r * x2.r - x1.i * x2.i;
	ans.i = x1.r * x2.i + x1.i * x2.r;
	return (ans);
}

static Rcomplex Conj(Rcomplex x) {
    Rcomplex y;
    y.r = x.r;
	y.i = -x.i;

    return y;
}

SEXP R_kztp(SEXP x, SEXP coeff, SEXP m, SEXP Rk, SEXP p, SEXP n, 
            SEXP rp1, SEXP rp2, SEXP cp1, SEXP cp2)
{
    SEXP ans;
    int N, M, L, T;
    int w;
    
    N = LENGTH(x);
    M = (INTEGER_VALUE(m)-1)*INTEGER_VALUE(Rk)+1;
    L = rint(round(M * REAL(p)[0]));
    T = floor((N-M)/L)+1;
    w = rint(REAL(n)[0]*INTEGER_VALUE(m));

    double *a = REAL(x);
    Rcomplex *c = calloc(w*T, sizeof(Rcomplex));
	if (c == NULL) {
	    error("Warning: Not enough memory.");
	    return NILSXP;
	}
    int rm1=MAX(rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(rp1)[0])),1);
    int rm2=round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(rp2)[0]);
    int cm1=MAX(rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(cp1)[0])),1);
    int cm2=rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(cp2)[0]));
    int delta_rm=rm2-rm1+1;
    int delta_cm=cm2-cm1+1;

    PROTECT( ans = allocMatrix(CPLXSXP, delta_rm, delta_cm));
    if (ans == NULL) {
        free(c);
        return NILSXP;
    }

	double *kc = REAL(coeff);
    for (int i=0; i<T; i++) {
        kzft(a+(i*L), 1, M, M, w, c+(i*w), kc);
    }

    Rcomplex *kztp=COMPLEX(ans);
    Rcomplex tmp;
    for (int i=0; i<delta_rm; i++) {
        for (int j=0; j<delta_cm; j++) {
            kztp[i+j*delta_rm].r = 0;
            kztp[i+j*delta_rm].i = 0;
            for (int k=0; k<T; k++) {
                tmp=cmult(c[(i+rm1-1)+k*w],c[(j+cm1-1)+k*w]);
                tmp=cmult(tmp,Conj(c[(i+j+rm1+cm1-1)+k*w]));
                kztp[i+j*delta_rm].r += tmp.r*M*M;
                kztp[i+j*delta_rm].i += tmp.i*M*M;
            }
        }
    }
         
    for (int i=0; i<delta_rm*delta_cm; i++) {
        kztp[i].r = kztp[i].r/T;
        kztp[i].i = kztp[i].i/T;
    }

    free(c);
    UNPROTECT(1);
    
    return ( ans );
}

SEXP R_kztp_test(SEXP x, SEXP m, SEXP Rk, SEXP p, SEXP n, 
            SEXP rp1, SEXP rp2, SEXP cp1, SEXP cp2)
{
    SEXP ans;
    int M, nc, nr;
    int w;
    Rcomplex *z;
    Rcomplex tmp;
    
    M = (INTEGER_VALUE(m)-1)*INTEGER_VALUE(Rk)+1;

	nc = ncols(x);
	nr = nrows(x);
	
    switch (TYPEOF(x)) {
    case INTSXP:
    case LGLSXP:
    case REALSXP:
	x = coerceVector(x, CPLXSXP);
	break;
    case CPLXSXP:
	if (NAMED(x)) x = duplicate(x);
	break;
    default:
	error("non-numeric argument");
    }
    PROTECT(x);
    
    Rcomplex *c = COMPLEX(x);
    
    int rm1=MAX(rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(rp1)[0])),1);
    int rm2=round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(rp2)[0]);
    int cm1=MAX(rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(cp1)[0])),1);
    int cm2=rint(round(INTEGER_VALUE(m)*INTEGER_VALUE(n)*REAL(cp2)[0]));
    int delta_rm=rm2-rm1+1;
    int delta_cm=cm2-cm1+1;

    PROTECT( ans = allocMatrix(CPLXSXP, delta_rm, delta_cm));
    if (ans == NULL) {
	    error("Warning: Not enough memory.");
        return NILSXP;
    }

    z=COMPLEX(ans);
    for (int i=0; i<delta_rm; i++) {
        for (int j=0; j<delta_cm; j++) {
            z[i+j*delta_rm].r = 0;
            z[i+j*delta_rm].i = 0;
            for (int k=0; k<nr; k++) {
                tmp=cmult(c[(i+rm1-1)+k*w],c[(j+cm1-1)+k*nc]);
                tmp=cmult(tmp,Conj(c[(i+j+rm1+cm1-1)+k*nc]));
                z[i+j*delta_rm].r += tmp.r*M*M;
                z[i+j*delta_rm].i += tmp.i*M*M;
            }
        }
    }
         
    for (int i=0; i<delta_rm*delta_cm; i++) {
        z[i].r = z[i].r/nr;
        z[i].i = z[i].i/nr;
    }

    UNPROTECT(2);
    return ( ans );
}
