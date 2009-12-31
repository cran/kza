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
    int K=INTEGER_VALUE(R_K);
    double c=REAL(R_c)[0];
    double *spg=REAL(R_spg);
    int N=nrow;
    double *pg=REAL(R_pg);

    variation(pg, N, K, c, spg);
    
    return R_spg;    
}
