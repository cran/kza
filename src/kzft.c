#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include "kzft.h"
#include "config.h"

#ifdef HAVE_FFTW
#include "fftw3.h"
#endif

SEXP check_fftw()
{
	SEXP ans;
	PROTECT (ans = NEW_LOGICAL(1));

	#ifdef FFTW
	LOGICAL_DATA(ans)[0] = 1;
	#else
	LOGICAL_DATA(ans)[0] = 0;
	#endif
	UNPROTECT(1);
	return ans;
}

#ifdef HAVE_FFTW
//for data that is C99 complex array
int kzfftw(double *data, int *t, int len, int nr, int nc, double *results)
{
  	fftw_plan plan;
  	fftw_complex *in, *out;
  	double *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  	int index_set_start=0;
	int index_set;
	p = results;
	for (int i=0; i<nr; i++) {
		int size=0;
		memset(in, 0, sizeof(fftw_complex) * len);
		index_set=index_set_start;
		for (int j = i; j < len+i && index_set<nr; j++) {
			if (t[index_set] > (j+1)) continue;
			in[j-i][0] = data[2*index_set];
			in[j-i][1] = data[2*index_set+1];
			index_set++;
			size++;
			if (j==i) {index_set_start++;}
		}

		fftw_execute(plan); /* repeat as needed */
  
		for (int j = 0; j < len; j++) {
			p[2*(i+j*nr)] = out[j][0] / size;
			p[2*(i+j*nr)+1] = out[j][1] / size;
		}
	}

  	fftw_destroy_plan(plan);
	fftw_free(in); fftw_free(out);
	
	return 0;
}

SEXP kzftw(SEXP x, SEXP index_set, SEXP M, SEXP ans)
{
	PROTECT(x = AS_NUMERIC(x));
	PROTECT(ans = AS_NUMERIC(ans));

	double *a = REAL(x);
	double *p = REAL(ans);

	int m=INTEGER_VALUE(M);

	SEXP dims = getAttrib(ans, R_DimSymbol); 
	int ncols = INTEGER(dims)[1]; 
	int nfaces = INTEGER(dims)[2]; 

	int *t = INTEGER(index_set);

	kzfftw(a, t, m, ncols, nfaces, p);

	UNPROTECT(2);
	return(ans);
}

int fftwz(Rcomplex *data, int *t, int len, int nr, int nc, Rcomplex *results)
{
  	fftw_plan plan_forward, plan_backward;
  	fftw_complex *in, *out;
	Rcomplex *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan_forward = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(len, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  	int index_set_start=0;
	int index_set;
	p = results;
	for (int i=0; i<nr; i++) {
		int size=0;
		memset(in, 0, sizeof(fftw_complex) * len);
		index_set=index_set_start;
		for (int j = i; j < len+i && index_set<nr; j++) {
			if (ISNAN(t[index_set])) continue;
			if (t[index_set] > (j+1)) continue;
			in[j-i][0] = data[index_set].r;
			in[j-i][1] = data[index_set].i;
			index_set++;
			size++;
			if (j==i) {index_set_start++;}
		}

		fftw_execute(plan_forward);
		for (int j = 0; j < len; j++) {
			p[i+j*nr].r = out[j][0] / size;
			p[i+j*nr].i = out[j][1] / size;
		}
	}
	
  	fftw_destroy_plan(plan_backward);
  	fftw_destroy_plan(plan_forward);
	fftw_free(in); fftw_free(out);
	
	return 0;
}

SEXP kzftwz(SEXP z, SEXP index_set, SEXP M, SEXP ans)
{
	PROTECT(z = AS_COMPLEX(z));
	PROTECT(ans = AS_COMPLEX(ans));

	SEXP dims = getAttrib(ans, R_DimSymbol); 
	int nrows = INTEGER(dims)[0];
	int ncols = INTEGER(dims)[1];

	Rcomplex *a = COMPLEX(z);
	Rcomplex *p = COMPLEX(ans);
	int *t = INTEGER(index_set);

	int m=INTEGER_VALUE(M);

	fftwz(a, t, m, nrows, ncols, p);

	UNPROTECT(2);
	return (ans);
}

int kzp_fftw_1k(double *data, int len, int nr, double *results)
{
  	fftw_plan plan_forward, plan_backward;
  	fftw_complex *in;
  	fftw_complex *out_forward, *out_backward;
	double *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out_forward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out_backward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan_forward = fftw_plan_dft_1d(len, in, out_forward, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(len, in, out_backward, FFTW_BACKWARD, FFTW_ESTIMATE);

	p = results;
	for (int i=0; i<nr-len; i++) {
		memset(in, 0, sizeof(fftw_complex) * len);
		for (int j = i; j < len+i; j++) {
			in[j-i][0] = data[j];
		}

		fftw_execute(plan_forward);
		fftw_execute(plan_backward);

		for (int j = 0; j < len; j++) {
			p[j] += fabs(out_forward[j][0]*out_backward[j][0] + out_forward[j][1]*out_backward[j][1])/(len*len*len);
		}
	}
	
  	fftw_destroy_plan(plan_backward);
  	fftw_destroy_plan(plan_forward);
	fftw_free(in); 
	fftw_free(out_forward); 
	fftw_free(out_backward);
	
	return 0;
}

int kzp_fftwz_1k(Rcomplex *data, int len, int nr, double *results)
{
  	fftw_plan plan_forward, plan_backward;
  	fftw_complex *in;
  	fftw_complex *out_forward, *out_backward;
	double *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out_forward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out_backward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan_forward = fftw_plan_dft_1d(len, in, out_forward, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(len, in, out_backward, FFTW_BACKWARD, FFTW_ESTIMATE);

	p = results;
	for (int i=0; i<nr-len; i++) {
		memset(in, 0, sizeof(fftw_complex) * len);
		for (int j = i; j < len+i; j++) {
			in[j-i][0] = data[j].r;
			in[j-i][1] = data[j].i;
		}

		fftw_execute(plan_forward);
		fftw_execute(plan_backward);

		for (int j = 0; j < len; j++) {
			p[j] += fabs(out_forward[j][0]*out_forward[j][0] - out_forward[j][1]*out_forward[j][1])/(len*len*len);
		}
	}
	
  	fftw_destroy_plan(plan_backward);
  	fftw_destroy_plan(plan_forward);
	fftw_free(in); 
	fftw_free(out_forward); 
	fftw_free(out_backward);
	
	return 0;
}

SEXP R_kzp_fftw(SEXP y, SEXP M, SEXP ans)
{
	Rcomplex *p=COMPLEX(y);
	double *a=REAL(ans);

	int n=LENGTH(y);
	int m=INTEGER_VALUE(M);
	
	kzp_fftwz_1k(p, m, n, a);

	return (ans);
}

SEXP R_kzp_fftw_1k(SEXP y, SEXP M, SEXP ans)
{
   	double *p=REAL(y);
	double *a=REAL(ans);

	int n=LENGTH(y);
	int m=INTEGER_VALUE(M);
	
	kzp_fftw_1k(p, m, n, a);

	return (ans);
}

int fftwz_1d(Rcomplex *data, int *t, int len, int nr, int freq, Rcomplex *results)
{
  	fftw_plan plan_forward, plan_backward;
  	fftw_complex *in, *out;
	Rcomplex *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan_forward = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(len, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  	int index_set_start=0;
	int index_set;
	p = results;
	for (int i=0; i<nr; i++) {
		int size=0;
		memset(in, 0, sizeof(fftw_complex) * len);
		index_set=index_set_start;
		for (int j = i; j < len+i && index_set<nr; j++) {
			if (ISNAN(t[index_set])) continue;
			if (t[index_set] > (j+1)) continue;
			in[j-i][0] = data[index_set].r;
			in[j-i][1] = data[index_set].i;
			index_set++;
			size++;
			if (j==i) {index_set_start++;}
		}

		fftw_execute(plan_forward);
		for (int j = 0, k=len-1; j < len; j++, k--) {
			p[i+j*nr].r = out[j][0] / size;
			p[i+j*nr].i = out[j][1] / size;
		}
	}
	
  	fftw_destroy_plan(plan_backward);
  	fftw_destroy_plan(plan_forward);
	fftw_free(in); fftw_free(out);
	
	return 0;
}

int fftwz_f(Rcomplex *data, int *t, int len, int nr, int freq, Rcomplex *results)
{
  	fftw_plan plan_forward, plan_backward;
  	fftw_complex *in, *out;
	Rcomplex *p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);

	plan_forward = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_backward = fftw_plan_dft_1d(len, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  	int index_set_start=0;
	int index_set;
	p = results;
	for (int i=0; i<nr-len; i++) {
		int size=0;
		memset(in, 0, sizeof(fftw_complex) * len);
		index_set=index_set_start;
		for (int j = i; j < len+i && index_set<nr; j++) {
			if (ISNAN(t[index_set])) continue;
			if (t[index_set] > (j+1)) continue;
			in[j-i][0] = data[index_set].r;
			in[j-i][1] = data[index_set].i;
			index_set++;
			size++;
			if (j==i) {index_set_start++;}
		}

		fftw_execute(plan_forward);

		p[i].r = out[freq][0] / size;
		p[i].i = out[freq][1] / size;
	}
	
  	fftw_destroy_plan(plan_backward);
  	fftw_destroy_plan(plan_forward);
	fftw_free(in); fftw_free(out);
	
	return 0;
}

SEXP R_kzftwzf_1d(SEXP z, SEXP index_set, SEXP F, SEXP M, SEXP ans)
{
	int n = LENGTH(z);

	PROTECT(z = AS_COMPLEX(z));
	PROTECT(ans = AS_COMPLEX(ans));

	Rcomplex *a = COMPLEX(z);
	Rcomplex *p = COMPLEX(ans);
	int *t = INTEGER(index_set);

	int m=INTEGER_VALUE(M);
	int f=INTEGER_VALUE(F);	

	fftwz_f(a, t, m, n, f, p);

	UNPROTECT(2);
	return (ans);
}

#endif //HAVE_FFTW
