#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include <R.h>

/* #define DEBBUG */
#ifdef DEBBUG
  int R_finite(double x) { return ( (x)==(x) ); }
  #define Calloc(b, t)  (t*) calloc(b,sizeof(t))
  #define Free free
  #define PRINT(x) { if ((x)==(x)) printf("%04.1f ",x); else printf("NaN "); }
#else
  #include <R.h>
  #include <Rinternals.h>
#endif

#define notNaN(x)   ((x)==(x))
#define isNaN(x)  (!((x)==(x)))
#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))
#define SQR(x) ((x)*(x))

/*============================================================================*/
/* The following macros were inspired by msum from                            */
/* http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090             */
/* Quote from it:                                                             */
/* "Full precision summation using multiple doubles for intermediate values   */
/* Rounded x+y stored in hi with the round-off stored in lo.  Together        */
/* hi+lo are exactly equal to x+y.  The loop applies hi/lo summation          */
/* to each partial so that the list of partial sums remains exact.            */
/* Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:   */
/* www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps"  */
/*============================================================================*/

/* SumErr - macro calculating error of the summing operation */
#define SumErr(a,b,ab) ((((a)>(b)) == ((a)>-(b))) ?  (b) - ((ab)-(a)) : (a) - ((ab)-(b)) )
/* SUM_1 - macro for calculating Sum+=x; Num+=n; Which is NaN aware and have minimal (single number) overflow error correction */
#define SUM_1(x,n, Sum, Err, Num)   if (R_finite(x)){ y=Sum; Err+=x; Sum+=Err; Num+=n; Err=SumErr(y,Err,Sum);  } 
#define mpartial 1024	


void SUM_N(double x, int n, double *partial, int *npartial, int *Num) 
{
  if (R_finite(x)){ 
    int j, i;
    double hi, lo, y;
    for (i=j=0; j<*npartial; j++) {
      y  = partial[j];
      hi = y + x;
      lo = SumErr(x,y,hi); 
      if (lo && i<mpartial) partial[i++] = lo;
      x = hi; 
    }
    partial[i] = x; 
    *npartial   = i+1;
    *Num+=n; 
  }
}

/*========================================================================================*/
/* Each iteration of an insertion sort removes an element from the input data, inserting  */
/*  it at the correct position in the already sorted list, until no elements are left     */
/* in the input. For unsorted data is much less efficient than the more advanced          */
/* algorithms such as Quicksort, Heapsort, or Merge sort, but it is very efficient on     */
/* data sets which are already substantially sorted (almost O(n))                         */
/* Referances:                                                                            */
/*   R. Sedgewick: Algorithms. Addison-Wesley (1988) (page 99)                            */
/*   http://en.wikipedia.org/wiki/Insertion_sort                                          */
/*   http://www.cs.ubc.ca/spider/harrison/Java/sorting-demo.html                          */
/* Input:                                                                                 */
/*   V    - data array we operate on remains unchanged (we assume that no number in idx   */
/*          array is longer than V                                                        */
/*   idx  - index numbers of elements in V to be partially sorted                         */
/*   nIdx - length of idx array                                                           */
/* Output:                                                                                */
/*   idx  - index numbers of sorted elements in V w                                       */
/*========================================================================================*/

void insertion_sort(const double *V, int *idx, const int nIdx) 
{ 
  int i, j, id;
  double v;
  for (i=1; i<nIdx; i++) {
    id = idx[i];
    v  = V[id];
    for(j=i; j>0; j--) {
      if (V[idx[j-1]]<v) break;
      idx[j] = idx[j-1];
    }
    idx[j] = id;
  }
}

/*==================================================================*/
/* Array Sum without round-off errors.                              */
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty double                                            */
/*   nIn  - size of In array                                        */
/* Output :                                                         */
/*   Out  - Array sum                                               */
/*==================================================================*/
void sum_exact(double *In, double *Out, const int *nIn)
{
  int i, j, n=*nIn, npartial=0, Num=0;
  double *in, x, partial[mpartial];
  in=In;
  for(i=0; i<n; i++, in++) {
    x = *in;
    SUM_N(x, 1, partial, &npartial, &Num);
  }
  *Out = partial[0];
  for(j=1; j<npartial; j++) *Out += partial[j];
}

/*===================================================================================*/
/* Mean function applied to (running) window. All additions performed using          */
/* addition algorithm which tracks and corrects addition round-off errors (see       */  
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps) */
/* Input :                                                                           */
/*   In   - array to run moving window over will remain unchanged                    */
/*   Out  - empty space for array to store the results                               */
/*   nIn  - size of arrays In and Out                                                */
/*   nWin - size of the moving window                                                */
/* Output :                                                                          */
/*   Out  - results of runing moving window over array In and collecting window mean */
/*===================================================================================*/
void runmean_exact(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation, full round-off correction*/
  int i, j, k2, n=*nIn, m=*nWin, npartial=0, Num=0;
  double *in, *out, partial[mpartial], Sum;
  double NaN = (0.0/0.0);

  k2 = m>>1;         /* right half of window size */
  in=In; out=Out; 
  /* step 1 - find mean of elements 0:(k2-1) */      
  for(i=0; i<k2; i++) {
    SUM_N(in[i], 1, partial, &npartial, &Num);
  }
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m; i++, out++) {
    SUM_N(in[i], 1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
  /* step 3: runsum of inner section. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m); but with round of error correction */
  for(i=m; i<n; i++, out++, in++) {
    SUM_N(in[m] , 1, partial, &npartial, &Num);
    SUM_N(-(*in),-1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++, out++, in++) {
    SUM_N(-(*in),-1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
}






/*==================================================================*/
/* 																	*/
/* Input :                                                          */
/*   outData  - empty double                                        */
/*   outDim  - output dimensions									*/
/*   inData - vector of data values									*/
/*   nIn - size of input vector										*/
/*   index - array of index locations								*/
/*   nDim - number of dimensions 									*/
/*	 inDim - array with size of each dimension of the index			*/
/*	 window - array with window filter dimensions					*/
/* Output :                                                         */
/*   out  - array of mean values                                    */
/*==================================================================*/
void kz_ma(double *outData, double *outDim,
	double *inData, const int *nIn, 
	double *index, const int *nInDim, double *inDim, 
	double *window)
{
	int nd;
	int n; // n is the size of input
	int i,j,k, zi;
	int offset, z_offset;
	int outSize;
	int x,y,z;
	int wl;
	
	double *data;
	double *df;
	double *w;
	double *out;
	double *od;
	double *id;
	double *cnts;
	double NaN = (0.0/0.0);
	
	out = outData;
	od = outDim;
	data = inData;
	n = *nIn;
	df = index;
	nd = *nInDim;
	id = inDim;
	
	w = window;
	
	for (i=0, outSize=1; i<nd; i++) {
		outSize *= outDim[i];
	}
    cnts = calloc(outSize, sizeof(double));
	if (cnts == NULL) { error("Unable to allocate space.\n"); return; }
	
	wl=0;
	for (i=0,x=0,y=0,z=0; i<n; i++) {
		y=df[i]-1;
		x=df[(int)(i+id[0])]-1;
		if (nd>2) z=df[(int)(i+id[0]*(id[1]-1))]-1;
		// given coordinates x,y,z find offset into output (zero based)
		
		for(zi=MAX(z-w[2],0); zi<=MIN(z+w[2],od[2]-1); zi++) {
			z_offset = MAX(zi*od[0]*od[1],0);
			offset = z_offset;
			for (j=MAX(x-w[0],0); j<=MIN(x+w[0],od[1]-1); j++) {
				wl = MIN(j*od[0]+y+w[1],(j+1)*od[0]-1) + z_offset;
				offset = z_offset;
				offset += MAX(j*od[0],0);  
				offset += MAX(y-w[1],0);
				for (k=offset; k<=wl; k++) {
					out[k] += data[i];
					cnts[k]++;
				}
			}
		}
	}
	
	for (i=0; i<outSize; i++) {out[i] = (cnts[i] ? out[i]/cnts[i] : NaN);}
	
	free(cnts);
}
