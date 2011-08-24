#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <memory.h>
#include <math.h>
#include <float.h>

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
	if (cnts == NULL) { printf("Unable to allocate space.\n"); return; }
	
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
