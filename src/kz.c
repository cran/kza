/*
    kz.c
    Copyright (C) 2011 Brian D. Close

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#include <R.h> 
#include <Rdefines.h>
#include "kz.h"

static double mavg2d(SEXP v, int row, int col, int w)
{
	double s=0.00;
	int i, j, z=0;
	int startrow, startcol, endrow, endcol;
	int nr=nrows(v);

    if (isMatrix(v)) {
        
        startrow = (row-w>0 ? row-w : 0);
        startcol = (col-w>0 ? col-w : 0);

        endrow = (row+w<nr ? row+w+1 : nrows(v));
        endcol = (col+w<ncols(v) ? col+w+1 : ncols(v));
        
	    for(i=startrow, z=0, s=0.00; i<endrow; i++) {
	        for(j=startcol; j<endcol; j++) {
	    	    if (R_FINITE(REAL(v)[i+nr*j])) { 
	    		    z++;
	    		    s += REAL(v)[i+nr*j];
	    	    }
	    	}
	    }
	} else error("Input is not a vector or Matrix.");
	if (z == 0) return R_NaN;
	return s/z;
}

static double mavg1d(SEXP v, int col, int w)
{
	double s=0.00;
	int i, z=0;
	int startcol, endcol;

	if (isVector(v)) {
        startcol = (col-w>0 ? col-w : 0);
        endcol = (col+w<LENGTH(v) ? col+w+1 : LENGTH(v));
        
	    for(i=startcol, z=0, s=0.00; i<endcol; i++) {
    	    if (R_FINITE(REAL(v)[i])) {
    		    z++;
    		    s += REAL(v)[i];
	    	}
	    }
	} else error("Input is not a vector or Matrix.");
	if (z == 0) return R_NaN;
	return s/z;
}

static double averaged(SEXP x, SEXP box_center, int width)
{
	SEXP dim;
	int ndim;
	int count=0;
	double s=0.00;
    SEXP index, box;
    int offset, offset_a;
    int i,j;
    
    dim = GET_DIM(x);
    ndim = LENGTH(dim);
    PROTECT(index = allocVector(INTSXP, LENGTH(dim)));
    PROTECT(box = allocMatrix(INTSXP, LENGTH(dim), 2));
    
    if (isArray(x)) {
        /* bounding box */
        for(i=0; i<ndim; i++) {
            /* box edges */
            INTEGER(box)[i] = (INTEGER(box_center)[i] - width>0 ? INTEGER(box_center)[i] - width : 0);
            INTEGER(box)[i+ndim] = (INTEGER(box_center)[i] + width <INTEGER(dim)[i] ? INTEGER(box_center)[i] + width : INTEGER(dim)[i]-1);
        }
        for(INTEGER(index)[0]=INTEGER(box)[0]; INTEGER(index)[0]<=INTEGER(box)[ndim]; INTEGER(index)[0]++) {
            for(INTEGER(index)[1]=INTEGER(box)[1]; INTEGER(index)[1]<=INTEGER(box)[1+ndim]; INTEGER(index)[1]++) {
              for(INTEGER(index)[2]=INTEGER(box)[2]; INTEGER(index)[2]<=INTEGER(box)[2+ndim]; INTEGER(index)[2]++) {
/*
        offset += INTEGER(index)[0];
        offset += INTEGER(dim)[0]*INTEGER(index)[1];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(index)[2];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(dim)[2]*INTEGER(index)[3]; 
*/
                    /* find offset into array */
                    for(i=0, offset=0; i<ndim; i++) {
                        for(j=0, offset_a=1; j<i; j++) {
                            offset_a *= INTEGER(dim)[j];
                        }
                        offset += offset_a * INTEGER(index)[i];
                    }
    	            if (R_FINITE(REAL(x)[offset])) {
    	                count++;
    		            s += REAL(x)[offset];
	                }
	            }
            }
        }
	} else error("Input is not a vector or Matrix.");
	
	UNPROTECT(2);
	if (count == 0) return R_NaN;
	return s/count;
}

SEXP kz1d(SEXP x, SEXP window, SEXP iterations)
{
	int p;
	int i, k;
	int m;
	SEXP ans, tmp, dim;
	int n;

	m = (2 * INTEGER_VALUE(window)) + 1; 
	p = (m-1)/2;
	
	dim = GET_DIM(x);
	n = LENGTH(x);
	
	PROTECT(tmp = allocVector(REALSXP, LENGTH(x)));
	PROTECT(ans = allocVector(REALSXP, LENGTH(x)));
    copyVector(tmp, x);

	for(k=0; k<INTEGER_VALUE(iterations); k++) {
        for(i=0;i<LENGTH(x);i++) {
            REAL(ans)[i] = mavg1d(tmp, i, p);
        }
        /* copyMatrix (destination, source, byrow) */
   	    copyVector(tmp, ans); 
	}
	UNPROTECT(2);
	return ans;
}

SEXP kz2d(SEXP x, SEXP window, SEXP iterations)
{
	int p;
	int i, j, k;
	int m;
	SEXP ans, tmp, dim;
	int nr, nc;

	m = (2 * INTEGER_VALUE(window)) + 1; 
	p = (m-1)/2;
	
	dim = GET_DIM(x);
	nr = INTEGER(dim)[0]; nc = INTEGER(dim)[1];
	
	PROTECT(ans = allocMatrix(REALSXP, nr, nc));
	PROTECT(tmp = allocMatrix(REALSXP, nr, nc));
    copyMatrix(tmp, x, 1);

	for(k=0; k<INTEGER_VALUE(iterations); k++) {
	    for(i=0;i<nr; i++) {
	        for(j=0;j<nc;j++) {
	            REAL(ans)[i+nr*j] = mavg2d(tmp, i, j, p);
	        }
        }
        /* copyMatrix (destination, source, byrow) */
   	    copyMatrix(tmp, ans, 1); 
	}
	UNPROTECT(2);
	return ans;
}

void copyArray(SEXP destination, SEXP source)
{
    int i;
    
    for(i=0; i<LENGTH(source);i++) {
        REAL(destination)[i] = REAL(source)[i];       
    }
}

SEXP kz3d(SEXP x, SEXP window, SEXP iterations)
{
	int p;
	int i, j, l;
	int m;
	SEXP ans, tmp, dim;
	SEXP index;
	int offset, offset_a;

	m = (2 * INTEGER_VALUE(window)) + 1; 
	p = (m-1)/2;
	
	dim = GET_DIM(x);
	
	PROTECT(index = allocVector(INTSXP, LENGTH(dim)));
	PROTECT(ans = allocArray(REALSXP, dim));
	PROTECT(tmp = allocArray(REALSXP, dim));
	copyArray(ans, x);

	for(l=0; l<INTEGER_VALUE(iterations); l++) {
   	    copyArray(tmp, ans); 
        for(INTEGER(index)[0]=0; INTEGER(index)[0]<INTEGER(dim)[0]; INTEGER(index)[0]++) {
            for(INTEGER(index)[1]=0; INTEGER(index)[1]<INTEGER(dim)[1]; INTEGER(index)[1]++) {
              for(INTEGER(index)[2]=0; INTEGER(index)[2]<INTEGER(dim)[2]; INTEGER(index)[2]++) {
/*
        offset += INTEGER(index)[0];
        offset += INTEGER(dim)[0]*INTEGER(index)[1];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(index)[2];
        offset += INTEGER(dim)[0]*INTEGER(dim)[1]*INTEGER(dim)[2]*INTEGER(index)[3]; 
*/
                    /* find offset into array */
                    for(i=0, offset=0; i<LENGTH(dim); i++) {
                        for(j=0, offset_a=1; j<i; j++) {
                            offset_a *= INTEGER(dim)[j];
                        }
                        offset += offset_a * INTEGER(index)[i];
                    }
                    REAL(ans)[offset] = averaged(tmp, index, p);
                }
            }
        }
	}
	UNPROTECT(3);
	return ans;
}

SEXP kz(SEXP x, SEXP window, SEXP iterations)
{
    SEXP ans=R_NilValue, dim;

    dim = getAttrib(x, R_DimSymbol);
    
    if (isArray(x) && LENGTH(dim) >= 3) {
        if (LENGTH(dim) > 3) {
            error("Too many dimensions -- not yet implemented, please contact the author for more info.");
            return(ans);
        }
        PROTECT(ans = kz3d(x, window, iterations)); 
        UNPROTECT(1);
    } else if (isMatrix(x)) {
        PROTECT(ans = kz2d(x, window, iterations)); 
        UNPROTECT(1);
    } else if (isVector(x)) {
        PROTECT(ans = kz1d(x, window, iterations)); 
        UNPROTECT(1);
    }
    return(ans);
}

// x is the space or independent variable or time variable
// y is the response variable to be averaged
// ans is (max(x)-min(x))*scale
SEXP kzs(SEXP ans, SEXP y, SEXP x, SEXP dx, SEXP q, SEXP iterations, SEXP minx, SEXP maxx)
{
	int l;

	// q must be odd integer
	l = INTEGER_VALUE(q); 
	
	// use x to define values of y to be averaged
	// initial
	int qi=INTEGER_VALUE(minx);
	for (int i=0, qj=qi+l; i<INTEGER_VALUE(maxx); i++) 
	{
		int s=0,size=0;
		for (int j=qi;j<qj;j++)
		{
			if (REAL(x)[j]<qj) 
			{
				s += REAL(y)[j]; 
				size++;
			}
			else break;
		}
		if (size) REAL(ans)[i] = s/size;
		qi += REAL(dx)[0];
	}
	
	return(ans);
}	
